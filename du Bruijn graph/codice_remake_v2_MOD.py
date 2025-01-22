import pickle
import gzip
import random
import plotly.express as px
import multiprocessing as mp
from collections import Counter, defaultdict
import time

from ipywidgets import interact


###############################################################################
# CLASSI PER NODI, ECC.
###############################################################################
class Node:


    def __init__(self, lab):
        self.label = lab
        self.indegree = 0
        self.outdegree = 0

    def __eq__(self, other):
        return self.label == other.label

    def incrementa_indegree(self, value=1):
        self.indegree += value

    def incrementa_outdegree(self, value=1):
        self.outdegree += value

###############################################################################
# FUNZIONI PER IL GRAFO (COSTRUZIONE NODI IN PARALLELO)
###############################################################################
def partial_nodes(_chunk):

    local_nodes = {}
    for key in _chunk:
        prefix = key[:-1]
        suffix = key[1:]

        if prefix not in local_nodes:
            local_nodes[prefix] = Node(prefix)
        local_nodes[prefix].incrementa_outdegree()

        if suffix not in local_nodes:
            local_nodes[suffix] = Node(suffix)
        local_nodes[suffix].incrementa_indegree()

    return local_nodes

def merge_nodes(results):

    merged_nodes = {}
    for partial_n in results:
        for key, node in partial_n.items():
            if key not in merged_nodes:
                merged_nodes[key] = node
            else:
                merged_nodes[key].incrementa_indegree(node.indegree)
                merged_nodes[key].incrementa_outdegree(node.outdegree)
    return merged_nodes

def get_nodes_parallel(_dict_kmer_count, num_processes=6):

    # Divide il dizionario in chunks
    keys = list(_dict_kmer_count.keys())
    chunk_size = len(keys) // (num_processes + 1) if (num_processes + 1) > 0 else len(keys)
    chunks = [keys[i:i + chunk_size] for i in range(0, len(keys), chunk_size)]

    # Crea un pool di processi
    with mp.Pool(processes=num_processes) as pool:
        # Applica la funzione su ciascun chunk
        results = pool.map(partial_nodes, chunks)

    # Unisci i risultati
    _nodes = merge_nodes(results)
    return _nodes

###############################################################################
# DISTRIBUZIONE K-MER (GRAFICO)
###############################################################################
def distribuzione_kmer(_dict_kmer_count):

    valori_numerici = list(_dict_kmer_count.values())
    contatore = Counter(valori_numerici)
    max_numerosita = max(valori_numerici)
    kmer_val = list(range(0, max_numerosita + 1))  # Da 0 al massimo valore di numerosità
    frequenze = [contatore.get(n, 0) for n in kmer_val]
    fig = px.bar(
        x=kmer_val,
        y=frequenze,
        title='Distribuzione delle Numerosità dei k-mer',
        labels={'x': 'Numerosità', 'y': 'Frequenza'},
        text=frequenze
    )
    # Aggiorno layout ed estetica
    fig.update_layout(
        xaxis=dict(tickmode='linear'),
        yaxis=dict(title='Frequenza'),
        title=dict(x=0.5, xanchor='center'),  # Centra il titolo
        uniformtext_minsize=8,
        uniformtext_mode='hide'
    )
    fig.update_traces(texttemplate='%{text}', textposition='outside')
    fig.show()

###############################################################################
# CONTA HUB
###############################################################################
def count_hub(_nodes):

    _candidates = [n_id for n_id, obj in _nodes.items() if (obj.indegree + obj.outdegree) > 2]
    return _candidates

###############################################################################
# PICK RANDOM HUB
###############################################################################
def pick_random_branching_node(_candidates, _filtered_candidates=None):

    if not _candidates:
        return None
    if _filtered_candidates is None:
        _filtered_candidates = []
    current_candidates = [item for item in _candidates if item not in _filtered_candidates]
    if not current_candidates:
        return None
    return random.choice(current_candidates)

###############################################################################
# STRUTTURE DI ADIACENZA (PER EVITARE SCANSIONI TOTALI)
###############################################################################
def build_adjacency_structures(dict_kmer_count):

    adjacency_right = defaultdict(list)
    adjacency_left = defaultdict(list)

    for kmer, coverage in dict_kmer_count.items():
        prefix = kmer[:-1]
        suffix = kmer[1:]
        adjacency_right[prefix].append((suffix, coverage))
        adjacency_left[suffix].append((prefix, coverage))

    return adjacency_right, adjacency_left

###############################################################################
# EXTEND_RIGHT
###############################################################################
def extend_right(node_id, adjacency_right, adjacency_left,_candidates, _filtered_candidates, min_coverage=1):

    contig_ext = []
    current_node = node_id

    while True:
        successors = adjacency_right[current_node]
        if not successors:
            # Nessuna estensione possibile
            break

        # Ordina i successori per coverage decrescente
        successors.sort(key=lambda x: x[1], reverse=True)
        best_next_node, best_cov = successors[0]

        # Se la copertura è troppo bassa, interrompi
        if best_cov < min_coverage:
            break

        # 1) Rimuovo dalla lista adjacency_right[current_node]
        adjacency_right[current_node] = [(s, c) for (s, c) in successors if s != best_next_node]
        # 2) Rimuovo da adjacency_left[best_next_node] il current_node
        adjacency_left[best_next_node] = [(p, c) for (p, c) in adjacency_left[best_next_node] if p != current_node]

        # Se best_next_node è un hub e non è ancora filtrato, aggiungilo
        if best_next_node in _candidates and best_next_node not in _filtered_candidates:
            _filtered_candidates.append(best_next_node)

        # Aggiungo al contig l'ultimo nucleotide di best_next_node
        contig_ext.append(best_next_node[-1])

        # Avanzo al nodo successivo
        current_node = best_next_node
        # print(f"Sto percorrendo l'arco a destra e aggiungengo il nodo {current_node}")

    return "".join(contig_ext), _filtered_candidates

###############################################################################
# EXTEND_LEFT
###############################################################################
def extend_left(node_id, adjacency_right, adjacency_left,_candidates, _filtered_candidates, min_coverage=1):

    contig_ext = []
    current_node = node_id

    while True:
        predecessors = adjacency_left[current_node]
        if not predecessors:
            break

        # Ordina i predecessori per coverage
        predecessors.sort(key=lambda x: x[1], reverse=True)
        best_prev_node, best_cov = predecessors[0]

        if best_cov < min_coverage:
            break

        adjacency_left[current_node] = [(p, c) for (p, c) in predecessors if p != best_prev_node]
        adjacency_right[best_prev_node] = [(s, c) for (s, c) in adjacency_right[best_prev_node] if s != current_node]

        # Se best_prev_node è un hub e non è ancora filtrato, aggiungilo
        if best_prev_node in _candidates and best_prev_node not in _filtered_candidates:
            _filtered_candidates.append(best_prev_node)

        # Aggiungo il primo nucleotide di best_prev_node in coda a contig_ext
        contig_ext.append(best_prev_node[0])

        current_node = best_prev_node

    # Ricorda di invertire la lista contig_ext costruita al contrario
    contig_ext.reverse()
    return "".join(contig_ext), _filtered_candidates

###############################################################################
# COSTRUZIONE DI UN CONTIG A PARTIRE DA UN HUB
###############################################################################
def build_local_contig(adjacency_right,adjacency_left,_candidates,_filtered_candidates,min_coverage=1):

    start_node = pick_random_branching_node(_candidates, _filtered_candidates)
    if not start_node:
        # Nessun nodo con (in+out) > 2 => potremmo restituire stringa vuota
        return "", _filtered_candidates

    # Estensione a sinistra
    contig_left, filt_candidates_left = extend_left(start_node,adjacency_right,adjacency_left,_candidates,_filtered_candidates,min_coverage=min_coverage)
    # Estensione a destra (partendo dal medesimo start_node)
    contig_right, new_filtered_candidates = extend_right(start_node,adjacency_right,adjacency_left,_candidates,filt_candidates_left,min_coverage=min_coverage)

    # Contig completo
    contig = contig_left + start_node + contig_right
    return contig, new_filtered_candidates

###############################################################################
# CALCOLO DELL'N50 SULLA BASE DELLE LUNGHEZZE DEI CONTIG PRODOTTI
###############################################################################
def calcola_n50(lunghezze_contig:list):
    try:
        # Ordina le lunghezze dei contig in ordine decrescente
        lunghezze_contig.sort(reverse=True)

        # Calcola la lunghezza totale dei contig
        lunghezza_totale = sum(lunghezze_contig)

        # Somma le lunghezze dei contig finché non si raggiunge almeno il 50% della lunghezza totale
        somma_parziale = 0
        for contig in lunghezze_contig:
            somma_parziale += contig
            if somma_parziale >= lunghezza_totale / 2:
                return contig  # Restituisce la lunghezza del contig che raggiunge o supera il 50%
    except IndexError:
        print("Lista di contig vuota o dati non validi...")

###############################################################################
# COSTRUZIONE DEI CONTIG SULLA BASE DELLE NECESSITA' DELL'UTENTE
###############################################################################
def iterative_contig_generation(_candidates, _start, _min_coverage=7):
    filtered_candidates = []
    counter = 1
    lenghts_contig = []

    # Chiede all'utente il tipo di iterazione
    while True:
        try:
            iter_choice = input("\nQuanto profonda vuoi sia la lettura del grafo?\n\t- '1' per ricoprire tutti i nodi hub del grafo"
                                " (operazione computazionalmente onerosa)\n\t- '2' per un numero di contig da te desiderato\n")
            if iter_choice not in ['1', '2']:
                raise ValueError
            break
        except ValueError:
            print("Attenzione a cosa digiti... Inserisci '1' o '2' per scegliere il tipo di iterazione.")

    # Se l'utente sceglie iterazione definita, chiede il numero massimo di iterazioni
    max_iterations = None
    if iter_choice == '2':
        while True:
            try:
                max_iterations = int(input("Inserisci il numero massimo di contig (intero positivo):\n"))
                if max_iterations <= 0:
                    raise ValueError
                break
            except ValueError:
                print("Devi inserire un valore intero positivo, non caratteri a casaccio!!!")

    # Esegui il loop in base alla scelta dell'utente
    if iter_choice == '1':
        # Iterazione basata sulla condizione while
        while len(filtered_candidates) < len(_candidates):
            # Costruzione strutture di adiacenza
            _adjacency_right, _adjacency_left = build_adjacency_structures(dict_kmer_count)
            contig, filtered_candidates = build_local_contig(_adjacency_right,_adjacency_left,_candidates,filtered_candidates,_min_coverage)
            print(f"Il {counter}° contig generato misura {len(contig)} basi")
            print(f"Hub filtrati finora: {len(filtered_candidates)} / {len(_candidates)}")
            lenghts_contig.append(len(contig))
            counter += 1
    else:
        # Iterazione fissa
        while counter <= max_iterations:
            # Costruzione strutture di adiacenza
            _adjacency_right, _adjacency_left = build_adjacency_structures(dict_kmer_count)
            contig, filtered_candidates = build_local_contig(_adjacency_right,_adjacency_left,_candidates,filtered_candidates,_min_coverage)
            print(f"Il {counter}° contig generato misura {len(contig)} basi")
            print(f"Hub filtrati finora: {len(filtered_candidates)} / {len(_candidates)}")
            lenghts_contig.append(len(contig))
            counter += 1

    print(f"L'N50 ottenuto dai contig generati è pari a {calcola_n50(lenghts_contig)}")
    end = time.time()
    print(f"Il processo ha richiesto {int((end - _start) // 3600)} ore e {int(((end - _start)/3600 - (end - _start)//3600)*60)} minuti")

def intestazione_progetto():
    print("\n{:^80}\n".format("Progetto: Grafo di de Bruijn"))
    print("Partendo da read paired end di 150pb, ottenute da un sequenziamento Illumina")
    print("condotto su genoma di Phocaeicola vulgatus, proveremo a realizzare dei contig.")
    print("Al termine di questa operazione restituiremo il valore N50 calcolato per")
    print("verificare la bont\u00e0 nella ricostruzione della sequenza.\n")
    print("I dati iniziali sono stati trattati con lo scopo di ottenere dei k-mer di")
    print("lunghezza pari a 17nt, calcolarne la rappresentativit\u00e0 nelle read e infine")
    print("costruire dei nodi prefix e suffix con cui percorrere il grafo.\n")

###############################################################################
# MAIN
###############################################################################
if __name__ == '__main__':
    start = time.time()
    intestazione_progetto()

    # Caricamento del dizionario {kmer -> coverage}
    with gzip.open('kmer_17_NEW_PhoeVul_num1filtered.pkl.gz', 'rb') as f:
        dict_kmer_count = pickle.load(f)

    print(f"Lavoreremo con {len(dict_kmer_count)} k-mer")
    distribuzione_kmer(dict_kmer_count)

    # Costruzione dei nodi del grafo (in parallelo)
    nodes = get_nodes_parallel(dict_kmer_count)
    print(f"Lavoreremo con {len(nodes)} nodi")

    # Troviamo i cosiddetti 'hub'
    candidates = count_hub(nodes)
    print(f"Gli hub nel grafo saranno {len(candidates)}")

    # Produco contig ed estrapolo l'N50 sulla base dei limiti che mi inpone l'utente
    iterative_contig_generation(candidates, start)

