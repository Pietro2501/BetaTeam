import pickle
import gzip
import random
import plotly.express as px
import multiprocessing as mp
from collections import Counter, defaultdict
import time

###############################################################################
# CLASSI PER NODI, ECC.
###############################################################################
class Node:
    """Classe associata al nodo, utile a definire in seguito il grafo di de Brujin"""

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
    """Costruisce nodi di de Bruijn (senza edges, qui si accumulano solo indegree/outdegree)."""
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
    """Unisce i dizionari di nodi dai risultati dei processi."""
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
    """Costruisce i nodi di de Bruijn usando multiprocessing."""
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
    """Grafico di distribuzione della numerosità dei k-mer."""
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
    """Ritorna la lista dei nodi con (indegree + outdegree) > 2."""
    _candidates = [n_id for n_id, obj in _nodes.items() if (obj.indegree + obj.outdegree) > 2]
    return _candidates

###############################################################################
# PICK RANDOM HUB
###############################################################################
def pick_random_branching_node(_candidates, _filtered_candidates=None):
    """Seleziona un nodo a caso tra i 'branching', escludendo quelli già filtrati."""
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
    """
    Costruisce due strutture di adiacenza:
      - adjacency_right[prefix] = lista di (suffix, coverage)
      - adjacency_left[suffix]  = lista di (prefix, coverage)
    così da evitare di scansionare tutto il dict ad ogni passo.
    """
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
def extend_right(node_id, adjacency_right, adjacency_left, dict_kmer_count,
                 _candidates, _filtered_candidates, min_coverage=1):
    """
    Estende il contig verso destra partendo da node_id,
    scegliendo il successore con coverage più alto dagli adjacency.
    """
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

        # kmer effettivo = primo nucleotide di current_node + best_next_node
        best_kmer = current_node[0] + best_next_node

        # "Consumo" il k-mer rimuovendolo dalle strutture:
        # 1) Rimuovo dal dict_kmer_count
        if best_kmer in dict_kmer_count:
            del dict_kmer_count[best_kmer]
        # 2) Rimuovo dalla lista adjacency_right[current_node]
        adjacency_right[current_node] = [(s, c) for (s, c) in successors if s != best_next_node]
        # 3) Rimuovo da adjacency_left[best_next_node] il current_node
        adjacency_left[best_next_node] = [
            (p, c) for (p, c) in adjacency_left[best_next_node] if p != current_node
        ]

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
def extend_left(node_id, adjacency_right, adjacency_left, dict_kmer_count,
                _candidates, _filtered_candidates, min_coverage=1):
    """
    Estende il contig a sinistra partendo da node_id,
    scegliendo il predecessore con coverage più alto dagli adjacency.
    """
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

        # kmer effettivo = best_prev_node + ultimo nucleotide di current_node
        best_kmer = best_prev_node + current_node[-1]

        # "Consumo" il k-mer rimuovendolo dalle strutture:
        if best_kmer in dict_kmer_count:
            del dict_kmer_count[best_kmer]

        adjacency_left[current_node] = [
            (p, c) for (p, c) in predecessors if p != best_prev_node
        ]
        adjacency_right[best_prev_node] = [
            (s, c) for (s, c) in adjacency_right[best_prev_node] if s != current_node
        ]

        # Se best_prev_node è un hub e non è ancora filtrato, aggiungilo
        if best_prev_node in _candidates and best_prev_node not in _filtered_candidates:
            _filtered_candidates.append(best_prev_node)

        # Aggiungo il primo nucleotide di best_prev_node in coda a contig_ext
        contig_ext.append(best_prev_node[0])

        current_node = best_prev_node
        # print(f"Sto percorrendo l'arco a destra e aggiungengo il nodo {current_node}")



    # Ricorda di invertire la lista contig_ext costruita al contrario
    contig_ext.reverse()
    return "".join(contig_ext), _filtered_candidates

###############################################################################
# COSTRUZIONE DI UN CONTIG A PARTIRE DA UN HUB
###############################################################################
def build_local_contig(dict_kmer_count,
                       adjacency_right,
                       adjacency_left,
                       _candidates,
                       _filtered_candidates,
                       min_coverage=1):
    """
    1) Sceglie un nodo 'branching' a caso,
    2) Estende a sinistra e a destra,
    3) Restituisce l'intero contig (stringa nucleotidica) e la lista aggiornata.
    """
    start_node = pick_random_branching_node(_candidates, _filtered_candidates)
    if not start_node:
        # Nessun nodo con (in+out) > 2 => potremmo restituire stringa vuota
        return "", _filtered_candidates

    # Estensione a sinistra
    contig_left, filt_candidates_left = extend_left(
        start_node,
        adjacency_right,
        adjacency_left,
        dict_kmer_count,
        _candidates,
        _filtered_candidates,
        min_coverage=min_coverage
    )

    # Estensione a destra (partendo dal medesimo start_node)
    contig_right, new_filtered_candidates = extend_right(
        start_node,
        adjacency_right,
        adjacency_left,
        dict_kmer_count,
        _candidates,
        filt_candidates_left,
        min_coverage=min_coverage
    )

    # Contig completo
    contig = contig_left + start_node + contig_right
    return contig, new_filtered_candidates

###############################################################################
# MAIN
###############################################################################
if __name__ == '__main__':
    start = time.time()

    # Caricamento del dizionario {kmer -> coverage}
    with gzip.open('kmer_17_NEW_PhoeVul_num1filtered.pkl.gz', 'rb') as f:
        dict_kmer_count = pickle.load(f)

    print(f"Lavoreremo con {len(dict_kmer_count)} k-mer")
    distribuzione_kmer(dict_kmer_count)

    # Costruzione dei nodi del grafo (in parallelo)
    nodes = get_nodes_parallel(dict_kmer_count)
    print(f"Lavoreremo con {len(nodes)} nodi")

    # Costruzione strutture di adiacenza
    adjacency_right, adjacency_left = build_adjacency_structures(dict_kmer_count)

    # Troviamo i cosiddetti 'hub'
    candidates = count_hub(nodes)
    print(f"Gli hub nel grafo saranno {len(candidates)}")

    filtered_candidates = []
    counter = 1

    # Finché non abbiamo "filtrato" tutti gli hub
    while len(filtered_candidates) < len(candidates):
        # Costruisce un contig a partire da un hub casuale
        contig, filtered_candidates = build_local_contig(
            dict_kmer_count,
            adjacency_right,
            adjacency_left,
            candidates,
            filtered_candidates,
            min_coverage=7
        )
        print(f"Il {counter}° contig generato misura {len(contig)} basi")
        print(f"Hub filtrati finora: {len(filtered_candidates)} / {len(candidates)}")
        counter += 1

    end = time.time()
    print(f"Il processo ha richiesto {(end - start) / 3600:.2f} ore")
