import pickle
import gzip
import random
import plotly.express as px
import multiprocessing as mp
from collections import Counter
import time

from pip._internal.resolution.resolvelib import candidates


class Node:
    """ Classe associata al nodo, utile a definire in seguito il grafo di de Brujin"""

    def __init__(self, lab):
        self.label = lab
        self.indegree= 0
        self.outdegree = 0

    def __eq__(self, other):
        return self.label == other.label

    def incrementa_indegree(self, value=1):
        self.indegree += value

    def incrementa_outdegree(self, value=1):
        self.outdegree += value

# class Edge:
#     """ Classe associata all'edge, utile a definire in seguito il grafo di de Brujin"""
#
#     def __init__(self, kmer, counter):
#         self.label = kmer
#         self.contatore = counter
#
#     def __eq__(self, other):
#         return self.label == other.label
#
#     def incrementa_contatore(self):
#         self.contatore += 1

with gzip.open('kmer_17_NEW_PhoeVul_num1filtered.pkl.gz', 'rb') as f:
    dict_kmer_count = pickle.load(f)

def partial_nodes(_chunk):
    """Costruisco il grafo di de Bruijn ottimizzando la ricerca con un dizionario"""
    # edges = {}
    local_nodes = {}

    for key in _chunk:
        prefix = key[:-1]
        suffix = key[1:]

        if prefix not in local_nodes:
            local_nodes[prefix] = Node(prefix)
            local_nodes[prefix].incrementa_outdegree()
        else:
            local_nodes[prefix].incrementa_outdegree()

        if suffix not in local_nodes:
            local_nodes[suffix] = Node(suffix)
            local_nodes[prefix].incrementa_indegree()
        else:
            local_nodes[prefix].incrementa_indegree()
    return local_nodes

def merge_nodes(results):
    """Unisce i dizionari di nodi dai risultati dei processi."""
    merged_nodes = {}
    for partial_nodes in results:
        for key, node in partial_nodes.items():
            if key not in merged_nodes:
                merged_nodes[key] = node
            else:
                merged_nodes[key].incrementa_indegree(node.indegree)
                merged_nodes[key].incrementa_outdegree(node.outdegree)
    return merged_nodes

def get_nodes_parallel(_dict_kmer_count, num_processes=6):
    """Costruisce il grafo di de Bruijn usando multiprocessing."""
    # Divide il dizionario in chunks
    keys = list(_dict_kmer_count.keys())
    chunk_size = len(keys) // (num_processes+1)
    chunks = [keys[i:i + chunk_size] for i in range(0, len(keys), chunk_size)]

    # Crea un pool di processi
    with mp.Pool(processes=num_processes) as pool:
        # Applica la funzione su ciascun chunk
        results = pool.map(partial_nodes, chunks)

    # Unisci i risultati
    _nodes = merge_nodes(results)
    return _nodes

def distribuzione_kmer(_dict_kmer_count):
    valori_numerici = list(_dict_kmer_count.values())
    contatore = Counter(valori_numerici)
    max_numerosita = max(valori_numerici)
    kmer_val = list(range(0, max_numerosita + 1)) # Da 0 al massimo valore di numerosità
    frequenze = [contatore.get(n, 0) for n in kmer_val]
    fig = px.bar(
        x=kmer_val,
        y=frequenze,
        title='Distribuzione delle Numerosità dei k-mer',
        labels={'x': 'Numerosità', 'y': 'Frequenza'},
        text=frequenze,  # Aggiunge le frequenze come testo sopra le barre
        # color=frequenze,  # Colora le barre in base alla frequenza
        # color_continuous_scale='Greens'  # Palette di colori
    )
    #Aggiornare il layout per migliorare l'estetica
    fig.update_layout(
        xaxis=dict(tickmode='linear'),
        yaxis=dict(title='Frequenza'),
        title=dict(x=0.5, xanchor='center'),  # Centra il titolo
        uniformtext_minsize=8,
        uniformtext_mode='hide'
    )
    #Aggiornare le tracce per posizionare il testo fuori dalle barre
    fig.update_traces(texttemplate='%{text}', textposition='outside')
    #Mostrare il grafico
    fig.show()


###############################################################################
# 1. Scelta di un nodo "branching" in modo casuale
###############################################################################
def pick_random_branching_node(_candidates, _filtered_candidates=None):
    """
    Seleziona un nodo a caso tra quelli con (indegree + outdegree) > 2.
    nodes: dict { node_id (str) : Node(...) }
    Restituisce: la chiave (node_id) di un nodo valido, o None se non ce n'è.
    """
    if not _filtered_candidates:
        random_node = random.choice(_candidates)
        if not _candidates:
            return None
    else:
        current_candidates = [item for item in _candidates if item not in _filtered_candidates]
        random_node = random.choice(current_candidates)
    return random_node

def count_hub(_nodes):
    _candidates = [n_id for n_id, obj in _nodes.items() if (obj.indegree + obj.outdegree) > 2]
    return _candidates

###############################################################################
# 2. Estensione a destra (extend_right)
###############################################################################
def extend_right(node_id, _dict_kmer_count, _candidates, _filtered_candidates, min_coverage=1):
    """
    Estende il contig a destra partendo da un (k-1)-mer (node_id),
    scegliendo il k-mer con coverage più alto tra quelli che iniziano con node_id.

    Parametri:
    - node_id: str, il (k-1)-mer di partenza
    - dict_kmer_count: dict { kmer (str) : coverage (int) }
    - k: lunghezza del k-mer
    - min_coverage: soglia minima di coverage, per evitare percorsi ciechi

    Restituisce:
    - Una stringa contenente i nucleotidi aggiunti (non include node_id iniziale).
    """
    contig_ext = []
    current_node = node_id

    while True:
        successors = []

        # Scansiona tutti i k-mer per trovare quelli che iniziano con current_node
        for kmer, coverage in dict_kmer_count.items():
            # Esempio: se current_node="ATG", cerchiamo kmer come "ATGC", "ATGA", ...
            if kmer.startswith(current_node):
                next_node = kmer[1:] # il (k-1)-mer successivo
                successors.append((next_node, coverage))

        if not successors:
            # Nessun k-mer inizia con current_node => estensione finita
            break

        # Ordina i successori per coverage in modo decrescente
        successors.sort(key=lambda x: x[1], reverse=True)

        # Prende il successore con coverage più alto
        best_next_node, best_cov = successors[0]
        del dict_kmer_count[current_node[0]+best_next_node] #FATTO ALLA FINE E DA VERIFICAREEEEEEEEEEEEEEEEEEE

        # Aggiungo il nodo nella lista dei nodi filtrati (per la ricerca degli hub) se questo è presente nella lista candidates (costruita a monte valutando indregree+outdegree > 2)
        if best_next_node in _candidates and best_next_node not in _filtered_candidates:
            _filtered_candidates.append(best_next_node)

        # Se la copertura è troppo bassa, interrompiamo
        if best_cov < min_coverage:
            break

        # Aggiungo al contig l'ultimo nucleotide di best_next_node
        # Esempio: current_node="ATG", best_next_node="TGC" => aggiungiamo "C"
        contig_ext.append(best_next_node[-1])

        # Aggiorno current_node
        current_node = best_next_node

    return "".join(contig_ext), _filtered_candidates


###############################################################################
# 3. Estensione a sinistra (extend_left)
###############################################################################
def extend_left(node_id, dict_kmer_count, _candidates, _filtered_candidates, min_coverage=1):
    """
    Estende il contig a sinistra partendo da un (k-1)-mer (node_id),
    scegliendo il k-mer con coverage più alto tra quelli che finiscono con node_id.

    Parametri:
    - node_id: str, il (k-1)-mer di partenza
    - dict_kmer_count: dict { kmer (str) : coverage (int) }
    - k: lunghezza del k-mer
    - min_coverage: soglia minima di coverage

    Restituisce:
    - Una stringa (al contrario!) dei nucleotidi aggiunti a sinistra,
      che poi andrà ribaltata prima di unirla al contig finale.
    """
    contig_ext = []
    current_node = node_id

    while True:
        predecessors = []

        # Scansiona tutti i k-mer per trovare quelli che finiscono con current_node
        for kmer, coverage in dict_kmer_count.items():
            # Esempio: if current_node="ATG", cerchiamo kmer come "CATG", "GATG", ...
            if kmer.endswith(current_node):
                # best_prev_node = kmer[:-1] (il (k-1)-mer precedente)
                best_prev_node = kmer[:-1]
                predecessors.append((best_prev_node, coverage))

        if not predecessors:
            break

        # Ordina i predecessori per coverage in modo decrescente
        predecessors.sort(key=lambda x: x[1], reverse=True)

        best_prev_node, best_cov = predecessors[0]
        del dict_kmer_count[best_prev_node+current_node[-1]] #FATTO ALLA FINE E DA VERIFICAREEEEEEEEEEEEEEEEEEE

        # Aggiungo il nodo nella lista dei nodi filtrati (per la ricerca degli hub) se questo è presente nella lista candidates (costruita a monte valutando indregree+outdegree > 2)
        if best_prev_node in _candidates and best_prev_node not in _filtered_candidates:
            _filtered_candidates.append(best_prev_node)

        if best_cov < min_coverage:
            break

        # Aggiungo il primo nucleotide di best_prev_node in coda a contig_ext
        # Esempio: current_node="TGC", best_prev_node="ATG" => aggiungiamo "A"
        contig_ext.append(best_prev_node[0])

        # Aggiorno current_node
        current_node = best_prev_node

    # L'espansione left l'abbiamo costruita "al contrario",
    # quindi invertiamo la lista contig_ext
    contig_ext.reverse()
    return "".join(contig_ext), _filtered_candidates


###############################################################################
# 4. Costruzione del contig "locale" da un nodo scelto
###############################################################################
def build_local_contig(_dict_kmer_count, _candidates, _filtered_candidates,  _min_coverage=1):
    """
    1) Sceglie un nodo 'branching' a caso,
    2) Estende a sinistra e a destra,
    3) Restituisce l'intero contig (stringa nucleotidica).

    Parametri:
    - dict_kmer_count: { kmer : coverage }
    - nodes: { node_id : Node(...) } con .indegree e .outdegree
    - k: lunghezza dei k-mer
    - min_coverage: soglia minima per considerare un ramo "credibile"

    Restituisce:
    - contig (str) risultante
    """
    start_node = pick_random_branching_node(_candidates)
    if not start_node:
        # Nessun nodo con (in+out) > 2 => potremmo restituire stringa vuota
        return ""

    # Estensione a sinistra
    contig_left, filt_candidates_left = extend_left(start_node, _dict_kmer_count, _candidates, _filtered_candidates, min_coverage=_min_coverage)

    # Estensione a destra
    contig_right, filtered_candidates_tot = extend_right(start_node, _dict_kmer_count, _candidates, filt_candidates_left, min_coverage=_min_coverage)

    # Il contig è: contig_left + start_node + contig_right
    #  - contig_left (nucleotidi aggiunti a sinistra)
    #  - start_node: è il (k-1)-mer di partenza
    #  - contig_right: nucleotidi aggiunti a destra
    contig = contig_left + start_node + contig_right

    return contig, filtered_candidates_tot

if __name__ == '__main__':
    start = time.time()
    print(f"Lavoreremo con {len(dict_kmer_count)} kmer")
    distribuzione_kmer(dict_kmer_count)
    nodes = get_nodes_parallel(dict_kmer_count)
    print(f"Lavoreremo con {len(nodes)} nodi")
    candidates = count_hub(nodes)
    print(f"Gli hub nel grafo saranno {len(candidates)}")
    filtered_candidates = []
    counter = 1
    while len(filtered_candidates) != len(candidates):
        contig, list_new_candidates = build_local_contig(dict_kmer_count, candidates, filtered_candidates, _min_coverage=7)
        filtered_candidates += list_new_candidates
        print(f"Il {counter}° contig generato misura {len(contig)} basi")
        counter += 1
    end = time.time()
    print(f"Il processo ha richiesto SOLO {(end - start)/60:.2f} minuti")