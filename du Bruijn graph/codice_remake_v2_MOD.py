import pickle
import gzip
import random
import plotly.express as px
from multiprocessing import Pool
from collections import Counter, defaultdict
import time


class Node:
    """
    Represents a node in a de Brujin graph.

    Each node has a label to identify it, along with indegree and outdegree
    attributes to track incoming and outgoing edges respectively.

    :param label: The label identifying the node.
    :type label: Any
    :param indegree: The count of incoming edges to the node.
    :type indegree: int
    :param outdegree: The count of outgoing edges from the node.
    :type outdegree: int
    """

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


    # CREAZIONE DEI NODI E PROCESSAMENTO PARALLELIZZATO #


def partial_nodes(_chunk):
    """
    Generate a dictionary of local nodes with their in-degrees and out-degrees
    based on the given chunk of keys.

    The function takes a sublist of _dict_kmer_count keys (_chunk).
    It processes each key to form 'prefix' and 'suffix' substrings by slicing the
    first and last characters, respectively. For each unique prefix and suffix,
    it creates a new node if it does not already exist in the local_nodes dictionary.
    The prefix increments its out-degree and the suffix increments its in-degree.

    :param _chunk: A list containing string keys to process.
    :type _chunk: list

    :return: A dictionary of processed nodes, where each key is a string
        and the value is a Node object with updated in-degrees and
        out-degrees.
    :rtype: dict[str, Node]
    """
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
    """
    Merges multiple partial node dictionaries into a single dictionary.
    This is used to combine results from different computations into a single node
    dictionary with updated degrees.

    :param results: A list of dictionaries, where every dictionary represents partial nodes.
        Each dictionary maps a key to a node object
    :type results: list[dict]
    :return: A single dictionary of nodes, where the keys are unique and the nodes have
        aggregated indegree and outdegree values.
    :rtype: dict
    """
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
    """
    Divides a dictionary of k-mers into chunks, processes them in parallel

    The function employs multiprocessing to speed up the computation by distributing
    chunks of keys from the dictionary across multiple worker processes.
    The final results are merged into a unified dataset.

    :param _dict_kmer_count: A dictionary where keys are k-mer strings and values
        represent associated counts or related data.
    :param num_processes: An optional integer specifying the number of processes
        to use for parallel processing. Defaults to 6.

    :return: A unified dataset of nodes derived from processing the dictionary of k-mers
        in parallel.
    :rtype: dict
    """
    # divido il dizionario in chunks
    keys = list(_dict_kmer_count.keys())
    chunk_size = len(keys) // num_processes
    chunks = [keys[i:i + chunk_size] for i in range(0, len(keys), chunk_size)]

    # creo un pool di processi
    with Pool(processes=num_processes) as pool:
        # applica la funzione su ciascun chunk
        results = pool.map(partial_nodes, chunks)

    # unisco i risultati
    _nodes = merge_nodes(results)
    return _nodes


    # GRAFICO DI DISTRIBUZIONE K-MER  #


def distribuzione_kmer(_dict_kmer_count):
    """
    Generates a histogram displaying the distribution of k-mer multiplicities from the given
    k-mer counts dictionary. The x-axis represents k-mer multiplicities, and the y-axis
    represents k-mer respective frequencies.

    :param _dict_kmer_count: A dictionary where keys are k-mers and values are their corresponding
        multiplicities (integer counts).
    :type _dict_kmer_count: dict
    :return: None. The function does not return a value; it displays frequency distribution.
    """
    valori_numerici = list(_dict_kmer_count.values())
    contatore = Counter(valori_numerici) # conto le frequenze dei kmer con stessa molteplicità (contatore = dict)
    max_numerosita = max(valori_numerici)
    kmer_val = list(range(0, max_numerosita + 1))  # da 0 al massimo valore di numerosità
    frequenze = [contatore.get(n, 0) for n in kmer_val] # lista delle frequenze per ogni valore di numerosità (se una certa numerosità non è presente, di dafult:0)
    fig = px.bar(x=kmer_val,y=frequenze,title='Distribuzione delle Numerosità dei k-mer',labels={'x': 'Numerosità', 'y': 'Frequenza'},text=frequenze)
    # aggiorno layout ed estetica
    fig.update_layout(xaxis=dict(tickmode='linear'),yaxis=dict(title='Frequenza'),title=dict(x=0.5, xanchor='center'),uniformtext_minsize=8,uniformtext_mode='hide')
    fig.update_traces(texttemplate='%{text}', textposition='outside') # serve per stampare il valore già associato tramite "text=frequenze" nel metodo "px.bar"
    fig.show()


    # SELEZIONE NODI HUB #


def hub_finding(_nodes):
    """
    Identifies and returns a list of hub nodes from the given node dictionary.

    :param _nodes: A dictionary where keys represent node identifiers and values 
        are objects that must have `indegree` and `outdegree` attributes.
    :return: A list containing node identifiers of all nodes whose total degree 
        (indegree + outdegree) exceeds 2.
    :rtype: 
        list
    """
    _candidates = [n_id for n_id, obj in _nodes.items() if (obj.indegree + obj.outdegree) > 2]
    return _candidates


    # SCELTA DI NODO HUB RANDOM #


def pick_random_branching_node(_candidates, _filtered_candidates=None):
    """
    Selects a random branching node from the provided candidates, excluding
    those in the filtered candidates list if it is provided.

    :param _candidates: A list of all candidate nodes from which a branching
        node will be randomly selected.
    :type _candidates: list
    :param _filtered_candidates: A list of candidates that should be excluded
        from selection. Default is None.
    :type _filtered_candidates: list or None
    :return: A randomly selected branching node from the remaining candidates
        or None if no eligible candidates exist after filtering.
    :rtype: object or None
    """
    if not _candidates:
        return None
    if _filtered_candidates is None:
        _filtered_candidates = []
    current_candidates = [item for item in _candidates if item not in _filtered_candidates]
    if not current_candidates:
        return None
    return random.choice(current_candidates)


    # CREAZIONE DI STRUTTURE DI PROSSIMITA' PER EVITARE SCANSIONI TOTALI #


def get_proximity_dict(dict_kmer_count):
    """
    Computes proximity dictionaries for k-mers based on their prefixes and suffixes.

    This function processes a dictionary of k-mers and their coverages
    to construct two proximity dictionaries based on k-mers prefix or suffix
    and their associated coverages.

    :param dict_kmer_count: A dictionary where the key is a k-mer (string) and the
                            value is its corresponding coverage (int).
    :return: A tuple consisting of two dictionaries:
             - The first dictionary maps a k-mer prefix to a list of tuples, where
               each tuple contains a suffix (string) and its coverage (int).
             - The second dictionary maps a k-mer suffix to a list of tuples, where
               each tuple contains a prefix (string) and its coverage (int).
    :rtype: Tuple[DefaultDict, DefaultDict]
    """
    dict_prox_right = defaultdict(list)
    dict_prox_left = defaultdict(list)

    for kmer, coverage in dict_kmer_count.items():
        prefix = kmer[:-1]
        suffix = kmer[1:]
        dict_prox_right[prefix].append((suffix, coverage))
        dict_prox_left[suffix].append((prefix, coverage))

    return dict_prox_right, dict_prox_left


    # FUNZIONE DI ESTENSIONE DEL CONTIG A DESTRA #


def extend_right(node_id, dict_prox_right, dict_prox_left,_candidates, _filtered_candidates, min_coverage):
    """
    Extends a contig starting from a given node by iteratively appending the most suitable
    successor node (based on coverage) to the sequence, until no suitable extension is
    possible or the coverage drops below a defined threshold. This function updates the
    proximity dictionaries and filters candidate nodes if certain conditions are met.

    :param node_id: The starting node for the extension.
    :type node_id: str
    :param dict_prox_right: A dictionary mapping a node to its list of right-hand-side
                            successors with their respective coverages.
    :type dict_prox_right: dict
    :param dict_prox_left: A dictionary mapping a node to its list of left-hand-side
                           predecessors with their respective coverages.
    :type dict_prox_left: dict
    :param _candidates: A list of candidate hub nodes.
    :type _candidates: list[str]
    :param _filtered_candidates: A list of nodes that have already been identified
                                  and marked as hubs during the extension process.
    :type _filtered_candidates: list[str]
    :param min_coverage: The minimum required coverage value for a successor node
                         to be considered for extension.
    :type min_coverage: int
    :return: A tuple containing the extended contig (as a string) and the updated
             list of filtered candidates.
    :rtype: tuple[str, list[str]]
    """

    contig_ext = []
    current_node = node_id

    while True:
        successors = dict_prox_right[current_node]
        if not successors:
            break

        # ordino in modo decrescente i successors sulla base della coverage
        successors.sort(key=lambda x: x[1], reverse=True)
        best_next_node, best_cov = successors[0]

        if best_cov < min_coverage:
            break

        # rimuovo dalla lista dict_prox_right[current_node] il suffisso associato al kmer percorso
        dict_prox_right[current_node] = [(suffix, coverage) for (suffix, coverage) in successors if suffix != best_next_node]
        # rimuovo da dict_prox_left[best_next_node] il prefisso associato al kmer percorso
        dict_prox_left[best_next_node] = [(prefix, coverage) for (prefix, coverage) in dict_prox_left[best_next_node] if prefix != current_node]

        # se best_next_node è un hub e non è ancora filtrato, lo aggiungo alla lista _filtered_candidates
        if best_next_node in _candidates and best_next_node not in _filtered_candidates:
            _filtered_candidates.append(best_next_node)

        # aggiungo al contig l'ultimo nucleotide di best_next_node
        contig_ext.append(best_next_node[-1])
        current_node = best_next_node

    return "".join(contig_ext), _filtered_candidates


    # FUNZIONE DI ESTENSIONE DEL CONTIG A SINISTRA #

def extend_left(node_id, dict_prox_right, dict_prox_left,_candidates, _filtered_candidates, min_coverage):
    """
    Extends a contig sequence to the left based on a directed graph structure. The function
    traverses the graph by selecting the best predecessor node iteratively, based on coverage,
    and stops when a suitable node is no longer found. It also manages candidate nodes and
    tracks nodes that are filtered during the process.

    :param node_id: The starting node ID from which the extension begins.
    :param dict_prox_right: A dictionary containing the mapping of nodes to their right
        neighboring nodes along with associated coverage.
    :param dict_prox_left: A dictionary containing the mapping of nodes to their left
        neighboring nodes along with associated coverage.
    :param _candidates: A list of candidate nodes that may need further processing.
    :param _filtered_candidates: A list of nodes that have been filtered out during
        the traversal.
    :param min_coverage: The minimum coverage required for a predecessor node to be considered
        as the next node in the extension.
    :return: A tuple containing the extended sequence as a string and the updated filtered
        candidates list.
    """
    contig_ext = []
    current_node = node_id

    while True:
        predecessors = dict_prox_left[current_node]
        if not predecessors:
            break

        # ordino in modo decrescente i successors sulla base della coverage
        predecessors.sort(key=lambda x: x[1], reverse=True)
        best_prev_node, best_cov = predecessors[0]

        if best_cov < min_coverage:
            break

        dict_prox_left[current_node] = [(prefix, coverage) for (prefix, coverage) in predecessors if prefix != best_prev_node]
        dict_prox_right[best_prev_node] = [(suffix, coverage) for (suffix, coverage) in dict_prox_right[best_prev_node] if suffix != current_node]

        # se best_next_node è un hub e non è ancora filtrato, lo aggiungo alla lista _filtered_candidates
        if best_prev_node in _candidates and best_prev_node not in _filtered_candidates:
            _filtered_candidates.append(best_prev_node)

        # aggiungo il primo nucleotide di best_prev_node a contig_ext
        contig_ext.append(best_prev_node[0])

        current_node = best_prev_node

    # inverto la lista contig_ext in quanto costruita al contrario
    contig_ext.reverse()
    return "".join(contig_ext), _filtered_candidates


    # COSTRUZIONE DI UN CONTIG A PARTIRE DA UN NODO HUB #

def build_local_contig(dict_prox_right,dict_prox_left,_candidates,_filtered_candidates,min_coverage):
    """
    Constructs a local contig by selecting a starting branching node and extending it in both left
    and right directions. A contig represents a sequence of nodes, which is constructed from
    proximal node connections.

    :param dict_prox_right: Right proximal nodes dictionary for connections.
    :type dict_prox_right: dict
    :param dict_prox_left: Left proximal nodes dictionary for connections.
    :type dict_prox_left: dict
    :param _candidates: list of candidate nodes available for extension.
    :type _candidates: list
    :param _filtered_candidates: list of already filtered candidates that will be further pruned.
    :type _filtered_candidates: list
    :param min_coverage: Minimum coverage threshold for extending the contig.
    :type min_coverage: int
    :return: A tuple containing the constructed contig string and the updated filtered
             candidates list after contig extension.
    :rtype: tuple
    """
    start_node = pick_random_branching_node(_candidates, _filtered_candidates)
    if not start_node:
        print("!!! Non ho trovato nessun nodo hub per questo contig!!!")
        return "", _filtered_candidates

    contig_left, filt_candidates_left = extend_left(start_node,dict_prox_right,dict_prox_left,_candidates,_filtered_candidates,min_coverage=min_coverage)
    contig_right, new_filtered_candidates = extend_right(start_node,dict_prox_right,dict_prox_left,_candidates,filt_candidates_left,min_coverage=min_coverage)

    # genero il contig completo
    contig = contig_left + start_node + contig_right
    return contig, new_filtered_candidates


    # CALCOLO DELL'N50 SULLA BASE DELLE LUNGHEZZE DEI CONTIG PRODOTTI #

def calcola_n50(lunghezze_contig:list):
    """
    Calculates the N50 value of a list of contig lengths. The N50 is a statistical
    measure of sequence assembly quality in genomics. It is defined as the smallest
    contig size at which the sum of the lengths of the contigs contains at least
    50% of the total assembly length.

    :param lunghezze_contig: A list of integers representing the lengths of the contigs.
    :return: An integer representing the N50 value, which is the contig length
        that meets or exceeds 50% of the total length of all contigs.
    """
    try:
        # ordino le lunghezze dei contig in ordine decrescente
        lunghezze_contig.sort(reverse=True)
        lunghezza_totale = sum(lunghezze_contig)

        # sommo le lunghezze dei contig finché non si raggiunge almeno il 50% della lunghezza totale
        somma_parziale = 0
        for contig in lunghezze_contig:
            somma_parziale += contig
            if somma_parziale >= lunghezza_totale / 2:
                return contig  # restituisco la lunghezza del contig che raggiunge o supera il 50%
    except IndexError:
        print("Lista di contig vuota o dati non validi...")


    # COSTRUZIONE DEI CONTIG SULLA BASE DELLE NECESSITA' DELL'UTENTE #


def iterative_contig_generation(dict_kmer_count, _candidates, _start, _min_coverage=7):
    """
    This function performs an iterative generation of contigs from a dictionary of k-mer counts. The user has 2 options:
    covering all candidate hubs in the graph or specifying a desired number of contigs to generate.
    The process involves proximity dictionary construction, filtered candidates tracking, and N50 calculation
    for the generated contigs, finally reporting time elapsed for the entire operation.

    :param dict_kmer_count: A dictionary containing k-mer counts used to generate proximity dictionaries and perform
                             local contig construction.
    :type dict_kmer_count: dict
    :param _candidates: A list of candidate nodes (hubs) in the graph representing potential starting points for contig
                        generation.
    :type _candidates: list
    :param _start: Start time stamp obtained before initiating the iterative generation process, used for calculating
                   time elapsed.
    :type _start: float
    :param _min_coverage: Minimum coverage required for a k-mer to be included in contig generation. Defaults to 7.
    :type _min_coverage: int
    :return: None. The results and metrics are printed directly to the console.
    :rtype: None
    """
    filtered_candidates = []
    counter = 1
    lenghts_contig = []

    while True:
        try:
            iter_choice = input("\nQuanto profonda vuoi sia la lettura del grafo?\n\t- '1' per ricoprire tutti i nodi hub del grafo"
                                " (operazione computazionalmente onerosa)\n\t- '2' per un numero di contig da te desiderato\n")
            if iter_choice not in ['1', '2']:
                raise ValueError
            break
        except ValueError:
            print("Attenzione a cosa digiti... Inserisci '1' o '2' per scegliere il tipo di iterazione.")

    # se l'utente sceglie '2', chiedo il numero massimo di iterazioni
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

    # eseguo il loop in base alla scelta dell'utente
    if iter_choice == '1':
        while len(filtered_candidates) < len(_candidates):
            _dict_prox_right, _dict_prox_left = get_proximity_dict(dict_kmer_count)
            contig, filtered_candidates = build_local_contig(_dict_prox_right,_dict_prox_left,_candidates,filtered_candidates,_min_coverage)
            print(f"Il {counter}° contig generato misura {len(contig)} basi")
            print(f"Hub filtrati finora: {len(filtered_candidates)} / {len(_candidates)}")
            lenghts_contig.append(len(contig))
            counter += 1
    else:
        while counter <= max_iterations:
            _dict_prox_right, _dict_prox_left = get_proximity_dict(dict_kmer_count)
            contig, filtered_candidates = build_local_contig(_dict_prox_right,_dict_prox_left,_candidates,filtered_candidates,_min_coverage)
            print(f"Il {counter}° contig generato misura {len(contig)} basi")
            print(f"Hub filtrati finora: {len(filtered_candidates)} / {len(_candidates)}")
            lenghts_contig.append(len(contig))
            counter += 1

    print(f"\nL'N50 ottenuto dai contig generati è pari a {calcola_n50(lenghts_contig)}")
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

def kmer_filter_coverage():
    """
    Prompts the user to provide the minimum coverage value for k-mers to be used
    in graph construction for contigs.

    :return: The minimum coverage value specified by the user.
    :rtype: int
    """
    while True:
        try:
            min_coverage_choice = int(input("\nInserisci il valore minimo di coverage che vorresti abbiano i kmer del grafo per la costruzione dei contig:\n"))
            if min_coverage_choice <= 0:
                raise ValueError("Inserisci un valore intero positivo")
            break
        except ValueError:
            print("Attenzione a cosa digiti... Inserisci un valore intero")
    return min_coverage_choice



    #MAIN#


def main():
    start = time.time()
    intestazione_progetto()

    with gzip.open('kmer_17_NEW_PhoeVul_num1filtered.pkl.gz', 'rb') as f:
        dict_kmer_count = pickle.load(f)
    print(f"Lavoreremo con {len(dict_kmer_count)} k-mer")

    distribuzione_kmer(dict_kmer_count)

    nodes = get_nodes_parallel(dict_kmer_count)
    print(f"Lavoreremo con {len(nodes)} nodi")

    candidates = hub_finding(nodes)
    print(f"Gli hub nel grafo saranno {len(candidates)}")

    iterative_contig_generation(dict_kmer_count, candidates, start, kmer_filter_coverage())

if __name__ == '__main__':
    main()

