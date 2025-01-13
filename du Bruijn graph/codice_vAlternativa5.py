from BioSequence import tools_Bio as t
import time

from trimming import dinamic_trimming


class Node:
    """
    Represents a node with a label and an internal counter.

    The Node class provides the ability to uniquely identify nodes
    through a label and includes an internal counter that can be
    incremented. Nodes are hashable and support equality comparison
    based on their label.

    :ivar label: The label identifying the node.
    :type label: Any
    :ivar contatore: An internal counter associated with the node.
    :type contatore: int
    """

    def __init__(self, lab):
        self.label = lab
        self.contatore = 0

    def __eq__(self, other):
        return self.label == other.label

    def incrementa_contatore(self):
        self.contatore += 1

class Edge:
    """
    Represents an edge in the graph structure using a k-1 mer tuple.

    This class is designed to model edges where the labels are constructed
    from a pair of k-1 mers and extend their sequence by overlapping the final
    character of the second k-1 mer. The equality of edges is determined by
    comparing their labels.

    :ivar label: Label of the edge formed by merging two k-1 mers.
    :type label: str
    """

    def __init__(self, km1mer_tuple):
        self.label = km1mer_tuple[0] + km1mer_tuple[1][-1:]

    def __eq__(self, other):
        return self.label == other.label

def parse_fastq(file_obj, offset=33) -> dict:
    """
    Parses a FASTQ file and retrieves its content in the form of a dictionary. Each record in the
    FASTQ file includes the accession line, sequence, and quality values. The function also
    calculates the decoded quality scores from the given ASCII quality values based on the
    provided offset. The function ensures coherence in the format of the accession and plus lines.

    :param file_obj: The file path to the FASTQ file to be parsed.
    :type file_obj: str
    :param offset: Quality score offset, which is expected to be either 33 or 64. Used to decode
        ASCII quality scores into numeric values. Default is 33.
    :type offset: int, optional
    :return: A dictionary representing the parsed FASTQ file. The keys are the accession lines, and
        the values contain a dictionary of sequence, ASCII quality values, and decoded quality scores.
        If the file is corrupted or the offset is invalid, returns None.
    :rtype: dict or None
    """
    with open(file_obj, "r") as file_obj:
        def check_offset(offset: int) -> bool:
            res = False
            if type(offset) == int and offset in [33, 64]:
                res = True
            return res

        def check_coherence(acc: str, plus: str) -> bool:
            """
            Checks the coherence of a FASTQ record based on the given accession line and
            plus line. A FASTQ record is considered coherent if the accession line starts
            with "@" and the plus line starts with "+".

            :param acc: The accession line from the FASTQ record.
            :type acc: str
            :param plus: The plus line from the FASTQ record.
            :type plus: str
            :return: A boolean indicating whether the FASTQ record is coherent.
            :rtype: bool
            """
            res = False
            if acc.startswith("@") and plus.startswith("+"):
                res = True
            return res

        if check_offset(offset):
            fastq_dict = {}
            line = file_obj.readline()
            while line:
                acc = line.strip()
                fastq_dict.setdefault(acc, {})
                # qui andrebbe inserito un try/expcet per verificare che effettivamente tutte le linee siano lette
                seq, plus, qual = file_obj.readline().strip(), file_obj.readline().strip(), file_obj.readline().strip()
                if not check_coherence(acc, plus):
                    print("Il file è corrotto!!!")
                    fastq_dict = None
                    break
                if seq and qual:
                    fastq_dict[acc]["seq"] = seq
                    fastq_dict[acc]["ASCII_qual"] = qual
                    fastq_dict[acc]["qual"] = [ord(q) - offset for q in fastq_dict[acc]["ASCII_qual"]]
                    # fastq_dict[acc]["Pe"] = [10 ** (q / -10) for q in fastq_dict[acc]["qual"]]
                else:
                    print("Il file è corrotto!!!")
                    fastq_dict = None
                    break
                line = file_obj.readline()
            return fastq_dict
        elif not check_offset(offset):
            print(f"l'offset indicato non è numerico o è un valore differente da 33 o 64!!!")


def exp_err(qual, is_ascii=False):
    """Calculate the total expected error of a sequence's quality scores.

    Args:
        qual (str or list): Quality scores, either as a string of ASCII-encoded characters
                            or a list of numeric values.
        is_ascii (bool): A flag indicating if the input is ASCII-encoded (True) or numeric (False).

    Returns:
        float: The total expected error, calculated as a sum of probabilities.

    Raises:
        ValueError: If input is not a string or a list of numeric values.
    """
    if is_ascii:
        try:
            ee = sum(map(lambda q: 10 ** (q / -10), map(lambda x: ord(x) - 33, qual)))
        except TypeError:
            raise ValueError("Input quality scores must be a string for ASCII-encoded values.")
    else:
        try:
            ee = sum(map(lambda q: 10 ** (q / -10), qual))
        except TypeError:
            raise ValueError("Input quality scores must be a list of numeric values.")
    return ee

def hard_trimming(dict_fastq: dict, treshold=20):
    """
    Perform hard trimming of FASTQ sequence data by removing low-quality bases from the
    beginning of sequences. This method updates the sequence (`seq`), quality scores
    (`qual`), and ASCII-encoded quality scores (`ASCII_qual`) for each key in the
    provided dictionary, starting from the first base that meets or exceeds the
    specified quality threshold.

    :param dict_fastq: Dictionary containing FASTQ data for processing. Each key in
        the dictionary corresponds to a record, and its value is a dictionary with
        the keys: `seq` (sequence string), `qual` (quality score list),
        and `ASCII_qual` (ASCII-encoded quality scores).
    :type dict_fastq: dict
    :param treshold: Quality score threshold to determine the starting position for
        trimming. Default is 20.
    :type treshold: int
    :return: Updated `dict_fastq` dictionary after hard trimming has been applied
        to all entries.
    :rtype: dict
    """

    for key in dict_fastq.keys():
        for value in range(len(dict_fastq[key]["qual"])):
            if dict_fastq[key]["qual"][value] > treshold:
                dict_fastq[key].update({"seq": dict_fastq[key]["seq"][value:]})
                dict_fastq[key].update({"qual": dict_fastq[key]["qual"][value:]})
                dict_fastq[key].update({"ASCII_qual": dict_fastq[key]["ASCII_qual"][value:]})
                #print(f"ho aggiornato i valori associati alla chiave {key}")
                break
    return dict_fastq

def dinamic_trimming(dict_fastq: dict, treshold=25, window=15):
    """
    Perform dynamic trimming of FASTQ sequence data by removing low-quality bases from the
    beginning of sequences. This method updates the sequence (`seq`), quality scores
    (`qual`), and ASCII-encoded quality scores (`ASCII_qual`) for each key in the
    provided dictionary, starting from the first base that meets or exceeds the
    specified quality threshold.

    :param dict_fastq: Dictionary containing FASTQ data for processing. Each key in
        the dictionary corresponds to a record, and its value is a dictionary with
        the keys: `seq` (sequence string), `qual` (quality score list),
        and `ASCII_qual` (ASCII-encoded quality scores).
    :type dict_fastq: dict
    :param treshold: Quality score threshold to determine the starting position for
        trimming. Default is 20.
    :type treshold: int
    :return: Updated `dict_fastq` dictionary after hard trimming has been applied
        to all entries.
    :rtype: dict

    """
    for key in dict_fastq.keys():
        quality = dict_fastq[key]["qual"]
        list_mean = quality[:window]
        if sum(list_mean)/len(list_mean) < treshold:
            for value in range(len(list_mean), len(quality)):
                if (sum(list_mean)+quality[value])/(len(list_mean)+1) < treshold:
                    list_mean.append(quality[value])
                else:
                    break
            dict_fastq[key].update({"seq": dict_fastq[key]["seq"][len(list_mean):]})
            dict_fastq[key].update({"qual": dict_fastq[key]["qual"][len(list_mean):]})
            dict_fastq[key].update({"ASCII_qual": dict_fastq[key]["ASCII_qual"][len(list_mean):]})
    return dict_fastq

def dict_filter(_dict_fastq):
    """
    Filters a dictionary of FASTQ entries by expected error threshold.

    This function takes a dictionary of FASTQ entries and evaluates the expected
    error based on the "qual" field within each entry. If the expected error
    exceeds a specified threshold (here, 3), the corresponding entry is removed
    from the dictionary.

    :param _dict_fastq: Dictionary containing FASTQ entries, where keys are
        identifiers and values are dictionaries with a "qual" field
        representing quality scores.
    :type _dict_fastq: dict

    :return: A filtered dictionary containing only the FASTQ entries
        where the expected error is 3 or less.
    :rtype: dict
    """
    list_keys = list(_dict_fastq.keys())
    for key in list_keys:
        expect_err = exp_err(_dict_fastq[key]["qual"])
        if expect_err > 3:
            del _dict_fastq[key]
    return _dict_fastq


file_f= t.extract_info(r"PhoeVulATCC8482_R1.fq.gz")
dict_fastq_f= parse_fastq(file_f)
keys = list(dict_fastq_f.keys())[:10]
# for acc in keys:
#     print(dict_fastq_f[acc]["seq"], len (dict_fastq_f[acc]["seq"]))
#     print(dict_fastq_f[acc]["qual"])
dict_filtered_f = dinamic_trimming(dict_fastq_f)
# for acc in keys:
#     print(dict_filtered_f[acc]["seq"], len (dict_filtered_f[acc]["seq"]))
#     print(dict_filtered_f[acc]["qual"])
seq_filtered_f = [dict_filtered_f[key]["seq"] for key in dict_filtered_f.keys() if len(dict_filtered_f[key]["seq"]) > 100]


file_r= t.extract_info(r"PhoeVulATCC8482_R2.fq.gz")
dict_fastq_r= parse_fastq(file_r)
keys_r = list(dict_fastq_r.keys())[:10]
# for acc in keys_r:
#     print(dict_fastq_r[acc]["seq"], len (dict_fastq_r[acc]["seq"]))
#     print(dict_fastq_r[acc]["qual"])
dict_filtered_r = dinamic_trimming(dict_fastq_r)
# for acc in keys_r:
#     print(dict_filtered_r[acc]["seq"], len (dict_filtered_r[acc]["seq"]))
#     print(dict_filtered_r[acc]["qual"])
seq_filtered_r = [dict_filtered_r[key]["seq"] for key in dict_filtered_r.keys() if len(dict_filtered_r[key]["seq"]) > 100]
print(len(seq_filtered_r), len(seq_filtered_f))

mucchio_selvaggio = seq_filtered_f + seq_filtered_r



def rev_comp(seqlist):
    """
    Generate the reverse complement for a list of DNA sequences and combine it with the
    original sequences.

    Given a list of DNA sequences, this function computes the reverse complement of
    each sequence by translating the nucleotides (A ↔ T, C ↔ G), and then reverses it.
    The resulting reverse complement sequences are combined with the original list
    of sequences and returned as a single list.

    :param seqlist: List of DNA sequences for which reverse complements are to
        be generated.
    :type seqlist: list[str]
    :return: Combined list containing original sequences and their respective
        reverse complements.
    :rtype: list[str]
    """
    table= str.maketrans('ACTG','TGAC')
    total_list = [seq.translate(table)[::-1] for seq in seqlist] + seqlist
    return total_list


# def read_reads(fname):
#     """ FUnzione di prova che genera una lista di sequenza e prende in input un file FASTA (per ora) """
#     f = open(fname, 'r')
#     lines = f.readlines()
#     f.close()
#     reads = []

#     for line in lines:
#         if line[0] != '>':
#             reads = reads + [line.rstrip()]

#     return reads

# def parser_fastq(filename):
#     from Bio import SeqIO
#     import gzip
#     with gzip.open(filename, "rt") as file:
#         listaseq= [str(seq_record.seq) for seq_record in SeqIO.parse(file,"fastq")]
#     return listaseq

# def construct_graph(reads, k):
#     """Costruisco il grafo di De Bruijn da set di letture di lunghezza k, contando le ripetizioni dei nodi """
#     edges = []
#     nodes = set()
#
#     for read in reads:
#         for i in range(len(read) - k + 1):
#             edges.append(Edge((read[i:i + k - 1], read[i + 1:i + k])))
#             n1 = Node(read[i:i + k - 1])
#             n2 = Node(read[i + 1:i + k])
#
#             for node in nodes:
#                 if node == n1:
#                     node.incrementa_contatore()
#                     break
#             else:
#                 nodes.add(n1)
#
#             for node in nodes:
#                 if node == n2:
#                     node.incrementa_contatore()
#                     break
#             else:
#                 nodes.add(n2)
#
#     return edges, nodes

def construct_graph(reads, k):
    """
    Constructs a graph from a list of reads using a k-mer approach. The function generates edges
    and nodes based on overlapping k-1 subsequences from each read. The edges and nodes are
    stored in dictionaries.

    :param reads: A list of strings representing DNA reads.
    :type reads: list[str]
    :param k: The size of the k-mers used for graph construction.
    :type k: int
    :return: A tuple containing a dictionary of edges and a dictionary of nodes. The edges
             dictionary maps edge labels to Edge objects, and the nodes dictionary maps node
             labels to Node objects.
    :rtype: tuple[dict[str, Edge], dict[str, Node]]
    """
    edges = {}
    nodes = {}

    for read in reads:
        for i in range(len(read) - k + 1):
            # Creazione e gestione degli archi
            edge_label = read[i:i + k - 1] + read[i + 1:i + k][-1]
            if edge_label not in edges:
                edges[edge_label] = Edge((read[i:i + k - 1], read[i + 1:i + k]))

            # Gestione dei nodi con dizionari
            n1_label = read[i:i + k - 1]
            n2_label = read[i + 1:i + k]

            if n1_label not in nodes:
                nodes[n1_label] = Node(n1_label)
            else:
                nodes[n1_label].incrementa_contatore()

            if n2_label not in nodes:
                nodes[n2_label] = Node(n2_label)
            else:
                nodes[n2_label].incrementa_contatore()

    return edges, nodes

def output_contigs(g):
    """
    Generate a contiguous sequence (contig) from a graph structure.

    The function attempts to generate a contig from a graph represented by
    nodes and edges. Starting from a node with outgoing edges but no
    incoming edges, it traverses the graph by following edges with criteria
    to select the next edge when multiple options are available.

    :param g:
        A tuple containing the graph data structure. The first element is a
        dictionary of edges, where each edge is represented by an object with
        a `label` attribute. The second element is a dictionary of nodes,
        each represented by an object with labels (`label`) and a `contatore`
        attribute that acts as a counter.
    :return:
        A string representing the generated contig by traversing the graph
        according to the described logic.
    :raises ValueError:
        If a starting node cannot be found, invalid edges exist, or a
        destination node is missing during traversal.
    """
    edges, nodes = g
    start = None
    for node in nodes.values():
        outcoming_edges = [edge for edge in edges.values() if
                           edge.label.startswith(node.label)]  # trovo i possibili archi per partono dal nodo corrente
        incoming_edges = [edge for edge in edges.values() if
                          edge.label.endswith(node.label)]  # trovo i possibili archi per arrivano nel nodo corrente

        # il nodo di partenza è quello che ha almeno un arco in uscita, ma nessun arco in ingresso
        if outcoming_edges and not incoming_edges:
            start = node
            break

    if start is None:
        raise ValueError("Nessun nodo di partenza trovato. Verifica la struttura del grafo.")

    contig = start.label
    current = start

    while True:
        outcoming_edges = [edge for edge in edges if
                           edge.label.startswith(current.label)]  # trovo gli archi che partono dal nodo corrente

        if len(outcoming_edges) == 0:
            break  # se non ci sono più archi esce dal ciclo
        elif len(outcoming_edges) == 1:
            next_edge = outcoming_edges[0]
        else:  # se ci sono più archi, scegli quello il cui nodo di arrivo abbia il contatore maggiore
            best_edge = None
            max_contatore = -1

            for edge in outcoming_edges:
                destination_label = edge.label[1:]
                # trova il nodo di destinazione per l'arco sulla base della label identificata
                destination_node = None
                for node in nodes.values():
                    if node.label == destination_label:
                        destination_node = node
                        break
                if destination_node.contatore > max_contatore:
                    best_edge = edge
                    max_contatore = destination_node.contatore

            if best_edge is None:
                raise ValueError("Nessun arco valido trovato con un nodo di destinazione valido.")

            next_edge = best_edge

        edges.remove(next_edge)  # rimuovi l'arco dalla lista degli archi
        contig += next_edge.label[-1]  # aggiungi l'ultimo carattere del k-mer al contig

        destination_label = next_edge.label[1:]  # definisco il nodo successivo per ricostruire il contig
        current = None
        for node in nodes.values():
            if node.label == destination_label:
                current = node
                break

        if current is None:
            raise ValueError(f"Errore: il nodo successivo con label {destination_label} non è stato trovato nel grafo.")

    return contig

#print(len(seq_filtered_f), len(seq_filtered_r), len(mucchio_selvaggio)) # mucchio_selvaggio = seq_filtered_f + seq_filtered_r
#print(len(rev_comp(mucchio_selvaggio)))
#prova = construct_graph(rev_comp(mucchio_selvaggio), 23)
#contig = output_contigs(prova)
#print(contig)