from BioSequence import tools_Bio as t

class Node:
    """ Classe associata al nodo, utile a definire in seguito il grafo di de Brujin"""

    def __init__(self, lab):
        self.label = lab
        self.contatore = 0

    def __hash__(self):
        return hash(self.label)

    def __eq__(self, other):
        return self.label == other.label

    def incrementa_contatore(self):
        self.contatore += 1

class Edge:
    """ Classe associata all'edge, utile a definire in seguito il grafo di de Brujin"""

    def __init__(self, km1mer_tuple):
        self.label = km1mer_tuple[0] + km1mer_tuple[1][-1:]

    def __eq__(self, other):
        return self.label == other.label

def parse_fastq(file_obj, offset=33) -> dict:
    """

    :param file_obj:
    :param offset:
    :return:
    """
    with open(file_obj, "r") as file_obj:
        def check_offset(offset: int) -> bool:
            res = False
            if type(offset) == int and offset in [33, 64]:
                res = True
            return res

        def check_coherence(acc: str, plus: str) -> bool:
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
    if is_ascii:
        ee = sum(map(lambda q: 10 ** (q / -10), map(lambda x: ord(x) - 33, qual)))
    else:
        ee = sum(map(lambda q: 10 ** (q / -10), qual))
    return ee

def hard_trimming(dict_fastq: dict, treshold=20):
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

def dict_filter(_dict_fastq, k):
    list_keys = list(_dict_fastq.keys())
    for key in list_keys:
        expect_err = exp_err(_dict_fastq[key]["qual"])
        if expect_err > k:
            del _dict_fastq[key]
    return _dict_fastq


file_f= t.extract_info(r"PhoeVulATCC8482_R1.fq.gz")
dict_fastq_f= parse_fastq(file_f)
keys = list(dict_fastq_f.keys())[:10]
# for acc in keys:
#     print(dict_fastq_f[acc]["seq"], len (dict_fastq_f[acc]["seq"]))
#     print(dict_fastq_f[acc]["qual"])
dict_filtered_f = dinamic_trimming(dict_fastq_f)
dict_double_filtered_f = dict_filter(dict_filtered_f, 3)
# for acc in keys:
#     print(dict_filtered_f[acc]["seq"], len (dict_filtered_f[acc]["seq"]))
#     print(dict_filtered_f[acc]["qual"])
seq_filtered_f = [dict_double_filtered_f[key]["seq"] for key in dict_double_filtered_f.keys() if len(dict_double_filtered_f[key]["seq"]) > 100]


file_r= t.extract_info(r"PhoeVulATCC8482_R2.fq.gz")
dict_fastq_r= parse_fastq(file_r)
keys_r = list(dict_fastq_r.keys())[:10]
# for acc in keys_r:
#     print(dict_fastq_r[acc]["seq"], len (dict_fastq_r[acc]["seq"]))
#     print(dict_fastq_r[acc]["qual"])
dict_filtered_r = dinamic_trimming(dict_fastq_r)
dict_double_filtered_r = dict_filter(dict_filtered_r, 4)
# for acc in keys_r:
#     print(dict_filtered_r[acc]["seq"], len (dict_filtered_r[acc]["seq"]))
#     print(dict_filtered_r[acc]["qual"])
seq_filtered_r = [dict_double_filtered_r[key]["seq"] for key in dict_double_filtered_r.keys() if len(dict_double_filtered_r[key]["seq"]) > 100]
print(len(seq_filtered_r), len(seq_filtered_f))

mucchio_selvaggio = seq_filtered_f + seq_filtered_r



def rev_comp(seqlist):
    table= str.maketrans('ACTG','TGAC')
    total_list = [seq.translate(table)[::-1] for seq in seqlist]
    return total_list

mucchio_selvaggio_comp_rev = rev_comp(mucchio_selvaggio)
delirio_totale = mucchio_selvaggio + mucchio_selvaggio_comp_rev

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
    """Costruisco il grafo di de Bruijn ottimizzando la ricerca con un dizionario"""
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
    """ Applica il percorso Euleriano per ricostruire la sequenza originaria """
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
        outcoming_edges = [edge for edge in edges.values() if
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

        del edges[next_edge.label]  # rimuovi l'arco dalla lista degli archi
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

# print(len(seq_filtered_f), len(seq_filtered_r), len(mucchio_selvaggio)) # mucchio_selvaggio = seq_filtered_f + seq_filtered_r
# print(len(rev_comp(mucchio_selvaggio)))
# prova = construct_graph(rev_comp(mucchio_selvaggio), 23)
# contig = output_contigs(prova)
# print(contig)

print(len(seq_filtered_f), len(seq_filtered_r), len(mucchio_selvaggio), len(mucchio_selvaggio_comp_rev), len(delirio_totale)) # mucchio_selvaggio = seq_filtered_f + seq_filtered_r
de_brujin_graph = construct_graph(delirio_totale, 29)
contig = output_contigs(de_brujin_graph)
print(len(contig))