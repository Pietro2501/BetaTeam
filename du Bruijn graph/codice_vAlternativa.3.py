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


def dict_filter(_dict_fastq):
    list_keys = list(_dict_fastq.keys())
    for key in list_keys:
        expect_err = exp_err(_dict_fastq[key]["qual"])
        if expect_err < 2:
            del _dict_fastq[key]
    return _dict_fastq


file = t.extract_info("PhoeVulATCC8482_R1.fq.gz")
dict_fastq = parse_fastq(file)

dict_filtered = dict_filter(dict_fastq)
seq_filtered = [dict_filtered[key]["seq"] for key in dict_filtered.keys()]


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

def construct_graph(reads, k):
    """Costruisco il grafo di De Bruijn da set di letture di lunghezza k, contando le ripetizioni dei nodi """
    edges = []
    nodes = set()

    for read in reads:
        for i in range(len(read) - k + 1):
            edges.append(Edge((read[i:i + k - 1], read[i + 1:i + k])))
            n1 = Node(read[i:i + k - 1])
            n2 = Node(read[i + 1:i + k])

            for node in nodes:
                if node == n1:
                    node.incrementa_contatore()
                    break
            else:
                nodes.add(n1)

            for node in nodes:
                if node == n2:
                    node.incrementa_contatore()
                    break
            else:
                nodes.add(n2)

    return edges, nodes


def output_contigs(g):
    """ Applica il percorso Euleriano per ricostruire la sequenza originaria """
    edges, nodes = g
    start = None
    for node in nodes:
        outcoming_edges = [edge for edge in edges if
                           edge.label.startswith(node.label)]  # trovo i possibili archi per partono dal nodo corrente
        incoming_edges = [edge for edge in edges if
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
                for node in nodes:
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
        for node in nodes:
            if node.label == destination_label:
                current = node
                break

        if current is None:
            raise ValueError(f"Errore: il nodo successivo con label {destination_label} non è stato trovato nel grafo.")

    return contig

len(dict_fastq)
len(seq_filtered)
g = construct_graph(seq_filtered, 25)
contig = output_contigs(g)
print(len(contig))