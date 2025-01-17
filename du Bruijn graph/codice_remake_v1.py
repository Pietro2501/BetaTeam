import pickle
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

class Node:
    """ Classe associata al nodo, utile a definire in seguito il grafo di de Brujin"""

    def __init__(self, lab):
        self.label = lab
        self.indegree= 0
        self.outdegree = 0

    def __eq__(self, other):
        return self.label == other.label

    def incrementa_indegree(self):
        self.indegree += 1

    def incrementa_outdegree(self):
        self.outdegree += 1

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

with gzip.open('kmer_diz_PhoeVul.pkl.gz', 'rb') as f:
    dict_kmer_count = pickle.load(f)

# print(len(reads), len(dict_kmer_count))
# keys = list(dict_kmer_count.keys())[:20]
# for key in keys:
#     print(f" {key}: {dict_kmer_count[key]}")

def get_nodes(_dict_kmer_count):
    """Costruisco il grafo di de Bruijn ottimizzando la ricerca con un dizionario"""
    # edges = {}
    nodes = {}

    for key in _dict_kmer_count.keys():
        prefix = key[:-1]
        suffix = key[1:]

        if prefix not in nodes:
            nodes[prefix] = Node(prefix)
            nodes[prefix].incrementa_outdegree()
        else:
            nodes[prefix].incrementa_outdegree()

        if suffix not in nodes:
            nodes[suffix] = Node(suffix)
            nodes[prefix].incrementa_indegree()
        else:
            nodes[prefix].incrementa_indegree()
    return nodes

def distribuzione_kmer(_dict_kmer_count):
    valori_numerici = list(_dict_kmer_count.values())
    contatore = Counter(valori_numerici)
    kmer_val = list(range(0, 16))  # Da 0 a 15
    frequenze = [contatore.get(n, 0) for n in kmer_val]
    plt.figure(figsize=(10, 6))
    plt.bar(kmer_val, frequenze, color='skyblue', edgecolor='black')
    plt.title('Distribuzione delle Numerosità', fontsize=16)
    plt.xlabel('Numerosità', fontsize=14)
    plt.ylabel('Frequenza', fontsize=14)

    plt.xticks(kmer_val, fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

distribuzione_kmer(dict_kmer_count)


def output_contigs(g):
    """ Applica il percorso Euleriano per ricostruire la sequenza originaria """
    edges, nodes = g
    start = None
    for node in nodes.values():
        outcoming_edges = [edge for edge in edges.values() if edge.label.startswith(node.label)]  # trovo i possibili archi per partono dal nodo corrente
        incoming_edges = [edge for edge in edges.values() if edge.label.endswith(node.label)]  # trovo i possibili archi per arrivano nel nodo corrente
        # il nodo di partenza è quello che ha almeno un arco in uscita, ma nessun arco in ingresso
        if outcoming_edges and not incoming_edges:
            start = node
            break
    if start is None:
        raise ValueError("Nessun nodo di partenza trovato. Verifica la struttura del grafo.")

    contig = start.label
    current = start

    while True:
        outcoming_edges = [edge for edge in edges.values() if edge.label.startswith(current.label)]  # trovo gli archi che partono dal nodo corrente

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