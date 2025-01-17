import pickle
import gzip

class Node:
    """ Classe associata al nodo, utile a definire in seguito il grafo di de Brujin"""

    def __init__(self, lab):
        self.label = lab
        self.indegree= 0
        self.outdegree = 0

    def __eq__(self, other):
        return self.label == other.label

class Edge:
    """ Classe associata all'edge, utile a definire in seguito il grafo di de Brujin"""

    def __init__(self, km1mer_tuple):
        self.label = km1mer_tuple[0] + km1mer_tuple[1][-1:]
        self.contatore = 0

    def __eq__(self, other):
        return self.label == other.label

    def incrementa_contatore(self):
        self.contatore += 1

with gzip.open('list_seq_PhoeVul_filtered.pkl.gz', 'rb') as f:
    reads = pickle.load(f)

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