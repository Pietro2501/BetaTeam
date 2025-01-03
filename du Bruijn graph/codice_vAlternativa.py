class Node:
    """Class Node to represent a vertex in the De Bruijn graph with frequency count."""

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
    def __init__(self, km1mer_tuple):
        self.label = km1mer_tuple[0] + km1mer_tuple[1][-1:]

    def __eq__(self, other):
        return self.label == other.label


def read_reads(fname):
    """ Read short reads in FASTA format. It is assumed that one line in the input file correspond to one read. """
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    reads = []

    for line in lines:
        if line[0] != '>':
            reads = reads + [line.rstrip()]

    return reads


def construct_graph(reads, k):
    """Construct De Bruijn graph from sets of short reads with k-length word, counting node repetitions."""
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
    """ Perform searching for Eulerian path in the graph to output genome assembly"""
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

fname = 'g200reads.fa'
reads = read_reads(fname)
g = construct_graph(reads, 5)
contig = output_contigs(g)
print(contig)
