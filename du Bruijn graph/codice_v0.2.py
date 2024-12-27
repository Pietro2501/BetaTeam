class Node:
    """ Class Node to represent a vertex in the de bruijn graph """
    def __init__(self, lab):
        self.label = lab
        self.indegree = 0
        self.outdegree = 0

class Edge:
    def __init__(self, lab):
        self.label = lab

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
    """ Construct de bruijn graph from sets of short reads with k length word"""
    edges = dict()
    vertices = dict()

    for read in reads:
        i = 0
        while i+k < len(read):
            v1 = read[i:i+k]
            v2 = read[i+1:i+k+1]
            if v1 in edges.keys():
                vertices[v1].outdegree += 1
                # for instance in edges[v1]:
                #     if instance.label != v2:
                edges[v1].append(Edge(v2))
            else:
                vertices[v1] = Node(v1)
                vertices[v1].outdegree += 1
                edges[v1] = []
                edges[v1].append(Edge(v2))
            if v2 in edges.keys():
                vertices[v2].indegree += 1
            else:
                vertices[v2] = Node(v2)
                vertices[v2].indegree += 1
                edges[v2] = []
            i += 1

    return (vertices, edges)

def output_contigs(g):
    """ Perform searching for Eulerian path in the graph to output genome assembly"""
    V = g[0]
    E = g[1]
    # Pick starting node (the vertex with zero in degree)
    start = list(V.keys())[0]
    for k in V.keys():
        if V[k].indegree < V[start].indegree:
            start = k

    contig = start
    current = start
    while len(E[current]) > 0:
        edges_current = {i.label for i in E[current]}
        if len(edges_current)> 1:
            best = E[current][0]
            counter = 0
            best_position = 0
            for instance in E[current][1:]:
                if V[instance.label].outdegree > V[best.label].outdegree:
                # if V[instance.label].indegree + V[instance.label].outdegree > V[best.label].indegree + V[best.label].outdegree:
                    best = instance
                    best_position = counter
                counter += 1
            del E[current][best_position]
            contig += best.label[-1]
            current = best.label
        else:
            next_edge = E[current][0]
            del E[current][0]
            contig += next_edge.label[-1]
            current = next_edge.label

    return contig

def print_graph(g):
    """ Print the information in the graph to be (somewhat) presentable """
    V = g[0]
    E = g[1]
    for k in V.keys():
        print (f"Node name: {V[k].label}. indegree: {V[k].indegree}. outdegree: {V[k].outdegree}")
        print ("Edges: ")
        edges=[]
        for e in E[k]:
            # print(e.label)
            edges.append(e.label)
        print(set(edges))