from random import random
from random import sample
from collections import deque
import numpy as np


class Graph(object):
    def __init__(self, graph_dict={}):
        """ initializes a graph object """
        self.__graph_dict = graph_dict

    def vertices(self):
        """ returns the vertices of a graph """
        return sorted(list(self.__graph_dict.keys()))

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def __len__(self):
        return len(self.edges())

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in
        self.__graph_dict, a key "vertex" with an empty
        list as a value is added to the dictionary.
        Otherwise nothing has to be done.
        """
        if vertex not in self.vertices():
            self.__graph_dict[vertex] = []

    def add_edge(self, edge):
        """ It the edge is not in self.__graph_dict,
        the vertices of the edge are added to each other keys
        The function assumes that edge is of type set, tuple or list;
        (no multiple edges)
        """
        if edge[1] not in self.__graph_dict[edge[0]]:
            self.__graph_dict[edge[0]].append(edge[1])
        if edge[0] not in self.__graph_dict[edge[1]]:
            self.__graph_dict[edge[1]].append(edge[0])

    def __generate_edges(self):
        """ A static method generating the edges of the
        graph "graph". Edges are represented as sets
        two vertices (no loop)
        """
        edges = []
        for vertex_in in self.vertices():
            for vertex_out in self.__graph_dict[vertex_in]:
                if vertex_in < vertex_out:
                    edges.append((vertex_in, vertex_out))
        return edges

    def vertex_degree(self):
        """ returns a list of sets containing the
        name and degree of each vertex of a graph """
        return [(vertex, len(self.__graph_dict[vertex])) for vertex in self.vertices()]

    def degree_sequence(self):
        """ returns as a non-increasing list sequence of the vertex degrees of a graph """
        return [degree for _, degree in sorted(self.vertex_degree(), key=lambda x: x[1], reverse=True)]

    def find_isolated_vertices(self):
        """ returns the list of zero-degree vertices of a graph """
        return [vertex for vertex, degree in self.vertex_degree() if degree == 0]

    def density(self):
        """ returns the density of a graph """
        nbr_edges, nbr_vertices = float(len(self.edges())), float(len(self.vertices()))
        return 2 * nbr_edges / (nbr_vertices * (nbr_vertices - 1))

    def adjacency_matrix(self):
        """ returns the ajacency matrix of a graph
         in the form of a numpy array"""
        edges = self.edges()
        n = len(self.vertices())
        adj_matrix = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j and (self.vertices()[i], self.vertices()[j]) in edges:
                    adj_matrix[i][j], adj_matrix[j][i] = 1, 1
        return np.array(adj_matrix)

    def global_clustering_coeff(self):
        """ returns the global clustering coefficient of a graph """
        adj_mtrx = self.adjacency_matrix()
        path_length_two = np.linalg.matrix_power(adj_mtrx, 2)
        closed_triple_mtrx = np.multiply(adj_mtrx, path_length_two)
        n = len(self.vertices())
        nbr_closed_triple, nbr_all_triple = 0.0, 0.0  # float because of division
        nbr_closed_triple += sum(closed_triple_mtrx[i][e] for i in range(n) for e in range(n) if i != e)
        nbr_all_triple += sum(path_length_two[i][e] for i in range(n) for e in range(n) if i != e)
        # instead of not computing the diagonal
        # we could have substract np.trace(path_length_two) from nbr_triple
        return nbr_closed_triple / nbr_all_triple if nbr_all_triple != 0 else 0  # avoid 0 division

    def shortest_path(self, a, b):
        """ returns the shortest path distance between two given vertices a, b"""
        queue = deque()
        distance = {a: 0}
        queue.append(a)
        while len(queue) > 0:
            current = queue.popleft()
            for vertex in self.__graph_dict[current]:
                if vertex == b:
                    return distance[current] + 1
                if vertex not in distance:
                    queue.append(vertex)
                    distance[vertex] = distance[current] + 1
        return float("inf")

    def connected_components(self):
        """ returns a list of sets composed of two elements: the vertices list and the size
        of each connected components of a graph """
        components = []
        set_vertices = set(self.vertices())
        queue = deque()
        while len(set_vertices) > 0:
            init = set_vertices.pop()
            visited = {init: True}
            queue.append(init)
            while len(queue) > 0:
                current = queue.popleft()
                for vertex in self.__graph_dict[current]:
                    if vertex not in visited:
                        queue.append(vertex)
                    visited[vertex] = True
            set_vertices -= set(visited.keys())
            components.append(list(visited.keys()))
        return zip(components, map(lambda e: len(e), components))

    def connected_component_elements(self):
        """ returns a list of the vertices list of each connected components of a graph """
        return map(lambda x: x[0], self.connected_components())

    def component_diameter(self, component):
        """ returns the diameter of a given connected component element of a graph"""
        diameter = 0
        for init in component:
            queue = deque()
            distance = {init: 0}
            queue.append(init)
            while len(queue) > 0:
                current = queue.popleft()
                for vertex in self.__graph_dict[current]:
                    if vertex not in distance:
                        queue.append(vertex)
                        distance[vertex] = distance[current] + 1
            diameter = max((diameter, max(distance.values())))

        return diameter

    def forest_diameters(self):
        """ returns the list of the diameter of each connected components of a graph """
        return [self.component_diameter(component) for component in self.connected_component_elements()]

    def biggest_component_diameter(self):
        """ returns the diameter of the biggest connected component of a graph """
        return self.component_diameter(max(self.connected_component_elements(), key=len))

    def component_spanning_tree(self, component):
        """ returns the spanning tree of a given connected component of a graph """
        spanning_tree = Graph({})
        queue = deque()
        spanning_tree.add_vertex(component.pop())
        queue.extend(spanning_tree.vertices())
        while len(queue) > 0:
            current = queue.popleft()
            for vertex in self.__graph_dict[current]:
                if vertex not in spanning_tree.vertices():
                    queue.append(vertex)
                    spanning_tree.add_vertex(vertex)
                    spanning_tree.add_edge((current, vertex))
        return spanning_tree

    def spanning_forest(self):
        """ returns the list of spanning trees of each connected component of a graph """
        return [self.component_spanning_tree(component) for component in self.connected_component_elements()]

    def biggest_component_spanning_tree(self):
        """ returns the spanning tree of a the biggest connected component of a graph """
        return self.component_spanning_tree(max(self.connected_component_elements(), key=lambda c: len(c)))

    def __str__(self):
        """ A better way for printing a graph """
        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
            res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res

nucleotides = ["A", "T", "G", "C"]
codontable = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'stop', 'TAG': 'stop',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'stop', 'TGG': 'W'}

codon_graph = Graph({})


def distance_str(str_1, str_2):
    assert len(str_1) == len(str_2)
    return len([0 for i, j in zip(str_1, str_2) if i != j])

for pos_1 in nucleotides:
    for pos_2 in nucleotides:
        for pos_3 in nucleotides:
            codon_graph.add_vertex(pos_1 + pos_2 + pos_3)

for codon_ref in codon_graph.vertices():
    for codon_alt in codon_graph.vertices():
        if distance_str(codon_ref, codon_alt) == 1:
            codon_graph.add_edge((codon_ref, codon_alt))

assert len(codon_graph.vertices()) == 64
assert len(codon_graph.edges()) == 64*9//2
assert len(codon_graph.biggest_component_spanning_tree()) == 64 - 1
assert codon_graph.biggest_component_diameter() == 3
assert codon_graph.degree_sequence() == [9] * 64

from pepdata import *
print(amino_acid.amino_acid_letters)
print(amino_acid.hydropathy("K"))
print(amino_acid.hydropathy)
print(amino_acid.volume)
print(amino_acid.polarity)
print(amino_acid.pK_side_chain)
print(amino_acid.prct_exposed_residues)
print(amino_acid.hydrophilicity)
print(amino_acid.accessible_surface_area)
print(amino_acid.local_flexibility)
print(amino_acid.accessible_surface_area_folded)
print(amino_acid.refractivity)
