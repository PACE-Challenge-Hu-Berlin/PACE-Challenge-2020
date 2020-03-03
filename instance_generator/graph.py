import numpy as np
import numpy.random
import trees


class Graph:
    """
    A representation of an undirected graph
    Attributes:
        vertices : int
            total amount of vertices in the graph
        edges : ndarray
            the adjacency matrix containing edges in both directions
        edge_amount : int
            total amount of edges in the graph
        tree_depth : int
            the longest path in the generated tree
        seed : int
            the seed of the random number generator
    """
    def __init__(self, vertices, generator=trees.random_tree, p_general=0.5, p_clique=1, seed=None):
        """
        Creates a graph build from a generated tree and determines the tree_depth

        Parameters
            vertices : int
                the amount of vertices
            generator : function (int, int) -> (trees.TreeNode, trees.TreeNode[])
                the generating function, that creates a tree with the
                given amount of vertices (and the given seed)
                and outputs the root and the whole tree (default: trees.random_tree)
            p_general : float
                the probability that an edge should be added for allowed pairs of vertices (default: 0.5)
            p_clique : float
                the probability that an edge is added on the path of the longest path (default: 1)
            seed : int
                the seed for the random number generator
        """

        self.vertices = vertices
        self.edges = np.zeros((vertices, vertices))
        self.edge_amount = 0
        self.seed = seed

        if self.seed is None:
            self.seed = np.random.randint(2**31)

        np.random.seed(self.seed)

        root, tree = generator(vertices, np.random.randint(2**31))

        # lets find paths to all nodes of the tree
        paths = [[] for i in range(vertices)]
        open_ = [(root, [root.i])]

        while len(open_) > 0:
            node, path = open_.pop()
            if paths[node.i]:
                continue
            paths[node.i] = path
            for node2 in node.children:
                open_.append((node2, path + [node2.i]))

        for i in range(vertices):
            if not paths[i]:
                print("Something went wrong while analyzing the tree")
                return

        # transfer tree into graph
        for i in range(vertices):
            for child in tree[i].children:
                self.add_edge(i, child.i)

        # clutter the tree alot
        for i in range(vertices):
            for j in range(i+1, vertices):
                # determine if they are on the same path:
                if i in paths[j] or j in paths[i]:
                    if np.random.random() < p_general:
                        self.add_edge(i, j)

        # form a clique (or just part of it) on the longest path
        self.tree_depth, child = root.calculate_tree_depth()

        longest_path = paths[child.i]
        for i in range(len(longest_path)):
            for j in range(i+1, len(longest_path)):
                self.remove_edge(i, j)
                if np.random.random() < max(p_general, p_clique):
                    self.add_edge(i, j)

    def remove_edge(self, i, j):
        if self.edges[i, j] != 0:
            self.edge_amount -= 1
            self.edges[i, j] = 0
            self.edges[j, i] = 0

    def add_edge(self, i, j):
        if self.edges[i, j] == 0:
            self.edge_amount += 1
            self.edges[i, j] = 1
            self.edges[j, i] = 1

    def export(self):
        """Exports the graph into a PACE compatible string"""
        result = "p tdp {} {}\n".format(self.vertices, self.edge_amount)

        for i in range(self.vertices):
            for j in range(i+1, self.vertices):
                if self.edges[i, j] != 0:
                    result += "{} {}\n".format(i+1, j+1)

        return result

    def permute(self):
        """Randomly permutes the labels of the vertices"""
        permutation = np.arange(self.vertices)
        np.random.shuffle(permutation)

        old_copy = self.edges.copy()

        for i in range(self.vertices):
            for j in range(self.vertices):
                self.edges[i, j] = old_copy[permutation[i], permutation[j]]


