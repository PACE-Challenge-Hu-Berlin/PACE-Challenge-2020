import numpy as np
import numpy.random


class TreeNode:
    """A node in a tree

    Attributes:
        children : TreeNode[]
            the children of this node
        i : int
            the label of this node
    """
    def __init__(self, i):
        self.children = []
        self.i = i

    def add_child(self, node):
        for c in self.children:
            if c == node:
                return
        self.children.append(node)

    def calculate_tree_depth(self):
        """
        Calculate the deepest children in the subtree rooted
        at this node and it's corresponding depth

        Returns
            (int, TreeNode):
                the depth and leaf of maximal depth
        """

        # do a depth first search for the deepest node
        open_ = [(self, 1)]
        done = []
        max_node, max_depth = None, None

        while len(open_) > 0:
            min_index = -1
            min_depth = 0
            for i in range(len(open_)):
                node2, depth2 = open_[i]
                if min_index == -1 or depth2 < min_depth:
                    min_depth = depth2
                    min_index = i

            node, depth = open_.pop(min_index)

            if node in done:
                continue
            done.append(node)
            if max_node is None or depth > max_depth:
                max_node, max_depth = node, depth
            for node2 in node.children:
                open_.append((node2, depth + 1))

        return max_depth, max_node


def random_tree(vertices, seed):
    """
    Generates a random tree by generating a random
    Pruefer sequence and calculating the corresponding
    tree
    """
    np.random.seed(seed)

    tree = [TreeNode(i) for i in range(vertices)]
    if vertices == 1:
        return tree[0], tree

    # build the random tree
    pruefer = np.random.random_integers(0, vertices - 1, vertices - 2)
    d = np.zeros(vertices)
    for i in pruefer:
        d[i] += 1

    for i in pruefer:
        for j in range(vertices):
            if d[j] == 0:
                tree[j].add_child(tree[i])
                tree[i].add_child(tree[j])
                d[i] -= 1
                d[j] -= 1
                break
    v, w = -1, -1
    for i in range(vertices):
        if d[i] == 0:
            if v == -1:
                v = i
            else:
                w = i

    tree[w].add_child(tree[v])
    tree[v].add_child(tree[w])

    # take a random root
    root = np.random.randint(vertices)

    return tree[root], tree


def path(vertices, seed):
    tree = [TreeNode(i) for i in range(vertices)]
    for i in range(vertices - 1):
        tree[i].add_child(tree[i+1])

    return tree[0], tree


def star(vertices, seed):
    tree = [TreeNode(i) for i in range(vertices)]
    for i in range(vertices - 1):
        tree[0].add_child(tree[i+1])

    return tree[0], tree
