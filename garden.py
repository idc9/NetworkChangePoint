import numpy as np
import networkx as nx

from utils import get_degree_counts
from sklearn.externals import joblib


class Tree(object):
    """
    Nodes are zero indexed
    """
    def __init__(self):
        """
        Initializes a tree with 1 vertices and no edges.

        For a tree of size n, the nodes are named: 0, 1, ..., n - 1.

        Attributes
        ----------
        edges (list): parent of node m is edges[m - 1]

        out_degrees (list): the out degree of node m is out_degrees[m]

        change_points (None, float, list): keeps track of change points
        """
        self.edges = []
        self.out_degrees = [0]

        self.change_points = None

    def __repr__(self):
        return 'Tree with %d nodes' % self.n_nodes

    def save(self, fname):
        joblib.dump(self.__dict__, fname)

    @classmethod
    def load(cls, fname):
        data = joblib.load(fname)
        tree = Tree()
        tree.__dict__ = data
        return tree

    @classmethod
    def from_edgelist(self, edges):
        tree = Tree()
        for parent in edges:
            tree.add_node(parent)

        return tree

    @property
    def n_nodes(self):
        return len(self.edges) + 1

    def out_degree_counts(self):
        """
        Returns a list of out degree counts.
        """
        return get_degree_counts(self.out_degrees)

    def add_node(self, parent):
        """
        Adds a new node to the tree by specifying its parent vertex
        """
        assert (0 <= parent) and (parent < self.n_nodes)
        self.edges.append(parent)  # attach new node to parent vertex

        self.out_degrees.append(0)  # new node has out degree 0
        self.out_degrees[parent] += 1  # parent's out degree increases

    def to_graph(self):
        """
        Returns the tree as a networkx object
        """

        G = nx.DiGraph()
        G.add_node(0)
        for i in range(1, self.n_nodes):
            G.add_node(i)
            G.add_edge(self.edges[i - 1], i)

        return G

    def tree_generator(self, limit=None):
        if limit is None:
            limit = self.n_nodes
        # else:
            # assert (0 <= limit) and (limit <= self.n_nodes)

        t = Tree()
        yield t
        # for k in range(1, self.n_nodes - 1):
        for k in range(limit - 1):
            t.add_node(self.edges[k])
            yield t

    def get_outdegree_counts_history(self, sizes=None):
        """
        Returns the empirical degree distribution at each point in the
        tree's history

        Parameters
        ----------
        sizes (None, list, int): the sizes at which to return the degree counts

        """

        if sizes is not None and type(sizes) == int:
            sizes = [sizes]

        for t in self.tree_generator():
            if (sizes is None) or (t.n_nodes in sizes):
                yield t.n_nodes, t.out_degree_counts()


class TreeGrower(object):

    def __init__(self, f, seed_tree=None):
        """
        Grows a single Generalized Attachment Process tree.

        Parameters
        ----------
        f (calllable): the attachment function

        seed_tree (Tree, None): the tree to initialize from. If None the will initalize tree
        with 2 nodes and 1 edge
        """

        if seed_tree is None:
            self.tree = Tree()
        else:
            self.tree = seed_tree

        self.f = f

        # node weight is the f(d) where d is the node's out degree
        self.node_weights = np.array([self.f(d) for d in self.tree.out_degrees])
        self.total_weight = sum(self.node_weights) + 0.0  # sum of node weights

    def grow_new_node(self):
        """
        Adds a new now to the tree by randomly selecting a parent node
        """

        # sample from 0,..., (number of nodes - 1) according
        # to weight distribution
        parent = np.random.choice(self.tree.n_nodes,
                                  p=self.node_weights/self.total_weight)

        self.add_node_and_update_weights(parent)

    def grow_many_nodes(self, num):
        for _ in range(num):
            self.grow_new_node()

    def grow_to_size(self, n):
        num_to_add = n - self.tree.n_nodes
        assert num_to_add > 0

        self.grow_many_nodes(num_to_add)

    def add_node_and_update_weights(self, parent):
        """
        Updates the tree and the node weights
        """

        old_weight = self.node_weights[parent]

        self.tree.add_node(parent)  # add node to tree

        new_weight = self.f(self.tree.out_degrees[parent])
        self.node_weights[parent] = new_weight  # update parent's node weight

        # update node weights and total weight
        self.node_weights = np.append(self.node_weights, self.f(0))
        self.total_weight += self.f(0) + new_weight - old_weight


def cp_arg_checker(n, funs, change_point_sizes):

    # should have one more function than change points
    assert len(funs) == len(change_point_sizes) + 1

    # make sure each funstion is a function
    for i in range(len(funs)):
        assert callable(funs[i])

    # make sure the first change point size is at least 3
    assert change_point_sizes[0] >= 3

    for i in range(len(change_point_sizes)):

        if i != 0:  # make sure change point sizes are in creasing
            assert change_point_sizes[i - 1] < change_point_sizes[i]

        # the change ponints should all happen stricly before n
        assert change_point_sizes[i] < n


class ChangePointTreeGrower(object):

    def __init__(self, n, funs, change_point_props):

        if type(change_point_props) in [int, float]:
            change_point_props = [change_point_props]

        self.n = n
        self.funs = funs
        self.change_point_props = change_point_props
        self.change_point_sizes = [int(n * p) for p in change_point_props]
        cp_arg_checker(n, funs,  self.change_point_sizes)

        self.tree = Tree()
        self.tree.change_points = self.change_point_sizes

    def grow_cp_tree(self):

        for i in range(len(self.funs)):

            if i == len(self.funs) - 1:  # final segment
                new_size = self.n
            else:
                new_size = self.change_point_sizes[i]

            # grow the tree using the new function
            tree_grower = TreeGrower(f=self.funs[i],
                                     seed_tree=self.tree)
            tree_grower.grow_to_size(new_size)

            # update the tree
            self.tree = tree_grower.tree

        return self.tree

    def get_cp_degree_counts(self):
        """
        Returns the empirical degree counts at the end of each change point
        """
        deg_counts = []
        sizes = self.change_point_sizes
        sizes.append(self.n)
        for s, dc in self.tree.get_outdegree_counts_history(sizes=sizes):
            deg_counts.append(dc)

        return {sizes[i]: deg_counts[i] for i in range(len(sizes))}


def get_attachment_fun(family='affine', alpha=1.0, p=.5):
    """

    Parameters
    ----------
    family (str): ['affine', 'polynomial', 'uniform']


    returns a function f
    def f(d):
        # d = out-edgree
        return -----

    """

    if family == 'affine':
        assert alpha > -1

        def f(d):
            """
            Parameters
            ----------
            d (int): out-degree of a given vertex

            Returns
            --------
            d + 1 + alpha
            """
            return d + 1 + alpha

        return f

    elif family == 'polynomial':
        assert alpha > -1

        def f(d):
            """
            Parameters
            ----------
            d (int): out-degree of a given vertex

            Returns
            --------
            (d + alpha)^e
            """
            return (d + 1 + alpha) ** p

        return f

    elif family == 'uniform':
        def f(d):
            """
            Parameters
            ----------
            d (int): out-degree of a given vertex

            Returns
            -------
            1
            """
            return 1

        return f

    else:
        raise ValueError('%s is not a valid function family' % family)
