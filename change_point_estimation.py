from scipy.sparse import dok_matrix, csr_matrix
import numpy as np
from collections import Counter


def compute_difference_statistics(tree, h='loglog'):
    """

    Parameters
    ----------
    tree (Tree object): computes the differnce statistc for each size of
        the tree

    h (callable, str, int): the size of the tree to start computing
    the differene statistics (this is n/h_n in the paper).
        if str then must be one of ['log', 'loglog', 'nto1overloglog']
        if calllable then anchor will be set to n/h(n)
        if int then anchor will be set to h


    Output
    ------
    diff_stats (d_n(m) from the paper)
    anchor = n/h(n)
    """
    n = tree.n_nodes

    anchor = int(get_anchor(h, n))

    diff_stats = np.array([0.0 for _ in range(n)])
    # leaf_diffs = np.array([0.0 for _ in range(n)])

    for t in tree.tree_generator():
        m = t.n_nodes - 1

        if m == anchor:
            out_degrees = t.out_degrees
            anchor_dd = sp_from_dict(Counter(out_degrees), n)/(anchor + 1)
            # anchor_leaf_prop = anchor_dd[1, 0]
            deg_dist = sp_from_dict(Counter(out_degrees), n)
            out_degrees = sp_deg_dist(out_degrees, n)

        elif m > anchor:
            # update degree distribution/outdegree trackers
            parent = tree.edges[m - 1]
            parent_degree = out_degrees[parent, 0]
            deg_dist -= sp_unit_vector(parent_degree, n)
            deg_dist += sp_unit_vector(parent_degree + 1, n)
            deg_dist += sp_unit_vector(0, n)
            # leaf_prop = deg_dist[1, 0]/(m + 1)
            out_degrees += sp_unit_vector(parent, n)

            # compute degree distribution difference
            diffs = abs(deg_dist/(m + 1) - anchor_dd)
            degs_non_zero = diffs.nonzero()[0]
            mult = sp_from_dict({k: .5 ** k for k in degs_non_zero}, n)
            diff_stats[m] = diffs.T.dot(mult)[0, 0]

            # leaf_diffs[m] = abs(leaf_prop - anchor_leaf_prop)

    return diff_stats, anchor


def estimate_change_point(diff_stats, h='loglog', b='nto1overloglog',
                          upcrossing='first'):
    """

    Parameters
    ----------
    diff_stats: the output of compute_difference_statistics

    h: the argument supplied to compute_difference_statistics to compute
    the anchor

    b: how to compute the threshold

    upcrossing (str): one of ['first', 'last']. Return either the first
    upcrossing of the threshold or the last upcrossing of the threshold.

    Output
    ------
    cp_estimate (int): the estimated change point

    threshold (float): the threshold value used to compute the change point
    """
    n = len(diff_stats)

    anchor = get_anchor(h, n)
    threshold = get_threshold(b, n)

    up_crossings = get_upcrossings(np.array(diff_stats), threshold)
    up_crossings = up_crossings[up_crossings > anchor]

    if(len(up_crossings) == 0):
        print('warning no change point detected')
        return None, threshold
    elif upcrossing == 'first':
        return int(min(up_crossings)), threshold

    elif upcrossing == 'last':
        return int(max(up_crossings)), threshold


def get_seq_fun(kind='nto1overloglog'):
    if kind == 'nto1overloglog':
        def f(n):
            return n ** (1/np.log(np.log(n)))

    elif kind == 'loglog':
        def f(n):
            return np.log(np.log(n))

    elif kind == 'log':
        def f(n):
            return np.log(n)

    elif kind == 'sqrtoverlog':
        def f(n):
            return np.sqrt(n)/np.log(n)

    else:
        raise ValueError('{} is improper kind'.format(kind))

    return f


def get_anchor(x, n):
    if type(x) == callable:
        return n/x(n)
    elif type(x) == str:
        f = get_seq_fun(x)
        return n/f(n)

    elif type(x) in [int, float]:
        return x

    else:
        raise ValueError('{} is improper argument'.format(x))


def get_threshold(x, n):
    if type(x) == callable:
        return 1/x(n)
    elif type(x) == str:
        f = get_seq_fun(x)
        return 1/f(n)

    elif type(x) in [int, float]:
        return x

    else:
        raise ValueError('{} is improper argument'.format(x))


def get_upcrossings(vals, threshold):
    a = np.array(vals) - threshold
    return np.where(np.diff(np.sign(a)) > 0)[0] + 1


def sp_from_dict(d, n):
    sp = dok_matrix((n, 1))
    for k in d.keys():
        sp[k] = d[k]
    return csr_matrix(sp)


def sp_deg_dist(od, n):
    out_degrees = dok_matrix((n, 1))
    out_degrees[0:len(od), 0] = csr_matrix(od).T
    return csr_matrix(out_degrees)


def sp_unit_vector(i, n):
    e = dok_matrix((n, 1))
    e[i] = 1
    return csr_matrix(e)
