from collections import Counter
import numpy as np


def get_degree_counts(degrees):
    """
    Given a list of degrees, returns the empirical degree counts.

    Parameters
    ----------
    degrees (list): degrees of each node

    Output
    ------
    degree_counts (list) where degree_counts[m] is the number of nodes
    with degree m where m can range between 0 and n - 1
    """

    deg_counter = Counter(degrees)

    # deg_counts = [0 for _ in range(max(deg_counter.keys()) + 1)]
    deg_counts = np.zeros(len(degrees))
    for d in deg_counter.keys():
        deg_counts[d] = deg_counter[d]

    return deg_counts


def get_pdf_pairs(s1, d1, s2, d2):
    """
    Returns the empirical PDFs for two observed degree counts
    Parameters
    ----------
    s1 (int): size of graph 1
    d1 (list): out degree counts of graph 1
    s2 (int): size of graph 1
    d2 (list): out degree counts of graph 1
    """

    p1_full = np.zeros(max(s1, s2))
    p1_full[0:len(d1)] = d1
    p1_full = p1_full/s1
    p2_full = np.zeros(max(s1, s2))
    p2_full[0:len(d2)] = d2
    p2_full = p2_full/s2

    return p1_full, p2_full
