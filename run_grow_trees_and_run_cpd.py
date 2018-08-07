from time import time
import os
import sys
from sklearn.externals.joblib import dump

from garden import ChangePointTreeGrower, get_attachment_fun
from sim_utils import get_trees_fname, get_diffstat_fname
from change_point_estimation import compute_difference_statistics


def main(params, i):
    data_dir = 'data/'
    print(params)
    treepath = os.path.join(data_dir, get_trees_fname(params))
    treepath += '_r{}'.format(i)
    print(treepath)

    start_time = time()
    tree = grow_tree(params)
    dump({'params': params, 'tree': tree}, treepath)
    print('growing tree took {} seconds'.format(time() - start_time))

    dpath = os.path.join(data_dir, get_diffstat_fname(params))
    dpath += '_r{}'.format(i)
    dout = compute_difference_statistics(tree, h=params['h'])
    dump({'params': params, 'dout': dout}, dpath)
    print('change point took {} seconds'.format(time() - start_time))


def grow_tree(params):
    funs = [get_attachment_fun(**params['f0']),
            get_attachment_fun(**params['f1'])]

    cp_tree_grower = ChangePointTreeGrower(n=params['n'], funs=funs,
                                           change_point_props=params['cp'])

    return cp_tree_grower.grow_cp_tree()


if __name__ == "__main__":
    params = {'n': int(2e2), 'cp': .5, 'h': 'loglog',
              'f0': {'family': 'affine', 'alpha': 1},
              'f1': {'family': 'polynomial', 'alpha': 1, 'p': .5}}

    i = sys.argv[1]
    start_time = time()
    main(params, i)
    print('simulation took {} seconds'.format(time() - start_time))
