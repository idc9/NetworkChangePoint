def dict2str(d):
    s = ''
    first = True
    for k in sorted(d.keys()):
        if not first:
            s += '_'
        first = False

        if type(d[k]) == dict:
            s += '{}_{}'.format(k, dict2str(d[k]))
        else:
            s += '{}{}'.format(k, d[k])

    return s


def get_trees_fname(params):
    return 'trees_{}'.format(dict2str(params))


def get_diffstat_fname(params):
    return 'diffstats_{}'.format(dict2str(params))
