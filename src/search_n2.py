import pandas as pd


def two_compatible_2tuples(I: list[tuple[int,int]], J: list[tuple[int,int]], dist:list[int]):
    j = 0
    lj = len(J)
    res = []
    for a in I:
        if j >= lj:
            break

        d = J[j][0] - a[1]
        while d < dist[0]:
            j += 1
            if j >= lj:
                break
            d = J[j][0] - a[1]

        jj = j
        while dist[0] <= d < dist[1]:
            res.append((a, J[jj]))
            jj += 1
            if jj >= lj:
                break
            d = J[jj][0] - a[1]
    return res


# reference implementation
# def _resolve(res: list[tuple[tuple[int,int],tuple[int, int],...]], resout: list[tuple[tuple[int,int],tuple[int, int]]]):
#     nres = []
#     for r in res:
#         ridx = [r[-1] == j[0] for j in resout]
#         if any(ridx):
#             for e, j in enumerate(ridx):
#                 if j:
#                     nres.append((*r, resout[e][-1]))
#     return nres

# reference implementation
# def n_compatible_2tuples(n_lists, dists):
#     """Find compatible 2-tuples between given lists
#     assumptions for tuple (a, b):
#     - all tuples within a list must be same width that is b-a is constant for all tuples in a list
#     - b > a for all tuples in a list
#     """
#     # initialize search by extracting first compatible tuples
#     res = two_compatible_2tuples(n_lists[0], n_lists[1], dists[0])
#
#     for i in range(2, len(n_lists)):
#         resout = two_compatible_2tuples(
#             sorted({i[-1] for i in res}),
#             n_lists[i],
#             dists[i-1]
#         )
#
#         res = _resolve(res, resout)
#     return res


def n_compatible_2tuples(n_lists, dists):
    """Find compatible 2-tuples between given lists
    assumptions for tuple (a, b):
    - all tuples within a list must be same width that is b-a is constant for all tuples in a list
    - b > a for all tuples in a list
    """
    # make resolve in pandas
    # ok, this is good, use this onwards and cleanup it
    # todo: cleanup, maybe rewrite the downstream code

    # initialize search by extracting first compatible tuples
    results = []
    res = two_compatible_2tuples(n_lists[0], n_lists[1], dists[0])
    si = sorted({i[-1] for i in res})
    results.append(res)
    for i in range(2, len(n_lists)):
        res = two_compatible_2tuples(
            si,
            n_lists[i],
            dists[i-1]
        )
        si = sorted({i[-1] for i in res})
        results.append(res)

    a = pd.DataFrame(results[0], columns=[0, 'x'])
    for i in range(1, len(results)):
        b = pd.DataFrame(results[i], columns=['x', 'y'])
        r = pd.merge(a, b, on='x')
        a = r.rename(columns={'x': i}).rename(columns={'y': 'x'})

    a.rename(columns={'x': len(results)-1}, inplace=True)
    return [a.iloc[i].to_list() for i in range(len(a))]


def wrapper(n_lists, dists):
    # sort input
    n_lists = [sorted(l) for l in n_lists]

    # check b-a=C and b>a
    for _l in n_lists:
        _c = {b-a for a, b in _l}
        if len(_c) > 1:
            raise ValueError(f'All tuples in one of the lists are not the same width (widths: {_c}).')
        if any(i <= 0 for i in _c):
            raise ValueError(f'Tuples (a,b) appears not to satisfy "b>a" assumption.')

    return n_compatible_2tuples(n_lists, dists)


if __name__ == "__main__":
    import itertools
    dist1 = (12, 18)
    dist2 = (10, 15)

    nsize = 20

    A = list((i, i+6) for i in range(5, 15+nsize))
    B = list((i, i+6) for i in range(11, 30+nsize))
    C = list((i, i+7) for i in range(25, 41+nsize))

    ref = [(a, b, c) for a, b, c in itertools.product(A, B, C) if b[0]-a[1] in range(*dist1) and c[0]-b[1] in range(*dist2)]

    test_lists = [A, B, C]
    test_dists = [dist1, dist2]

    res = wrapper(test_lists, test_dists)

    print('ref:', ref)
    print('res:', res)
    assert res == ref

    #ab = compatible_tuples(A, B, dist1)
    #bc = compatible_tuples(B, C, dist2)
    #assert ab == ref1
    #assert bc == ref2

