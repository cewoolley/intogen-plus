import operator

import numpy as np

# Cython version
from intogen_combination.schulze_strongest_path_cython import strongest_path

# Python version
#from intogen_combination.schulze_strongest_path import strongest_path


def combination_ranking(ballot_dict, weights):

    # INIT
    all_candidates = list(set.union(*[set(ballot_dict[voter].keys()) for voter in ballot_dict]))
    all_candidates_to_idx = {k: i for i, k in enumerate(all_candidates)}
    size = len(all_candidates)

    pref = np.zeros(size**2, dtype=np.float64)  #[0.0] * (size**2)
    spath = np.zeros(size**2, dtype=np.float64)  #[0.0] * (size**2)

    # PREPARE
    if weights is None:
        weights = {}
        for voter in ballot_dict:
            weights[voter] = 1.
        weights = dict(weights)

    # A voter prefers j over i if j is in its ballot and i is either not in the ballot or ranked worse.
    # The weight of each voter is added in the same order as looping over voters and pairs.
    for voter in ballot_dict:
        d = ballot_dict[voter]
        if len(d) == 0:
            continue
        ranks = np.full(size, np.inf)
        ranks[[all_candidates_to_idx[j] for j in d]] = [d[j] for j in d]
        prefers = ranks[:, None] < ranks[None, :]  # rows: j (preferred), columns: i
        pref += prefers.ravel() * weights[voter]

    # STRONGEST PATH
    strongest_path(size, pref, spath)

    # COMBINATION RANKING
    s = spath.reshape(size, size)
    scores = (s < s.T).sum(axis=1)
    scores_dict = {}
    for i in range(size):
        scores_dict[all_candidates[i]] = int(scores[i])
    sorted_scores = sorted(scores_dict.items(), key=operator.itemgetter(1), reverse=True)

    ranking = {}
    prev_score = None
    prev_rank = None
    counter = 1
    while len(sorted_scores) > 0:
        c = sorted_scores.pop()
        if prev_score is None:
            ranking[c[0]] = 1
            prev_score = c[1]
            prev_rank = 1
        elif prev_score == c[1]:
            ranking[c[0]] = prev_rank
        elif prev_score < c[1]:
            ranking[c[0]] = counter
            prev_score = c[1]
            prev_rank = counter
        counter += 1

    return ranking
