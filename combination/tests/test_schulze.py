import operator
import random

import numpy as np

from intogen_combination.schulze_election import combination_ranking
from intogen_combination.schulze_strongest_path import strongest_path as python_strongest_path
from intogen_combination.schulze_strongest_path_cython import strongest_path


def legacy_combination_ranking(ballot_dict, weights):
    """Loop-based implementation (before vectorisation)"""
    all_candidates = list(set.union(*[set(ballot_dict[voter].keys()) for voter in ballot_dict]))
    idx = {k: i for i, k in enumerate(all_candidates)}
    size = len(all_candidates)
    pref = np.zeros(size ** 2)
    spath = np.zeros(size ** 2)
    for voter in ballot_dict:
        d = ballot_dict[voter]
        for i in all_candidates:
            if i not in d.keys():
                for j in d:
                    pref[idx[j] * size + idx[i]] += weights[voter]
            else:
                for j in d:
                    if d[j] < d[i]:
                        pref[idx[j] * size + idx[i]] += weights[voter]
    python_strongest_path(size, pref, spath)
    scores = {}
    for i in range(size):
        scores[all_candidates[i]] = sum(1 for j in range(size) if spath[i * size + j] < spath[j * size + i])
    sorted_scores = sorted(scores.items(), key=operator.itemgetter(1), reverse=True)
    ranking, prev_score, prev_rank, counter = {}, None, None, 1
    while sorted_scores:
        c = sorted_scores.pop()
        if prev_score is None:
            ranking[c[0]], prev_score, prev_rank = 1, c[1], 1
        elif prev_score == c[1]:
            ranking[c[0]] = prev_rank
        elif prev_score < c[1]:
            ranking[c[0]], prev_score, prev_rank = counter, c[1], counter
        counter += 1
    return ranking


def test_cython_strongest_path_is_identical_to_python():
    for n in [1, 2, 5, 40]:
        pref = np.random.default_rng(n).choice([0, 0.05, 0.1, 0.35, 0.7], n * n).astype(float)
        a, b = np.zeros(n * n), np.zeros(n * n)
        python_strongest_path(n, pref, a)
        strongest_path(n, pref, b)
        assert np.array_equal(a, b)


def test_combination_ranking_is_identical_to_legacy():
    rng = random.Random(3)
    for _ in range(25):
        genes = [f'G{i}' for i in range(rng.randint(3, 80))]
        ballots = {}
        for m in range(rng.randint(1, 9)):
            if rng.random() < 0.1:
                ballots[f'm{m}'] = {}
                continue
            selected = rng.sample(genes, rng.randint(1, min(40, len(genes))))
            ballots[f'm{m}'] = {g: rng.randint(1, 45) for g in selected}   # ties allowed
        if all(len(b) == 0 for b in ballots.values()):
            continue
        weights = {m: rng.choice([0.0, 0.05, 0.1, 0.15000000000000002, 0.3]) for m in ballots}
        assert combination_ranking(ballots, weights) == legacy_combination_ranking(ballots, weights)
