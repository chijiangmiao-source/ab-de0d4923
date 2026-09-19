"""Brute-force cross-check on small random pedigrees + edge-case tests."""
from __future__ import annotations

import random

import numpy as np
import pytest

from app.model import parse_case
from app.phasing import (
    INF,
    MendelError,
    _build_allowed,
    _is_heterozygous,
    solve,
)

_RANK = {None: 0, "A": 1, "C": 2, "G": 3, "T": 4}


def tuple_key(a, fi, mi, w, n):
    """Independent lexicographic key (plain tuples) for cross-checking."""
    fh = a.father_haps[fi]
    mh = a.mother_haps[mi]
    key = [_RANK[fh[0]], _RANK[fh[1]], _RANK[mh[0]], _RANK[mh[1]]]
    for c in range(n):
        key.append((w >> (2 * c)) & 1)
        key.append((w >> (2 * c + 1)) & 1)
    return tuple(key)

ALLELES = ["A", "C", "G", "T", None, None]  # weight missing a bit


def rand_gt(rng):
    a = rng.choice(ALLELES)
    b = rng.choice(ALLELES)
    return (a, b)


def gt_str(gt):
    return "." if gt is None else gt


def make_payload(rng, m, n):
    markers = [f"mk{i}" for i in range(m)]
    father = [rand_gt(rng) for _ in range(m)]
    mother = [rand_gt(rng) for _ in range(m)]
    children = [[rand_gt(rng) for _ in range(m)] for _ in range(n)]
    return {
        "markers": markers,
        "father": {"genotypes": [[gt_str(x) for x in g] for g in father]},
        "mother": {"genotypes": [[gt_str(x) for x in g] for g in mother]},
        "children": [
            {"name": f"c{c}", "genotypes": [[gt_str(x) for x in g] for g in children[c]]}
            for c in range(n)
        ],
    }


def feasible_case(payload):
    case = parse_case(payload)
    n = len(case.children_names)
    word_idx = np.arange(1 << (2 * n), dtype=np.int64)
    per_site = []
    for site in case.sites:
        a = _build_allowed(site, n, word_idx)
        if a is None:
            return None
        opts = [(w, fi, mi) for w, pl in a.pairs().items() for fi, mi in pl]
        per_site.append((a, opts))
    return case, per_site


def brute_optimal_paths(case, per_site):
    """Enumerate ALL full paths (small cases only), return optimal ones."""
    m = len(case.sites)
    best = INF
    optimal = []

    def rec(si, prev_w, cost, path):
        nonlocal best, optimal
        if si == m:
            if cost < best:
                best = cost
                optimal = [path[:]]
            elif cost == best:
                optimal.append(path[:])
            return
        # prune: remaining edges each cost >= 0 (no admissible pruning needed)
        for w, fi, mi in per_site[si][1]:
            add = 0 if prev_w is None else (w ^ prev_w).bit_count()
            path.append((w, fi, mi))
            rec(si + 1, w, cost + add, path)
            path.pop()

    rec(0, None, 0, [])
    return best, optimal


def canonicalize_path(case, per_site, path):
    """Apply the same global anchor convention used by the solver."""
    m = len(case.sites)
    n = len(case.children_names)
    pat_mask = sum(1 << (2 * c) for c in range(n))
    mat_mask = pat_mask << 1

    def flip(which):
        for si in range(m):
            a = per_site[si][0]
            w, fi, mi = path[si]
            if which == "F":
                fi = a.flip_fi[fi]
                w ^= pat_mask
            else:
                mi = a.flip_mi[mi]
                w ^= mat_mask
            path[si] = (w, fi, mi)

    fa = next((si for si, s in enumerate(case.sites) if _is_heterozygous(s.father)), None)
    mo = next((si for si, s in enumerate(case.sites) if _is_heterozygous(s.mother)), None)
    if fa is not None:
        w, fi, mi = path[fa]
        fh = per_site[fa][0].father_haps[fi]
        if fh[0] > fh[1]:
            flip("F")
    if mo is not None:
        w, fi, mi = path[mo]
        mh = per_site[mo][0].mother_haps[mi]
        if mh[0] > mh[1]:
            flip("M")

    # Parents lacking any heterozygous marker keep their label-flip symmetry:
    # no pseudo-anchor is imposed, so symmetric optimal paths remain in the
    # optimal set and show up in the ambiguity masks.  The canonical response
    # is simply the lexicographically smallest such path (path-level order).
    return tuple(path)


def engine_result_tuples(result, case, per_site):
    """Turn the engine result back into (w, fi, mi) per site for comparison."""
    out = []
    for si, site_r in enumerate(result.sites):
        a = per_site[si][0]
        fi = a.father_haps.index(tuple(site_r.father_hap))
        mi = a.mother_haps.index(tuple(site_r.mother_hap))
        w = 0
        for c, (pb, mb) in enumerate(site_r.transmissions):
            w |= pb << (2 * c)
            w |= mb << (2 * c + 1)
        out.append((w, fi, mi))
    return tuple(out)


@pytest.mark.parametrize("seed", range(300))
def test_brute_force_small(seed):
    rng = random.Random(seed)
    m = rng.randint(2, 4)
    n = rng.randint(1, 2)
    payload = make_payload(rng, m, n)
    built = feasible_case(payload)
    case = parse_case(payload)
    if built is None:
        # engine must report the earliest infeasible marker
        earliest = None
        nn = len(case.children_names)
        widx = np.arange(1 << (2 * nn), dtype=np.int64)
        for si, site in enumerate(case.sites):
            if _build_allowed(site, nn, widx) is None:
                earliest = si
                break
        with pytest.raises(MendelError) as exc:
            solve(case)
        assert exc.value.marker_index == earliest
        return

    _case, per_site = built
    best, paths = brute_optimal_paths(case, per_site)
    if best >= INF:
        return

    can_paths = {canonicalize_path(case, per_site, list(p)) for p in paths}
    result = solve(case)
    assert result.min_crossovers == best

    chosen = engine_result_tuples(result, case, per_site)
    assert chosen in can_paths

    # canonical solution = lexicographically smallest canonical optimal path
    def path_key(p):
        return tuple(
            tuple_key(per_site[si][0], fi, mi, w, n)
            for si, (w, fi, mi) in enumerate(p)
        )

    expected = min(can_paths, key=path_key)
    assert chosen == expected

    # ambiguity masks: field varies across the canonical optimal set
    for si in range(m):
        a = per_site[si][0]
        w, fi, mi = chosen[si]
        f_vals = {a.father_haps[p[si][1]] for p in can_paths}
        m_vals = {a.mother_haps[p[si][2]] for p in can_paths}
        sr = result.sites[si]
        assert sr.father_fixed == (
            len({v[0] for v in f_vals}) == 1,
            len({v[1] for v in f_vals}) == 1,
        )
        assert sr.mother_fixed == (
            len({v[0] for v in m_vals}) == 1,
            len({v[1] for v in m_vals}) == 1,
        )
        for c in range(n):
            p_vals = {(p[si][0] >> (2 * c)) & 1 for p in can_paths}
            q_vals = {(p[si][0] >> (2 * c + 1)) & 1 for p in can_paths}
            assert sr.child_fixed[c] == (len(p_vals) == 1, len(q_vals) == 1)


def test_simple_family_anchor():
    payload = {
        "markers": ["m1", "m2"],
        "father": {"genotypes": ["AC", "TT"]},
        "mother": {"genotypes": ["GG", "AC"]},
        "children": [{"name": "k", "genotypes": ["AG", "TC"]}],
    }
    case = parse_case(payload)
    r = solve(case)
    assert r.min_crossovers == 0
    # father anchored at first het m1: hap0 = A
    assert r.sites[0].father_hap == ("A", "C")
    # child at m1: paternal A -> bit 0, maternal G -> either, fixed to hap0
    assert r.sites[0].transmissions[0][0] == 0
    # mother anchored at m2: hap0 = A... "A" < "C"
    assert r.sites[1].mother_hap == ("A", "C")
    # child at m2 is TC = paternal T, maternal C(bit1): mother bit must be 1
    assert r.sites[1].transmissions[0][1] == 1
    # with 0 crossovers the m1 maternal bit equals m2 maternal bit -> 1
    assert r.sites[0].transmissions[0][1] == 1


def test_determinism_and_anchors_large():
    rng = random.Random(42)
    n = 6
    m = 1500
    markers = [f"snp{i:04d}" for i in range(m)]
    # generate a consistent trio-ish pedigree from latent parental haplotypes
    f_h0 = [rng.choice("AC") for _ in range(m)]
    f_h1 = [rng.choice("AC") for _ in range(m)]
    m_h0 = [rng.choice("GT") for _ in range(m)]
    m_h1 = [rng.choice("GT") for _ in range(m)]
    father = [f_h0[i] + f_h1[i] for i in range(m)]
    mother = [m_h0[i] + m_h1[i] for i in range(m)]
    children = []
    for c in range(n):
        cp = rng.randbytes(m)
        cm = rng.randbytes(m)
        g = []
        for i in range(m):
            pa = f_h0[i] if cp[i] % 2 == 0 else f_h1[i]
            ma = m_h0[i] if cm[i] % 2 == 0 else m_h1[i]
            g.append(pa + ma)
        children.append({"name": f"ch{c}", "genotypes": g})
    payload = {
        "markers": markers,
        "father": {"genotypes": father},
        "mother": {"genotypes": mother},
        "children": children,
    }
    case = parse_case(payload)
    r1 = solve(case)
    r2 = solve(case)
    assert r1.min_crossovers == r2.min_crossovers
    for s1, s2 in zip(r1.sites, r2.sites):
        assert s1 == s2
    # first paternal het marker must be anchored with smaller allele on hap0
    for si, site in enumerate(case.sites):
        if _is_heterozygous(site.father):
            a, b = site.father
            assert r1.sites[si].father_hap == (min(a, b), max(a, b))
            break
