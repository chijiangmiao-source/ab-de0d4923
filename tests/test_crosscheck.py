"""Cross-validation of the solver against independent oracles.

* ``enumerate_optima`` brute-forces the *entire* optimum set on tiny cases
  and derives optimum, canonical key and fixed masks from first principles.
* ``direct_optimum`` recomputes the optimum with a plain O(16^K) transition
  DP, validating the hypercube distance transform on medium cases.
* ``check_response`` re-verifies any solution against the raw genotypes.
"""

import random

import pytest

from app.errors import MendelianInconsistency
from app.solver import solve
from tests.oracle import check_response, direct_optimum, enumerate_optima

CALLED = ["0/0", "0/1", "1/1"]
WITH_MISSING = ["0/0", "0/1", "1/1", "./.", "0/.", "1/."]


def random_case(rng, n_markers, n_children, alphabet):
    markers = [f"m{i}" for i in range(n_markers)]
    father = [rng.choice(alphabet) for _ in range(n_markers)]
    mother = [rng.choice(alphabet) for _ in range(n_markers)]
    children = [
        [rng.choice(alphabet) for _ in range(n_markers)] for _ in range(n_children)
    ]
    child_ids = [f"c{k}" for k in range(n_children)]
    return markers, father, mother, child_ids, children


def flatten_canonical(resp):
    key = []
    for entry in resp["solution"]:
        key += entry["father"]["alleles"]
        key += entry["mother"]["alleles"]
        for child in entry["children"]:
            key += [child["paternal"], child["maternal"]]
    return tuple(key)


def assert_matches_oracle(markers, father, mother, child_ids, children):
    oracle = enumerate_optima(markers, father, mother, children)
    try:
        resp = solve(markers, father, mother, child_ids, children)
    except MendelianInconsistency as exc:
        assert oracle.get("infeasible_marker") == exc.marker_index
        return
    assert "infeasible_marker" not in oracle, "solver feasible but oracle is not"
    assert resp["min_crossovers"] == oracle["cost"]
    assert flatten_canonical(resp) == oracle["key"]
    for i, entry in enumerate(resp["solution"]):
        assert entry["father"]["fixed"] == oracle["father_fixed"][i], f"marker {i}"
        assert entry["mother"]["fixed"] == oracle["mother_fixed"][i], f"marker {i}"
        for k, child in enumerate(entry["children"]):
            p_fixed, m_fixed = oracle["child_fixed"][i][k]
            assert child["paternal_fixed"] == p_fixed, f"marker {i} child {k}"
            assert child["maternal_fixed"] == m_fixed, f"marker {i} child {k}"
    check_response(markers, father, mother, child_ids, children, resp)


@pytest.mark.parametrize(
    "n_markers,n_children,alphabet,cases",
    [
        (2, 1, CALLED, 15),
        (2, 2, CALLED, 12),
        (2, 3, CALLED, 8),
        (3, 1, CALLED, 15),
        (3, 2, CALLED, 4),
        (4, 1, CALLED, 8),
        (2, 1, WITH_MISSING, 15),
        (2, 2, WITH_MISSING, 8),
        (3, 1, WITH_MISSING, 8),
    ],
)
def test_against_enumerative_oracle(n_markers, n_children, alphabet, cases):
    rng = random.Random(20260919 + n_markers * 100 + n_children)
    for _ in range(cases):
        case = random_case(rng, n_markers, n_children, alphabet)
        assert_matches_oracle(*case)


@pytest.mark.parametrize(
    "n_markers,n_children,missing_rate,cases",
    [
        (2, 1, 0.0, 12),
        (2, 2, 0.0, 10),
        (3, 1, 0.0, 10),
        (3, 2, 0.0, 4),
        (4, 1, 0.0, 6),
        (2, 3, 0.0, 4),
        (2, 1, 0.25, 12),
        (3, 1, 0.25, 8),
        (2, 2, 0.25, 8),
        (3, 2, 0.15, 3),
    ],
)
def test_feasible_against_enumerative_oracle(n_markers, n_children, missing_rate, cases):
    """Guaranteed-feasible simulated families: exercise canonical + masks."""
    rng = random.Random(31337 + n_markers * 10 + n_children)
    for _ in range(cases):
        markers, father, mother, child_ids, children, _ = simulate_family(
            rng, n_markers, n_children, switch_rate=0.15, missing_rate=missing_rate
        )
        assert_matches_oracle(markers, father, mother, child_ids, children)


@pytest.mark.parametrize("n_markers,n_children,cases", [(25, 3, 6), (12, 4, 4), (6, 5, 2)])
def test_optimum_against_direct_dp(n_markers, n_children, cases):
    rng = random.Random(777 + n_markers * 10 + n_children)
    for _ in range(cases):
        markers, father, mother, child_ids, children = random_case(
            rng, n_markers, n_children, WITH_MISSING
        )
        expected = direct_optimum(father, mother, children)
        try:
            resp = solve(markers, father, mother, child_ids, children)
        except MendelianInconsistency:
            assert expected is None
            continue
        assert expected is not None
        assert resp["min_crossovers"] == expected
        check_response(markers, father, mother, child_ids, children, resp)


def simulate_family(rng, n_markers, n_children, switch_rate=0.02, missing_rate=0.0):
    """Generate a guaranteed-feasible family by simulating transmission.

    Returns the solver inputs plus the number of switches in the true
    transmission (an upper bound on the optimum).
    """

    def random_haplotypes():
        return [[rng.randint(0, 1) for _ in range(n_markers)] for _ in range(2)]

    father_haps, mother_haps = random_haplotypes(), random_haplotypes()

    def call(a, b):
        if missing_rate and rng.random() < missing_rate:
            if rng.random() < 0.5:
                return "./."
            return f"{rng.choice([a, b])}/."
        return f"{a}/{a}" if a == b else "0/1"

    father = [call(father_haps[0][i], father_haps[1][i]) for i in range(n_markers)]
    mother = [call(mother_haps[0][i], mother_haps[1][i]) for i in range(n_markers)]
    children, true_switches = [], 0
    for _ in range(n_children):
        p, m = rng.randint(0, 1), rng.randint(0, 1)
        gts = []
        for i in range(n_markers):
            if i:
                if rng.random() < switch_rate:
                    p ^= 1
                    true_switches += 1
                if rng.random() < switch_rate:
                    m ^= 1
                    true_switches += 1
            gts.append(call(father_haps[p][i], mother_haps[m][i]))
        children.append(gts)
    markers = [f"m{i}" for i in range(n_markers)]
    child_ids = [f"c{k}" for k in range(n_children)]
    return markers, father, mother, child_ids, children, true_switches


def test_large_instance_smoke():
    """Worst-case size: 1500 markers x 6 children."""
    rng = random.Random(42)
    markers, father, mother, child_ids, children, bound = simulate_family(
        rng, 1500, 6, switch_rate=0.02
    )
    resp = solve(markers, father, mother, child_ids, children)
    assert resp["marker_count"] == 1500
    assert resp["child_count"] == 6
    assert resp["min_crossovers"] <= bound  # optimum never exceeds the truth
    check_response(markers, father, mother, child_ids, children, resp)
    # deterministic across repeated runs
    assert solve(markers, father, mother, child_ids, children) == resp


def test_large_instance_with_missing():
    rng = random.Random(1337)
    markers, father, mother, child_ids, children, bound = simulate_family(
        rng, 1500, 6, switch_rate=0.03, missing_rate=0.15
    )
    resp = solve(markers, father, mother, child_ids, children)
    assert resp["min_crossovers"] <= bound
    check_response(markers, father, mother, child_ids, children, resp)
