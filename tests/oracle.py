"""Independent verification helpers for the test-suite.

These helpers re-implement the problem *definition* directly (enumeration,
quadratic DP, plain solution checking) without any of the solver's clever
machinery (distance transform, greedy canonical extraction).  They are only
used by tests; the service itself never enumerates optimal solutions.
"""

from itertools import product

from app.solver import INF, first_het_index, parse_genotype


def enumerate_optima(marker_ids, father, mother, children):
    """Exhaustively enumerate all feasible assignments (tiny cases only).

    Returns ``{"infeasible_marker": i}`` if marker ``i`` is the earliest
    locally infeasible marker, otherwise a dict with the optimum cost, the
    canonical (lexicographically smallest) key, the canonical assignment and
    per-item fixed flags derived from the projection of *all* optima.
    """
    n_markers = len(father)
    n_children = len(children)
    father_first = first_het_index(father)
    mother_first = first_het_index(mother)

    per_marker = []
    for i in range(n_markers):
        fpairs = sorted(parse_genotype(father[i]))
        mpairs = sorted(parse_genotype(mother[i]))
        if father_first == i:
            fpairs = [p for p in fpairs if p[0] == 0]
        if mother_first == i:
            mpairs = [p for p in mpairs if p[0] == 0]
        child_opts = [parse_genotype(children[k][i]) for k in range(n_children)]
        entries = []
        for fp in fpairs:
            for mp in mpairs:
                allowed = [
                    [
                        (p, m)
                        for p in (0, 1)
                        for m in (0, 1)
                        if (fp[p], mp[m]) in child_opts[k]
                    ]
                    for k in range(n_children)
                ]
                if any(len(options) == 0 for options in allowed):
                    continue
                for choice in product(*allowed):
                    bits = tuple(bit for pair in choice for bit in pair)
                    entries.append((fp, mp, bits))
        if not entries:
            return {"infeasible_marker": i}
        per_marker.append(entries)

    best = None
    best_key = None
    best_choice = None
    projections = None

    def fresh_projections():
        return [
            {"fp": set(), "mp": set(), "bits": [set() for _ in range(2 * n_children)]}
            for _ in range(n_markers)
        ]

    def absorb(choice):
        for i, (fp, mp, bits) in enumerate(choice):
            proj = projections[i]
            proj["fp"].add(fp)
            proj["mp"].add(mp)
            for j, bit in enumerate(bits):
                proj["bits"][j].add(bit)

    def rec(i, cost, prev_bits, key, choice):
        nonlocal best, best_key, best_choice, projections
        if best is not None and cost > best:
            return
        if i == n_markers:
            if best is None or cost < best:
                best = cost
                best_key = key
                best_choice = choice
                projections = fresh_projections()
                absorb(choice)
            elif cost == best:
                if key < best_key:
                    best_key = key
                    best_choice = choice
                absorb(choice)
            return
        for fp, mp, bits in per_marker[i]:
            step = 0 if i == 0 else sum(a != b for a, b in zip(prev_bits, bits))
            rec(i + 1, cost + step, bits, key + (fp + mp + bits), choice + [(fp, mp, bits)])

    rec(0, 0, None, (), [])

    father_fixed = [
        [len({fp[j] for fp in projections[i]["fp"]}) == 1 for j in (0, 1)]
        for i in range(n_markers)
    ]
    mother_fixed = [
        [len({mp[j] for mp in projections[i]["mp"]}) == 1 for j in (0, 1)]
        for i in range(n_markers)
    ]
    child_fixed = [
        [
            tuple(len(projections[i]["bits"][2 * k + side]) == 1 for side in (0, 1))
            for k in range(n_children)
        ]
        for i in range(n_markers)
    ]
    return {
        "cost": best,
        "key": best_key,
        "choice": best_choice,
        "father_fixed": father_fixed,
        "mother_fixed": mother_fixed,
        "child_fixed": child_fixed,
    }


def direct_optimum(father, mother, children):
    """Optimum cost via a plain O(4^K * 4^K) transition DP (no transform)."""
    n_markers = len(father)
    n_children = len(children)
    n_states = 1 << (2 * n_children)
    father_first = first_het_index(father)
    mother_first = first_het_index(mother)

    def bits(state, k):
        sh = 2 * (n_children - 1 - k)
        return (state >> (sh + 1)) & 1, (state >> sh) & 1

    feasible = []
    for i in range(n_markers):
        fpairs = sorted(parse_genotype(father[i]))
        mpairs = sorted(parse_genotype(mother[i]))
        if father_first == i:
            fpairs = [p for p in fpairs if p[0] == 0]
        if mother_first == i:
            mpairs = [p for p in mpairs if p[0] == 0]
        child_opts = [parse_genotype(children[k][i]) for k in range(n_children)]
        row = [False] * n_states
        for state in range(n_states):
            for fp in fpairs:
                for mp in mpairs:
                    if all(
                        (fp[bits(state, k)[0]], mp[bits(state, k)[1]]) in child_opts[k]
                        for k in range(n_children)
                    ):
                        row[state] = True
                        break
                    # continue to next mp/fp otherwise
                if row[state]:
                    break
        feasible.append(row)
        if not any(row):
            return None  # infeasible; solver raises instead

    costs = [0 if feasible[0][s] else INF for s in range(n_states)]
    for i in range(1, n_markers):
        nxt = []
        for t in range(n_states):
            best = min(
                costs[s] + bin(s ^ t).count("1") for s in range(n_states)
            )
            nxt.append(best if feasible[i][t] else INF)
        costs = nxt
    return min(costs)


def check_response(marker_ids, father, mother, child_ids, children, resp):
    """Validate a service response against the raw problem definition.

    Recomputes the crossover count of the reported canonical solution and
    verifies Mendelian consistency, the orientation rule and the mask shape.
    Raises AssertionError on any violation.
    """
    n_markers = len(marker_ids)
    n_children = len(children)
    assert resp["status"] == "ok"
    assert resp["marker_count"] == n_markers
    assert resp["child_count"] == n_children
    assert len(resp["solution"]) == n_markers

    father_first = first_het_index(father)
    mother_first = first_het_index(mother)
    switches = 0
    prev = None
    for i, entry in enumerate(resp["solution"]):
        assert entry["index"] == i
        assert entry["id"] == marker_ids[i]
        fp = tuple(entry["father"]["alleles"])
        mp = tuple(entry["mother"]["alleles"])
        assert len(fp) == 2 and len(mp) == 2
        assert len(entry["father"]["fixed"]) == 2
        assert len(entry["mother"]["fixed"]) == 2
        # parental haplotypes consistent with the observed genotypes
        assert fp in parse_genotype(father[i]), (i, "father", fp)
        assert mp in parse_genotype(mother[i]), (i, "mother", mp)
        # orientation rule at the first heterozygous marker
        if father_first == i:
            assert fp[0] == 0, (i, "father orientation", fp)
        if mother_first == i:
            assert mp[0] == 0, (i, "mother orientation", mp)
        assert len(entry["children"]) == n_children
        state = []
        for k, child_out in enumerate(entry["children"]):
            assert child_out["index"] == k
            assert child_out["id"] == child_ids[k]
            p, m = child_out["paternal"], child_out["maternal"]
            assert p in (0, 1) and m in (0, 1)
            inherited = (fp[p], mp[m])
            assert inherited in parse_genotype(children[k][i]), (i, k, inherited)
            state.extend((p, m))
        if prev is not None:
            switches += sum(a != b for a, b in zip(prev, state))
        prev = state
    assert switches == resp["min_crossovers"], (switches, resp["min_crossovers"])
    return True
