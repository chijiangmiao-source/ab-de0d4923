"""Minimum-recombination pedigree phasing solver.

Model
-----
Markers are processed in the submitted order.  At every marker the solver
chooses

* an ordered allele pair for each parent -- the two haplotype alleles, and
* for every child a pair of transmission bits: which paternal and which
  maternal haplotype was inherited.

A marker assignment is *feasible* when it is consistent with every observed
(possibly partially missing) genotype.  The objective is the total number of
transmission-chain switches between adjacent markers, summed over all
children and both parental sides.

Because the parents' per-marker phasing carries no cost, the problem
decomposes into a chain DP over the joint child-transmission state (two bits
per child, i.e. at most 4^6 = 4096 states).  Transition cost between two
states is their Hamming distance, so the min-plus product is evaluated with
an O(2K * 4^K) hypercube distance transform instead of an O(16^K) sweep.

Nothing here enumerates the optimal-solution set:

* the optimum comes from the forward/backward DP,
* the canonical solution (lexicographically smallest optimal solution in
  marker -> father -> mother -> children bit order) is extracted greedily,
* the fixed/ambiguous flags are derived from the sets of states and parent
  phases that lie on *some* optimal path, again via the forward/backward
  tables.

Determinism: every choice is resolved by explicit minima over totally
ordered domains, so identical inputs always produce identical outputs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from app.errors import MendelianInconsistency

# Comfortably larger than any attainable cost (<= 2 * 6 * 1499 switches),
# small enough that INF + accumulated relaxations never overflows int32.
INF = 100_000_000

# Ordered-pair set of a definitely heterozygous call ("0/1" or "1/0").
HET = frozenset({(0, 1), (1, 0)})


def parse_genotype(gt: str) -> frozenset:
    """Map a validated genotype string to its possible ordered allele pairs.

    Accepted forms (validated upstream): ``0/0``, ``0/1``, ``1/1``, ``./.``
    and half-missing calls such as ``0/.`` or ``./1``.  ``|`` is accepted as
    a separator alias.  Alleles are biallelic (0/1); ``.`` is a missing allele.
    """
    sep = "/" if "/" in gt else "|"
    left, right = gt.split(sep)
    called = [c for c in (left, right) if c != "."]
    if len(called) == 2:
        if left == right:
            allele = int(left)
            return frozenset({(allele, allele)})
        return HET
    if len(called) == 1:
        allele = int(called[0])
        other = 1 - allele
        return frozenset({(allele, allele), (allele, other), (other, allele)})
    return frozenset({(0, 0), (0, 1), (1, 0), (1, 1)})


def first_het_index(genotypes) -> "int | None":
    """Index of the first definitely heterozygous marker, or ``None``.

    Only calls with both alleles observed and different (``0/1``/``1/0``)
    count; partially missing calls (``0/.``) are not definitely heterozygous.
    """
    for i, gt in enumerate(genotypes):
        if parse_genotype(gt) == HET:
            return i
    return None


@dataclass(frozen=True)
class MarkerPlan:
    """Precomputed per-marker choice structure.

    fpairs / mpairs: ordered allele pairs the father / mother may take at
        this marker (after applying the orientation fix), sorted.
    combos: all ``(father_pair, mother_pair)`` combinations, father-major.
    combo_index: ``(father_pair, mother_pair) -> position in combos``.
    trans_masks: for every combo, for every child, a 4-bit mask of the
        compatible transmissions (bit ``t`` = transmission ``2*p + m``).
    """

    fpairs: tuple
    mpairs: tuple
    combos: tuple
    combo_index: dict
    trans_masks: tuple


def build_marker_plan(father_gt, mother_gt, child_gts, fix_father, fix_mother):
    """Enumerate the local choices at one marker (never the global solution)."""
    fpairs = sorted(parse_genotype(father_gt))
    mpairs = sorted(parse_genotype(mother_gt))
    if fix_father:  # first heterozygous marker: haplotype 0 carries allele 0
        fpairs = [pair for pair in fpairs if pair[0] == 0]
    if fix_mother:
        mpairs = [pair for pair in mpairs if pair[0] == 0]
    child_options = [parse_genotype(gt) for gt in child_gts]
    combos = tuple((fp, mp) for fp in fpairs for mp in mpairs)
    trans_masks = []
    for fp, mp in combos:
        per_child = []
        for options in child_options:
            mask = 0
            for trans in range(4):
                paternal, maternal = trans >> 1, trans & 1
                if (fp[paternal], mp[maternal]) in options:
                    mask |= 1 << trans
            per_child.append(mask)
        trans_masks.append(tuple(per_child))
    return MarkerPlan(
        fpairs=tuple(fpairs),
        mpairs=tuple(mpairs),
        combos=combos,
        combo_index={combo: ci for ci, combo in enumerate(combos)},
        trans_masks=tuple(trans_masks),
    )


class Solver:
    """Chain DP over joint child-transmission states for one request."""

    def __init__(self, n_children: int):
        self.n_children = n_children
        self.n_bits = 2 * n_children
        self.n_states = 1 << self.n_bits
        idx = np.arange(self.n_states, dtype=np.int64)
        self.state_index = idx
        # Big-endian layout: child k occupies bits 2*(K-1-k)+1 (paternal) and
        # 2*(K-1-k) (maternal), so numeric state order equals the canonical
        # lexicographic order of [p1, m1, p2, m2, ...].
        self.shifts = [2 * (n_children - 1 - k) for k in range(n_children)]
        self.trans = [(idx >> sh) & 3 for sh in self.shifts]
        self.partner = [idx ^ (1 << bit) for bit in range(self.n_bits)]
        popcount = np.zeros(self.n_states, dtype=np.int64)
        for value in range(1, self.n_states):
            popcount[value] = popcount[value >> 1] + (value & 1)
        self.popcount = popcount
        self._plan_cache = {}

    # -- per-marker feasibility structure ---------------------------------

    def plan_for(self, father_gt, mother_gt, child_gts, fix_father, fix_mother):
        key = (father_gt, mother_gt, child_gts, fix_father, fix_mother)
        hit = self._plan_cache.get(key)
        if hit is None:
            plan = build_marker_plan(
                father_gt, mother_gt, child_gts, fix_father, fix_mother
            )
            hit = (plan, self._compatibility(plan))
            self._plan_cache[key] = hit
        return hit

    def _compatibility(self, plan: MarkerPlan) -> np.ndarray:
        """Bitmask of compatible parent-phase combos for every state."""
        comp = np.zeros(self.n_states, dtype=np.uint16)
        for ci, per_child in enumerate(plan.trans_masks):
            ok = np.ones(self.n_states, dtype=bool)
            for k in range(self.n_children):
                ok &= (np.right_shift(per_child[k], self.trans[k]) & 1).astype(bool)
            comp |= ok.astype(np.uint16) << ci
        return comp

    # -- min-plus Hamming distance transform -------------------------------

    def relax(self, x: np.ndarray) -> np.ndarray:
        """``out[s'] = min_s x[s] + Hamming(s, s')`` on the 2K-bit hypercube.

        Processing one bit at a time relaxes exactly the paths whose differing
        bits form a subset of the bits handled so far, so after all bits the
        result is the exact distance transform.  Runs in O(2K * 4^K).
        """
        for partner in self.partner:
            x = np.minimum(x, x[partner] + 1)
        return x


def _combo_bit(comp: np.ndarray, ci: int) -> np.ndarray:
    return ((comp >> ci) & 1).astype(bool)


def solve(marker_ids, father, mother, child_ids, children):
    """Run the full analysis and return a JSON-serialisable result dict.

    Raises :class:`MendelianInconsistency` with the earliest marker that has
    no Mendelian-consistent local assignment.  (Markers are independent with
    respect to feasibility: parental phasing and child transmissions are free
    per marker, so a globally consistent explanation exists iff every marker
    is locally feasible.)
    """
    n_markers = len(marker_ids)
    n_children = len(children)
    solver = Solver(n_children)
    per_marker_children = tuple(zip(*children))
    father_first = first_het_index(father)
    mother_first = first_het_index(mother)

    plans, comps = [], []
    for i in range(n_markers):
        plan, comp = solver.plan_for(
            father[i],
            mother[i],
            per_marker_children[i],
            father_first == i,
            mother_first == i,
        )
        if not comp.any():
            context = {
                "father": father[i],
                "mother": mother[i],
                "children": [children[k][i] for k in range(n_children)],
            }
            raise MendelianInconsistency(i, marker_ids[i], context)
        plans.append(plan)
        comps.append(comp)

    # Forward pass: forward[i][s] = min cost of markers 0..i ending in state s.
    forward = [np.where(comps[0] != 0, 0, INF).astype(np.int32)]
    for i in range(1, n_markers):
        relaxed = solver.relax(forward[-1])
        forward.append(np.where(comps[i] != 0, relaxed, INF).astype(np.int32))
    optimum = int(forward[-1].min())

    # Backward pass: backward[i][s] = min cost of markers i+1..end from state s.
    backward = [None] * n_markers
    backward[-1] = np.zeros(solver.n_states, dtype=np.int32)
    for i in range(n_markers - 2, -1, -1):
        nxt = np.where(comps[i + 1] != 0, backward[i + 1], INF)
        backward[i] = solver.relax(nxt).astype(np.int32)

    solution = []
    all_fixed = True
    prefix_cost = 0
    prev_state = -1
    for i in range(n_markers):
        plan, comp = plans[i], comps[i]
        n_combos = len(plan.combos)

        # States lying on at least one globally optimal path.
        on_optimum = (forward[i].astype(np.int64) + backward[i]) == optimum

        # -- ambiguity mask: values realised by some optimal solution ------
        used = [
            ci
            for ci in range(n_combos)
            if np.any(on_optimum & _combo_bit(comp, ci))
        ]
        father_fixed = [
            len({plan.combos[ci][0][j] for ci in used}) == 1 for j in (0, 1)
        ]
        mother_fixed = [
            len({plan.combos[ci][1][j] for ci in used}) == 1 for j in (0, 1)
        ]
        child_fixed = []
        for k in range(n_children):
            sh = solver.shifts[k]
            p_bits = ((solver.state_index >> (sh + 1)) & 1).astype(bool)
            m_bits = ((solver.state_index >> sh) & 1).astype(bool)
            p_fixed = not (np.any(on_optimum & p_bits) and np.any(on_optimum & ~p_bits))
            m_fixed = not (np.any(on_optimum & m_bits) and np.any(on_optimum & ~m_bits))
            child_fixed.append((p_fixed, m_fixed))

        # -- greedy canonical choice: smallest (fp, mp, state) such that the
        # chosen prefix still extends to a globally optimal solution --------
        if i == 0:
            reach = np.zeros(solver.n_states, dtype=np.int64)
        else:
            reach = prefix_cost + solver.popcount[solver.state_index ^ prev_state]
        extendable = (comp != 0) & ((reach + backward[i].astype(np.int64)) == optimum)
        usable = [
            ci
            for ci in range(n_combos)
            if np.any(extendable & _combo_bit(comp, ci))
        ]
        fp = min(plan.combos[ci][0] for ci in usable)
        mp = min(plan.combos[ci][1] for ci in usable if plan.combos[ci][0] == fp)
        star = plan.combo_index[(fp, mp)]
        state = int(np.argmax(extendable & _combo_bit(comp, star)))
        if i > 0:
            prefix_cost += int(solver.popcount[state ^ prev_state])
        prev_state = state

        children_out = []
        for k in range(n_children):
            sh = solver.shifts[k]
            p_fixed, m_fixed = child_fixed[k]
            children_out.append(
                {
                    "index": k,
                    "id": child_ids[k],
                    "paternal": (state >> (sh + 1)) & 1,
                    "maternal": (state >> sh) & 1,
                    "paternal_fixed": p_fixed,
                    "maternal_fixed": m_fixed,
                }
            )
        all_fixed = (
            all_fixed
            and all(father_fixed)
            and all(mother_fixed)
            and all(flag for pair in child_fixed for flag in pair)
        )
        solution.append(
            {
                "index": i,
                "id": marker_ids[i],
                "father": {"alleles": [fp[0], fp[1]], "fixed": father_fixed},
                "mother": {"alleles": [mp[0], mp[1]], "fixed": mother_fixed},
                "children": children_out,
            }
        )

    return {
        "status": "ok",
        "marker_count": n_markers,
        "child_count": n_children,
        "min_crossovers": optimum,
        "all_fixed": all_fixed,
        "solution": solution,
    }
