"""Minimum-crossover family phasing via min-sum DP over transmission states.

No complete solution is ever enumerated.

State model
-----------
At marker i the joint word W_i packs two bits per child c:
  bit 2c     -> which paternal homolog (0/1) child c inherited
  bit 2c + 1 -> which maternal homolog (0/1) child c inherited
A switch of such a bit on an edge is one crossover in that meiosis, hence the
edge cost is ``popcount(W_{i-1} xor W_i)`` (K = 2 * #children <= 12 bits).

At every marker each parent's unphased genotype admits one or two ordered
haplotype pairs (at most 4 parental orientation combinations).  For a fixed
combination the children constrain their two-bit blocks independently, so
feasibility is a vectorised boolean predicate over all 2**K words; the word
cost does not depend on the orientation combination, which is needed only for
the lexicographic decision and the emitted haplotypes.

Symmetry handling
-----------------
Swapping the two homolog labels of one parent at *every* marker (and the
corresponding inherited bits of every child) leaves every cost invariant.
The canonical 0-direction is anchored at each parent's first truly
heterozygous marker (haplotype 0 carries the lexicographically smaller allele
there).  A parent with no heterozygous marker cannot be physically labelled;
that symmetry is intentionally left in the feasibility space so the
ambiguity masks report the unobservable transmission bits, while the greedy
canonical reconstruction still emits one unique lexicographic representative.

DP machinery
------------
``min_x cost[x] + popcount(x xor y)`` on the K-cube is a Hamming distance
transform in O(K * 2**K).  Forward/backward layers certify which local values
participate in a globally optimal path; the canonical optimum is recovered by
greedy lexicographic minimisation under that certificate, and the ambiguity
masks derive from every certified local value -- without enumerating any
complete solution.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import Case, SiteInput

INF = 1_000_000_000

_RANK = {None: 0, "A": 1, "C": 2, "G": 3, "T": 4}


@dataclass
class SiteAllowed:
    father_haps: tuple[tuple[str | None, str | None], ...]
    mother_haps: tuple[tuple[str | None, str | None], ...]
    combos: tuple[tuple[int, int], ...]          # (fi, mi) list
    combo_feasible: np.ndarray                   # bool[J, 2**K]
    allowed: np.ndarray                          # bool[2**K], OR over combos
    flip_fi: dict[int, int]
    flip_mi: dict[int, int]
    _pairs: dict[int, list[tuple[int, int]]] | None = None

    def pairs(self) -> dict[int, list[tuple[int, int]]]:
        """Lazy expansion (used only by the small-case brute-force tests)."""
        if self._pairs is None:
            out: dict[int, list[tuple[int, int]]] = {}
            for j, (fi, mi) in enumerate(self.combos):
                for w in np.flatnonzero(self.combo_feasible[j]):
                    out.setdefault(int(w), []).append((fi, mi))
            self._pairs = out
        return self._pairs


@dataclass
class SiteResult:
    marker: str
    father_hap: tuple[str | None, str | None]
    mother_hap: tuple[str | None, str | None]
    transmissions: tuple[tuple[int, int], ...]
    father_fixed: tuple[bool, bool]
    mother_fixed: tuple[bool, bool]
    child_fixed: tuple[tuple[bool, bool], ...]


@dataclass
class PhaseResult:
    min_crossovers: int
    sites: tuple[SiteResult, ...]
    children_names: tuple[str, ...]
    child_crossovers: tuple[int, ...]


class MendelError(Exception):
    def __init__(self, marker_index: int, marker: str):
        self.marker_index = marker_index
        self.marker = marker
        super().__init__(f"no Mendelian-consistent interpretation at marker {marker!r}")


# ---------------------------------------------------------------------------
# Per-site feasibility (fully vectorised)
# ---------------------------------------------------------------------------

def _hap_options(gt: tuple[str | None, str | None]) -> list[tuple[str | None, str | None]]:
    """Distinct ordered haplotype pairs for one unphased diploid genotype."""
    a, b = gt
    if a is not None and b is not None:
        return [(a, a)] if a == b else [(a, b), (b, a)]
    if a is None and b is None:
        return [(None, None)]
    x = b if a is None else a
    return [(x, None), (None, x)]


def _compatible(a: str | None, b: str | None) -> bool:
    return a is None or b is None or a == b


def _child_requirements(cgt: tuple[str | None, str | None]):
    if cgt[0] is None and cgt[1] is None:
        return ((None, None),)
    return tuple({(cgt[0], cgt[1]), (cgt[1], cgt[0])})


def _is_heterozygous(gt: tuple[str | None, str | None]) -> bool:
    return gt[0] is not None and gt[1] is not None and gt[0] != gt[1]


def _build_allowed(
    site: SiteInput, n_children: int, word_idx: np.ndarray
) -> SiteAllowed | None:
    f_opts = _hap_options(site.father)
    m_opts = _hap_options(site.mother)
    reqs = [_child_requirements(cgt) for cgt in site.children]
    combos = tuple(
        (fi, mi) for fi in range(len(f_opts)) for mi in range(len(m_opts))
    )
    words = word_idx.shape[0]
    combo_feasible = np.zeros((len(combos), words), dtype=bool)

    for j, (fi, mi) in enumerate(combos):
        fh = f_opts[fi]
        mh = m_opts[mi]
        feas = np.ones(words, dtype=bool)
        for c in range(n_children):
            allow = np.zeros(4, dtype=bool)  # indexed by pb + 2*mb
            for pb in (0, 1):
                for mb in (0, 1):
                    fa, ma = fh[pb], mh[mb]
                    for ca, cb in reqs[c]:
                        if _compatible(ca, fa) and _compatible(cb, ma):
                            allow[pb + 2 * mb] = True
                            break
            block = (word_idx >> (2 * c)) & 3
            feas &= allow[block]
        combo_feasible[j] = feas

    allowed = combo_feasible.any(axis=0)
    if not allowed.any():
        return None

    flip_fi = {0: 1, 1: 0} if len(f_opts) == 2 else {0: 0}
    flip_mi = {0: 1, 1: 0} if len(m_opts) == 2 else {0: 0}
    return SiteAllowed(
        tuple(f_opts), tuple(m_opts), combos, combo_feasible, allowed,
        flip_fi, flip_mi,
    )


# Hamming distance transform ------------------------------------------------

def _cube_index_sets(k: int, words: int):
    base = np.arange(words, dtype=np.int64)
    sets = []
    for b in range(k):
        s = 1 << b
        lo = base[(base & s) == 0]
        sets.append((lo, lo + s))
    return sets


def _hamming_transform(cost: np.ndarray, index_sets) -> np.ndarray:
    """out[y] = min_x (cost[x] + popcount(x xor y))."""
    out = cost.copy()
    for lo, hi in index_sets:
        a = out[lo]
        b = out[hi]
        out[lo] = np.minimum(a, b + 1)
        out[hi] = np.minimum(b, a + 1)
    return out


# Lexicographic keys --------------------------------------------------------

def _anchor_tuple(gt: tuple[str | None, str | None]) -> tuple[str, str]:
    a, b = gt
    return (a, b) if a < b else (b, a)  # type: ignore[return-value]


def _combo_base_key(fh, mh) -> int:
    """Key contribution of the four parental allele digits."""
    key = _RANK[fh[0]]
    key = key * 5 + _RANK[fh[1]]
    key = key * 5 + _RANK[mh[0]]
    key = key * 5 + _RANK[mh[1]]
    return key


def _word_child_keys(words: int, n: int) -> np.ndarray:
    """Per-word lexicographic value of all child transmission bits."""
    w = np.arange(words, dtype=np.int64)
    key = np.zeros(words, dtype=np.int64)
    for c in range(n):
        pb = (w >> (2 * c)) & 1
        mb = (w >> (2 * c + 1)) & 1
        key = key * 25 + pb * 5 + mb
    return key


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------

def solve(case: Case) -> PhaseResult:
    m = len(case.sites)
    n = len(case.children_names)
    k = 2 * n
    words = 1 << k
    word_idx = np.arange(words, dtype=np.int64)
    index_sets = _cube_index_sets(k, words)
    child_keys = _word_child_keys(words, n)
    child_scale = 25 ** n  # base key sits above all child digits

    allowed_sites: list[SiteAllowed] = []
    for si, site in enumerate(case.sites):
        a = _build_allowed(site, n, word_idx)
        if a is None:
            raise MendelError(si, site.marker)
        allowed_sites.append(a)

    # -- anchors: first truly heterozygous marker per parent ---------------
    f_anchor_site = next(
        (si for si, s in enumerate(case.sites) if _is_heterozygous(s.father)), None
    )
    m_anchor_site = next(
        (si for si, s in enumerate(case.sites) if _is_heterozygous(s.mother)), None
    )
    f_anchor = (
        _anchor_tuple(case.sites[f_anchor_site].father) if f_anchor_site is not None else None
    )
    m_anchor = (
        _anchor_tuple(case.sites[m_anchor_site].mother) if m_anchor_site is not None else None
    )

    # Active orientation combos per site after imposing the real anchors.
    # Parents without any heterozygous marker keep full label-flip symmetry
    # here on purpose: it must remain visible in the ambiguity masks; the
    # greedy reconstruction still picks one unique canonical representative.
    active: list[np.ndarray] = []
    allowed_masks: list[np.ndarray] = []
    for si, a in enumerate(allowed_sites):
        act = np.ones(len(a.combos), dtype=bool)
        for j, (fi, mi) in enumerate(a.combos):
            if f_anchor_site == si and a.father_haps[fi] != f_anchor:
                act[j] = False
            if m_anchor_site == si and a.mother_haps[mi] != m_anchor:
                act[j] = False
        active.append(act)
        mask = a.combo_feasible[act].any(axis=0)
        if not mask.any():  # defensive; anchor symmetry guarantees feasibility
            raise MendelError(si, case.sites[si].marker)
        allowed_masks.append(mask)

    # ---- forward ----
    fwd: list[np.ndarray] = []
    cur = np.full(words, INF, dtype=np.int32)
    cur[allowed_masks[0]] = 0
    fwd.append(cur)
    for si in range(1, m):
        moved = _hamming_transform(fwd[-1], index_sets)
        fwd.append(np.where(allowed_masks[si], moved, INF).astype(np.int32))

    best = int(fwd[-1].min())
    if best >= INF:
        raise MendelError(m - 1, case.sites[-1].marker)

    # ---- backward ----
    bwd: list[np.ndarray] = [np.empty(words, dtype=np.int32) for _ in range(m)]
    last = np.full(words, INF, dtype=np.int32)
    last[allowed_masks[-1]] = 0
    bwd[-1] = last
    for si in range(m - 2, -1, -1):
        moved = _hamming_transform(bwd[si + 1], index_sets)
        bwd[si] = np.where(allowed_masks[si], moved, INF).astype(np.int32)

    cert = [
        (fwd[si].astype(np.int64) + bwd[si].astype(np.int64)) == best
        for si in range(m)
    ]

    # ---- greedy canonical reconstruction ----
    chosen: list[tuple[int, int, int]] = []
    prefix_cost = 0
    prev_word: int | None = None
    for si in range(m):
        a = allowed_sites[si]
        if prev_word is None:
            reachable = cert[si]
        else:
            edge = _popcount_xor(word_idx, prev_word)
            reachable = cert[si] & (
                prefix_cost + edge + bwd[si].astype(np.int64) == best
            )
        best_key = None
        best_pick: tuple[int, int, int] | None = None
        for j, (fi, mi) in enumerate(a.combos):
            if not active[si][j]:
                continue
            cand = reachable & a.combo_feasible[j]
            if not cand.any():
                continue
            base = _combo_base_key(a.father_haps[fi], a.mother_haps[mi]) * child_scale
            keys = base + child_keys
            keys = np.where(cand, keys, np.iinfo(np.int64).max)
            w = int(np.argmin(keys))
            key = int(keys[w])
            if best_key is None or key < best_key:
                best_key = key
                best_pick = (w, fi, mi)
        if best_pick is None:  # defensive; mathematically unreachable
            raise MendelError(si, case.sites[si].marker)
        chosen.append(best_pick)
        w = best_pick[0]
        prefix_cost += 0 if prev_word is None else (w ^ prev_word).bit_count()
        prev_word = w

    # ---- emit canonical solution + fixedness masks ----
    results: list[SiteResult] = []
    for si, site in enumerate(case.sites):
        a = allowed_sites[si]
        w, fi, mi = chosen[si]

        f_idx = set()
        m_idx = set()
        for j, (fj, mj) in enumerate(a.combos):
            if active[si][j] and (cert[si] & a.combo_feasible[j]).any():
                f_idx.add(fj)
                m_idx.add(mj)
        father_fixed = (
            len({a.father_haps[j][0] for j in f_idx}) == 1,
            len({a.father_haps[j][1] for j in f_idx}) == 1,
        )
        mother_fixed = (
            len({a.mother_haps[j][0] for j in m_idx}) == 1,
            len({a.mother_haps[j][1] for j in m_idx}) == 1,
        )

        child_fixed: list[tuple[bool, bool]] = []
        trans: list[tuple[int, int]] = []
        for c in range(n):
            trans.append(((w >> (2 * c)) & 1, (w >> (2 * c + 1)) & 1))
            has0p = bool((cert[si] & (((word_idx >> (2 * c)) & 1) == 0)).any())
            has1p = bool((cert[si] & (((word_idx >> (2 * c)) & 1) == 1)).any())
            has0m = bool((cert[si] & (((word_idx >> (2 * c + 1)) & 1) == 0)).any())
            has1m = bool((cert[si] & (((word_idx >> (2 * c + 1)) & 1) == 1)).any())
            child_fixed.append((not (has0p and has1p), not (has0m and has1m)))

        results.append(
            SiteResult(
                marker=site.marker,
                father_hap=a.father_haps[fi],
                mother_hap=a.mother_haps[mi],
                transmissions=tuple(trans),
                father_fixed=father_fixed,
                mother_fixed=mother_fixed,
                child_fixed=tuple(child_fixed),
            )
        )

    chosen_words = np.array([c[0] for c in chosen], dtype=np.int64)
    child_crossovers: list[int] = []
    for c in range(n):
        pb = (chosen_words >> (2 * c)) & 1
        mb = (chosen_words >> (2 * c + 1)) & 1
        child_crossovers.append(int(np.count_nonzero(np.diff(pb)) + np.count_nonzero(np.diff(mb))))

    return PhaseResult(
        best,
        tuple(results),
        case.children_names,
        tuple(child_crossovers),
    )


def _popcount_xor(word_idx: np.ndarray, prev: int) -> np.ndarray:
    x = word_idx ^ prev
    pc = np.zeros_like(x)
    while x.any():
        pc += x & 1
        x >>= 1
    return pc
