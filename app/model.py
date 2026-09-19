"""Family phasing domain model and validation (422 errors)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


class PayloadError(Exception):
    """Structurally invalid request. ``path`` is a JSON-pointer-ish location."""

    def __init__(self, message: str, path: str):
        super().__init__(message)
        self.message = message
        self.path = path


_MISSING = {".", "-", ""}
_ALLELES = {"A", "C", "G", "T"}
_MISSING_TOKEN = "."

MAX_MARKERS = 1500
MAX_CHILDREN = 6


@dataclass(frozen=True)
class SiteInput:
    """One marker column."""

    marker: str
    father: tuple[str | None, str | None]
    mother: tuple[str | None, str | None]
    children: tuple[tuple[str | None, str | None], ...]


@dataclass(frozen=True)
class Case:
    sites: tuple[SiteInput, ...]
    children_names: tuple[str, ...]
    marker_names: tuple[str, ...]


def _err(message: str, path: str) -> PayloadError:
    return PayloadError(message, path)


def _parse_genotype(value: Any, path: str, who: str) -> tuple[str | None, str | None]:
    if isinstance(value, str):
        token = value
    elif isinstance(value, list) and len(value) == 2 and all(isinstance(x, str) for x in value):
        token = value[0] + value[1]
    else:
        raise _err(f"{who} genotype must be a string (e.g. 'AT', '..') or a 2-char array", path)

    if len(token) != 2:
        raise _err(f"{who} genotype must have exactly 2 alleles, got {token!r}", path)

    out: list[str | None] = []
    for ch in token:
        up = ch.upper()
        if up in _MISSING:
            out.append(None)
        elif up in _ALLELES:
            out.append(up)
        else:
            raise _err(f"illegal allele {ch!r}; expected one of A/C/G/T or missing '.'", path)
    return (out[0], out[1])


def parse_case(payload: Any) -> Case:
    if not isinstance(payload, dict):
        raise _err("request body must be a JSON object", "$")

    markers = payload.get("markers")
    if not isinstance(markers, list):
        raise _err("'markers' must be a list", "$.markers")
    if not 2 <= len(markers) <= MAX_MARKERS:
        raise _err(
            f"between 2 and {MAX_MARKERS} ordered markers are required, "
            f"got {len(markers)}",
            "$.markers",
        )

    father = payload.get("father")
    if not isinstance(father, dict):
        raise _err("'father' must be an object", "$.father")
    mother = payload.get("mother")
    if not isinstance(mother, dict):
        raise _err("'mother' must be an object", "$.mother")

    children = payload.get("children")
    if not isinstance(children, list):
        raise _err("'children' must be a list", "$.children")
    if not 1 <= len(children) <= MAX_CHILDREN:
        raise _err(
            f"between 1 and {MAX_CHILDREN} children are required, got {len(children)}",
            "$.children",
        )

    children_names: list[str] = []
    for ci, child in enumerate(children):
        cpath = f"$.children[{ci}]"
        if not isinstance(child, dict):
            raise _err("child must be an object", cpath)
        name = child.get("name", f"child_{ci + 1}")
        if not isinstance(name, str) or not name:
            raise _err("child 'name' must be a non-empty string", f"{cpath}.name")
        if name in children_names:
            raise _err(f"duplicate child name {name!r}", f"{cpath}.name")
        children_names.append(name)
        if "genotypes" not in child:
            raise _err("child requires 'genotypes'", cpath)

    father_g = father.get("genotypes")
    if not isinstance(father_g, list):
        raise _err("'genotypes' must be a list", "$.father.genotypes")
    mother_g = mother.get("genotypes")
    if not isinstance(mother_g, list):
        raise _err("'genotypes' must be a list", "$.mother.genotypes")

    seen_markers: set[str] = set()
    marker_names: list[str] = []
    sites: list[SiteInput] = []

    for mi, marker in enumerate(markers):
        mpath = f"$.markers[{mi}]"
        if isinstance(marker, str):
            mname = marker
        elif isinstance(marker, dict) and isinstance(marker.get("name"), str):
            mname = marker["name"]
        else:
            raise _err("marker must be a string or an object with string 'name'", mpath)
        if mname in seen_markers:
            raise _err(f"duplicate marker name {mname!r}", mpath)
        seen_markers.add(mname)
        marker_names.append(mname)

        if mi >= len(father_g):
            raise _err(f"missing father genotype at marker {mname!r}", f"$.father.genotypes[{mi}]")
        if mi >= len(mother_g):
            raise _err(f"missing mother genotype at marker {mname!r}", f"$.mother.genotypes[{mi}]")
        fa = _parse_genotype(father_g[mi], f"$.father.genotypes[{mi}]", "father")
        mo = _parse_genotype(mother_g[mi], f"$.mother.genotypes[{mi}]", "mother")

        child_gts: list[tuple[str | None, str | None]] = []
        for ci, child in enumerate(children):
            glist = child["genotypes"]
            cpath = f"$.children[{ci}].genotypes"
            if not isinstance(glist, list):
                raise _err("'genotypes' must be a list", cpath)
            if mi >= len(glist):
                raise _err(
                    f"missing {children_names[ci]!r} genotype at marker {mname!r}",
                    f"{cpath}[{mi}]",
                )
            child_gts.append(_parse_genotype(glist[mi], f"{cpath}[{mi}]", "child"))
        sites.append(SiteInput(marker=mname, father=fa, mother=mo, children=tuple(child_gts)))

    extra_father = len(father_g) - len(markers)
    extra_mother = len(mother_g) - len(markers)
    if extra_father:
        raise _err(
            f"father has {extra_father} extra genotype(s); length must equal number of markers",
            "$.father.genotypes",
        )
    if extra_mother:
        raise _err(
            f"mother has {extra_mother} extra genotype(s); length must equal number of markers",
            "$.mother.genotypes",
        )
    for ci, child in enumerate(children):
        extra = len(child["genotypes"]) - len(markers)
        if extra:
            raise _err(
                f"child {children_names[ci]!r} has {extra} extra genotype(s); "
                "length must equal number of markers",
                f"$.children[{ci}].genotypes",
            )

    return Case(tuple(sites), tuple(children_names), tuple(marker_names))


MISSING_TOKEN = _MISSING_TOKEN
