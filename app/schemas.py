"""Request and response schemas for the phasing API."""

from typing import Annotated, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, StringConstraints

# Biallelic, possibly (partially) missing genotype call: "0/0", "0/1", "1/1",
# "./.", "0/.", "./1", ...  "|" is accepted as a separator alias.
GenotypeStr = Annotated[str, StringConstraints(pattern=r"^[01.][/|][01.]$")]

MarkerId = Annotated[str, StringConstraints(min_length=1, max_length=128)]
ChildId = Annotated[str, StringConstraints(min_length=1, max_length=128)]

MIN_MARKERS = 2
MAX_MARKERS = 1500
MIN_CHILDREN = 1
MAX_CHILDREN = 6


class ChildInput(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: Optional[ChildId] = None
    genotypes: list[GenotypeStr] = Field(
        min_length=MIN_MARKERS, max_length=MAX_MARKERS
    )


class PhaseRequest(BaseModel):
    """One phasing job: ordered markers plus unphased family genotypes."""

    model_config = ConfigDict(extra="forbid")

    markers: list[MarkerId] = Field(
        min_length=MIN_MARKERS,
        max_length=MAX_MARKERS,
        description="Marker ids in map order (2..1500).",
    )
    father: list[GenotypeStr] = Field(
        min_length=MIN_MARKERS, max_length=MAX_MARKERS
    )
    mother: list[GenotypeStr] = Field(
        min_length=MIN_MARKERS, max_length=MAX_MARKERS
    )
    children: list[ChildInput] = Field(min_length=MIN_CHILDREN, max_length=MAX_CHILDREN)


class ParentPhase(BaseModel):
    """Ordered haplotype alleles [haplotype0, haplotype1] at one marker."""

    alleles: list[int]
    fixed: list[bool]


class ChildPhase(BaseModel):
    index: int
    id: str
    paternal: int
    maternal: int
    paternal_fixed: bool
    maternal_fixed: bool


class MarkerPhase(BaseModel):
    index: int
    id: str
    father: ParentPhase
    mother: ParentPhase
    children: list[ChildPhase]


class PhaseResponse(BaseModel):
    status: Literal["ok"] = "ok"
    marker_count: int
    child_count: int
    min_crossovers: int
    all_fixed: bool
    solution: list[MarkerPhase]
