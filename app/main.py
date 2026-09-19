"""HTTP API for the pedigree phasing service (pure backend, no frontend)."""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from app import __version__
from app.errors import MendelianInconsistency, StructureError
from app.schemas import PhaseRequest, PhaseResponse
from app.solver import solve

app = FastAPI(
    title="Pedigree Phasing Service",
    version=__version__,
    description=(
        "Minimum-recombination phasing for nuclear families. "
        "Deterministic, fully offline, no frontend."
    ),
)


@app.exception_handler(RequestValidationError)
async def request_validation_handler(request: Request, exc: RequestValidationError):
    """Re-emit Pydantic errors as a stable, JSON-safe, locatable 422 body."""
    detail = [
        {
            "loc": list(err.get("loc", ())),
            "msg": err.get("msg", ""),
            "type": err.get("type", ""),
        }
        for err in exc.errors()
    ]
    return JSONResponse(status_code=422, content={"detail": detail})


@app.exception_handler(StructureError)
async def structure_handler(request: Request, exc: StructureError):
    return JSONResponse(status_code=422, content={"detail": exc.errors})


@app.exception_handler(MendelianInconsistency)
async def mendelian_handler(request: Request, exc: MendelianInconsistency):
    return JSONResponse(
        status_code=409,
        content={
            "detail": {
                "error": "mendelian_inconsistency",
                "message": (
                    "no Mendelian-consistent explanation exists; earliest "
                    f"unsolvable marker: index {exc.marker_index} "
                    f"(id {exc.marker_id!r})"
                ),
                "marker_index": exc.marker_index,
                "marker_id": exc.marker_id,
                "genotypes": exc.genotypes,
            }
        },
    )


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/api/v1/phase", response_model=PhaseResponse)
def phase(request: PhaseRequest):
    _validate_lengths(request)
    child_ids = [
        child.id if child.id is not None else f"child{k + 1}"
        for k, child in enumerate(request.children)
    ]
    return solve(
        marker_ids=request.markers,
        father=request.father,
        mother=request.mother,
        child_ids=child_ids,
        children=[child.genotypes for child in request.children],
    )


def _validate_lengths(request: PhaseRequest) -> None:
    """Cross-field check: every genotype list must match the marker count."""
    expected = len(request.markers)
    errors = []
    for field, values in (("father", request.father), ("mother", request.mother)):
        if len(values) != expected:
            errors.append(
                {
                    "loc": ["body", field],
                    "msg": (
                        f"expected {expected} genotypes (one per marker), "
                        f"got {len(values)}"
                    ),
                    "type": "length_mismatch",
                }
            )
    for k, child in enumerate(request.children):
        if len(child.genotypes) != expected:
            errors.append(
                {
                    "loc": ["body", "children", k, "genotypes"],
                    "msg": (
                        f"expected {expected} genotypes (one per marker), "
                        f"got {len(child.genotypes)}"
                    ),
                    "type": "length_mismatch",
                }
            )
    if errors:
        raise StructureError(errors)
