"""HTTP layer: validation (422), phasing endpoint, health check."""
from __future__ import annotations

import json

from fastapi import FastAPI, Request, Response
from fastapi.responses import JSONResponse

from .model import PayloadError, parse_case
from .phasing import MendelError, solve

app = FastAPI(title="Family Phasing Service", version="1.0.0")


def _canonical_json(payload) -> bytes:
    # Deterministic serialization: identical inputs must produce byte-identical
    # responses (sorted keys, fixed separators, no host-dependent whitespace).
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")


class CanonicalResponse(Response):
    media_type = "application/json"

    def render(self, content) -> bytes:
        return _canonical_json(content)


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.exception_handler(PayloadError)
async def payload_error_handler(_request: Request, exc: PayloadError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "error": {
                "type": "invalid_structure",
                "message": exc.message,
                "path": exc.path,
            }
        },
    )


@app.post("/phase")
async def phase(request: Request) -> Response:
    try:
        payload = await request.json()
    except Exception:
        return JSONResponse(
            status_code=422,
            content={
                "error": {
                    "type": "invalid_structure",
                    "message": "request body must be valid JSON",
                    "path": "$",
                }
            },
        )

    case = parse_case(payload)
    try:
        result = solve(case)
    except MendelError as exc:
        return JSONResponse(
            status_code=422,
            content={
                "error": {
                    "type": "mendel_inconsistent",
                    "message": (
                        "no Mendelian-consistent interpretation exists at the "
                        f"earliest unsolvable marker {exc.marker!r}"
                    ),
                    "marker_index": exc.marker_index,
                    "marker": exc.marker,
                    "path": f"$.markers[{exc.marker_index}]",
                }
            },
        )

    markers_out = []
    for site in result.sites:
        markers_out.append(
            {
                "marker": site.marker,
                "father": {
                    "haplotype": [
                        "." if a is None else a for a in site.father_hap
                    ],
                    "ambiguous": [not site.father_fixed[0], not site.father_fixed[1]],
                },
                "mother": {
                    "haplotype": [
                        "." if a is None else a for a in site.mother_hap
                    ],
                    "ambiguous": [not site.mother_fixed[0], not site.mother_fixed[1]],
                },
                "children": [
                    {
                        "name": result.children_names[c],
                        "transmitted_haplotype": [pb, mb],
                        "ambiguous": [
                            not site.child_fixed[c][0],
                            not site.child_fixed[c][1],
                        ],
                    }
                    for c, (pb, mb) in enumerate(site.transmissions)
                ],
            }
        )

    body = {
        "min_crossovers": result.min_crossovers,
        "child_crossovers": {
            result.children_names[c]: result.child_crossovers[c]
            for c in range(len(result.children_names))
        },
        "solution": markers_out,
    }
    return CanonicalResponse(content=body)
