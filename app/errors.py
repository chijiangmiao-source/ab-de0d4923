"""Domain exceptions shared by the solver and the HTTP layer."""


class StructureError(Exception):
    """Cross-field structural validation failure (reported as HTTP 422).

    ``errors`` is a list of ``{"loc": [...], "msg": str, "type": str}`` items
    following the shape of FastAPI/Pydantic validation errors so clients can
    locate every offending field.
    """

    def __init__(self, errors):
        self.errors = list(errors)
        super().__init__(str(self.errors))


class MendelianInconsistency(Exception):
    """Raised when a marker admits no Mendelian-consistent assignment.

    Carries the 0-based index and the id of the *earliest* marker for which
    no choice of parental phasing and child transmissions can explain the
    observed genotypes (reported as HTTP 409), plus the observed genotypes
    at that marker for context.
    """

    def __init__(self, marker_index: int, marker_id: str, genotypes=None):
        self.marker_index = marker_index
        self.marker_id = marker_id
        self.genotypes = genotypes
        super().__init__(
            f"no Mendelian-consistent assignment for marker "
            f"{marker_index} (id {marker_id!r})"
        )
