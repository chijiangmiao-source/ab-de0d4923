"""Unit tests for the core solver primitives and hand-verifiable cases."""

import pytest

from app.errors import MendelianInconsistency
from app.solver import first_het_index, parse_genotype, solve


class TestParseGenotype:
    def test_homozygous(self):
        assert parse_genotype("0/0") == {(0, 0)}
        assert parse_genotype("1/1") == {(1, 1)}

    def test_heterozygous_is_unordered(self):
        assert parse_genotype("0/1") == {(0, 1), (1, 0)}
        assert parse_genotype("1/0") == {(0, 1), (1, 0)}
        assert parse_genotype("0|1") == {(0, 1), (1, 0)}

    def test_fully_missing(self):
        assert parse_genotype("./.") == {(0, 0), (0, 1), (1, 0), (1, 1)}

    def test_half_missing(self):
        assert parse_genotype("0/.") == {(0, 0), (0, 1), (1, 0)}
        assert parse_genotype("./1") == {(1, 1), (0, 1), (1, 0)}


class TestFirstHet:
    def test_basic(self):
        assert first_het_index(["0/0", "0/1", "1/1"]) == 1
        assert first_het_index(["1/0"]) == 0

    def test_missing_is_not_het(self):
        assert first_het_index(["./.", "0/.", "0/1"]) == 2

    def test_none(self):
        assert first_het_index(["0/0", "./.", "1/1"]) is None


class TestHandComputedExample:
    """Fully hand-verified three-marker family (see README worked example)."""

    def setup_method(self):
        self.result = solve(
            ["rs01", "rs02", "rs03"],
            ["0/1", "0/1", "0/0"],
            ["0/0", "0/1", "0/0"],
            ["proband"],
            [["0/0", "0/1", "0/0"]],
        )

    def test_optimum(self):
        assert self.result["min_crossovers"] == 0

    def test_canonical_solution(self):
        sol = self.result["solution"]
        assert [m["father"]["alleles"] for m in sol] == [[0, 1], [1, 0], [0, 0]]
        assert [m["mother"]["alleles"] for m in sol] == [[0, 0], [0, 1], [0, 0]]
        for entry in sol:
            (child,) = entry["children"]
            assert (child["paternal"], child["maternal"]) == (0, 0)

    def test_fixed_mask(self):
        sol = self.result["solution"]
        assert sol[0]["father"]["fixed"] == [True, True]
        assert sol[0]["children"][0]["paternal_fixed"] is True
        assert sol[0]["children"][0]["maternal_fixed"] is False
        # two optimal phasings of the father coexist at marker 2
        assert sol[1]["father"]["fixed"] == [False, False]
        assert sol[1]["mother"]["fixed"] == [True, True]
        assert sol[2]["children"][0]["maternal_fixed"] is False
        assert self.result["all_fixed"] is False


class TestOrientation:
    def test_first_het_forces_zero_on_haplotype_0(self):
        result = solve(
            ["m1", "m2"],
            ["1/0", "0/0"],  # first het marker: orientation fixed to (0, 1)
            ["0/0", "0/0"],
            ["c1"],
            [["0/0", "0/0"]],
        )
        assert result["solution"][0]["father"]["alleles"] == [0, 1]
        assert result["solution"][0]["father"]["fixed"] == [True, True]

    def test_parent_without_het_site(self):
        result = solve(
            ["m1", "m2"],
            ["0/0", "0/0"],
            ["0/1", "0/0"],
            ["c1"],
            [["0/0", "0/0"]],
        )
        assert result["solution"][0]["father"]["alleles"] == [0, 0]
        assert result["solution"][0]["father"]["fixed"] == [True, True]


class TestAllMissing:
    def test_everything_free(self):
        result = solve(
            ["m1", "m2", "m3"],
            ["./.", "./.", "./."],
            ["./.", "./.", "./."],
            ["c1", "c2"],
            [["./.", "./.", "./."], ["./.", "./.", "./."]],
        )
        assert result["min_crossovers"] == 0
        for entry in result["solution"]:
            assert entry["father"]["alleles"] == [0, 0]  # canonical minimum
            assert entry["father"]["fixed"] == [False, False]
            assert entry["mother"]["fixed"] == [False, False]
            for child in entry["children"]:
                assert (child["paternal"], child["maternal"]) == (0, 0)
                assert child["paternal_fixed"] is False
                assert child["maternal_fixed"] is False
        assert result["all_fixed"] is False


class TestMendelianInconsistency:
    def test_earliest_marker_reported(self):
        with pytest.raises(MendelianInconsistency) as excinfo:
            solve(
                ["a", "b", "c"],
                ["0/0", "0/0", "0/0"],
                ["0/0", "0/0", "0/0"],
                ["c1"],
                [["0/0", "1/1", "0/0"]],  # marker 1 impossible
            )
        assert excinfo.value.marker_index == 1
        assert excinfo.value.marker_id == "b"

    def test_first_marker_reported(self):
        with pytest.raises(MendelianInconsistency) as excinfo:
            solve(["a", "b"], ["0/0", "0/0"], ["0/0", "0/0"], ["c1"], [["1/1", "0/0"]])
        assert excinfo.value.marker_index == 0

    def test_later_markers_do_not_mask_earliest(self):
        # both markers 0 and 2 are inconsistent; earliest must win
        with pytest.raises(MendelianInconsistency) as excinfo:
            solve(
                ["a", "b", "c"],
                ["0/0", "0/0", "0/0"],
                ["0/0", "0/0", "0/0"],
                ["c1"],
                [["1/1", "0/0", "1/1"]],
            )
        assert excinfo.value.marker_index == 0


class TestDeterminism:
    def test_repeated_solve_identical(self):
        args = (
            ["m1", "m2", "m3", "m4"],
            ["0/1", "./.", "0/1", "0/0"],
            ["0/0", "0/1", "0/.", "0/1"],
            ["c1", "c2"],
            [["0/1", "0/0", "0/1", "0/0"], ["0/0", "0/1", "./.", "0/1"]],
        )
        assert solve(*args) == solve(*args)
