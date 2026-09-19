"""HTTP API tests: happy path, locatable 422s, 409 Mendelian, determinism."""

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

EXAMPLE = {
    "markers": ["rs01", "rs02", "rs03"],
    "father": ["0/1", "0/1", "0/0"],
    "mother": ["0/0", "0/1", "0/0"],
    "children": [{"id": "proband", "genotypes": ["0/0", "0/1", "0/0"]}],
}


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json() == {"status": "ok"}


def test_happy_path_full_body():
    resp = client.post("/api/v1/phase", json=EXAMPLE)
    assert resp.status_code == 200
    assert resp.json() == {
        "status": "ok",
        "marker_count": 3,
        "child_count": 1,
        "min_crossovers": 0,
        "all_fixed": False,
        "solution": [
            {
                "index": 0,
                "id": "rs01",
                "father": {"alleles": [0, 1], "fixed": [True, True]},
                "mother": {"alleles": [0, 0], "fixed": [True, True]},
                "children": [
                    {
                        "index": 0,
                        "id": "proband",
                        "paternal": 0,
                        "maternal": 0,
                        "paternal_fixed": True,
                        "maternal_fixed": False,
                    }
                ],
            },
            {
                "index": 1,
                "id": "rs02",
                "father": {"alleles": [1, 0], "fixed": [False, False]},
                "mother": {"alleles": [0, 1], "fixed": [True, True]},
                "children": [
                    {
                        "index": 0,
                        "id": "proband",
                        "paternal": 0,
                        "maternal": 0,
                        "paternal_fixed": True,
                        "maternal_fixed": False,
                    }
                ],
            },
            {
                "index": 2,
                "id": "rs03",
                "father": {"alleles": [0, 0], "fixed": [True, True]},
                "mother": {"alleles": [0, 0], "fixed": [True, True]},
                "children": [
                    {
                        "index": 0,
                        "id": "proband",
                        "paternal": 0,
                        "maternal": 0,
                        "paternal_fixed": True,
                        "maternal_fixed": False,
                    }
                ],
            },
        ],
    }


def test_repeated_calls_byte_identical():
    other = {
        "markers": ["a", "b"],
        "father": ["0/1", "0/1"],
        "mother": ["0/1", "0/0"],
        "children": [{"genotypes": ["0/1", "0/0"]}],
    }
    first = client.post("/api/v1/phase", json=EXAMPLE)
    client.post("/api/v1/phase", json=other)  # interleave a different request
    second = client.post("/api/v1/phase", json=EXAMPLE)
    assert first.status_code == second.status_code == 200
    assert first.content == second.content


def test_default_child_ids():
    payload = {
        "markers": ["a", "b"],
        "father": ["0/1", "0/0"],
        "mother": ["0/0", "0/1"],
        "children": [{"genotypes": ["0/0", "0/0"]}, {"genotypes": ["0/1", "0/0"]}],
    }
    resp = client.post("/api/v1/phase", json=payload)
    assert resp.status_code == 200
    ids = [c["id"] for c in resp.json()["solution"][0]["children"]]
    assert ids == ["child1", "child2"]


def test_pipe_separator_and_half_missing_accepted():
    payload = {
        "markers": ["a", "b"],
        "father": ["0|1", "0/."],
        "mother": ["./.", "0/0"],
        "children": [{"genotypes": ["0/1", "0/."]}],
    }
    assert client.post("/api/v1/phase", json=payload).status_code == 200


class TestUnprocessableStructure:
    def test_bad_genotype_located(self):
        payload = {
            "markers": ["a", "b", "c"],
            "father": ["0/1", "0/2", "0/0"],
            "mother": ["0/0", "0/0", "0/0"],
            "children": [{"genotypes": ["0/0", "0/0", "0/0"]}],
        }
        resp = client.post("/api/v1/phase", json=payload)
        assert resp.status_code == 422
        detail = resp.json()["detail"]
        assert any(item["loc"] == ["body", "father", 1] for item in detail)

    def test_length_mismatch_located(self):
        payload = {
            "markers": ["a", "b"],
            "father": ["0/1", "0/0", "0/0"],
            "mother": ["0/0", "0/0"],
            "children": [{"genotypes": ["0/0", "0/0", "0/1"]}],
        }
        resp = client.post("/api/v1/phase", json=payload)
        assert resp.status_code == 422
        locs = [item["loc"] for item in resp.json()["detail"]]
        assert ["body", "father"] in locs
        assert ["body", "children", 0, "genotypes"] in locs
        assert ["body", "mother"] not in locs

    def test_too_few_markers(self):
        payload = {
            "markers": ["a"],
            "father": ["0/1"],
            "mother": ["0/0"],
            "children": [{"genotypes": ["0/0"]}],
        }
        resp = client.post("/api/v1/phase", json=payload)
        assert resp.status_code == 422
        assert any(item["loc"][:2] == ["body", "markers"] for item in resp.json()["detail"])

    def test_child_count_bounds(self):
        base = {
            "markers": ["a", "b"],
            "father": ["0/1", "0/0"],
            "mother": ["0/0", "0/1"],
        }
        resp = client.post("/api/v1/phase", json={**base, "children": []})
        assert resp.status_code == 422
        seven = [{"genotypes": ["0/0", "0/0"]} for _ in range(7)]
        resp = client.post("/api/v1/phase", json={**base, "children": seven})
        assert resp.status_code == 422

    def test_extra_field_rejected(self):
        payload = {**EXAMPLE, "phased": True}
        resp = client.post("/api/v1/phase", json=payload)
        assert resp.status_code == 422

    def test_empty_marker_id_rejected(self):
        payload = {**EXAMPLE, "markers": ["rs01", "", "rs03"]}
        resp = client.post("/api/v1/phase", json=payload)
        assert resp.status_code == 422
        assert any(
            item["loc"] == ["body", "markers", 1] for item in resp.json()["detail"]
        )


class TestMendelianConflict:
    def test_earliest_unsolvable_marker(self):
        payload = {
            "markers": ["m1", "m2", "m3"],
            "father": ["0/0", "0/0", "0/0"],
            "mother": ["0/0", "0/0", "0/0"],
            "children": [{"genotypes": ["0/0", "1/1", "0/0"]}],
        }
        resp = client.post("/api/v1/phase", json=payload)
        assert resp.status_code == 409
        detail = resp.json()["detail"]
        assert detail["error"] == "mendelian_inconsistency"
        assert detail["marker_index"] == 1
        assert detail["marker_id"] == "m2"
        assert detail["genotypes"]["father"] == "0/0"
        assert detail["genotypes"]["children"] == ["1/1"]

    def test_conflict_response_is_deterministic(self):
        payload = {
            "markers": ["m1", "m2"],
            "father": ["0/0", "0/0"],
            "mother": ["0/0", "0/0"],
            "children": [{"genotypes": ["1/1", "0/0"]}],
        }
        first = client.post("/api/v1/phase", json=payload)
        second = client.post("/api/v1/phase", json=payload)
        assert first.status_code == second.status_code == 409
        assert first.content == second.content
