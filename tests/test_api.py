"""End-to-end API tests: 200, 422 locators, Mendel error, byte determinism."""
from __future__ import annotations

import json

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def post(payload):
    return client.post("/phase", content=json.dumps(payload),
                       headers={"content-type": "application/json"})


def test_health():
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json() == {"status": "ok"}


def test_basic_ok_shape():
    payload = {
        "markers": ["m1", "m2"],
        "father": {"genotypes": ["AC", "TT"]},
        "mother": {"genotypes": ["GG", "AC"]},
        "children": [{"name": "k", "genotypes": ["AG", "TC"]}],
    }
    r = post(payload)
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["min_crossovers"] == 0
    assert len(body["solution"]) == 2
    assert set(body["child_crossovers"]) == {"k"}


def test_byte_identical_repeat():
    payload = {
        "markers": ["b", "a"],
        "father": {"genotypes": ["AC", "AC"]},
        "mother": {"genotypes": ["GT", "GT"]},
        "children": [
            {"name": "x", "genotypes": ["AG", "CT"]},
            {"name": "y", "genotypes": [".T", "AG"]},
        ],
    }
    r1 = post(payload)
    r2 = post(payload)
    assert r1.status_code == 200
    assert r1.content == r2.content


def test_422_too_few_children():
    payload = {
        "markers": ["m1", "m2"],
        "father": {"genotypes": ["AC", "AC"]},
        "mother": {"genotypes": ["GT", "GT"]},
        "children": [],
    }
    r = post(payload)
    assert r.status_code == 422
    err = r.json()["error"]
    assert err["path"] == "$.children"


def test_422_bad_allele_location():
    payload = {
        "markers": ["m1", "m2"],
        "father": {"genotypes": ["AC", "ZZ"]},
        "mother": {"genotypes": ["GG", "AC"]},
        "children": [{"name": "k", "genotypes": ["AG", "TC"]}],
    }
    r = post(payload)
    assert r.status_code == 422
    err = r.json()["error"]
    assert err["path"] == "$.father.genotypes[1]"


def test_422_length_mismatch():
    payload = {
        "markers": ["m1", "m2"],
        "father": {"genotypes": ["AC"]},
        "mother": {"genotypes": ["GG", "AC"]},
        "children": [{"name": "k", "genotypes": ["AG", "TC"]}],
    }
    r = post(payload)
    assert r.status_code == 422
    assert "$.father.genotypes[1]" in r.json()["error"]["path"]


def test_422_too_few_markers():
    payload = {
        "markers": ["only"],
        "father": {"genotypes": ["AC"]},
        "mother": {"genotypes": ["GT"]},
        "children": [{"name": "k", "genotypes": ["AG"]}],
    }
    r = post(payload)
    assert r.status_code == 422
    assert r.json()["error"]["path"] == "$.markers"


def test_422_too_many_markers():
    payload = {
        "markers": [f"m{i}" for i in range(1501)],
        "father": {"genotypes": ["AC"] * 1501},
        "mother": {"genotypes": ["GT"] * 1501},
        "children": [{"name": "k", "genotypes": ["AG"] * 1501}],
    }
    r = post(payload)
    assert r.status_code == 422
    assert r.json()["error"]["path"] == "$.markers"


def test_422_bad_json():
    r = client.post("/phase", content="{not json",
                    headers={"content-type": "application/json"})
    assert r.status_code == 422


def test_mendel_inconsistent_earliest_marker():
    # father AA x mother CC cannot produce a G child at m2
    payload = {
        "markers": ["m1", "m2"],
        "father": {"genotypes": ["AC", "AA"]},
        "mother": {"genotypes": ["GG", "CC"]},
        "children": [{"name": "k", "genotypes": ["AG", "GG"]}],
    }
    r = post(payload)
    assert r.status_code == 422
    err = r.json()["error"]
    assert err["type"] == "mendel_inconsistent"
    assert err["marker_index"] == 1
    assert err["marker"] == "m2"


def test_missing_genotypes_still_solvable():
    payload = {
        "markers": ["m1", "m2"],
        "father": {"genotypes": ["..", "AC"]},
        "mother": {"genotypes": ["..", "GT"]},
        "children": [{"name": "k", "genotypes": ["..", "AG"]}],
    }
    r = post(payload)
    assert r.status_code == 200, r.text
    assert r.json()["min_crossovers"] == 0
