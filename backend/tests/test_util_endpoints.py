"""Tests for /api/util/mol-image + /api/util/mol-pdb.

These endpoints are RDKit-only — no upstream services involved — so the
tests run hermetically without any respx mock. The point is to lock in:
    * happy path returns the right content-type + non-trivial payload
    * bad SMILES doesn't 500 the table (mol-image returns 200 blank PNG,
      mol-pdb returns 404 — different conventions because the table can
      tolerate a placeholder PNG, while the NGL viewer must know to fall
      back to receptor-only rendering)
    * size + width clamping rejects abusive inputs without server load
"""

from __future__ import annotations

import pytest


def test_mol_image_valid_smiles(client):
    """A 24-char drug-like SMILES should render as a non-trivial PNG."""
    r = client.get("/api/util/mol-image", params={"smiles": "CCO", "w": 120, "h": 80})
    assert r.status_code == 200
    assert r.headers["content-type"] == "image/png"
    # Header bytes should be the PNG magic number — protects against
    # an accidental JSON / text response slipping past content-type.
    assert r.content.startswith(b"\x89PNG")
    # Sanity: ethanol PNG > 200 B; a 1x1 sentinel is ~68 B.
    assert len(r.content) > 200


def test_mol_image_invalid_smiles_returns_blank_png(client):
    """Garbage SMILES must NOT raise 500 — return the 1x1 transparent PNG
    so a 200-row table doesn't break on one bad entry."""
    r = client.get("/api/util/mol-image", params={"smiles": "$%^not_a_smiles"})
    assert r.status_code == 200
    assert r.headers["content-type"] == "image/png"
    assert r.content.startswith(b"\x89PNG")
    # The sentinel blank is ~68 B; a real ethanol render is 1.5+ KB.
    assert len(r.content) < 200


@pytest.mark.parametrize("w,h", [(23, 80), (601, 80), (80, 23), (80, 601)])
def test_mol_image_rejects_out_of_range_dimensions(client, w, h):
    """w/h are clamped to [24, 600] in the FastAPI Query() validator —
    out-of-range values must 422 instead of being silently coerced."""
    r = client.get("/api/util/mol-image", params={"smiles": "CCO", "w": w, "h": h})
    assert r.status_code == 422


def test_mol_image_cache_header_set(client):
    """We rely on the browser cache for hot scroll perf — make sure the
    long-lived Cache-Control header doesn't regress."""
    r = client.get("/api/util/mol-image", params={"smiles": "CCO"})
    assert "max-age" in r.headers.get("cache-control", "")


def test_mol_pdb_valid_smiles(client):
    """3D embedding for a common drug-like SMILES (celecoxib subgraph)."""
    smi = "Cc1ccc(cc1)c2cc(nn2c3ccc(cc3)S(N)(=O)=O)C(F)(F)F"
    r = client.get("/api/util/mol-pdb", params={"smiles": smi})
    assert r.status_code == 200
    assert "x-pdb" in r.headers["content-type"]
    body = r.text
    # PDB block sanity: starts with a record type, has 3D coords, ends with END.
    assert body.startswith(("HETATM", "ATOM", "REMARK", "COMPND", "MODEL"))
    assert "END" in body
    # 3D embedding produces ≥3 atoms even for tiny molecules.
    assert body.count("HETATM") + body.count("ATOM") >= 3


def test_mol_pdb_invalid_smiles_returns_404(client):
    """Unlike mol-image which serves a blank sentinel, mol-pdb returns
    404 so the NGL viewer knows to switch to receptor-only rendering."""
    r = client.get("/api/util/mol-pdb", params={"smiles": "$%^bad"})
    assert r.status_code == 404
    assert "Could not generate" in r.json()["detail"]


def test_mol_pdb_empty_smiles_rejected(client):
    """Empty SMILES hits the FastAPI Query(min_length=1) validator."""
    r = client.get("/api/util/mol-pdb", params={"smiles": ""})
    assert r.status_code == 422


def test_mol_image_caching_same_smiles_returns_identical_bytes(client):
    """Per-process lru_cache makes the second call return the same PNG
    bytes as the first. Side-effect: regression guard against accidental
    introduction of nondeterminism in the RDKit pipeline (random coords,
    timestamps, etc.)."""
    r1 = client.get("/api/util/mol-image", params={"smiles": "CCO"})
    r2 = client.get("/api/util/mol-image", params={"smiles": "CCO"})
    assert r1.status_code == r2.status_code == 200
    assert r1.content == r2.content
