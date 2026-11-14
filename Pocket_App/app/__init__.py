"""Pocket_App — pocket detection & 3D pharmacophore extraction service.

Phase 3 (SBDD) component of the AI Pharma platform. Exposes two HTTP
endpoints over the production discovery pipeline:

    POST /detect          — locate a druggable pocket on a protein
    POST /pharmacophore   — extract a 3D pharmacophore from a co-crystal ligand

Used by the backend ``production_pipeline`` to constrain REINVENT 4
generation and pre-filter candidates before docking.
"""
