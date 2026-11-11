# Data and model availability — preliminary policy

**This candidate does not redistribute raw DUD-E files, model weights, PDB structures, proprietary checkpoints, patient information, or original source data by default.** It is a safety-conscious code candidate, not a complete independent data replication archive.

- EGFR DUD-E data source: target `egfr`, files `actives_final.ism` and `decoys_final.ism`; local Phase 10 matched panel 50 active + 50 decoy; confirm the license and obtain from an authorized upstream archive before redistribution.
- ALK5 DUD-E source: target `tgfr1`, 133 clustered active entries and 8,500 decoy entries in previously downloaded files. Known observed SHA-256 digests: active `6a1df16bfa5477c932057e574146c2d58c5854768525feb900aa8042c0930829`; decoy `503596358d8d72ed29f22595c8d1c5ab702768f4a30e9d1502a43c47130d8440`. These identify observed files, not ownership or distribution rights.
- ALK5 matched frozen panel SHA-256 `42461c5b2677e5b4a2df8ac6ba6c81bfaec2b319f261709c88dd591a84f36d64`.
- Crystal structures: EGFR 1M17; ALK5 1VJY and 3HMM. Provide PDB identifiers, ligand IDs, chain selections, preparation and grid settings. Vina binding scores are not experimental binding affinities.
- AIPHARMA generated score tables and evidence reports may be publishable **after** checking rights and review. Full reproduction also needs panel identities/SMILES, native binaries and receptor preparation artifacts.
- Provide exact ADMET-AI model version and weights/checkpoint provenance, and REINVENT checkpoint availability separately.

**Publication gap:** a legally redistributable, curated data+receptor bundle or documented exact acquisition and preprocessing recipe is still required. A source-only ZIP cannot solve dataset reproducibility on its own.
