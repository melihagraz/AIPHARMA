# Reproducibility and source release

Start with [`reproducibility/README.md`](reproducibility/README.md), then
[`INSTALL_RELEASE_TR.md`](INSTALL_RELEASE_TR.md) for WSL commands.

**Release status:** v0.1.0-rc1 is a **local source candidate**, not a GitHub release,
not a software DOI and not an approved open-source distribution.

The `scripts/aipharma_release.py` tool inventories current local experiments,
records precise installed Python and native Vina versions, executes EAVE tests,
and builds a fail-closed source ZIP with a file hash manifest. It does **not**
reproduce model results from a clean public machine or distribute DUD-E data.

See [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md),
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md), and
[`LICENSE_PENDING.md`](LICENSE_PENDING.md).
