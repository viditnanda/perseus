# Releasing Perseus

Everything is automated by `.github/workflows/release.yml` once two one-time steps are done.

## One-time setup

1. **PyPI.** Create an account at https://pypi.org (or log in). Under *Your account → Publishing → Add a new
   pending publisher*, enter: PyPI project name `perseus-tda`, owner `viditnanda`, repository `perseus`,
   workflow `release.yml`, environment `pypi`. No password or token is stored anywhere.
2. **Zenodo (DOI).** Log in at https://zenodo.org with GitHub, open *GitHub* in the account menu, and switch
   on `viditnanda/perseus`. Every GitHub release then gets a DOI; `.zenodo.json` and `CITATION.cff` supply
   the metadata.

## Each release

1. Update the version in `python/pyproject.toml`, `python/perseus/__init__.py` and `CITATION.cff`.
2. Tag and push: `git tag v6.0.0 && git push origin v6.0.0`.
3. The workflow builds the Linux, macOS (Apple silicon and Intel) and Windows executables and attaches them
   to a GitHub release, then builds the Python wheels (Linux, macOS) and source package and uploads them to
   PyPI. The web page links to `releases/latest`, so it needs no change.

Windows Python wheels are not built: the code relies on 128-bit integers, which Microsoft's compiler lacks.
Windows users can use the executable, or the Python package under WSL.

## Rerunning a release

If the release workflow fails after the tag exists, fix the workflow on `main`, then rerun it for the same
tag without moving the tag (which Zenodo has already archived): push a branch `release/<tag>`, e.g.
`git push origin main:release/v6.0.0`, or use "Run workflow" on the release workflow in the Actions tab.
The code built is always the tag's; only the workflow comes from the branch.
