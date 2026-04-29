# projectmd/v2/ — refactored copies of the main analysis scripts

This folder holds **refactored copies** of selected scripts from `projectmd/`.
The originals are kept untouched; nothing here ever modifies files in
`projectmd/` or in the legacy `output/` tree.

## Goals of the refactor

1. Remove copy-pasted helper functions (`bh_adjust`, `np.fromstring`,
   GEO loading boilerplate, Seurat `setwd + source('workspace_paths.R')`
   prelude, trajectory scoring, ...) by routing them through
   `projectmd/v2/helpers/`.
2. Fix the deprecated `np.fromstring` calls (NumPy 2.0 removes this).
3. Standardise `workspace_paths` usage:
   * Python: import only `output_path`-style functions from
     `helpers.helpers_common`. No more `OUTPUT_DIR / "x"` patterns.
   * R: a one-liner `init_workspace()` instead of the 7-line
     `setwd + source` prelude.
4. Write every output under `output/v2/...`, mirroring the existing
   `output/figXX/` folder layout. This keeps v2 results side-by-side
   with the originals so they can be diffed before any switch-over.

## Scripts kept untouched (HPC / long-running)

These are **not** rewritten because they require either an HPC SLURM
allocation or many CPU-hours to re-run, and we should not change them
without an explicit go-ahead from the user:

* `fig02_infercnv_run_gte009.R`
* `fig02_infercnv_genelevel_panels.R`
* `fig06_pyscenic_run.py`
* `fig06_scenic_ctx_aucell_pipeline.sh`
* every script under `projectmd/hpc/`

## How dependencies are resolved

Each v2 script writes outputs to `output/v2/figXX/`.
When a v2 script depends on another script's output it uses
`dep_file("figXX", "filename.csv")`, which looks first in
`output/v2/figXX/` and falls back to the legacy `output/figXX/`.

This means each v2 script can run independently against the existing
v1 outputs **and** can be chained with its v2 siblings.

## Running

Run from the repository root, e.g.

```bash
cd /home/zlcmc/wslproject
conda run -n epn2 python projectmd/v2/fig04_gillen_scrna_validation.py
```

Outputs land in `output/v2/fig4_scrna/` instead of `output/fig4_scrna/`.
