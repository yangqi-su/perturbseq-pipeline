# perturbseq-pipeline: fix list

Code review of `src/perturbseq_pipeline/`, started 2026-10-03, consolidated 2026-10-07, status updated 2026-10-22.

**Line numbers.** Sections 1-12 were written against `a08071c` (upstream `1c48f9d` plus the knockdown filter), before the reformat; those line numbers are stale and must be looked up again when an item is started. Sections 13-14, and entries marked "(reformatted)", refer to `feat/msigdb-gene-sets` (`cb7679d`). Section 15 refers to `feat/guide-calling-fishash-pgmm` (`3d4376d`).

**Status (2026-10-22).** Upstream `main` is `195fb8d`.
- Merged (ours): #7 knockdown filter (`3164a97`); #8 reformat (`78e8510`); #9 MSigDB gene sets + gseapy enrichment (`238f3ef`); #16 knockdown non-carrier baseline (`195fb8d`, 8.7).
- Merged (others): #10 four P0s, including 2.1 and 11.2 (`c2e3501`); #11 nine P1s, including 11.17 (`bef4abb`); #12-#15 high-MOI mode (membership calling, membership-aware statistics, per-cell scores on membership, joint membership regression; `dd9a34f`..`9e26f8b`). Entries marked fixed by #10/#11 are taken from the PR descriptions, not re-checked here. Items in sections 2-14 were not re-checked against the high-MOI code.
- Open: #17 guide calling, fishash + Poisson-Gaussian mixture (`feat/guide-calling-fishash-pgmm`, `3d4376d`; section 15).
- Worktrees: `data_perturb_v2/worktrees/pp-guide-calling` (#17), `data_perturb_v2/worktrees/pp-kd-non-carriers` (#16, merged, can be removed).
- Package size: 42.0k lines originally → 21.7k on `feat/msigdb-gene-sets`; not re-measured after #10-#17. Remaining size work: sections 13-14 (~4-5k more).

How to read the evidence column:

- **verified**: confirmed by reading the code.
- **reported**: found by a code-reading pass but not re-checked by hand.
- **estimate**: reasoned from the code. Nothing here was run, timed or memory-profiled.

---

## 1. Size and duplication

| # | Problem | Evidence | Fix |
|---|---|---|---|
| 1.1 | 42k lines, of which only ~42% carried real code; the rest was formatting (one argument per line, bracket-only lines, a blank line after nearly every statement). | verified | **Done (PR #8, `df394f7`):** ruff, line length 120, `skip-magic-trailing-comma = true` (without it ruff keeps every exploded call and the code grows), plus a one-off pass dropping blank lines inside functions and `# ----` banner lines. 49.0k → 36.8k lines including tests; syntax tree identical in all 49 files; test outcomes identical. Ruff settings are in `pyproject.toml`. |
| 1.5 | `gene_sets.py` hard-coded 74 gene sets as Python lists (8.3k lines once formatted), and they were **wrong**: no Hallmark set matched any MSigDB release (6.0, 7.0, 2024.1; median Jaccard 0.15, e.g. KRAS_SIGNALING_DN shares 2 of 112 genes), 6 sets exceeded Hallmark's 200-gene maximum, `HALLMARK_TNF_SIGNALING_VIA_NFKB` should be `TNFA`; Reactome/GO BP/KEGG were 24 hand-picked partial sets. Hallmark is a default source, so default annotations used Hallmark names on wrong genes. | verified | **Done (`feat/msigdb-gene-sets`, `cb7679d`; PR after #8):** full MSigDB collections through gseapy (`Msigdb.get_gmt`), release `modules.program_enrichment.msigdb_version: 2026.1`; new `input.species` (human/mouse) picks native mouse collections (replaces `program_enrichment.species` and the symbol-capitalisation step); MSigDB has no mouse KEGG, so mouse + `kegg` without a GMT fails at validation. QC mito/ribo/hb matching is already case-insensitive, so species does not set the prefixes. ORA, BH, GMT reading and the dot plot moved to gseapy (7.5). New dependency: `gseapy` (installed in `py313_torch211`). A failed download only logged a warning (11.17, fixed in #11). **Merged as #9.** |
| 1.2 | `run_pipeline` is one function of ~800 lines (reformatted; `cli.py`), ~90 of them record-keeping (13 A2). | verified | Split into one function per stage, with `run_pipeline` calling them in order. |
| 1.3 | 8 standard/large function pairs, chosen by `cfg.use_large_mode` (auto: ≥1M cells or ≥5,000 targets). In most pairs the large version handles every case and the standard version is the same thing without caps. See section 2. | verified | Merge each pair (section 2). Estimated saving ~700-1,000 formatted lines (the 8 standard versions total 670). |
| 1.4 | 333 config settings across 29 dataclasses; several do nothing (section 6). | verified | Delete dead settings and their checks (13 A5) and make the remaining caps explicit (section 2). |

## 2. Standard vs large mode: merging, and bugs where the modes disagree

The two modes differ in four ways. Each one points to a shared building block:

| Kind of difference | Where | Shared building block |
|---|---|---|
| Reading one gene for a subset of cells. Standard reads the whole column and then masks it; large reads only the needed rows. | perturbation, ps_score | `gene_values(layer, gene, cell_idx)`. It already exists twice: `perturbation.py:235 _gene_values_at_indices` and `ps_score.py:284 _extract_gene_values`. |
| Per-group sums. Standard makes the selected genes dense; large processes chunks of 256 genes and keeps them sparse. | modules effect matrix (`modules.py:861` vs `:1513`) | Loop over gene chunks and compute `indicator @ block` to get sum, sum of squares and linear sum. Standard is the same loop with one chunk. |
| Counting cells. Standard sums boolean masks; large builds a target × cluster count table once. | enrichment (`enrichment.py:535` vs `:1231`) | Build one crosstab `counts[target, stratum, cluster]`, then delete the standard version. |
| Filtering. Standard filters step by step; large builds one keep-mask and slices once. | qc (`qc.py:353` vs `:527`) | Build one mask, slice once, then recount per-gene detection on the kept cells. |

**The real switch is a set of sampling caps, not "large vs standard".** These are the places where the modes differ scientifically:

| Cap | Where | Today |
|---|---|---|
| 100k control cells | perturbation `perturbation.py:105` | hard-coded constant |
| 50k control cells | PS `ps_score.py:139` | hard-coded constant |
| 200k cells for marker-gene selection | modules `_select_genes_large` | `scaling.marker_max_cells` |
| 100 permutations for the omnibus test | enrichment `enrichment.py:107` | constant |

Fix: make each cap an explicit setting, where `None` means "use everything". Today, a dataset that crosses the 1M-cell threshold switches to subsampled statistics without telling the user.

The guide step is already merged this way: `top_two_guides` (`guides.py:1006`) has a single sparse/dense switch, and `tests/test_scaling_consistency.py:194` tests it. The large guide path can call `top_two_guides(X, force_sparse=True)`.

**Bugs where results differ between modes:**

| # | Bug | Evidence |
|---|---|---|
| 2.1 | **Enrichment, large mode, "other targets" comparison.** The reference cluster counts (`c_ref`) include only targets that pass `min_cells_per_target` (`enrichment.py:1463`, `:1500`). The reference total (`n_ref = total_targeting - n_target`, `:1356`, `:1698`) includes every targeting cell. Cells of small targets therefore end up in "not in cluster" for every cluster, which biases odds ratios and p-values toward **enrichment** (corrected 2026-10-06: `c_ref` is too small and `d` too large, so OR = (a·d)/(b·c) is inflated; re-checked on the reformatted code at `enrichment.py:586`, `:612-623`, `:684-691` of `feat/msigdb-gene-sets`). Standard mode is correct (`_reference_mask`, `:200`). | verified; **fixed in #10** (LARGE now matches STANDARD) |
| 2.2 | **QC, large mode, gene filter.** Genes are kept using `var["n_cells_by_counts"]`, which was computed *before* cell filtering. Standard mode recounts on the surviving cells, so large mode keeps genes that fall below `min_cells_per_gene` after QC. Large mode also reuses QC columns already in the input file, which may predate prefiltering (`qc.py:242-256`). | verified |
| 2.3 | Enrichment composition %: standard divides by all of a target's cells, including those in dropped small clusters (`enrichment.py:813`, `:830`). Large divides by kept clusters only (`:1476-1492`). | reported |
| 2.4 | Enrichment large mode returns early before resolving controls when there are no testable targets (`:1341`), and skips the "primary control unavailable" warning. | reported |
| 2.5 | PS large mode promises an LDA on a bounded subset (`ps_score.py:77-79`, `lda_large_max_cells` `config.py:829`), but it only skips the LDA (`:1956-1964`). | reported |
| 2.6 | Precision differences: perturbation large uses float32 vs float64 in standard. PS large casts scores to float32 (`ps_score.py:1018`, `:1777`). Cells sitting exactly on a threshold can be classified differently. | verified (perturbation) / reported (PS) |
| 2.7 | Modules effect matrix edge cases: large mode clips negative variances to 0, standard does not. Control variance is two-pass in standard and one-pass in large. With a single control cell, standard gives p=1 and large computes p. | reported |
| 2.8 | Modules `_select_genes_standard` writes `rank_genes_groups` into the caller's `expr.uns` (`modules.py:306`). There is a dead `try/finally: pass` at `:709-729`. | reported |
| 2.9 | Each stage decides its mode from different inputs. QC and guides look only at cell count (`qc.py:892`, `guides.py:1145`); PS also looks at the number of perturbations (`ps_score.py:2350`). One run can therefore mix modes. | reported |

**First step before merging anything:** write a test that forces every stage through both modes on small synthetic data and compares the output tables. Why first:
- There is no reference today: no test runs the same data through both modes (11.23), so a merge cannot show it reproduces standard mode instead of carrying over a large-mode bug.
- On day one it should fail on exactly 2.1 and 2.2, which confirms the known bugs and shows whether there are others.
- After 2.1/2.2 are fixed it must pass, and every merge step keeps it passing; afterwards it becomes "one mode, same answers as before" and also guards the data-layout change (12.3).
- After the merge the standard enrichment path (~215 lines) disappears (7.8).

## 3. Memory

Notation: C = one full sparse copy of the cells × genes matrix (nnz × 8-12 bytes). The docs estimate KOLF (2.66M cells) at 5-15 GB per copy and 135-155 GB peak for the run (`docs/PIPELINE_DEVELOPMENT.md:575`, `:582`). `docs/slurm.md:58` asks for 180-200 GB for about 2.5M cells.

| # | Problem | Evidence | Fix |
|---|---|---|---|
| 3.1 | **The h5ad is always loaded fully into memory** (`io.py:315 sc.read_h5ad(path)`, no backed or lazy option). | verified | Section 5 covers on-disk options. |
| 3.2 | **At least 2C held for the whole run**: `counts = X.copy()` at load (`io.py:478`), then `X = counts.copy()` and normalization (`cluster.py:149`). Standard mode also copies `lognorm` (`cluster.py:174`), giving 3C. | reported | Accept 2C (counts plus normalized) as the minimum. Drop the extra copy in standard mode. |
| 3.3 | **Memory peaks at QC because of copies.** `all_cells = expr.copy()` is made by default to write the unfiltered h5ad (`cli.py:1009`; `write_unfiltered_h5ad=True`, `config.py:1179`). In standard mode `_write_unfiltered_object` copies again to merge guides before writing (`cli.py:568`). Pair mode adds `_tmp = expr.copy()` just to compute a few QC columns (`cli.py:1026`). Each filter step slices into a new copy while the old one is still alive. Total: 4C, or about 6C in standard mode. The `cli.py:1009` copy happens **in large mode too**, which undoes the no-copy large branch inside `_write_unfiltered_object`. | verified (1009, 568) / reported (1026, large-mode note) | Write the unfiltered file without copying (attach the guide obsm, write, detach). Compute pair-mode QC columns without copying. Make the unfiltered output opt-in. |
| 3.4 | The HVG subset for PCA, `expr[:, hvg_mask].copy()` (`cluster.py:392`), copies X, counts and lognorm separately. The 200k-cell HVG sample copies X and lognorm separately even though they are the same object. | reported | Slice only the layer that PCA needs. |
| 3.5 | Standard PCA uses `scale(zero_center=True)`, which makes the matrix dense: cells × 3000 genes is about 30-60 GB (`docs/PIPELINE_DEVELOPMENT.md:580`). | reported | Use `zero_center=False` with a sparse PCA everywhere, as large mode already does. |
| 3.6 | Modules: `sc.tl.score_genes` (`modules.py:2844`) reads all genes on every call, so there is one temporary full copy per program (scanpy `tools/_score_genes.py:216,231`). Standard mode also builds a dense cells × genes float64 block (`modules.py:827-852`) and copies `expr` if the obs column isn't categorical (`:316`). | reported | Compute program scores from the gene subset directly. |
| 3.7 | PS standard mode copies the whole lognorm layer (`ps_score.py:348`). | reported | Fixed by 7.4 (new pertps). |
| 3.8 | lochNESS builds a k=300 neighbour graph (`lochness.py:314`, `config.py:850`): about 15-25 GB at 2.5M cells. | estimate | Check whether k=300 is needed, or reuse the stage-4 graph. |
| 3.9 | Output: in large mode X and lognorm are the same object but are written twice; the file is gzip-compressed (`io.py:1008`). | reported | Write each matrix once. Consider lzf or no compression for speed. |
| 3.10 | `_collect` is only `gc.collect()` (`cli.py:289`, `cluster.py:77`) and frees nothing that is still referenced. | reported | Delete references explicitly where copies are made. |

**Which way each stage reads the matrix.** The matrix is stored by rows (CSR) everywhere, so reading one gene scans every stored value.

| Reads one gene across all cells (expensive on CSR) | Reads a subset of cells (cheap) | Uses the whole matrix | Never reads expression |
|---|---|---|---|
| perturbation standard (one scan per target, about 2k scans, `perturbation.py:699`); PS standard (`ps_score.py:690`); modules large (chunks of 256 genes, `modules.py:~1650`); knockdown filter (all targets in one pass, `knockdown_filter.py:93`); QC gene filter (`qc.py:809`); HVG subset (`cluster.py:392`); split-features load (`io.py:457`) | perturbation large (`perturbation.py:257`); PS large (`ps_score.py:441`); QC cell filter (`qc.py:752`); HVG sample; `assigned_only` subset | QC metrics; normalize/log1p; `score_genes` | neighbours, UMAP, Leiden; lochNESS; distances (PCA only); enrichment; meta table |

Fix 3.11: make one column-oriented (CSC) copy of just the ~2k target genes and use it for stage 5, PS and the knockdown filter. **Superseded (2026-10-06) by section 12:** these stages need only each cell's own-target value and control statistics, which a row-wise pass provides without a CSC copy.

## 4. Speed

Rough ranking of compute cost at about 2.5M cells, 2k targets (estimate):

1. **Stage 11, distance between targets** (off by default). It computes all K(K-1)/2 pairs, and each pair runs `cdist` three times on up to 2,000 cells (`distance.py:789`, `:889-929`). At K=2,000 that is about 2M pairs and about 10^15 operations.
   - About two thirds of the work is repeated: each target's within-group term d(x,x') is recomputed in every one of its K-1 pairs.
   - Every pair's cells are sent separately to loky worker processes, about 800 KB per pair (`cpu_parallel_backend="loky"`, `config.py:1088`). That could add up to terabytes of data copying (estimate).
   - Fix: compute each within-group term once, share the embedding with the workers instead of sending it per pair, or use `rsc.ptg.Distance.pairwise` on GPU (section 5).
2. **Stage 10, distance permutation tests.** Per target: one (n+m)² distance matrix (up to 7,000², about 390 MB per worker), then 1000 permutations. Each permutation copies an n×n block with `D[np.ix_]`, about 32 MB at n=2000.
   - Cost grows with the number of targets times the cap squared. It is **not** O(N²) in total cells, because of the caps.
3. lochNESS k=300 neighbour search.
4. **Neighbours, UMAP, Leiden** on CPU (`cluster.py:592-620`). This is the heaviest step among those on by default.
5. Repeated reads of one gene across all cells (section 3 table).
6. gzip write of about 2C.

| # | Problem | Evidence | Fix |
|---|---|---|---|
| 4.1 | Within-group distances recomputed for every pair in stage 11. | verified (`compute_energy_distance` per pair) | Compute once per target. |
| 4.2 | Per-pair data sent to worker processes in stage 11. | verified (task tuples carry arrays) / estimate (size) | Shared memory (`SharedArrayBuffer` already exists in `data_access.py`), or threads. |
| 4.3 | 1000 permutations with BH correction over thousands of targets: the smallest possible p is about 1e-3, which caps how many discoveries are possible. | estimate | Check against the FDR level used. Consider an asymptotic or GPU test. |

## 5. GPU and on-disk processing

| # | Problem | Evidence | Fix |
|---|---|---|---|
| 5.1 | **The GPU is barely used.** The only real GPU computation is the CuPy correlation in modules (`modules.py:2270-2293`). | verified | — |
| 5.2 | **The clustering GPU decision is ignored.** `resolve_stage_backend("clustering")` is computed and logged, then neighbours, UMAP and Leiden always run as CPU scanpy (`cluster.py:540`). `rapids-singlecell` is listed as an optional dependency but only checked for availability, never imported (`compute.py:326`, `:422`). The `distance_space` GPU decision is ignored as well (`distance.py:889-921`). | verified | Wire up rapids-singlecell behind the existing decision (table below). |
| 5.3 | GPU is tested only with mocks (`tests/test_compute.py:91-154`). | reported | Add one real GPU smoke test, run on Slurm. |
| 5.4 | The GPU Spearman ranking uses `argsort` (ties broken arbitrarily; reformatted `modules.py:623-626`); CPU uses pandas (tied ranks averaged). The two can disagree when there are ties. | verified (2026-10-06 sweep) | Use average ranks on GPU too, or drop the GPU branch (13 A3). |

**rapids-singlecell 0.15.1 is already installed in `py313_torch211`.** Mapping:

| Pipeline step | rapids-singlecell | Expected gain (from general knowledge, not measured here) |
|---|---|---|
| neighbours / UMAP / Leiden | `rsc.pp.neighbors`, `rsc.tl.umap`, `rsc.tl.leiden` | largest |
| lochNESS k=300 graph | `rsc.pp.neighbors` | large |
| normalize, log1p, HVG, scale, regress_out, PCA | `rsc.pp.*` | moderate to large |
| Harmony (currently harmonypy, `cluster.py:785`) | `rsc.pp.harmony_integrate` | large |
| QC metrics, filters, Scrublet | `rsc.pp.calculate_qc_metrics`, `filter_cells`/`filter_genes`, `scrublet` | moderate |
| modules: markers, program scores | `rsc.tl.rank_genes_groups`, `rsc.tl.score_genes` | moderate |
| stage 11 distance matrix | `rsc.ptg.Distance("edistance").pairwise`: GPU, computes each within-group term once, supports multiple GPUs. It excludes self-distances, so values shift slightly. | large |
| stage 10 permutation test | not available in rsc (energy distance only, no permutations) | needs a custom CuPy port |

Pattern: `rsc.get.anndata_to_GPU` → `rsc.*` calls → `rsc.get.anndata_to_CPU`, placed behind `resolve_stage_backend`. The matrix must fit in GPU memory: at about 2.5M cells one copy is 5-15 GB, which fits on an H200 or A100.

**On disk / Dask: not usable in the current state.**

| # | Problem | Evidence |
|---|---|---|
| 5.5 | Nothing reads lazily or from disk: no Dask, no `backed=`, no `read_lazy` (every read at `io.py:108/215/315/341/1161/1205/1246` loads fully). `StorageConfig` and `data_access.resolve_storage_mode` / `log_storage_diagnostics` are never called. | verified |
| 5.6 | Code that would break, or silently load everything into memory, if X were a Dask array: 21 `issparse` checks (`cluster.py:396/443/478` send Dask data down the dense path, whose PCA needs `dask_ml`, which is not installed); `.toarray()` 14 times; `csr_matrix(` and `csc_matrix(` 30 times; `np.asarray(` 105 times; 21 in-place assignments to X or layers; full `.copy()` calls. | reported |
| 5.7 | scanpy 1.12.1 supports Dask in normalize_total, log1p, calculate_qc_metrics, HVG (seurat/seurat_v3), scale, and PCA on sparse Dask (covariance method). It has no Dask support in regress_out, neighbors, umap, leiden or score_genes. Installed: anndata 0.12.16 (`read_lazy`, `read_elem_lazy`), dask 2026.1.1. | reported |

Order of work: remove copies (3.3-3.6) → CSC copy of target genes (3.11) → GPU for stage 4 and stage 11 → Dask only if the 2C baseline still doesn't fit. Stages that could realistically run on Dask: QC, normalize, HVG, PCA.

## 6. Settings that do nothing

| # | Setting | Evidence |
|---|---|---|
| 6.1 | `storage.*` (`mode`, `backed_threshold_cells`, `prefer_backed_h5ad`, `keep_embeddings_in_memory`; `config.py:1110-1141`). Only `shared_worker_arrays` is read (`distance.py:632`). | verified |
| 6.2 | `compute.backend: gpu` for clustering and distance_space (5.2). | verified |
| 6.3 | `distance.primary_metric`: read but unused; stage 10 always computes energy distance (`distance.py:470`). | reported |
| 6.4 | `ps_score.lda_large_max_cells` (2.5). | reported |
| 6.5 | Never read anywhere: `input.guide_table_guide_column` (`config.py:186`; the loader at `io.py:618-640` reads only the cell, gene and count columns), `ps_score.lda_large_stratified` (`:833`), `lochness.store_all_target_scores` (`:873`), `storage.keep_embeddings_in_memory` (`:1137`), `qc.doublets.min_gene_counts` (`:1347`), `guides.fastq.restrict_to_gex_barcodes` (`:1413`; the restriction is always applied, `guide_counting.py:23`). | reported (field-by-field search of 343 settings) |
| 6.6 | Checked for valid values but never used: `lochness.target_chunk_size` (`:869`), `lochness.max_targets_in_obs` (`:877`), `modules.program_enrichment.method` (`:704`; ORA is always used). | reported |
| 6.7 | `program_enrichment.*` settings are read through `getattr(cfg, name, default)`, with the defaults hard-coded a second time (`gene_sets.py:1531-1539`, `plots.py:8214-8216`). | reported |
| 6.8 | `cluster.scale_max_value=10` only applies in standard mode; large mode never clips (`cluster.py:442-462`). | reported |

Removal of these settings and their validation blocks is item 13 A5; 6.2 goes with 13 A3.

Good news: `config/default.yaml` matches the dataclass defaults exactly (333 settings compared, no differences). Every shipped example config parses and passes validation.

## 7. Replacing code with existing libraries

**pertpy.** pertpy is never imported. One module is a real duplicate (distance) plus two small test functions, about 400-450 removable lines (~1% of the package).

| Pipeline part | pertpy | Relationship | Verdict |
|---|---|---|---|
| `distance.py` stage 10 (energy distance + permutation test) | `Distance`, `DistanceTest` | Same statistic | Keep the pipeline's version. |
| `distance.py` stage 11 (pair matrix) | `Distance.pairwise` | Same idea; pertpy excludes self-distances | Use `rsc.ptg` (GPU) instead. |
| `gene_sets.py` ORA test | `Enrichment.hypergeometric` | Same test | Keep. pertpy expects scanpy marker output, and the gene sets bundled with the pipeline would still be needed. |
| `perturbation.py` target-gene test | `WilcoxonTest` | Overlaps only on Mann-Whitney | Keep. |
| guide assignment, enrichment, knockdown filter, lochNESS, PS, modules | `GuideAssignment`, `Sccoda`/`Milo`, `Mixscape`, `PseudobulkSpace` | Different methods or questions | Keep. |

Reasons not to import pertpy:

- It pulls in a JAX stack (ott-jax, numpyro, funsor, flax), and `import pertpy` loads every tool, so there is no way to import only the distance code.
- It needs Python 3.12 or newer.
- Its permutation RNG is unseeded (`_distance_tests.py:264`).
- `DistanceTest` builds a dense matrix for every group (about 390 GB for 1,000 targets at the pipeline's caps; estimate).
- Its default p-value formula, multiple-testing correction (Holm-Šidák) and MMD kernel (linear) all differ from the pipeline's.

**Other replacements:**

| # | Replacement | Gain |
|---|---|---|
| 7.1 | Modules pseudobulk: `sc.get.aggregate(..., func=["sum","mean","var"])` (scanpy, already a dependency). | The standard builder (458 lines, about 70 of real logic) shrinks to about 25 lines. |
| 7.2 | Hand-written BH → `statsmodels.stats.multitest.multipletests(method="fdr_bh")` (statsmodels is already a dependency). The `gene_sets.py` copy is gone on `feat/msigdb-gene-sets` (gseapy does it). Remaining: `perturbation.benjamini_hochberg` (reformatted `perturbation.py:208`), imported by enrichment, distance, modules and `pair_guide_report`. Mask NaN first (`scipy.stats.false_discovery_control` raises on NaN). | One tested implementation, ~17 lines. |
| 7.3 | Stage 11 pair matrix → `rsc.ptg.Distance("edistance").pairwise` (section 5). | GPU speed, no repeated within-group work. |
| 7.4 | **PS stage → new `pertps.run_ps_score_exact_fast`** (local `data_perturb_v2/PS_python`, branch `import-pertscore-py`). One call scores all targets, streaming 8,192-cell row chunks from memory or from disk; uses every control cell; selects each target's signature genes with a Welch t-test and fits a ridge model; supports multilabel `GENE1+GENE2` input (pair-guide mode). | Replaces both PS paths of `ps_score.py` (888 lines → ~100): the full lognorm copy of standard mode (3.7), the per-target loop with a 50k-control cap and float32 scores of large mode (2.5, 2.6), so PS no longer depends on dataset size. |
| 7.5 | **Done (`feat/msigdb-gene-sets`):** program ORA → `gseapy.enrich` (p-values identical to the old code on 130 terms), FDR from gseapy, plot → `gseapy.dotplot` (terms picked by p-value first so the labelling term is always shown; no plot when nothing is significant, where the old plot drew 2 non-significant terms; dot size = fraction of set in program). **FDR family:** gseapy 1.3.1 skips sets sharing 0 genes with the program (`stats.calc_pvalues`, `if x < 1: continue`) and runs BH on the rest (`enrichr.py` ~707), undocumented; `min_overlap` (2) is applied afterwards. Old code counted only sets with ≥2 shared genes. Principled family would be all size-eligible sets: on 3,805 eligible sets and a planted 50-gene program, 31 / 52 / 105 significant terms for all / ≥1 / ≥2. **Decided 2026-10-06: keep gseapy's ≥1 family** (documented in `docs/methods.md`). | ~100-150 lines removed; p-values unchanged |
| 7.6 | **Done (`feat/msigdb-gene-sets`):** `parse_gmt` → `gseapy.read_gmt`. | ~25 lines |
| 7.7 | Enrichment omnibus "permutation" p (reformatted `enrichment.py:200-218`) is really a parametric multinomial simulation from global cluster frequencies (fixes row totals only). `chi2_contingency(table, method=PermutationMethod()/MonteCarloMethod())` (scipy 1.17 installed) fixes both margins: a slightly different null. | ~15-20 lines; p changes slightly |

**Library sweep (2026-10-06, read-only, reformatted line numbers).** Already library-backed: 10x reading, QC metrics, filters, scrublet, KS/MWU, Fisher, CMH, normalisation/HVG/PCA/neighbors/UMAP/Leiden, hierarchical clustering.

| # | Candidate | Replacement | Same? | Saves | Dependency |
|---|---|---|---|---|---|
| 7.8 | Standard enrichment path (`enrichment.py:283-497`) | keep only the large path (pandas crosstab) after fixing 2.1; part of the section 2 merge | same after the fix | ~215 | none |
| 7.9 | `Config.validate` (`config.py:1179-1598`, 420 lines, 130 raises) | pydantic dataclasses with `strict=True` (keep the two silent rewrites at `:1241`, `:1243` by hand); or read `Literal` type hints at runtime (~15 lines) | same pass/fail, different messages | ~200 (pydantic) / ~50 | pydantic 2.x (installed, not declared) |
| 7.10 | Top-two guides per cell ×4 (`guides.py:214-375`, `dual_guides.py:239-253`, `guide_qc.py:184-202`) | one numba kernel (no package does this) | same except tie order | ~110 | none |
| 7.11 | CPU detection, `limit_blas_threads` (`compute.py:50-93`, `:466-485`) | `os.process_cpu_count()`, `joblib.parallel_config(inner_max_num_threads=…)` | same if Slurm pins CPUs (unverified) | ~60 | none |
| 7.12 | `_run_harmony` (`cluster.py:405-459`) | `sc.external.pp.harmony_integrate` | same call; orientation check differs (harmonypy not installed to verify) | ~45 | none |
| 7.13 | Two 10x readers (`io.py:977-1072`) | one reader choosing by `path.is_dir()` | same | ~35 | none |
| 7.14 | `compute_pcoa_coordinates` (`distance.py:565-609`) | `sklearn.manifold.ClassicalMDS(metric="precomputed")` | same up to axis sign; still cap at positive eigenvalues | ~30 | none |
| 7.15 | Welch t ×2 (`modules.py:367`, `:541`), `_nb_logpmf` (knockdown filter), `_mad_bounds` (`qc.py:544`) | `ttest_ind_from_stats(equal_var=False)`, `nbinom.logpmf`, `median_abs_deviation(scale="normal")` | same | ~20 | none |
| 7.16 | Gzip FASTQ reading (`guide_counting.py:133-176`) | xopen | same | ~40 | xopen (skip: too little gain) |

Keep hand-written: single-guide assignment (crispat's rule differs and pulls torch/pyro; #17 adds crispat's Poisson-Gaussian mixture and fishash as numpy ports instead of dependencies, section 15), knockdown count model, FASTQ scaffold matching, lochNESS kernels, MMD, Haldane OR, stratified subsampling, energy-distance permutation shortcut (dcor rebuilds the full matrix per permutation, scales the statistic by nm/(n+m) and counts `>` instead of `>=`). The unused `energy_distance_from_cdist` is 13 A1.

What 7.4 (PS → pertps) needs before it can go in:

- The pipeline imports the legacy `PerturbAnalyzer`. It still exists on lab `main` (`pertps_project/pertps/analyzer.py`); on `import-pertscore-py` it survives only under `examples/legacy_ps_python/`. It is not installed in `py313_torch211`, so 6 PS tests are skipped (11.25).
- The pipeline's `ps` extra points at `weili-lab/PS_python#subdirectory=pertps_project` (legacy layout). The new code is on the fork with 3 uncommitted files (+315 lines, checked 2026-10-06), so it needs committing, merging and pinning first.
- The new pertps has no LDA, so the PS LDA/UMAP figure is either dropped or kept separately.
- Pass `layer="counts"` and match `target_sum` to the pipeline's normalization, because pertps normalizes from counts itself.
- Control cells get no score. The standard-mode quadrant columns for control cells, and possibly the LDA, need adapting.
- Scores will differ from the legacy `PerturbAnalyzer`. Compare the two on one dataset.

## 8. Knockdown filter (PR weili-lab/perturbseq-pipeline#7, merged; non-carrier baseline #16, merged)

Updated 2026-10-06; 8.7 added 2026-10-22. PR #7 was merged into upstream `main` on 2026-10-05 (`3164a97`) with all four commits: `a08071c` (filter), `02442f4` (mean replaces median), `4670ee9` (count model), `5b20baf` (exact fit instead of EM). Everything in 8.6 is follow-up work for a new PR. Benchmark scripts in `scripts/kd_filter_benchmark/`, results in `results/kd_filter_benchmark/` (both in `data_perturb_v2/`, outside the pipeline repo).

- 8.1 The filter uses the target's own mRNA, so it is valid for CRISPRi/knockdown screens. For CRISPR knockout, the mRNA often persists, and a method that looks at the cell's overall phenotype (e.g. `pertpy.tl.Mixscape`) is the right tool. Add a sentence to the PR and README saying this. **Open.**
- 8.2 `kd_keep` drops cells shown to fail (`escaper`, `failed_group`) **and** cells that can't be tested because the target is unmeasured or barely expressed (`not_measured`, `low_control_expression`). Cells in groups too small to test (`non_testable`) are kept. Intentional for now; note it in the docs. **Open.**
- 8.3 **Median group test failed under dropout. Fixed in `02442f4`.**
  - The old step 1 took the *median* of per-cell ratios. For a target detected in fewer than about 50% of controls, most target cells are zero from dropout alone, so the median ratio was 0 and the group passed with no knockdown; step 2 then marked every detected cell an escaper.
  - Mean of ratios equals ratio of means here (the control mean is constant within a group), so the fix was median → mean of `kd_ratio`. Status `failed_median` became `failed_group`, column `passed_median` became `passed_group`.
  - The proposed second detection gate (50-80% of controls) for the per-cell step was **rejected**: it would leave most targets with no per-cell filtering. The count model (8.4) handles low expression instead.

### 8.4 `method: count_model` (added in `4670ee9`, fitting changed in `5b20baf`)

`mean_ratio` stays the default. The count model works on raw counts (`layers['counts']`).

- **Known per cell**: raw target count xᵢ, size factor sᵢ (library / median library), control mean μ and dispersion φ of the gene in the cell's context. Expected unperturbed count mᵢ = μ·sᵢ.
- **Model**: each cell is an escaper, NB(mᵢ, φ), with probability π, or knocked down, NB(ρ·mᵢ, φ). ρ and π are estimated once per group.
- **Dispersion** is shared across genes, DESeq2-style: gene-wise moment estimate, trend φ(μ) = a + b/μ fitted robustly on log scale, shrinkage weighted by a sampling variance measured from split halves in 10 equal-count expression bins, prior variance floor 0.25. With 30 controls this cut the median log error from 0.458 to 0.355 (unit test).
- **Fit** (`5b20baf`): for each of 100 ρ values (0.01-1.00) the log-likelihood is concave in π and peaks where π equals the mean P(escaper), found by bisection; the best (ρ, π) pair wins. This replaced EM, which hit its 200-iteration cap in ~31% of simulated groups (flat likelihood at low expression or weak knockdown) and stopped short, underestimating π from its start at 0.1. The exact fit is never worse than EM in likelihood; it is 2-10× slower per group (median 1.5 ms at 30 cells, 70 ms at 2,000), which is still small.
- **Per cell**: pᵢ = P(escaper | xᵢ) = π·f_E / (π·f_E + (1−π)·f_K), i.e. posterior odds = prior odds × likelihood ratio. π is an empirical-Bayes prior estimated from the same group; uncertainty in ρ and π is not propagated into pᵢ.
- **Group passes** when ρ < `max_rho` (0.3) **and** π < `max_escaper_fraction` (0.5). The second condition is needed because a group with no knockdown can be fitted as "a few cells near zero plus mostly escapers", giving a tiny ρ.
- **Cell is an escaper** when pᵢ ≥ `min_escaper_prob` (0.9), equivalent to a likelihood ratio ≥ 9(1−π)/π (36× at π = 0.2, 171× at π = 0.05).
- `kd_ratio` and `group_mean_ratio` are still written but do not drive any count-model decision; `max_mean_ratio` and `max_cell_ratio` have no effect under it.
- **Raw-count check**: the run stops if the first 10,000 cells of `layers['counts']` are not non-negative integers. Applied only to the count model, because the pipeline officially supports normalized-only input (`input.normalized_layer: X`), which a pipeline-wide check would break (see 11.9).
- **Bug found and fixed during benchmarking**: a lane with 10 controls gave collapsed quantile bins, an empty bin, NaN prior variance and NaN φ for the whole lane; groups then passed at ρ = 0.01. Fixed with equal-count rank bins, shrinkage written in weight form, and an error on non-finite likelihoods.

### 8.5 Benchmark results

**Simulation** (NB counts with known escapers; μ 0.2-50, ρ 0.1/0.3/0.5/1.0, 0/10/30% escapers, 30/100/300 cells, φ 0.1 and 0.5, 50 repeats). Run with EM and before the escaper-fraction rule; not rerun.

| φ = 0.1 | count_model | mean_ratio |
|---|---|---|
| ρ = 0.1 with 30% escapers: pass rate | 93% | 11% |
| ρ = 0.1: knocked-down cells wrongly called escapers | ≤ 0.4% | ~3%, up to 11% at μ ≈ 1.5 |
| ρ = 0.1: precision of escaper calls | 98% | 69-79% |
| ρ = 0.1, 10% escapers: recall of escapers | 47% (4% at μ = 0.2, 96% at μ = 25) | 69% |
| ρ = 1 (no knockdown): pass rate | ~1% (4% at μ = 0.2) | ≤ 0.3% |
| ρ = 0.5 (should fail): pass rate | 6-7% (17-19% at μ ≤ 0.4) | 1-3% |
| Threshold-free ranking of escapers vs knockdowns (AUROC) | 0.865 | 0.858 |

**Nadig HepG2 CRISPRi** (145k cells, 9.6k genes; Slurm). Numbers are for the current code (exact fit, escaper-fraction rule); "per lane" uses `context_key: gem_group`.

| | count_model | mean_ratio |
|---|---|---|
| Real targets passing (one baseline / per lane) | 1,515 / 1,495 | 1,426 / 1,400 |
| Escapers called, share of cells in passing targets | 3.2% | 6.4% |
| Passing targets with zero escapers (one baseline) | 758 of 1,515 | 164 of 1,426 |
| Pseudo-target false passes (one baseline / per lane), of 1,029 | 3 (0.29%) / 4 (0.39%) | 2 (0.19%) / 3 (0.29%) |
| Same, before the escaper-fraction rule | 17 (1.65%) / 26 (2.53%) | — |
| Strong-effect targets: escapers closer to controls than knockdowns are (significant) | 97-98% (84%) | 91-92% (62-63%) |
| Targets evaluable under both (154-156): closer (significant) | 80-82% (51-53%) | 83% (55-60%) |

- Pseudo-targets: half the control cells relabelled as fake targets with real gene names and real group sizes, balanced across expression fifths, 30 draws. Any pass is false. Every false pass, for both methods, is in the lowest expression fifth (control mean 0.1-0.4 counts per cell).
- Escaper-vs-control test: energy distance in a PCA-50 embedding built with all target genes removed, permutation p-values, targets with ≥ 5 escapers and ≥ 5 knockdowns.
- The extra count-model passes are mostly strong knockdowns with escapers, which the mean-ratio group test fails because escapers raise the group mean.
- On the targets both methods can evaluate, the two methods' escaper calls are about equally control-like; mean ratio is slightly ahead on significance there.
- EM → exact fit on Nadig: 8 (one baseline) and 16 (per lane) targets changed verdict, almost all at a cutoff (ρ 0.29 vs 0.30, π 0.48-0.50 vs 0.50-0.51); escaper calls rose ~9% (3,196 → 3,469). Pre-change outputs are in `results/kd_filter_benchmark/nadig/em_baseline/`.

### 8.6 Open decisions

1. **`max_escaper_fraction` cutoff (0.5).** Two remaining false passes (RPP14, HERC1) have π = 0.38-0.45, so ~0.35 would remove them. But strong knockdowns sit right at 0.5 after the exact fit (PRMT5, UBA1, WDR74, CLOCK, DYNC1H1, SMG1: ρ 0.01-0.16, π 0.50-0.61) and now fail. Decide whether a group with about half escapers should pass (escapers marked) or fail; needs a sweep of cutoffs on Nadig.
2. **Default method**: keep `mean_ratio` or switch to `count_model`.
3. **`max_mean_ratio`** 0.3 vs ~0.5.
4. **Detection gate** (`min_pct_expressing_control` 10%): whether it should still apply to the count model.
5. **`min_control_cells`** = 10: enough for the per-context control fit? Since #16 it counts the baseline cells, i.e. the target's non-carriers in high-MOI mode.
6. **Weak-knockdown leniency at very low expression** (simulation: ρ = 0.5 passes 17-19% at μ ≤ 0.4): not rechecked with the escaper-fraction rule or the exact fit.
7. **PR #7 description** was out of date when merged (median, old names, no count model). The README and module docstring are current, so this only matters as a record; optionally add a comment on the merged PR.
8. **Scaling** to 10M cells is held for the pipeline-wide redesign (12.6).

### 8.7 Non-carrier baseline in high-MOI mode (#16, merged 2026-10-19, `c80647f`)

- **Why:** in high-MOI mode most cells carry several targeting guides, so NTC-only cells are few and carry fewer guides than the cells they are compared with. With fewer than `min_control_cells` of them in a context, the NTC baseline left every target `non_testable`; otherwise the baseline came from a small, systematically different group.
- **What:** new `knockdown_filter.control`: `ntc` | `non_carriers` | `null`. `non_carriers` (high-MOI only) takes each target's baseline from the targeting cells of the same context that do not carry it (`MembershipIndex.other_mask(target)`), the same "other" control the high-MOI statistics stages use. `null` (default) resolves to `non_carriers` in high-MOI mode and `ntc` otherwise, so single-guide and dual-guide-pair runs are unchanged. `control: ntc` restores the old high-MOI behaviour.
- **mean_ratio:** `kd_ratio` = target expression / mean in that target's non-carriers of the context.
- **count_model:** the dispersion trend over genes is fitted once per context on all targeting cells; each target's μ and gene-wise φ come from its own non-carriers and are shrunk toward that trend. `_control_fit` was split into `_dispersion_prior` and `_shrink_phi`; the `ntc` path gives the same numbers as before (unit test).
- Output layout unchanged (`obs['kd_ratio']`, `kd_status`, `kd_keep`, `obsm['kd_*_membership']`); there is no column recording which baseline was used (left out to keep the legacy golden tables unchanged).
- **Open:**
  1. Not benchmarked on a real high-MOI dataset (non_carriers vs ntc); tests are synthetic.
  2. Non-carriers carry other targets. If one of those lowers this target's gene (same pathway, or a CRISPRi guide near this gene's promoter), the baseline drops and the knockdown is underestimated. Excluding cells that carry known co-regulating targets would need a target × gene effect matrix (stage 5 output) first.

## 9. Packaging (Nextflow / nf-core)

Not now. The analysis is one in-memory AnnData passed through 14 dependent stages, so Nextflow would either wrap a single Python process or force an h5ad write and re-read between stages. nf-core compliance (template, modules, nf-test, containers, MultiQC, review) is a large ongoing commitment.

Do instead: ship an Apptainer/Docker image of the current command-line tool. Optionally, use Nextflow only for the per-sample front end (alignment and FASTQ guide counting per lane), which is genuinely parallel.

## 10. Scorpus / Lance

Not a speed win for this pipeline:

- Scorpus stores one row per cell, which is ideal for random minibatches in training.
- The pipeline's repeated reads are by gene, which is the worst case for that layout.
- Stage 4 needs the whole matrix, and stages 9-11 use only the PCA and neighbour graph.
- The published benchmarks compare against a naive on-disk h5ad, not against loading the file into memory.
- Conversion takes about 2 h per 3M cells plus a manual schema step. Counts must be integers, guide matrices are not stored, and results cannot be written back.

Possible later uses: memory-bounded streaming of pseudobulk or target-vs-control sums (scorpus `pp/de.py` already has a one-pass t-test), rerunning a subset of an already-converted screen, and keeping one copy of the data shared with pertTF training.

---

## 11. Other problems (second sweep)

### 11a. Scientific and statistical problems (most consequential first)

| # | Problem | Evidence | Fix |
|---|---|---|---|
| 11.1 | **The distance test cannot call targets significant in large screens.** p = (1+k)/(1+B) with B=1000 (`distance.py:441`, `config.py:904`), so the smallest possible p is about 1e-3. BH then runs over all targets (`:700`). At FDR 0.05, nothing is called until about T/50 targets sit at that floor, e.g. about 200 for a 10k-target screen. `distance_significant` then reflects how many permutations were run, not biology. | verified (arithmetic) | **Decided by simulation (2026-10-03): use the spectral null with k=50.** E = aᵀMa with M = −HDH. Null ≈ c·Σλᵢzᵢ², where c = (1/n+1/m)/(N−1). Take the top 50 eigenvalues of M with `scipy.sparse.linalg.eigsh`, approximate the rest by a gamma matched to the leftover trace and sum of squares, and get p with Imhof (clip at about 1e-12). Rejected: both gamma approximations are 5-7x too optimistic at p≤1e-3 and 100x+ too optimistic below 1e-4. Generalized Pareto: 11 of 90 alternative cases gave p=0, otherwise errors up to ±1.3 in log10. Always report the energy distance next to the FDR. See the simulation results below. |
| 11.2 | **Distance with the "other" control includes each target's own cells in its reference.** One control sample is drawn for all targets (`distance.py:579-599`), and "other" means all targeting cells (`perturbation.py:164`). Stages 5, 6 and modules remove the target's own cells. This shrinks the distance. It is reached when `perturbation.primary_control="other"` or through the automatic fallback when there are too few NTC cells. In theory, a fraction f of the target inside the reference scales the energy distance by exactly (1−f)²: about 4% smaller at f=2% (a 50-target screen), 19% at f=10%. Identical cells appearing on both sides shrink it further. | verified | **Agreed (2026-10-03):** for each target, the "other" reference excludes that target's cells, in every stage. **Done in #10** for the distance stage (the last stage that kept them); in high-MOI mode "other" is membership-aware (#13). |
| 11.3 | **"Hit" / "effective knockdown" means different things in different places.** Stage 5, the report and the meta table use KS FDR < 0.05 plus log2FC < **0.0** (`perturbation.py:1164-1177`, `config.py:593`), so any significant decrease counts, even 5%. The knockdown filter requires about 70% knockdown. PS uses quadrants; distance uses FDR on energy distance. When the KS column is missing, `meta.py:137-139` falls back to the one-sided Mann-Whitney FDR. The KS test is two-sided (`perturbation.py:381`). | verified | Define one knockdown criterion with a minimum effect size, and use it everywhere. |
| 11.4 | **Lane/batch is handled differently in each stage, and is ignored by default in the testing stages.** Stage 5 pools NTC cells across lanes. Stage 6 uses CMH only when `stratify_by` is set (default None, `config.py:680`), and its omnibus test is never stratified. Distance samples proportionally within target and within control separately but never matches lane composition between them. Distance uses uncorrected `X_pca`, lochNESS automatically prefers `X_pca_harmony` (`lochness.py:268-272`), and HVG selection ignores `batch_key`. | verified | Decide on one lane policy (stratify or match lanes) and one embedding for all stages. |
| 11.5 | Stage 6 with CMH takes the effect size and the direction from different analyses. The odds ratio is the CMH pooled one (`enrichment.py:1010-1012`), but `direction` comes from pooled unstratified percentages (`:1035-1039`), so under Simpson's paradox a positive log OR can be labelled "depleted". A non-finite CMH odds ratio silently falls back to the unstratified one while keeping the CMH p-value. `significant` ignores `low_power` (`:1914-1925`). | reported | Take direction from the CMH odds ratio, and flag the fallback. |
| 11.6 | The default primary control differs by stage: NTC for stages 5, 7 and 8, "other" for stage 6 (`config.py:583, 668, 744`). Distance follows the *configured* control, not the one stage 5 ended up using after a fallback. | reported | Use one control setting for the whole run. |
| 11.7 | Pseudocounts and NaN handling in log2FC/BH differ. Stage 5 uses 0.01 (`perturbation.py:98`) while modules use 1e-9 (`modules.py:99`), so lowly expressed genes get extreme log2FC in modules. `target_sum=None` makes the units depend on the dataset. Modules set NaN p-values to 1 before BH (enlarging the set of tests), while stage 5 drops them. Both hand-written BH implementations are correct, but they treat inf differently. | reported | Use one log2FC helper and one BH (statsmodels). |
| 11.8 | A precomputed `lognorm` layer skips normalization (`cluster.py:111-127`), but every stage applies `expm1`, i.e. assumes natural-log `log1p`. A log2 input would be silently wrong. | reported (effect inferred) | Check or document the expected log base. |
| 11.9 | **Raw counts are assumed without a check.** `_ensure_counts_layer` copies X into `layers['counts']` with no check that X holds integer counts (`io.py:475-478`). A normalized X silently becomes "counts", and the copy duplicates the matrix (3.2). | verified | **Agreed (2026-10-03):** check that the counts layer is non-negative integers (on `.data` for sparse) and stop with an error otherwise. **Partly done (`4670ee9`):** the count model checks the first 10,000 cells; a pipeline-wide check would break the supported normalized-only input (`normalized_layer: X`), so that route needs a decision first (8.4). |
| 11.10 | lochNESS uses a fixed cutoff that depends on target frequency and has no null model. The score is local/overall − 1 (`lochness.py:398-410`): for a target making up 1e-4 of cells, one same-target neighbour out of 300 gives a score of about 32, so a fixed `enrichment_cut=0.5` over-calls rare targets. Standard and large modes report different statistics (`pct_cells_enriched` over all cells vs `pct_own_cells_enriched`), and standard runs lack the `lochness_pct_enriched` and `lochness_median` meta columns (`meta.py:194-197`). | reported | Add a permutation or binomial null per target frequency, and use the same statistic in both modes. |
| 11.11 | Minimum-cell thresholds differ by stage (5:10, 6:10, modules:20, KD filter:30, distance:30, lochNESS:10, PS:10). The meta table is the union of the target sets, so rows have structurally different NaNs. | reported | Use one target-eligibility rule. |
| 11.12 | Distance results for a target change when the target list changes. Per-target seeds use the target's position in the sorted list (`distance.py:482`; distance space uses `seed + i*43`, `:860`). There are also two seed sources (`run.seed=0`, `distance.random_seed=123`). | reported | Seed from the target name only, and from `run.seed`. |
| 11.13 | `ntc_patterns` includes `safe[-_. ]?harbor` (`config.py:423`): safe-targeting guides, which cut DNA, are pooled with true non-targeting guides. This matters for CRISPR knockout. | reported | Keep safe-harbor guides as a separate control class. |
| 11.14 | Different embeddings per stage: when `batch_key` is set, distance and distance space default to `X_pca` while lochNESS prefers `X_pca_harmony`. `get_embedding` swaps in a different representation when one is missing, but DistanceResults records the *requested* name (`data_access.py:139-153`, `distance.py:575,693,723`). | reported | Record the representation actually used, and use one representation for all stages. |
| 11.15 | Large-mode effect size in enrichment uses the global targeting composition, including the focal target, for "other" (`enrichment.py:2009-2018`). This is acknowledged only in a comment. | reported | Merging the two modes (section 2) fixes this. |
| 11.39 | **Pair-guide mode reruns the stage 5 target-knockdown test on the same cells with different rules**, and the report shows both tables (different FDR and hit call for the identical test). Differences (reformatted `pair_guide_report.py` = PGR): no 1%-of-controls expression filter (PGR:1060 vs `perturbation.py:316-326`); BH within lane × control × stratum (PGR:1125-1127) vs across all targets per control (`perturbation.py:545-546`); pair mode tests `target_symbol` from the reference (PGR:1033-1049), stage 5 only exact gene names; stage 5 caps controls at 100k, pair mode uses all; pair mode always runs both controls and hard-codes NTC as primary (PGR:1084, 1150), ignoring `primary_control`; "detected guide" is `>= min_umi` (PGR:296) vs `> detection_threshold` (`dual_guides.py:402`). | reported (2026-10-06 review) | Take the pooled rows from stage 5; keep only the per-lane/per-stratum/per-pair analyses; one FDR family and one expression filter; decide whether the `target_symbol` lookup moves into stage 5 (policy, with 11.3/11.6). |
| 11.40 | In pair mode the main `guide_representation` table lists every individual guide with 0 cells, because `OBS_GUIDE` holds "A|C" pair labels (`dual_guides.py:466, 578`; `guides.py:649-664`). The `guide_top_vs_second` figure draws the single-guide dominance line, which pair mode does not use. | reported | Skip these two outputs in pair mode (PGR has the pair versions). |

**Simulation behind 11.1** (scripts in `data_perturb_v2/scripts/distance_null_sim/`, results in `data_perturb_v2/results/distance_null_sim/`; Slurm jobs 20879016, 20879019, 20879020).

Setup: synthetic 30-dimensional embeddings. G1 isotropic Gaussian, G2 anisotropic (variance ∝ 1/PC), G3 a 5-cluster mixture with t(3) noise. Group sizes (n, m) = (30, 500), (200, 1000), (1000, 3000). TRUTH = 1e5-1e7 permutations.

- **Null mean:** c·tr(M) is exact (empirical/formula = 0.998-1.001).
- **Null variance:** the Gaussian-label formula 2c²‖M‖²_F is too large by up to 19% at n=30 (empirical/formula ratio 0.81-0.995), so the spectral method is conservative for small, isotropic targets. Possible refinement: rescale to the exact permutation variance.
- **Calibration** (fraction of null p ≤ 0.001; R = 20000 for the small sizes): spectral 0.0003-0.0010, pipeline and GPD 0.0007-0.0017, gamma 0.0049-0.0070 for G2/G3.
- **Tail accuracy** (log10 of method p / TRUTH p, for TRUTH below 1e-4): spectral median +0.06 (range −0.30 to +0.74), gamma median −2.7 to −3.6, pipeline median +2.0 (floored).
- **k does not matter:** k=50, k=200 and the full spectrum give the same p (|Δlog10| ≤ 9e-5).
- **Runtime per target, one core, at the caps (2000, 5000):** current pipeline function 96.5 s. Setup (distances plus centring) 0.8 s + eigsh k=50 2.2 s + Imhof 5 ms, about 3 s in total. Batching the existing B=1000 permutations into one matrix multiply would take 1.9 s instead of 96.5 s.
- The simulation's statistic matches `distance_test_permutation` exactly (relative difference ≤ 9e-14).

**Real-data check (Nadig HepG2 CRISPRi, `perturb/petrub_datasets/h5ad/Nadig_hepg2_raw_singlecell_01.h5ad`; results in `data_perturb_v2/results/distance_real_nadig/`; Slurm jobs 20879120, 20879122, 20879127, 20879128, 20879129).**

Setup: pipeline embedding (normalize, log1p, 3000 HVGs, scale, PCA 50); caps 2000 target cells and 5000 control cells; 4,976 NTC cells; 1,818 targets with ≥30 cells (median 45 cells).

- **Where real effects sit:** 1,036 of 1,818 targets (57%) have spectral p ≤ 1e-3, and 486 are ≤ 1e-12. The current pipeline would tie all 1,036 at 0.001. **Correction to the framing above:** BH still calls targets at the floor once more than T/50 sit there (here about 37), so the floor blocks discoveries only in screens with few strong hits. In every screen it destroys the ranking among strong hits.
- **Accuracy (final, all truth jobs finished):** 58 of 60 chosen targets have a reliable permutation truth (≥10 exceedances; B up to 1e9 on an H200; truth p from 0.99 down to 3.0e-8). log10(spectral k=50 / truth): median −0.015, range −0.107 to +0.067, i.e. spectral p is 0.78-1.17x the truth.
  - By band: p > 1e-2 median −0.003; 1e-4 to 1e-2 median −0.016; 1e-6 to 1e-4 median −0.027; 1e-9 to 1e-6 (4 targets) median −0.072. Spectral is slightly anti-conservative (about 15%) below about 1e-6.
  - Deepest targets at B=1e9: XRCC6 truth 2.1e-7 vs spectral 1.8e-7; PPP2R1A 9.2e-8 vs 1.0e-7; GNPAT 3.0e-8 vs 2.3e-8.
  - FBL and TFAM had 0 exceedances at 1e8 (truth < 1e-8); spectral gives its 1e-12 floor for both, which is consistent.
  - Pipeline B=1000 vs truth: median +0.37 in log10, up to +4.5.
  - The variance formula is accurate here (empirical/formula 0.98-1.06).
- **Eigenvalues:** the top eigenvalue is about 10% of the trace. The top 50 hold 55% of the trace but 99.4% of the sum of squares; the top 200 hold 65% and 99.8%. k=200 vs k=50 differs by at most 8.5e-6 in log10. The full spectrum vs k=50 differs by at most 1e-4 relative, apart from the 1e-12 numerical floor.
- **Runtime, per target on one core with 24 running side by side:** building the kernel 4.8 s, eigsh k=50 5.7 s, compared with the current B=1000 test at 3.7 s. Targets are small here, so permutations are cheap and spectral costs about 2.8x more. The synthetic "30x faster" applies only near the 2000-cell cap. All 1,818 targets took 13 minutes on 24 cores. Possible speed-up: the control-control block is shared across targets.
- **How the costs work: permutation vs spectral.**
  - **Distances are computed once per target**, not per permutation. Pool target and control cells (N = n + m) and build the N×N Euclidean distance matrix D on the PCA embedding, at a cost of N²·d, at most 7000² × 50 at the caps.
  - **Each permutation re-sums D under new labels.** E = 2·mean(D_XY) − mean(D_XX) − mean(D_YY), and relabelling changes which pairs count as within or between. Done naively that is N² per permutation.
  - **The pipeline's trick** (`distance_test_permutation`): the row sums of D and the grand total S_total are computed once. Then R_X = Σ_{i∈X} rowsum_i = S_XX + S_XY, so S_XY = R_X − S_XX and S_YY = S_total − 2R_X + S_XX. Only the smaller group's within block is summed fresh: **min(n,m)² + min(n,m) per permutation** (`use_x = n <= m`).
  - **So the current test's cost grows with target size:** min(n,m)² × B. The cap `distance.max_cells_per_target = 2000` (`config.py:900`; controls capped at `max_control_cells = 5000`) bounds it at 2000² = 4M entries per permutation, i.e. 4e9 for B=1000. Each permutation also copies the block with `D[np.ix_(px, px)]`, about 32 MB at n=2000.
  - Measured on one core: n≈45 (typical Nadig) 3.7 s per target; n=2000 (cap) 96.5 s per target.
  - **Spectral cost does not depend on n:** building D (N²·d) plus the top-50 eigenvalues of the N×N matrix, about 3-10 s per target at N≈5000-7000. It is therefore slower than the current B=1000 test for small targets (2.8x on Nadig) and faster near the cap (about 30x).
  - **At equal accuracy** (permutations need about 10/p), spectral wins everywhere below about 1e-4. Measured permutation truth on an H200: 1e6 permutations 1-12 s; 1e8 permutations 145-390 s; 1e9 permutations 25-33 min per target; p ≤ 1e-12 is not feasible. Spectral takes about 10 s on one CPU core for any p.
  - **Note on the truth runs:** they used the all-pairs route (E = aᵀMa, batched as one matrix multiply on the GPU, N² per permutation) instead of the n² trick. For small targets the trick is about 25,000x less arithmetic per permutation, although it means scattered memory reads, so the real speed-up on a GPU would be smaller. The truths are exact either way; only the runtime was larger than necessary.
  - Speed-up for spectral on screens with mostly small targets: the control-control block of D (about 5000×5000) is the same for every target, so compute it once.
- **What the p-values can and cannot rank:**
  - Spectral p is verified only down to about 1e-5. Imhof bottoms out at 1e-12, and 486 of 1,818 Nadig targets (27%) are tied there.
  - p combines effect size and cell count: Spearman(−log10 p, E) = 0.88 overall but only 0.54 among the 814 targets with p < 1e-6, while Spearman with n_cells is 0.22-0.30.
  - At extreme p, small departures from interchangeability (lanes, guides) dominate.
  - **Decision:** use p only for the FDR call. Rank by a **bias-corrected energy distance**: subtract the null mean c·tr(M), or use the version without self-distances. The current V-statistic has expected value c·tr(M) > 0 under the null, larger for small targets. Report it with a bootstrap CI and n_cells.
- **Dataset notes for running the pipeline on it:** `var_names` are Ensembl IDs while targets are symbols (`var["gene_name"]`); it is a dual-guide library (`sgID_AB`); it has 56 lanes (`gem_group`).

### 11b. Errors that are swallowed and automatic fallbacks

52 `except` clauses: 36 catch every exception (25 return a default or `pass`, 11 log or collect the error), and only 1 re-raises. Most are harmless checks of the environment or version strings. The ones that matter:

| # | Problem | Evidence |
|---|---|---|
| 11.16 | **`run_parallel` reruns everything serially on any error** (`compute.py:608-616`), including a worker killed for running out of memory. Out of memory therefore turns into twice the runtime and probably another out-of-memory failure in the main process. Not recorded in the status table. Six call sites: `distance.py:659` (stage 10), `:917` (stage 11), `enrichment.py:1050`, `:1809`, `perturbation.py:742`, `:1041`. Likely cause of the memory problems (inferred, needs a measurement): the stage-5 worker functions capture the full expression layer (`_eval_target_std`) or the whole AnnData (`_eval_target_large`), and the default `loky` backend uses separate processes, so joblib must write out or send that data to the workers. Stage 11 sends each pair's cells separately. **Fix:** remove the catch-all and let errors stop the run; give workers small per-target arrays (e.g. from the CSC target-gene copy, 3.11) or use threads; set `n_jobs` explicitly for memory-heavy stages. Related (2026-10-06 sweep): `limit_blas_threads` (reformatted `compute.py:477`) has `except Exception: yield`, which turns a worker's real error into "generator didn't stop after throw()". | verified |
| 11.17 | If program enrichment fails, the code only logs a warning (`modules.py`, `except Exception` around `run_program_enrichment`). Programs keep the "unannotated" placeholder, and the stage is still marked COMPLETED, so a crash looks like "no significant terms". On `feat/msigdb-gene-sets` this also covers a failed MSigDB download. | verified (reformatted); **fixed in #11** (run warning, module status, report text; `modules.program_enrichment.require` aborts instead) |
| 11.18 | A custom GMT file that failed to load was dropped with only a log warning. **On `feat/msigdb-gene-sets` the per-file catch is gone**; a failure now falls under 11.17. | verified |
| 11.19 | If Scrublet crashes, scores are set to NaN (`doublets.py:103`), but the warning says "automatic threshold failed … (scores retained)", which is not what happened. | reported |
| 11.20 | Fallbacks that change the method without being visible: enrichment large-mode control switch (no warning, `enrichment.py:1391-1396`). The distance control ladder can end with control groups smaller than `min_cells`, even 1 cell, logged at info level only (`distance.py:554-567`). A missing `stratify_by` column in distance is skipped silently (`:592`). Modules "falling back to HVGs" uses all genes when there is no `highly_variable` column. Large-mode QC reuses precomputed QC metrics from the input file (`qc.py:241-258`). | reported |
| 11.21 | Plot fallbacks: a failed heatmap ordering silently uses input order (`_order_by_similarity`), and ARI/NMI are silently left out (`plot_module_concordance`, `except Exception: pass`). See 11h. | reported |

### 11c. Tests

- 11.22 253 tests in 15 files. About 41 call sites run the full pipeline (slow).
- 11.23 **No test runs the full pipeline in large mode.** `test_scaling_consistency.py` compares some private standard/large functions on a 400-cell × 80-gene object. Never compared: QC `_filter_large`, PS large, modules `_select_genes_large`, the large plot sampling, and the large branch of `_write_unfiltered_object`. (These are where bugs 2.2, 2.5 and 2.8 live.)
- 11.24 Modules never imported by any test: `basic_qc`, `report`, `report_markdown`, `qc_plots`, `pair_guide_report`, `provenance`. `plots.py` (8.9k lines) is imported for one helper.
- 11.25 Skipped in `py313_torch211`: 6 PS tests (pertps is not installed) and 2 harmony tests. Errors: the 14 `test_basic_qc.py` tests, because `openpyxl` (declared) is not installed; `gseapy` (declared on `feat/msigdb-gene-sets`) is not installed either. An `importorskip("statsmodels")` is pointless because statsmodels is a core dependency. With numba installed, the pure-Python fallbacks in `lochness.py:77` and `guides.py:591/915` are never exercised.
- 11.26 GPU is tested only with mocks (5.3).

### 11d. Dependencies

- 11.27 Imported but not declared: `numba` (optional, with a silent fallback), `torch` and `pynvml` (GPU detection), `sklearn` (`plots.py:8762`, in the `ps` extra only).
- 11.28 Declared but unused: `leidenalg` (Leiden runs with `flavor="igraph"`, `cluster.py:638`; verified), the entire `gpu` extra (`cuml` and `rapids-singlecell` are never imported), `umap-learn` in the `ps` extra.
- 11.29 The `ps` extra points at `weili-lab/PS_python#subdirectory=pertps_project` (legacy layout with `PerturbAnalyzer`, still on lab `main`). The new pertps keeps it only as a legacy example (7.4).

### 11e. Docs vs code

- 11.30 **Documented but not implemented:** storage/backed mode (README:743, `docs/methods.md:622-638` claims stages 4-11 use `data_access` accessors, which they don't; `docs/slurm.md:52`) and GPU acceleration (README:186, 742; `docs/methods.md:587-611`; `docs/slurm.md:53`).
- 11.31 The `peak_rss_mb` profiling column is the memory in use at the end of each stage, not the peak (`compute.py:715-719`; reformatted `:615`). Any memory measurement must take this into account; `resource.getrusage` gives the true peak. Goes away with 13 A2.
- 11.32 The README numbers stages 1-11 while the code logs "Stage N/14". The README's list of optional stages leaves out `knockdown_filter`, and its repository layout (README:771) leaves out 11 modules.

### 11f. Maintainability

- 11.33 Largest functions (reformatted): `run_pipeline` 800, `pair_perturbation` 499, `pair_guide_qc` 499, `run_basic_qc` 405, `assign_guide_pairs` 331, `Config.validate` 320, `build_report` 315, `write_markdown_report` 242, `_test_cluster_enrichment_large` 236, `resolve_stage_backend` 205.
- 11.34 Duplicated helpers: 2 BH implementations (one left on `feat/msigdb-gene-sets`, 7.2); 2 identical `_expression_layer` functions plus `data_access.get_layer`, `modules._dense_layer` and 11 inline layer fallbacks; 2 identical gene-value extractors; 3 NTC detectors (`guides.py:306`, `:348`, `guide_design.py:84`); the same function names in different modules (`plot_guide_qc`, `_write_table`, `_top_two`, `attach_scores`, …).
- 11.35 94 imports inside functions (26 in `cli.py`), which hide two import cycles (io→guides, guides→dual_guides).
- 11.36 Module-level side effects: `matplotlib.use("Agg")` (`plots.py:43`, `qc_plots.py:23`), `sns.set_theme` (`plots.py:75`), and `write_h5ad` setting a process-wide anndata option (`io.py:1005`).

### 11g. Outputs and record-keeping files

- 11.37 `lognorm` is probably stored twice on disk in both modes (as X and as a layer). Guide counts are stored three ways (obsm, `_guide_barcodes.txt` on by default, optional `_guides.h5ad`). `processed_all_cells.h5ad` is on by default. gzip is hard-coded (`io.py:969`), and files over 50 MB are moved (`relocate_if_large`).
- 11.38 (→ 13 A2) Overlapping record-keeping files written from pipeline code: `logs/resolved_config.yaml`, `run_manifest.json` (which embeds the config again), `logs/module_status.json` (rewritten after every stage), `tables/module_status.csv` (written twice), `compute_profile.csv`, `figure_manifest.csv`, a column-name mapping CSV plus a copy in `uns`, and basic_qc's own `provenance.json` plus `uns['provenance']` from a second module (`provenance.py`) that overlaps `run_manifest.py`. Two of these are enough: the resolved config and the command line.

### 11h. Figure bugs (2026-10-06 plotting review; reformatted line numbers)

- 11.41 **`_plot_ps_quadrants`** (`plots.py:1420`) colours quadrants with `c=codes, cmap="tab10"` and no `vmin`/`vmax`, so a quadrant's colour depends on which quadrants the target has, and there is no quadrant legend. **verified**
- 11.42 `_plot_program_umaps` draws `RdBu_r` without centring at 0, so white is not zero.
- 11.43 `plot_distance_overview` (`:2591`) labels every bar "FDR < threshold" when the `significant` column is missing.
- 11.44 `plot_perturbation_atlas` (`:2283`) z-scores inside the plot and fills missing values with 0, showing "no data" as average.
- 11.45 Statistics that exist only inside figures: ARI/NMI only in a title (11.21), phenocopy similarity (`:1053`), PS-vs-strength recomputed (`:1338`, already at `cli.py:808`), module correlation recomputed (`:2019`), Pearson r in titles (`:830`, `:1347`), qc_plots contingency table; `plot_guide_detection_sensitivity` mixes per-sample points with across-sample mean lines.
- 11.46 Pair mode draws several figures twice (main + pair versions: guide QC, clustering by lane, volcano, per-target ECDFs); QC distributions are drawn four ways; three style systems and two class palettes coexist, so figures in one report look different.

## 12. Data layout: fixed matrix, mask columns, planned passes

Added 2026-10-06. Everything in this section is **estimate** (reasoned from the code, nothing run or timed) unless marked verified.

### 12.1 Why the matrix is written out now (verified)

A default run writes the expression matrix twice and makes extra in-memory copies to do it:

1. **All-cells checkpoint** (`write_unfiltered_h5ad`, on by default). `all_cells = expr.copy()` of every loaded cell, QC metrics added, full h5ad written before any filtering (`cli.py:1006-1022`), so that QC-failed cells are not lost.
2. **Final processed h5ad** (stage 13, `cli.py:~2440`). QC-filtered cells and genes, with `counts`, `lognorm` and `X`, all results in `obs`/`obsm`, guides merged, gzip (`io.py:1008`). In large mode `X` and `lognorm` are one object written twice (3.9).
3. **Guide h5ad** (smaller).

None of this is needed for the analysis. It follows from the scanpy pattern "filter by slicing in memory, then write the processed file": once cells are sliced away, the checkpoint is the only way to keep them. Normalized layers can be recomputed from counts and library size. The only real need to write the matrix is a one-time conversion when the input is not already a usable on-disk store (e.g. per-lane 10x outputs to concatenate).

### 12.2 What each stage needs from the matrix

From the review and function names, not re-traced line by line.

| Stage | Needs | Kind |
|---|---|---|
| 2 QC | library size, genes detected, % mito per cell | per-cell values (row-local) |
| 4 HVG | per-gene mean and variance | per-gene totals |
| 4 PCA → neighbours, UMAP, Leiden; lochNESS, distances, enrichment | cells × ~3k HVGs only | small cells × HVG matrix |
| 5 perturbation strength | each target's own gene in its cells vs controls (log2fc, KS, MWU) | own-target value per cell + target genes in (sampled) controls |
| Knockdown filter | as stage 5, plus control sums per context | as stage 5 + group × gene sums |
| 7 modules: effect matrix | mean per target × gene | group × gene sums |
| 7 modules: program scores | per-cell scores over program genes | per-cell, but programs known only after the effect matrix |
| 7 markers | per-cluster statistics | cluster × gene sums; clusters known only after stage 4 |
| 8 PS | % control expressing, DE genes in target cells | own-target values + control sample; DE genes known only after stage 5 |

Stage 5, the knockdown filter and PS compute the same things separately (own-target values, control statistics for target genes); modules and the knockdown filter both need group sums.

### 12.3 Design A: fixed matrix + mask columns

- **Raw store**: counts written once (row-chunked CSR, h5ad or zarr) and never modified. It is the all-cells object, so the checkpoint write disappears.
- **One cell table** (`obs`, every loaded cell): each step adds columns instead of slicing (`qc_pass`, `qc_fail_reason`, `guide_assignment`, `is_singlet`, `kd_status`, `kd_keep`, `cluster` with NaN for cells not embedded, scores, and a derived `in_analysis`).
- **One gene table** (`var`): `qc_gene_pass`, `highly_variable`, ...
- **Derived results**: embeddings for all cells (NaN rows where not embedded), group × gene summary tables, per-target tables.
- **Normalized values** are never stored; recomputed per chunk from counts and the library-size column.
- **Deliverable**: a small `results.h5ad` with `X=None` (AnnData allows this) holding full `obs`/`var`/`obsm`, plus the path to the raw store. Opt-in export of one self-contained h5ad at the end for anyone who needs it.
- Middle steps receive an explicit mask:

| Step | Cells / genes used |
|---|---|
| QC metrics | all cells → `qc_pass` |
| Guide assignment | guide matrix only |
| Summary pass | `qc_pass` cells, `qc_gene_pass` genes |
| Knockdown filter | `qc_pass` + targeting/control labels (on the summaries) |
| HVG → PCA → clustering, lochNESS, distances | `qc_pass` & (singlets if `assigned_only`) |
| Stage 5, modules, PS | same, minus `kd_keep == False` if the mask is applied |

### 12.4 Design B: planned passes from disk

Each summary in 12.2 is row-local, so it can be built from one sequential pass over row chunks. That is the access CSR on disk is good at, so **keep CSR**; no CSC copy is needed (this replaces fix 3.11).

1. **Pass 1, QC**: per-cell metrics → `qc_pass`. Possibly already in `obs` from Cell Ranger.
2. **Pass 2, summaries** (after QC and guide assignment; normalize each chunk on the fly):
   - per-gene mean and variance for HVG;
   - Σx, Σx², number of non-zero cells per (group, context) × gene: targets, and controls per context and split half. One sparse product `Wᵀ @ X_chunk` per chunk;
   - each targeting cell's own-target value (raw and normalized): binary search in its row;
   - target genes in a sample of controls, as a small CSC matrix, for KS/MWU;
   - the size factor needs the median library size, known only at the end: accumulate Σx·library and Σlibrary² and divide afterwards, so no second pass.
3. **Pass 3, extraction**: cells × HVG matrix → PCA and everything downstream; optionally program or DE genes.
4. **Pass 4 (optional)**: cluster × gene sums for markers; anything needing clusters or DE genes.

After pass 3 the full matrix is needed only for the optional export. Memory: per-cell vectors, the HVG matrix, and group × gene sums. The largest piece is the group × gene sums: 10k targets × 36k genes ≈ 1.4 GB per statistic in float32, reducible to HVGs ∪ targets. At 10M cells this is maybe 10-30% of one full copy, against 2-6 copies today. Sums, recombined across chunks, parallelize over chunks (Dask or multiprocessing).

Tests that need individual values (KS, MWU, Wilcoxon) are why the own-target values and a control sample are kept; Wilcoxon markers need pass 4 or a switch to a sum-based test.

### 12.5 Pros and cons vs the current setup

**Pros**

| | Current | A: masks, no rewrite | B: + streaming |
|---|---|---|---|
| Peak memory | 2-6 copies (KOLF 2.66M cells: 135-155 GB per the docs) | ~1 copy | fraction of one copy |
| Disk written per run | checkpoint (~1 copy) + final h5ad (counts + lognorm + duplicated X), gzip | tables and embeddings only | same as A |
| QC-failed cells | kept only via the checkpoint | always kept (the raw store) | same as A |
| Re-running one step | rerun the pipeline, rewrite everything | rewrite a few columns | same as A |
| Which cells a step used | implicit, visible only in `cli.py` | explicit mask column per step | same as A |
| Size ceiling | fits in memory several times | fits in memory once | disk and time, not memory |
| Repo guidelines (immutable data, metadata separate from matrix) | no | yes | yes |

**Cons**

| | A: masks, no rewrite | B: + streaming |
|---|---|---|
| Refactor size | moderate: QC and `assigned_only` stop slicing, every step takes a mask, stage 13 changes | large: every expression-reading step rewritten around summaries |
| Silently wrong results | **new risk**: a step that forgets its mask includes filtered cells and looks normal. Needs one shared row-reading function and a cell-count test per step | same, plus the summary maths must be checked against current outputs |
| Fit with scanpy | good (in-memory subsets) | partial: no Dask in neighbours, UMAP, Leiden, score_genes, so those stay on the in-memory HVG matrix; stage 5, modules, PS move to own summary code |
| What users receive | no longer one self-contained h5ad: opt-in export, README and output docs change | same |
| Dependence on the raw store | results point to it; moving, deleting or editing it orphans them. Fixed location plus an ID or checksum check | same |
| Library size definition | must be fixed once (all genes vs QC-passing genes); today it depends on what is in memory | same |
| Reading cost | normalized values recomputed on use (cheap) | runtime set by I/O on the shared filesystem, not measured |
| Small datasets | negligible cost | overhead for no gain (Nadig, 145k cells, runs in minutes today) |
| Tests | moderate update | substantial: before/after equality per consumer |

### 12.6 Knockdown filter at 10M cells (optimisation held for this redesign)

The fit scales linearly with targeting cells (about 100 ρ × group size × 30 bisection steps per group, independent across targets). The code around it does not:

| Step | Cost at 10M cells, 36k genes, 10k targets, C contexts |
|---|---|
| `_linear_target_columns`, `_raw_target_columns` | one full pass per layer, plus two float64 CSC copies of perhaps 20-40% of all stored values (~100+ GB each) |
| `counts[ctrl_rows]`, `.multiply`, split halves | control rows copied about 3-4 times |
| Per-gene loop (`targets == gene`, `contexts == ctx`, `x[ctrl]`, `.toarray()`) | a pass over all N cells per target and per context: ~10¹¹ string comparisons, ×C with contexts (~2×10¹³ for 200 lanes) |
| `row["_cells"]` | a length-N boolean mask (10 MB) kept per (target, context) row: ~100 GB pooled, terabytes with lanes |

- The last two would stop a run. Fix: group cells once by (target, context) and keep integer indices.
- CSC vs CSR: CSC speeds the column slices, is neutral for row sums, slows the control-row selection, and converting 10M × 36k on disk is a full transpose. The real fix is not slicing columns: the filter needs only each cell's own-target value and control sums, both row-local (12.4).
- Fit-side savings (not the bottleneck): the gammaln terms cancel in log f_E − log f_K and in the comparison across ρ, so the fit needs none; coarse-to-fine ρ grid; ~20 bisection steps; fit all groups at once with `np.add.reduceat` or on GPU.

### 12.7 Recommendation

- **A at any scale**: removes most memory and disk waste and makes provenance explicit. Guard against forgotten masks with one row-reading function and one cell-count test.
- **B only when datasets stop fitting in memory** (beyond KOLF-scale on a 128-200 GB node). Before committing, time one pass over a KOLF-sized store on this filesystem.
- Build A behind a single function, "rows R of the counts, normalized", that reads from memory now and streams from disk later, so consumers do not change when B arrives.
- Prerequisites: the control/context policies (11.3, 11.4, 11.6, 11.14), a single library-size definition, and the standard/large merge (section 2), otherwise every consumer is ported twice.

---

## 13. Codebase size, group A: unused machinery (potential fixes, not started)

Measured 2026-10-06 on `feat/msigdb-gene-sets` (`cb7679d`, reformatted code; line numbers refer to that branch). Package there: 21.7k lines = ~15.9k code, 2.4k docstrings, 1.2k comments, 2.2k blank. Verbatim duplication is small: pylint `duplicate-code` (≥5 lines, ignoring docstrings/comments/imports) finds only two blocks (`perturbation.py:160-183` vs `ps_score.py:220-240`; `lochness.py:119-132` vs `ps_score.py:165-181`). The remaining size is structural, so the savings below come from removing things, not deduplicating. None of A changes a statistic.

| # | What | Where (lines) | Removable | What changes for users | Tests affected |
|---|---|---|---|---|---|
| A1 | **Unused functions** (vulture ≥60% confidence, each checked by grep for callers): `cluster_composition` (`cluster.py:490`), `energy_distance_from_cdist` (`distance.py:109`), `is_pair_mode` (`dual_guides.py:174`), `_parse_target_genes_unique` (`guides.py:171`), `_gene_vector` (`perturbation.py:145`), `_obs_values` (`plots.py:258`), the `top_hits`/`top_targets`/`top_effects` result methods (`enrichment.py:124`, `lochness.py:124`, `ps_score.py:170`, `perturbation.py:114`); in `data_access.py` everything except `get_embedding`, `get_target_indices_map`, `SharedArrayBuffer` (`is_backed`, `resolve_storage_mode`, `log_storage_diagnostics`, `get_layer`, `get_expression_vector`, `get_obs_column`, `get_control_indices`; the storage layer the docs describe, 11.30); in `compute.py` `get_gpu_memory_info`, `estimate_gpu_memory_safe`, `limit_blas_threads` (only tests call them). | ~260 lines in these functions | ~260 + their tests | nothing | parts of `test_storage_data_access.py` (160) and `test_compute.py` (303) |
| A2 | **Record-keeping files** (11.38): `run_manifest.py` (305: `ModuleStatusTracker`, `build_run_manifest`, `manifest_summary_rows`), `provenance.py` (136: git/package/Slurm info, `collect`, `uns_safe`, `write_json`; `collect` is called from 7 modules to put a provenance dict in `uns`), `ComputeProfiler`/`stage_profile` (`compute.py`, ~95), ~90 bookkeeping lines in `run_pipeline` (`cli.py:348-430`, `:969-1114`), `figure_manifest.csv` writer. Keep only `logs/resolved_config.yaml` (`cli.py:331`) and the command line (log it once at start). | ~600 | Gone: `logs/run_manifest.json`, `logs/module_status.json`, `tables/module_status.csv`, `tables/compute_profile.csv`, `tables/figure_manifest.csv`, basic_qc `provenance.json`, `uns['provenance']`; the HTML report loses its module-status table and provenance rows (`report.py:172`, `cli.py:1070`). A stage that fails still stops the run with its error. | `test_compute.py`, `test_unified_pipeline.py`, `test_basic_qc.py` assertions on these files |
| A3 | **CPU/GPU decision layer** (6.2): `resolve_stage_backend` (`compute.py:236`, 205 lines), `ComputeDecision`, `log_compute_decision`, `resolve_worker_count`, `is_gpu_available`, `is_package_available`, called in 7 stages (`cluster.py:320`, distance, enrichment, lochness, modules, perturbation, ps_score). The only GPU code it switches on is one cupy branch in `modules.py:617`. Replace with `compute.n_jobs` passed to `run_parallel`, plus the modules cupy branch on an explicit flag (or drop it). | ~300 | Settings removed: `compute.backend`, `gpu_device`, `gpu_min_cells`, `gpu_min_dense_elements`, `gpu_memory_fraction`, `log_backend_decisions`, and probably the per-stage `*_n_jobs` (check). **Overlaps with the GPU plan (suggested order 5)**: if rapids-singlecell comes in, a much smaller switch (one flag) is enough. | `test_compute.py` |
| A4 | **Markdown copy of the HTML report**: `report_markdown.py` (308; `write_markdown_report` 242), called at `cli.py:1074` when `output.report_markdown_name` is set (default `report.md`). | ~310 + setting | No `report.md` (HTML report unchanged). **Open: does anyone use `report.md`?** | none found by grep |
| A5 | **Dead settings and their checks** (section 6; 6.1, 6.4-6.8) and the matching blocks in `Config.validate` (320 lines, `config.py`). | ~100-200 | Config keys removed; old configs that set them fail with the unknown-key error (`config.py` `_from_dict`). | config tests |

Total A: ~1.4k lines (A1-A5), plus the related tests. Order inside a PR: A1 and A4 first (no output change), then A2 and A3 (outputs and settings change). Watch: `run_parallel`'s serial rerun on any error (11.16) sits in the same file as A2/A3 and should be fixed at the same time.

## 14. Codebase size: other groups (reviews 2026-10-06, read-only, reformatted line numbers)

Savings are reviewers' estimates from reading; nothing implemented.

### 14.1 `pair_guide_report.py` (1,462 lines; only when `guides.assignment_mode` is `dual_guide_pair`/`pair`, `cli.py:466`)

| Function | Verdict | Removable |
|---|---|---|
| `qc_before_after` (cli:542) | Duplicates `qc_steps` (`qc.py:205-249`), `qc_summary` (`qc.py:386`) and `plot_qc` (`plots.py:331-403`), plus an extra full `expr.copy()` and QC pass in `cli.py:475-498`. Its per-lane cascade ignores `min_counts_per_cell`/`max_pct_hb` and counts genes before the gene prefilter, so its counts can disagree with `qc_steps`. | ~130 |
| `pair_guide_qc` (cli:545) | Mostly pair-specific (scaffold/slot/feature-role diagnostics, construct counts). Duplicates: `pair_assignment_status_per_lane` (also written at cli:538), status counts (= `dual_guides.pair_assignment_per_lane`), class counts (= `guides.per_lane_assignment`); `guide_gex_barcode_overlap` is circular. | ~90 |
| `single_guide_diagnostic` | Off by default (`config.py:312`). | ~25 |
| `clustering_pair_figures` (cli:624) | Lane/Leiden UMAPs and lane composition repeat `plot_clustering`; pair-status panels are pair-specific. | ~60 |
| `pair_perturbation` (cli:674) | Per-lane/per-stratum/per-pair tests and the support matrix are pair-specific; the pooled rows retest stage 5 (11.39); volcano and ECDF overview repeat plots.py. | ~80 |
| Repetition + dead branches | ECDF loop ×4, crosstab+ALL ×5, same filter ×6, `LANE_KEY` checks that never fire, unused names. | ~60-80 |

Total ~350-450. Whatever is removed must also leave the report lists (`report.py:147-153`, `report_markdown.py:155-288`).

### 14.2 Plotting (`plots.py` 2,634, `qc_plots.py` 604)

About 700 lines (~22%) replaceable: ~10 hand-made `imshow` heatmaps → `sns.heatmap`/`sns.clustermap` (`plot_perturbation_atlas` ~70, `plot_enrichment` ~55, `plot_lochness` ~60); UMAPs → `sc.pl.embedding` on a throwaway AnnData (it writes `uns[...colors]` into the object it gets); one ranked-bar helper (7 copies); grey-background highlight UMAP (6 copies); stacked % bars (5 copies); qc_plots histograms → `sns.histplot(hue=…)`. Keep: alluvial, volcanoes, gseapy dot plot. Fix the figure bugs (11h) first.

### 14.3 Size summary

| Group | Where | Lines | Changes results? |
|---|---|---|---|
| A unused machinery | 13 | ~1.4k | no (some outputs and settings go) |
| Plots | 14.2 | ~700 | figures only |
| Standard/large merge | section 2, 7.8 | ~700-1,000 | yes: fixes 2.1/2.2 |
| `pair_guide_report.py` | 14.1, 11.39 | ~400 | yes: one FDR family |
| PS → pertps | 7.4 | ~800 | yes |
| Library swaps | 7.2, 7.7, 7.9-7.15 | ~500 | mostly no |

Together ~4-5k lines, which would take the package from 21.7k to ~17k.

## 15. Guide calling: fishash and Poisson-Gaussian mixture (#17, open)

Added 2026-10-22. Branch `feat/guide-calling-fishash-pgmm` (`3d4376d`, on `195fb8d`), worktree `data_perturb_v2/worktrees/pp-guide-calling`.

**What.** Two per-guide callers: `guides.calling_method` (single-guide mode) takes `dominance` (default, unchanged), `fishash`, `poisson_gauss`; `high_moi.method` takes `threshold`, `knee`, `fishash`, `poisson_gauss`. Each decides per cell × guide whether the cell carries the guide (`guide_calling.py`). Single-guide: 1 call → assigned to the called guide (not necessarily the top-UMI guide), 2+ → `ambiguous`, 0 → `unassigned`. High-MOI: the calls are the membership matrix, so they feed the membership statistics and the non-carrier knockdown baseline (8.7). Report, QC labels and warnings, and plots describe the selected method; the dominance line is drawn only for `dominance`.

- **fishash** (Kamm, Yeung & Forrest, BMC Bioinformatics 2026): one-sided hypergeometric test of each non-zero count against the count expected from the cell's and the guide's totals (so it accounts for guide-library size); Guo-Sarkar FDR (Bonferroni within a cell, BH across cells); refit that masks called entries and re-estimates the background by rank-1 Poisson matrix completion. Numpy port of the R package; the hypergeometric tail is a port of R's `phyper` (~200× faster than `scipy.stats.hypergeom.logsf`, which also underflows). Config `fishash`: `fdr` 0.05, `min_count` 2, `refit` 10 (0 = off).
- **poisson_gauss** (crispat `fit_PGMM`, with GWT's settings): per guide, log2 of the non-zero counts as a continuous Poisson (background) + Normal (carriers) mixture; cutoff = smallest count with P(Normal) > 0.5. crispat's procedure written in numpy: AutoDelta MAP with its priors, softmax weight logits, its ten data-independent starts (copied from pyro 1.9.1 / torch 2.11), Adam lr 0.01 betas (0.8, 0.99), ≤ 500 steps, same early stopping. The fit stops before converging (the posterior has no finite maximum as the Normal weight → 0), so matching crispat needs its exact procedure, not just its model. Fitted on distinct counts weighted by frequency (exact). Config `poisson_gauss`: `min_count` 3, `min_cells` 4 (GWT: UMI floor, guides in > 3 cells).

**Verification** (`tests/test_guide_calling.py`; reference files in `tests/golden/`, rebuilt by `tests/make_fishash_reference.py` with R fishash and `tests/make_poisson_gauss_reference.py` with crispat):
- fishash: identical calls to R on 2 synthetic scenarios × refit 0/10; log p within ~1e-6.
- poisson_gauss: loss equals pyro's at crispat's fitted parameters (rtol 1e-5); gradient check; cutoffs equal native crispat's on all 16 reference guides, including seed-sensitive ones. From identical starts parameters agree to ~1e-4 (pyro is float32).
- End-to-end pipeline runs for both modes × both methods. With #16 merged: 99 passed, 2 skipped; `test_golden_file_is_current` fails identically on `main` (local package versions).
- Second opinion (OpenCode, GPT-6 Astra): found the weight parameterisation bug (one logit instead of pyro's two softmax logits, i.e. half the effective step on the weights), fixed before the PR. Report: `data_perturb_v2/copilot/second-opinion/20261015-0202-opencode-numpy-crispat-pgmm-implementation.md`.
- Speed (Slurm 20960324, one thread, per guide, 2k-90k non-zero cells): crispat 1.9-7.9 s, numpy port ~0.016 s. That benchmark ran before the two fixes; its cutoff agreement (17/24) describes the old code.

**Known properties (not bugs):**
- poisson_gauss uses raw counts and ignores the cell's guide-library size (as crispat and GWT do); fishash accounts for it.
- crispat's cutoff depends on its seeds for about 1 guide in 5; the port reproduces crispat's seeds 0-9, so it inherits that.
- Guides with background only get cutoff 2, which is why GWT's UMI floor of 3 matters.
- Remaining small numerical differences from crispat (float64 vs float32, no weight clipping, cutoff compared in log space, ∞ instead of NaN when no count crosses 0.5) can only flip a cutoff that sits on the 0.5 boundary; none did on the reference guides.

**Open:**
1. No comparison on real data: dominance vs fishash vs poisson_gauss on one lane (low-MOI) and one high-MOI dataset.
2. Both callers run on the whole guide matrix; GWT fits poisson_gauss per lane. Per-lane fitting (or running fishash per lane) is not implemented.
3. Per-guide cutoffs and fitted parameters are only logged in summary, not written as a table.
4. Other crispat methods (Gauss, binomial, Poisson/NB "SCEPTRE-like", Beta mixtures, ratio/max/quantile/UMI rules) deliberately not added: the simple rules and Beta mixtures overlap `dominance`; the Poisson/NB/binomial models add guide-library size and batch, which fishash covers except batch. If a covariate model is wanted, port the Poisson one with a batch term, or use SCEPTRE itself.
5. Noted while writing #17, not fixed: `report_markdown.py` also printed the single-guide rule in high-MOI reports (upstream; re-check against current `main`).

## Suggested order

0. **Done / in flight:** reformat (1.1, #8, merged); MSigDB gene sets + gseapy enrichment (1.5, 7.5, 7.6; #9, merged); 2.1 and 11.2 (#10, merged); 11.17 (#11, merged); knockdown non-carrier baseline (8.7, #16, merged); guide calling (15, #17, open).
1. **Correctness first:** write the test that runs both modes on synthetic data (section 2), then fix 2.2 and 11.16 (2.1 and 11.2 were fixed in #10). Fix the figure bugs in 11h (small). Settle the knockdown-filter decisions in 8.6 in a follow-up PR (PR #7 is merged; 8.3 is fixed).
2. **Decide one policy each** for controls, lanes, hit criterion and embedding (11.3, 11.4, 11.6, 11.14), then implement it everywhere, including pair mode (11.39).
3. **Merge the standard/large pairs** (section 2), so the data-layout change is made once.
4. **Data layout A** (12.3): mask columns instead of slicing, no matrix rewrite, one row-reading function. This removes most copies (3.2-3.4) and replaces the CSC copy of target genes (3.11). Port the knockdown filter first (it has tests and the Nadig benchmark to compare against), then stage 5, PS and modules.
5. **Speed and GPU:** use rapids-singlecell for stage 4, lochNESS and stage 11, and fix stage 11's repeated work and per-pair data copying.
6. **Simplify / size:** unused machinery (13; A1 and A4 can go any time, they change no output), replace PS with the new pertps (7.4), library swaps (7.2, 7.7, 7.9-7.15, scanpy pseudobulk 7.1), `pair_guide_report.py` (14.1), plots (14.2). Summary in 14.3.
7. **Data layout B** (12.4), only when datasets stop fitting in memory; time one pass over a KOLF-sized store first.
8. **Docs:** remove the claims about storage and GPU until they are implemented.
