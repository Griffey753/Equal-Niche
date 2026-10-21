# Equal-Niche

MATLAB implementation of **Equal-Niche: Niche-Uniform First-Parent Sampling
for Multimodal Multiobjective Optimization**, with processed comparison,
ablation, sampling-audit and sensitivity results.

This is a curated publication package, not a copy of the research workspace.
The MATLAB algorithm files are byte-identical to the current experimental
implementation. No third-party framework, competitor implementation or
benchmark/reference-set files are redistributed here.

## Requirements and quick start

Tested environment: MATLAB R2026a and PlatEMO v4.9. The implementation uses
`pdist2` from the Statistics and Machine Learning Toolbox. Install PlatEMO
separately from its [official repository](https://github.com/BIMK/PlatEMO)
and select release v4.9 rather than an unpinned latest version.

After adding the PlatEMO v4.9 directory and its subdirectories to MATLAB's
path, run from this repository:

```matlab
addpath(fullfile(pwd,'matlab'));
addpath(fullfile(pwd,'tests'));
test_equal_niche_release;
run_equal_niche_example;
```

The example uses PlatEMO's built-in MMF1 solely as a runnable demonstration.
It is not a replacement for the paper's shared benchmark wrapper and common
initial populations. It does not overwrite published experiment results.

## Source map

| File | Responsibility |
| --- | --- |
| `matlab/EqualNicheAdapter.m` | Main loop, population selection and exact-deduplicated archive update |
| `matlab/+equalniche/temporaryNiches.m` | Deterministic exclusive temporary-niche construction |
| `matlab/+equalniche/reproduce.m` | Parent sampling, SBX and polynomial mutation |
| `matlab/+equalniche/fullNNSL.m` | Shared local-dominance ranking used by population and archive |
| `matlab/+equalniche/selectByFronts.m` | Complete-front priority and critical-front truncation |
| `matlab/+equalniche/normalizedDistance.m` | Bound-normalized decision-space distance |

The contribution under study is reproductive allocation: a member of niche
`i` has Parent 1 probability `1/(K*n_i)`. Parent 2 is population-uniform.
The individual-uniform ablation and interpolated sampler are exposed through
`parent1Mode`; the default is `EQUAL_NICHE`. Other optional controls are the
shared neighborhood radius, center order and observation-only sampling trace.

The local-dominance ranking and associated survival components are not claimed
as new contributions. See Q. Deng et al., *A niching-based nondominated sorting
for multimodal multiobjective optimization with local Pareto fronts*, Applied
Soft Computing 179 (2025), 113223,
[doi:10.1016/j.asoc.2025.113223](https://doi.org/10.1016/j.asoc.2025.113223).
For the platform, see Y. Tian et al., *PlatEMO: A MATLAB Platform for
Evolutionary Multi-Objective Optimization*, IEEE Computational Intelligence
Magazine 12(4), 73–87 (2017),
[doi:10.1109/MCI.2017.2742868](https://doi.org/10.1109/MCI.2017.2742868).

## Processed data and experimental scope

- `data/comparison`: 18 problems x 5 methods x 30 seeds, per-run indicators,
  problem summaries, average ranks and statistical comparisons.
- `data/ablation`: 18 problems x 2 Parent 1 samplers x 30 paired seeds.
- `data/auxiliary`: interpolation strength, shared-radius and center-order
  comparisons, plus generation-level multinomial sampling calibration.
- `protocols`: preserved machine-readable experiment settings. Paths identify
  locations in the original research project, not files supplied in this package.
- `file_manifest.json`: file sizes and SHA-256 digests.

Metrics are dimensionless. Lower IGDX and IGDF and higher PSP are preferred.
Problem, method/variant/condition and seed columns identify observations.
Columns labelled STANDARD describe the final-population evaluation reported
in the paper; NATIVE columns, where present in preserved data, describe archive
sensitivity outputs and must not be substituted into the main comparison.
Hashes identify the original run and initialization artifacts. `artifact_path`
values are project-relative provenance references; they do not imply that the
raw artifact is included. Local machine paths are removed without changing
numerical cells.

This package supplies algorithm source and processed data, **not yet the full
raw-data and five-solver reproduction bundle**. Raw MATLAB run artifacts,
shared initial populations, the common problem wrapper, benchmark reference
sets, comparator adapters and the Python analysis/plotting pipeline are not
included. Exact regeneration of all paper tables and figures requires those
additional materials. Third-party distribution rights must be settled before
adding benchmark or competitor sources. All optimization runs underlying these
processed results use MATLAB; orchestration, metrics, statistics and plotting
in the research pipeline also use Python.

## Rights and release status

No open-source license is assigned by this staging package. Public availability
is not a blanket grant to reuse code or third-party methods; the authors must
confirm an appropriate software/data license before a formal versioned release.
PlatEMO, benchmark definitions and comparator sources retain their respective
owners' terms. No repository publication or release DOI is asserted here.

Funding: Joint Funds of the National Natural Science Foundation of China,
Grant No. U24A20267.
