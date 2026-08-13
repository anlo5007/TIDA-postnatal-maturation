# `functional_connectivity_analysis.m`


## What it does

For every pair of segmented ROIs within the same arcuate hemisphere, the script
computes the normalised cross-covariance of their calcium traces. Connection
strength is the peak with the largest absolute correlation coefficient, keeping
the sign at that peak so that the direction of the interaction is preserved. One
trace is then shifted by the corresponding lag and the Pearson correlation of the
aligned pair is computed as a second measure.

Both measures are tested against a null distribution obtained by resampling the
two traces **independently**. Independent resampling preserves the temporal
structure of each trace while destroying any relationship between them, so the
resulting distribution is the distribution of the measure under the hypothesis of
no functional coupling. A pair is called connected when its coefficient falls
outside the configured percentiles of that distribution.

Significant pairs are drawn as edges on a map of the ROI centroids, coloured by
coefficient, optionally over the cell map image of the slice.

> **On terminology.** These bounds are percentiles of a null distribution, that
> is critical values for a test of no cross-correlation. They are not confidence
> intervals of the estimated coefficient, although an earlier version of this
> code named them `CI`.

## Pipeline

1. import the trace `.csv` and its companion `-props.csv`
2. keep the ROIs whose status is selected, typically `accepted`
3. choose the time section to analyse
4. detrend, by 4th-order polynomial fit (after Rahmati et al., PLoS Comput Biol
   2016) or by mean subtraction
5. cross-covariance of every pair, and the autocorrelation of every ROI
6. feature extraction, lag alignment, linear correlation
7. significance against the surrogate null distribution
8. centroid distances, connectivity map, per-pair figures
9. tables of all pairs, significant pairs and a summary

## Input data

Two files in the same folder, as exported by Inscopix Data Processing Software:

| File | Layout |
|---|---|
| `<name>.csv` | row 1 ROI names, row 2 ROI status, rows 3 onwards numeric; column 1 time in seconds starting at 0, then one column per ROI |
| `<name>-props.csv` | one row per ROI, with at least `Name`, `CentroidX`, `CentroidY`, `ColorR`, `ColorG`, `ColorB` |

Optionally a cell map image in the same folder (`cfg.cellmap_suffix`, default
`*map.png`), drawn behind the connectivity map.

The trace file layout is parsed explicitly rather than left to the automatic
detection of `readtable`, whose behaviour has changed across MATLAB releases.
Sections of the recording are detected from jumps in the timestamps of column 1.

## How to run

1. Edit the `CONFIGURATION` section.
2. Run the file. A dialog asks for the calcium trace `.csv`.
3. If `cfg.roi_status` is empty, a list dialog offers the ROI statuses present in
   the file. If `cfg.time_section_mode` is `BlockAnalysis`, a second list dialog
   offers the detected sections.

Set `cfg.n_iterations` to about 10 for a first pass. The script prints a runtime
estimate after the first pair.

## Configuration

| Setting | Meaning |
|---|---|
| `start_dir` | Folder the file dialog opens in |
| `roi_status` | ROI status to keep, `""` asks |
| `props_suffix` | Suffix of the companion props file |
| `time_section_mode` | `BlockAnalysis`, `DualCursor`, `SingleCursor` or `full` |
| `time_window_s` | Window length for `SingleCursor` |
| `detrend_method` | `polynomial`, `mean` or `none` |
| `detrend_order` | Polynomial order |
| `metric` | `AbsMax` or `MaxClosestToZero` |
| `run_significance_test` | Run the surrogate test |
| `surrogate_method` | `block_bootstrap` or `circular_shift` |
| `n_iterations` | Number of surrogates per pair |
| `block_size_samples` | Block length, `[]` selects `round(n^(1/5))` |
| `null_percentiles` | Lower and upper percentile of the null distribution |
| `multiple_comparison` | `none` or `fdr` |
| `fdr_alpha` | Level for Benjamini–Hochberg |
| `random_seed` | Integer for a repeatable run, `[]` to leave the RNG alone |
| `colour_scale` | `AbsMax`, `MinMax` or `FullScale` |
| `show_cellmap`, `cellmap_suffix` | Cell map background of the map |
| `save_files`, `save_pair_figures`, `figure_format` | Output control |

### Choosing the block size

The block must be long enough to contain the autocorrelation of the calcium
signal, which is dominated by the decay of the indicator: blocks shorter than
that leave the surrogates less autocorrelated than the data, and the test becomes
too permissive. It must also be short enough for the analysed section to be
divided into enough blocks; the script warns below 25. Block length is set in
**samples**, so it must be reconverted whenever the acquisition rate or the
temporal downsampling changes.

### Multiple comparisons

The number of tests grows as the square of the number of ROIs: 20 ROIs give 190
pairs, so at a nominal 5% about 10 pairs are expected to be called significant by
chance, which inflates the proportion of connected pairs. `cfg.multiple_comparison`
set to `fdr` applies Benjamini–Hochberg to the per-pair p-values instead of the
percentile criterion. The default is `none`, which reproduces the published
analysis.

## Output

One folder, `Correlation_results`, beside the input file:

```
correlation_results.xlsx    sheets all_pairs, significant_pairs, summary, parameters
crosscorrelograms.csv       lag and the correlogram of every pair
autocorrelograms.csv        lag and the autocorrelogram of every ROI
figures/                    raw traces, detrending, correlograms,
                            connectivity map, one figure per significant pair
```

`all_pairs` holds, per pair: the two ROI names, the maximum cross-correlation
coefficient and its lag, the linear correlation after alignment, the null bounds
and p-value for both measures, the significance verdict, and the distance between
centroids in pixels. The results table is also returned to the workspace as
`results`.

Node degree (Fig 6D) is obtained by counting, per neuron, its occurrences in
`significant_pairs`.

## Runtime

Cost is the number of pairs times the number of iterations, each iteration
building two surrogates and one cross-covariance. Surrogates are drawn afresh for
every pair, as in the original analysis; precomputing one surrogate set per ROI
would be considerably faster but would make the null distributions of different
pairs share draws, so it is deliberately not done.

## Known limitations

- The linear correlation after alignment is computed on the overlapping part of
  the two traces only, so a large lag shortens the sample it is based on.
- `MaxClosestToZero` is not recommended for noisy traces: a spurious peak near
  zero lag is then selected in preference to the real one.
- The map places nodes at ROI centroids in pixel coordinates; it carries no scale
  bar, which is added at figure assembly.
