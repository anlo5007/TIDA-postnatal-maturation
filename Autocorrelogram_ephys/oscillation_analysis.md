# `oscillation_analysis.m`

Classification of TIDA neurons as Regular or Non-regular oscillators, and
measurement of oscillation strength and frequency, from whole-cell current-clamp
recordings. One self-contained MATLAB script:
configuration, pipeline and all functions in a single file.


## What the code does

A membrane-potential trace is mean-subtracted and its normalised
autocorrelation is computed. The most prominent local maximum of the positive
half of the autocorrelogram gives two measures: its height is the
autocorrelation coefficient, used as an index of oscillation strength, and the
reciprocal of its lag is the oscillation frequency. The dominant frequency is
also estimated independently, from the amplitude spectrum of the
autocorrelogram.

Each coefficient is then tested against a null distribution built by block
bootstrap **of the same analysed window**: contiguous blocks of the trace are
resampled with replacement, which preserves the structure of the signal within
a block while destroying it across blocks. The resulting distribution is a null
for rhythmicity at timescales longer than the block. The reported p-value is
two-sided — the proportion of surrogate coefficients whose absolute value is at
least as large as the empirical one.

Everything runs in a single pass per recording:

1. import one channel of the `.abf`
2. lowpass filtering and downsampling, with the same settings for every file
3. selection of the analysed window
4. autocorrelogram, with optional manual marking of a different peak
5. amplitude spectrum of the autocorrelogram
6. block-bootstrap test

Because the coefficient and its p-value are computed on the same window, the
two always refer to the same stretch of data.

## Requirements

- MATLAB R2020a or newer
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- `abfload.m` by Harald Hentschke, which reads Axon Binary Files and is **not
  redistributed here**:
  <https://www.mathworks.com/matlabcentral/fileexchange/6190>
- EMF figure export requires Windows. On macOS and Linux the script warns and
  saves figures as PDF instead.

## How to run

1. Put `abfload.m` on the MATLAB path.
2. Open `oscillation_analysis.m` and edit the `CONFIGURATION` section.
   Nothing below that section needs to be changed.
3. Run the file. A dialog asks for the folder holding the recordings.

Before a long run, set `n_bootstrap` to about 10 and check that the windows,
figures and table are what you expect. The script prints a runtime estimate for
each recording after its first bootstrap iteration.

## Choosing the analysed window

| `window_mode` | Behaviour |
|---|---|
| `click_start` | The trace is plotted; one click sets the start of a window of `window_size_s` |
| `click_range` | The trace is plotted; two clicks set the start and the end |
| `fixed` | No plotting and no clicking; one dialog collects the start time of every recording, each window lasting `window_size_s` |

In the two click modes the user is asked first, per recording, and may decline,
in which case the whole recording is analysed.

`click_start` is the recommended default because it keeps the window length
constant across recordings. Autocorrelation coefficients and bootstrap null
distributions both depend on how much data went into them, so windows of very
different lengths are not directly comparable; the script warns when the
analysed durations differ by more than 10%.

The limits actually used are written to the results table, so any run can be
repeated exactly by switching to `fixed` and supplying those start times.

## Configuration

| Setting | Meaning |
|---|---|
| `channel` | Channel of the `.abf` to analyse |
| `search_subfolders` | Search the selected folder recursively |
| `start_dir` | Folder the selection dialog opens in |
| `lowpass_hz` | Passband in Hz, `[]` to skip filtering |
| `downsample_factor` | Keep 1 sample every n, `1` to skip |
| `window_mode` | `'click_start'`, `'click_range'` or `'fixed'` |
| `window_size_s` | Window length for `click_start` and `fixed` |
| `allow_manual_peak` | Offer manual marking of a different peak |
| `spectrum_xlim_hz` | x-axis limit of the spectrum figure only |
| `run_bootstrap` | Run the statistical test |
| `n_bootstrap` | Number of bootstrap repetitions |
| `ci_quantiles` | Lower and upper quantile of the null distribution |
| `block_size_s` | Bootstrap block length |
| `min_blocks` | Warn when the window yields fewer blocks than this |
| `random_seed` | Integer for a repeatable run, `[]` to leave the RNG untouched |
| `progress_interval` | How often progress is printed |
| `save_files` | Write results and figures to disk |
| `figure_format` | `'emf'` or `'pdf'` |

Preprocessing is applied identically to every recording rather than tuned file
by file, so it is fully described by the `parameters` sheet of the output
workbook and needs no further documentation in the methods.

## Choosing the block length

Two constraints pull in opposite directions, and the script warns when either
is violated:

- The block must be **long enough** to contain the correlation imposed by the
  lowpass filter, otherwise the surrogates are less autocorrelated than they
  should be and the test becomes too permissive. A block of at least about
  `4 / lowpass_hz` seconds is a reasonable rule of thumb.
- The block must be **short enough** for the window to be divided into enough
  blocks; fewer than `min_blocks` (default 25) gives a coarse null
  distribution.

With the defaults (2 Hz lowpass, 2 s blocks, 100 s windows) the window yields
50 blocks and the block is 4× the filter timescale. It is worth reporting a
sensitivity check at a longer block.

## The manual peak override

When `allow_manual_peak` is true, the user may mark a peak other than the
automatically detected one. The manual value is reported **alongside** the
automatic one and never replaces it, and the bootstrap always tests the
automatic coefficient.

This is deliberate. The surrogates are scored with the automatic criterion, so
testing a hand-picked empirical value against them would compare two different
statistics. Searching the surrogates for their most prominent peak, exactly as
is done for the recording, also means the null accounts for the search over
lags, which makes the test conservative in that respect.

## Output

One folder, `Autocorrelogram_results`, inside the selected folder:

- `autocorrelogram_results.xlsx` — sheet `results` with one row per recording
  (preprocessing actually applied, window limits, automatic and manual
  coefficients, lags and frequencies, spectral peak frequency, block count,
  confidence bounds and p-value) and sheet `parameters` with every setting, the
  MATLAB version and the run date
- `figures/` — per recording: preprocessing, analysed window, autocorrelogram,
  amplitude spectrum, and the bootstrap null distribution with the empirical
  value marked
- `Analysed_windows.pdf` — overview of every analysed window

The results table is also returned to the workspace as `results`.

## Reproducibility

Preprocessing is fixed in the configuration and recorded in the output.
The window depends on user clicks unless `window_mode` is `'fixed'`, and the
limits used are written to the results table, so a click-based run can be
replayed deterministically in `'fixed'` mode.

Setting `random_seed` to an integer makes the bootstrap draws repeatable.
Leaving it empty uses the default MATLAB stream, and p-values then vary between
runs on the order of `1/n_bootstrap`.

## Runtime

The bootstrap runs one autocorrelation per repetition, so cost scales with the
number of samples in the window and with `n_bootstrap`. Downsampling is the
main lever: a factor of 10 cuts the sample count tenfold and shortens the run
accordingly. The script prints its own estimate after the first iteration of
each recording, which is more reliable than any figure quoted here.

## File overview

| File | Contents |
|---|---|
| `oscillation_analysis.m` | Everything: configuration, the pipeline, and all local functions |
| `docs/oscillation_analysis.md` | This file |

Within the script, functions are grouped under banner comments: pipeline, data
import, preprocessing, window selection, analysis, figures, output, and
helpers. Each carries a standard MATLAB help block, so the documentation reads
directly from the source.
