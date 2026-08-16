%% ========================================================================
%  CORRELATION_NETWORK_ANALYSIS
%  Pairwise cross-correlation and functional connectivity map of calcium
%  imaging traces exported from Inscopix (.csv).
%
%  ------------------------------------------------------------------------
%  SUMMARY
%  ------------------------------------------------------------------------
%  Traces of the accepted ROIs are detrended and, for every pair of ROIs,
%  the normalised cross-covariance is computed. A feature is extracted from
%  it, either the coefficient of largest absolute value or the peak closest
%  to zero lag. One signal is then shifted by the corresponding lag and the
%  Pearson correlation of the aligned pair is computed.
%
%  Both measures are tested against a null distribution obtained by
%  resampling the two signals independently, which preserves the temporal
%  structure of each trace while destroying any relationship between them. A
%  pair is called significant when its measure falls outside the null
%  bounds. Significant pairs are drawn as edges on a map of the ROI
%  centroids, coloured by the correlation coefficient.
%
%  NOTE ON TERMINOLOGY. The bounds are percentiles of a null distribution,
%  that is critical values for a test of no cross-correlation. They are not
%  confidence intervals of the estimated coefficient, although the earlier
%  version of this code called them CI.
%
%  ------------------------------------------------------------------------
%  PIPELINE
%  ------------------------------------------------------------------------
%    1. import the trace .csv and its companion -props.csv
%    2. keep the ROIs whose status is selected (typically "accepted")
%    3. choose the time section to analyse (cfg.time_section_mode)
%    4. detrend (4th-order polynomial, after Rahmati et al. 2016, or mean
%       subtraction)
%    5. cross-covariance of every pair, plus the autocorrelation of each ROI
%    6. feature extraction, lag alignment and linear correlation per pair
%    7. significance against a surrogate null distribution
%    8. distance between centroids, connectivity map, per-pair figures
%    9. tables of all, significant and non-significant pairs
%
%  ------------------------------------------------------------------------
%  QUICK START
%  ------------------------------------------------------------------------
%   1. Edit the CONFIGURATION section below.
%   2. Run this file. A dialog asks for the calcium trace .csv.
%
%  Before a long run, set cfg.n_iterations to about 10 to check the sections,
%  figures and tables. The script prints a runtime estimate after the first
%  pair.
%
%  ------------------------------------------------------------------------
%  REQUIREMENTS
%  ------------------------------------------------------------------------
%   MATLAB R2020a or newer.
%   Signal Processing Toolbox   - xcov, findpeaks
%   Statistics and Machine Learning Toolbox - bootstrp, prctile, pdist2
%   Image Processing Toolbox    - imread/imshow, only when cfg.show_cellmap
%                                 is true
%
%  ------------------------------------------------------------------------
%  INPUT DATA
%  ------------------------------------------------------------------------
%   Two files in the same folder, as exported by Inscopix:
%     <name>.csv        row 1 ROI names, row 2 ROI status, rows 3 and below
%                       numeric, column 1 time in s starting at 0 and one
%                       column per ROI
%     <name>-props.csv  one row per ROI, with at least the variables Name,
%                       CentroidX, CentroidY, ColorR, ColorG and ColorB
%   Optionally a cell map image in the same folder (cfg.cellmap_suffix),
%   used as the background of the connectivity map.
%
%   The layout of the trace file is read explicitly rather than left to the
%   automatic detection of readtable, whose behaviour has changed across
%   MATLAB releases.
%
%  ------------------------------------------------------------------------
%  OUTPUT
%  ------------------------------------------------------------------------
%   One folder, Correlation_results, containing
%     correlation_results.xlsx  sheets 'all_pairs', 'significant_pairs',
%                               'summary' and 'parameters'
%     crosscorrelograms.csv     lag and the correlogram of every pair
%     autocorrelograms.csv      lag and the autocorrelogram of every ROI
%     figures/                  raw traces, detrending, correlograms,
%                               connectivity map, and one figure per
%                               significant pair
%
%  ------------------------------------------------------------------------
%  REPRODUCIBILITY
%  ------------------------------------------------------------------------
%   Set cfg.random_seed to an integer to make the surrogate draws
%   repeatable. Every setting, the MATLAB version and the run date are
%   written to the 'parameters' sheet.
%
%  ------------------------------------------------------------------------
%  AUTHORS, CITATION, LICENSE
%   Author:   Andrea Locarno, Stockholm University, ORCID: 0000-0003-0640-1510
%   Contact:  andrea.locarno@dbb.su.se; andrealocarno91@gmail.com
%   Citation: https://doi.org/10.64898/2026.07.16.737908
%   License:  GPLv3
%  ========================================================================

clearvars
close all
clc


%% ========================== CONFIGURATION ==============================
% Every user-modifiable setting lives here. Nothing below this section
% needs to be edited to run the analysis.

cfg = struct();

% --- input -------------------------------------------------------------
cfg.start_dir      = 'Y:\Ca_imaging_analysis\Tida_developing';   % folder the file dialog opens in
cfg.roi_status     = "accepted";   % "" asks which status to keep, or e.g. "accepted"
cfg.props_suffix   = '-props.csv';

% --- time section -------------------------------------------------------
% 'BlockAnalysis' pick one of the sections detected in the recording
% 'DualCursor'    click the start and the end on the trace figure
% 'SingleCursor'  click the start, take cfg.time_window_s from there
% 'full'          analyse the whole recording, no dialog
cfg.time_section_mode = 'BlockAnalysis';
cfg.time_window_s     = 100;

% --- detrending ---------------------------------------------------------
cfg.detrend_method = 'polynomial'; % 'polynomial', 'mean' or 'none'
cfg.detrend_order  = 4;            % polynomial order

% --- cross-correlation feature -----------------------------------------
% "AbsMax"           coefficient of largest absolute value
% "MaxClosestToZero" peak closest to zero lag, unreliable on noisy traces
cfg.metric = "AbsMax";

% --- significance test --------------------------------------------------
cfg.run_significance_test = true;
cfg.surrogate_method      = 'block_bootstrap'; % or 'circular_shift'
cfg.n_iterations          = 1000;
cfg.block_size_samples    = 200;  % block bootstrap only, [] = round(n^(1/5))
cfg.null_percentiles      = [2.5 97.5]; % [lower upper]
cfg.multiple_comparison   = 'none';     % 'none' or 'fdr'
cfg.fdr_alpha             = 0.05;       % used when multiple_comparison is 'fdr'
cfg.random_seed           = [];         % [] = do not touch the RNG

% --- connectivity map ---------------------------------------------------
cfg.colour_scale    = "AbsMax";   % "AbsMax", "MinMax" or "FullScale"
cfg.show_cellmap    = true;       % draw the map image behind the network
cfg.cellmap_suffix  = '*map.png';

% --- output -------------------------------------------------------------
cfg.save_files        = true;
cfg.save_pair_figures = true;  % one figure per significant pair
cfg.figure_format     = 'pdf';


%% ============================== RUN ====================================

cfg = validate_config(cfg);
check_dependencies(cfg);

results = run_analysis(cfg);


%% ========================================================================
%  PIPELINE
%  ========================================================================

function results = run_analysis(cfg)
%RUN_ANALYSIS  Cross-correlation network analysis of one calcium recording.
%
%   results = RUN_ANALYSIS(cfg) asks for a calcium trace .csv, analyses
%   every pair of selected ROIs and returns a table with one row per pair.
%   Tables and figures are written to disk when cfg.save_files is true.
%
%   See also SHIFTLINEARCORR, CORRELATION_MAP.

    % --- import -------------------------------------------------------
    [file_name, folder_path] = uigetfile(fullfile(cfg.start_dir, '*.csv'), ...
                                         'Select calcium csv trace');
    if isequal(file_name, 0)
        error('correlation:noFileSelected', 'No file was selected.');
    end

    dataset = import_recording(folder_path, file_name, cfg);
    fprintf('%s: %d ROIs, %d samples at %.2f Hz\n', file_name, ...
            numel(dataset.roi_names), numel(dataset.time_s), 1/dataset.dt_s);

    keep = select_roi_by_status(dataset.roi_status, cfg.roi_status);
    dataset.traces     = dataset.traces(:, keep);
    dataset.roi_names  = dataset.roi_names(keep);
    dataset.roi_status = dataset.roi_status(keep);
    dataset.props      = dataset.props(keep, :);
    n_roi              = numel(dataset.roi_names);
    fprintf('You have selected a total of %d neurons\n', n_roi);

    if n_roi < 2
        error('correlation:tooFewRois', ...
              'At least 2 ROIs are needed, %d were selected.', n_roi);
    end

    output_dir = prepare_output_folder(folder_path, cfg);
    set_random_seed(cfg);

    % --- time section --------------------------------------------------
    traces_fig = plot_raw_traces(dataset);
    [section_traces, time_s] = time_selection(dataset.traces, dataset.time_s, ...
                                              dataset.breakpoints, dataset.section_names, ...
                                              traces_fig, 'Mode', cfg.time_section_mode, ...
                                              'TimeWindow', cfg.time_window_s);
    save_figure(traces_fig, fullfile(output_dir, 'figures', 'raw_traces'), cfg);
    close(traces_fig);
    fprintf('Analysing %.2f s (%d samples)\n', time_s(end), numel(time_s));

    % --- detrending ----------------------------------------------------
    detrended    = detrend_traces(section_traces, cfg);
    detrend_fig  = plot_detrending(time_s, section_traces, detrended, dataset.roi_names);
    save_figure(detrend_fig, fullfile(output_dir, 'figures', 'detrending'), cfg);
    close(detrend_fig);

    % --- correlograms ---------------------------------------------------
    combinations = nchoosek(1:n_roi, 2);
    n_pairs      = size(combinations, 1);
    fprintf('%d ROIs give %d pairs\n', n_roi, n_pairs);

    corr = compute_correlograms(detrended, time_s, combinations);

    save_figure(plot_correlogram_grid(corr.lag_s, corr.cross, ...
                                      pair_labels(dataset.roi_names, combinations), ...
                                      'Cross-correlograms'), ...
                fullfile(output_dir, 'figures', 'crosscorrelograms'), cfg);
    save_figure(plot_correlogram_grid(corr.lag_s, corr.auto, dataset.roi_names, ...
                                      'Autocorrelograms'), ...
                fullfile(output_dir, 'figures', 'autocorrelograms'), cfg);
    close all

    % --- pairwise analysis ----------------------------------------------
    pair_results = analyse_pairs(detrended, time_s, combinations, cfg);

    % --- significance ----------------------------------------------------
    [is_significant, q_values] = apply_multiple_comparison(pair_results, cfg);

    % --- distance between centroids --------------------------------------
    centroids = [dataset.props.CentroidX, dataset.props.CentroidY];
    distances = sqrt(sum((centroids(combinations(:,1), :) - ...
                          centroids(combinations(:,2), :)).^2, 2));

    % --- results table ----------------------------------------------------
    results = build_results_table(dataset, combinations, pair_results, ...
                                  is_significant, q_values, distances);

    n_significant = sum(results.significant_max_corr);
    fprintf('\n%d of %d pairs significant (%.1f%%)\n', ...
            n_significant, n_pairs, 100*n_significant/n_pairs);

    % --- connectivity map --------------------------------------------------
    map_fig = correlation_map(results.max_corr_coeff, combinations, ...
                              results.significant_max_corr, dataset.props, ...
                              'ScaleColorMap', cfg.colour_scale, ...
                              'CellMap', cfg.show_cellmap, ...
                              'Path', folder_path, 'Suffix', cfg.cellmap_suffix);
    save_figure(map_fig, fullfile(output_dir, 'figures', 'connectivity_map'), cfg);
    close(map_fig);

    % --- saving ------------------------------------------------------------
    if cfg.save_files
        if cfg.save_pair_figures
            save_significant_pair_figures(corr, detrended, time_s, combinations, ...
                                          results, dataset.roi_names, output_dir, cfg);
        end
        write_outputs(output_dir, results, corr, combinations, dataset, cfg);
        fprintf('Results written to "%s".\n', output_dir);
    end

    fprintf('Done.\n');
end


%% ========================================================================
%  DATA IMPORT
%  ========================================================================

function dataset = import_recording(folder_path, file_name, cfg)
%IMPORT_RECORDING  Read an Inscopix trace .csv and its companion props file.
%
%   dataset = IMPORT_RECORDING(folder_path, file_name, cfg) returns a struct
%   with fields
%     traces         n-by-m double, one column per ROI
%     time_s         n-by-1 uniform time base rebuilt from the sampling
%                    interval, starting at 0
%     raw_time_s     n-by-1 timestamps as written in the file
%     dt_s           sampling interval
%     roi_names      1-by-m string
%     roi_status     1-by-m string, e.g. "accepted" or "rejected"
%     props          m-row table read from <name>-props.csv
%     breakpoints    indices at which the timestamps jump, marking the
%                    stitching between experimental sections
%     section_names  one name per section
%     name           file name without extension
%
%   The file layout is read explicitly: row 1 holds the ROI names, row 2
%   their status, and the numeric block starts at row 3 with time in column
%   1. Relying on readtable to detect this instead makes the import depend
%   on the MATLAB release.
%
%   See also SELECT_ROI_BY_STATUS.

    file_path = fullfile(folder_path, file_name);
    raw = readcell(file_path);

    if size(raw, 1) < 3 || size(raw, 2) < 2
        error('correlation:badTraceFile', ...
              '"%s" does not have the expected layout (names, status, data).', file_path);
    end

    dataset.roi_names  = string(raw(1, 2:end));
    dataset.roi_status = string(raw(2, 2:end));

    numeric_block = raw(3:end, :);
    is_missing = cellfun(@(v) ~isnumeric(v) || isempty(v) || (isscalar(v) && ismissing(v)), ...
                         numeric_block);
    numeric_block(is_missing) = {NaN};
    numeric_block = cell2mat(numeric_block);

    dataset.raw_time_s = numeric_block(:, 1);
    dataset.traces     = numeric_block(:, 2:end);

    if numel(dataset.roi_names) ~= size(dataset.traces, 2)
        error('correlation:inconsistentHeader', ...
              'The header lists %d ROIs but the file holds %d trace columns.', ...
              numel(dataset.roi_names), size(dataset.traces, 2));
    end

    % Uniform time base. The original code used linspace(0, n*dt, n), whose
    % spacing is dt*n/(n-1); (0:n-1)*dt is the intended time base and differs
    % from it by one sample over the whole recording.
    dataset.dt_s   = dataset.raw_time_s(2) - dataset.raw_time_s(1);
    n_samples      = size(dataset.traces, 1);
    dataset.time_s = (0:n_samples-1)' * dataset.dt_s;

    % Sections are separated by a jump in the recorded timestamps
    dataset.breakpoints = find(round(diff(dataset.raw_time_s), 3) ~= round(dataset.dt_s, 3));
    n_sections = numel(dataset.breakpoints) + 1;
    dataset.section_names = arrayfun(@(k) sprintf('Section %d', k), (1:n_sections)', ...
                                     'UniformOutput', false);
    for i_section = 1:numel(dataset.breakpoints)
        fprintf('Your section %d ends at %.2f s\n', ...
                i_section, dataset.time_s(dataset.breakpoints(i_section)));
    end

    [~, base_name] = fileparts(file_name);
    dataset.name  = base_name;
    dataset.props = import_props(folder_path, base_name, size(dataset.traces, 2), cfg);
end


function props = import_props(folder_path, base_name, n_roi, cfg)
%IMPORT_PROPS  Read the companion props file and check the columns needed.
    props_path = fullfile(folder_path, [base_name cfg.props_suffix]);
    if ~isfile(props_path)
        error('correlation:noPropsFile', ...
              'The props file "%s" was not found.', props_path);
    end

    props = readtable(props_path);

    required = {'Name', 'CentroidX', 'CentroidY', 'ColorR', 'ColorG', 'ColorB'};
    missing_columns = required(~ismember(required, props.Properties.VariableNames));
    if ~isempty(missing_columns)
        error('correlation:badPropsFile', ...
              '"%s" is missing the column(s): %s.', props_path, strjoin(missing_columns, ', '));
    end

    if height(props) ~= n_roi
        error('correlation:propsMismatch', ...
              'The props file describes %d ROIs but the trace file holds %d.', ...
              height(props), n_roi);
    end
end


function keep = select_roi_by_status(roi_status, requested_status)
%SELECT_ROI_BY_STATUS  Logical mask of the ROIs whose status is kept.
%
%   With an empty requested_status the available statuses are offered in a
%   list dialog; cancelling keeps every ROI.
    available = unique(roi_status, 'stable');

    if strlength(requested_status) > 0
        keep = ismember(roi_status, requested_status);
        if ~any(keep)
            error('correlation:noRoiWithStatus', ...
                  'No ROI has status "%s". Available: %s.', ...
                  requested_status, strjoin(cellstr(available), ', '));
        end
        return
    end

    [selection, was_chosen] = listdlg('PromptString', 'Select the neuron you want to analyze based on the tag', ...
                                      'ListString', cellstr(available), 'ListSize', [160 160]);

    if ~was_chosen
        warning('correlation:noStatusSelected', ...
                'No tag has been selected, all the neurons will proceed to the analysis pipeline.');
        keep = true(1, numel(roi_status));
        return
    end

    keep = ismember(roi_status, available(selection));
end


%% ========================================================================
%  TIME SECTION SELECTION
%  ========================================================================

function [selected_data, time_updated] = time_selection(data, mytime, breakpoints, section_names, myfigure, options)
%TIME_SELECTION  Select a portion of time-series data.
%
%   [selected_data, time_updated] = TIME_SELECTION(data, mytime,
%   breakpoints, section_names, myfigure) returns the part of data to
%   analyse and its time base restarted at 0, according to options.Mode:
%
%     "BlockAnalysis" pick one of the sections delimited by breakpoints
%     "DualCursor"    click the start and the end on myfigure
%     "SingleCursor"  click the start, take options.TimeWindow from there
%     "full"          the whole recording, without any dialog
%
%   data is n-by-m with one column per signal, mytime is n-by-1.
%   selected_data is a numeric matrix, not a cell array as in the earlier
%   version of this function.
%
%   See also DETREND_TRACES.
    arguments
        data                (:,:) {mustBeNumeric}
        mytime              (:,1) {mustBeNumeric}
        breakpoints         (:,1) {mustBeNumeric}
        section_names       (:,1) cell
        myfigure
        options.Mode        {mustBeMember(options.Mode, ...
                            ["BlockAnalysis", "DualCursor", "SingleCursor", "full"])} = "BlockAnalysis"
        options.TimeWindow  (1,1) {mustBeNumeric, mustBePositive} = 100
    end

    section_limits = [0; breakpoints; numel(mytime)];

    switch options.Mode
        case "full"
            index_start = 1;
            index_end   = numel(mytime);

        case "BlockAnalysis"
            [i_section, was_chosen] = listdlg('PromptString', 'Select the Time section', ...
                                              'ListString', section_names, ...
                                              'ListSize', [160 160], 'SelectionMode', 'single');
            if ~was_chosen
                error('correlation:noSectionSelected', 'No time-section selected.');
            end
            index_start = section_limits(i_section) + 1;
            index_end   = section_limits(i_section + 1);

        case "DualCursor"
            axes_handles = findall(myfigure.Children, 'Type', 'axes');
            first_click  = read_click(mytime, 'start', axes_handles);
            second_click = read_click(mytime, 'end',   axes_handles);
            clicked      = sort([first_click, second_click]); % tolerate reverse order
            waitfor(msgbox(sprintf('You selected %.2f s', diff(clicked))));

            [~, index_start] = min(abs(mytime - clicked(1)));
            [~, index_end]   = min(abs(mytime - clicked(2)));

        case "SingleCursor"
            axes_handles = findall(myfigure.Children, 'Type', 'axes');
            first_click  = read_click(mytime, 'start', axes_handles);

            [~, index_start] = min(abs(mytime - first_click));
            window_samples   = round(options.TimeWindow / (mytime(2) - mytime(1)));
            index_end        = index_start + window_samples - 1;

            if index_end > numel(mytime)
                warning('correlation:windowTruncated', ...
                        'Selected end point is outside the time limit, using trace end instead.');
                index_end = numel(mytime);
            end
            draw_marker(axes_handles, mytime(index_end));
    end

    if index_end - index_start < 1
        error('correlation:emptySection', 'The selected section contains fewer than 2 samples.');
    end

    selected_data = data(index_start:index_end, :);
    time_updated  = mytime(index_start:index_end) - mytime(index_start);
end


function time_point = read_click(mytime, which_end, axes_handles)
%READ_CLICK  Read one click, clamp it to the trace and mark it on every axes.
    point = ginput(1);
    if isempty(point)
        error('correlation:noClick', 'No point was selected on the trace.');
    end
    time_point = point(1,1);

    if time_point < min(mytime)
        waitfor(warndlg(sprintf('The %s point is outside the trace limits. Using trace start.', which_end)));
        time_point = min(mytime);
    elseif time_point > max(mytime)
        waitfor(warndlg(sprintf('The %s point is outside the trace limits. Using trace end.', which_end)));
        time_point = max(mytime);
    end

    draw_marker(axes_handles, time_point);
end


function draw_marker(axes_handles, time_point)
%DRAW_MARKER  Vertical line at a time point, on every axes of the figure.
    for i_axes = 1:numel(axes_handles)
        xline(axes_handles(i_axes), time_point, 'k--', 'LineWidth', 1.5);
    end
end


%% ========================================================================
%  PREPROCESSING
%  ========================================================================

function detrended = detrend_traces(traces, cfg)
%DETREND_TRACES  Remove the slow trend of every trace.
%
%   detrended = DETREND_TRACES(traces, cfg) subtracts, column by column,
%   either a polynomial fit of order cfg.detrend_order (after Rahmati et
%   al., PLoS Comput Biol 2016) or the mean of the trace. Removing the trend
%   matters here because a shared drift between two ROIs would otherwise
%   inflate their cross-correlation at every lag.
%
%   See also DETREND.
    arguments
        traces (:,:) {mustBeNumeric, mustBeNonempty}
        cfg    (1,1) struct
    end

    switch cfg.detrend_method
        case 'polynomial'
            detrended = detrend(traces, cfg.detrend_order);
        case 'mean'
            detrended = traces - mean(traces, 1);
        case 'none'
            detrended = traces;
    end
end


%% ========================================================================
%  CROSS-CORRELATION
%  ========================================================================

function corr = compute_correlograms(traces, time_s, combinations)
%COMPUTE_CORRELOGRAMS  Cross- and autocorrelograms of every ROI pair.
%
%   corr = COMPUTE_CORRELOGRAMS(traces, time_s, combinations) returns a
%   struct with fields
%     lag_s  (2n-1)-by-1 lag axis in s
%     cross  (2n-1)-by-p normalised cross-covariance, one column per pair
%     auto   (2n-1)-by-m normalised autocovariance, one column per ROI
%
%   XCOV is used rather than XCORR so that the plotted and exported
%   correlograms are computed exactly like the coefficient that is tested:
%   FEATURE_EXTRACTION also works on the cross-covariance, and with XCORR
%   the marked peak could sit slightly off the plotted curve.
%
%   See also XCOV, FEATURE_EXTRACTION.
    arguments
        traces       (:,:) {mustBeNumeric, mustBeNonempty}
        time_s       (:,1) {mustBeNumeric}
        combinations (:,2) {mustBeNumeric}
    end

    corr.lag_s = lag_axis(time_s);
    corr.cross = zeros(numel(corr.lag_s), size(combinations,1));
    corr.auto  = zeros(numel(corr.lag_s), size(traces,2));

    for i_pair = 1:size(combinations,1)
        corr.cross(:, i_pair) = xcov(traces(:, combinations(i_pair,1)), ...
                                     traces(:, combinations(i_pair,2)), 'normalized');
    end

    for i_roi = 1:size(traces,2)
        corr.auto(:, i_roi) = xcov(traces(:, i_roi), 'normalized');
    end
end


function lag_s = lag_axis(time_s)
%LAG_AXIS  Symmetric lag axis matching the output of XCOV, in seconds.
    time_from_zero = time_s - time_s(1);
    lag_s = sort([-time_from_zero; time_from_zero(2:end)]);
end


function [maxcorrcoeff, lag_max_s, shift_samples] = feature_extraction(x1, x2, lag, options)
%FEATURE_EXTRACTION  Extract a feature from the cross-covariance of two signals.
%
%   [maxcorrcoeff, lag_max_s, shift_samples] = FEATURE_EXTRACTION(x1, x2,
%   lag) computes the normalised cross-covariance of x1 and x2 and returns
%   the coefficient of largest absolute value, the lag at which it occurs
%   and the corresponding shift in samples relative to zero lag.
%
%   Name-value argument
%     'MyMetric'  "AbsMax" (default) takes the largest absolute
%                 coefficient. "MaxClosestToZero" takes the peak closest to
%                 zero lag, which is not recommended when the
%                 cross-covariance is noisy, since a spurious peak near zero
%                 is then selected in preference to the real one.
%
%   lag_max_s is expressed in the units of lag, that is seconds, not in
%   samples as the earlier help text stated.
%
%   See also XCOV, FINDPEAKS, SHIFTLINEARCORR.
    arguments
        x1  (:,1) {mustBeNumeric}
        x2  (:,1) {mustBeNumeric}
        lag (:,1) {mustBeNumeric}
        options.MyMetric {mustBeMember(options.MyMetric, ["AbsMax", "MaxClosestToZero"])} = "AbsMax"
    end

    cross_covariance = xcov(x1, x2, 'normalized');

    if numel(cross_covariance) ~= numel(lag)
        error('correlation:lagMismatch', ...
              'The lag vector (%d) does not match the cross-covariance (%d).', ...
              numel(lag), numel(cross_covariance));
    end

    switch options.MyMetric
        case "AbsMax"
            [~, i_selected] = max(abs(cross_covariance));

        case "MaxClosestToZero"
            [~, peak_lags] = findpeaks(abs(cross_covariance), lag);
            if isempty(peak_lags)
                error('correlation:noPeak', ...
                      'No peak was found in the cross-covariance.');
            end
            [~, i_closest] = min(abs(peak_lags));
            [~, i_selected] = min(abs(lag - peak_lags(i_closest)));
    end

    maxcorrcoeff  = cross_covariance(i_selected);
    lag_max_s     = lag(i_selected);
    shift_samples = i_selected - numel(x1); % zero lag sits at index numel(x1)
end


function pair = shiftlinearcorr(x1, x2, t, opts)
%SHIFTLINEARCORR  Cross-correlation feature, lag alignment and linear correlation.
%
%   pair = SHIFTLINEARCORR(x1, x2, t) extracts a feature from the
%   cross-covariance of x1 and x2 (see FEATURE_EXTRACTION), shifts one
%   signal by the corresponding lag, and computes the Pearson correlation of
%   the overlapping part of the aligned pair.
%
%   With opts.RunTest true, both measures are tested against a null
%   distribution: each surrogate is built by resampling x1 and x2
%   independently, which keeps the temporal structure of each trace and
%   destroys any relationship between them. A measure is significant when it
%   falls outside opts.NullPercentiles of that distribution. The p-value is
%   the proportion of surrogates at least as extreme in absolute value.
%
%   These bounds are critical values of a null distribution, not confidence
%   intervals of the estimate.
%
%   pair is a struct with fields
%     max_corr_coeff, lag_s, shift_samples
%     linear_corr_coeff
%     null_lower_max, null_upper_max, p_max, outside_max
%     null_lower_linear, null_upper_linear, p_linear, outside_linear
%   where outside_* is +1 above the upper bound, -1 below the lower bound
%   and 0 inside. All test fields are NaN when opts.RunTest is false.
%
%   See also FEATURE_EXTRACTION, BLOCKBOOTSTRP, SCRAMBLING.
    arguments
        x1 (:,1) {mustBeNumeric}
        x2 (:,1) {mustBeNumeric}
        t  (:,1) {mustBeNumeric}
        opts.Metric          {mustBeMember(opts.Metric, ["AbsMax", "MaxClosestToZero"])} = "AbsMax"
        opts.RunTest         (1,1) logical = false
        opts.NullPercentiles (1,2) {mustBeNumeric} = [2.5 97.5]
        opts.SurrogateMethod {mustBeMember(opts.SurrogateMethod, ...
                             ["block_bootstrap", "circular_shift"])} = "block_bootstrap"
        opts.BlockSize       (1,1) {mustBeInteger, mustBePositive} = 200
        opts.Iterations      (1,1) {mustBeInteger, mustBePositive} = 1000
    end

    if numel(x1) ~= numel(x2)
        error('correlation:unequalLength', 'x1 and x2 must have the same length.');
    end

    lag = lag_axis(t);

    [pair.max_corr_coeff, pair.lag_s, pair.shift_samples] = ...
        feature_extraction(x1, x2, lag, 'MyMetric', opts.Metric);
    pair.linear_corr_coeff = aligned_linear_correlation(x1, x2, pair.shift_samples);

    pair.null_lower_max    = NaN;   pair.null_upper_max    = NaN;
    pair.p_max             = NaN;   pair.outside_max       = NaN;
    pair.null_lower_linear = NaN;   pair.null_upper_linear = NaN;
    pair.p_linear          = NaN;   pair.outside_linear    = NaN;

    if ~opts.RunTest
        return % the earlier version left these outputs unassigned here
    end

    null_max    = nan(opts.Iterations, 1);
    null_linear = nan(opts.Iterations, 1);

    for i_iteration = 1:opts.Iterations
        surrogate_1 = make_surrogate(x1, opts);
        surrogate_2 = make_surrogate(x2, opts);

        [null_max(i_iteration), ~, surrogate_shift] = ...
            feature_extraction(surrogate_1, surrogate_2, lag, 'MyMetric', opts.Metric);
        null_linear(i_iteration) = ...
            aligned_linear_correlation(surrogate_1, surrogate_2, surrogate_shift);
    end

    [pair.null_lower_max, pair.null_upper_max, pair.p_max, pair.outside_max] = ...
        summarise_null(null_max, pair.max_corr_coeff, opts.NullPercentiles);
    [pair.null_lower_linear, pair.null_upper_linear, pair.p_linear, pair.outside_linear] = ...
        summarise_null(null_linear, pair.linear_corr_coeff, opts.NullPercentiles);
end


function linear_coeff = aligned_linear_correlation(x1, x2, shift_samples)
%ALIGNED_LINEAR_CORRELATION  Pearson correlation of two signals aligned at a lag.
%
%   The two signals are shifted against each other by shift_samples and the
%   correlation is computed on the overlapping part only.
    if shift_samples < 0
        x1_aligned = x1(1:end + shift_samples);
        x2_aligned = x2(-shift_samples + 1:end);
    elseif shift_samples > 0
        x1_aligned = x1(shift_samples + 1:end);
        x2_aligned = x2(1:end - shift_samples);
    else
        x1_aligned = x1;
        x2_aligned = x2;
    end

    if numel(x1_aligned) < 3
        linear_coeff = NaN; % the lag consumed almost the whole window
        return
    end

    coefficient_matrix = corrcoef(x1_aligned, x2_aligned);
    linear_coeff = coefficient_matrix(1,2);
end


function [lower_bound, upper_bound, p_value, outside] = summarise_null(null_values, observed, percentiles)
%SUMMARISE_NULL  Bounds, two-sided p-value and verdict for one statistic.
    lower_bound = prctile(null_values, min(percentiles));
    upper_bound = prctile(null_values, max(percentiles));
    p_value     = sum(abs(null_values) >= abs(observed)) / numel(null_values);

    if observed > upper_bound
        outside = 1;
    elseif observed < lower_bound
        outside = -1;
    else
        outside = 0;
    end
end


function surrogate = make_surrogate(x, opts)
%MAKE_SURROGATE  One surrogate trace, by block bootstrap or circular shift.
    switch opts.SurrogateMethod
        case "block_bootstrap"
            surrogate = blockbootstrp(x, opts.BlockSize, 1);
        case "circular_shift"
            surrogate = scrambling(x);
    end
end


function [btstrppd, numblocks] = blockbootstrp(data, blocksize, iterations, verbose)
%BLOCKBOOTSTRP  Block bootstrap resampling of time series data.
%
%   btstrppd = BLOCKBOOTSTRP(data, blocksize, iterations) divides the input
%   into contiguous blocks of blocksize samples and resamples at the block
%   level rather than at the sample level, which preserves the structure of
%   the signal within a block while destroying it across blocks.
%
%   If the length of data is not divisible by blocksize, a random subset of
%   samples equal to the remainder is removed before blocking and each
%   surrogate is linearly interpolated back onto the original sample
%   positions, so that all surrogates keep the input length.
%
%   btstrppd is length(data)-by-iterations. numblocks is the number of
%   blocks actually used.
%
%   This function is duplicated, identically, in autocorrelogram_analysis.m
%   so that each script stays self-contained. Keep the two copies in step.
%
%   See also BOOTSTRP, SCRAMBLING.
    arguments
        data       (:,1) {mustBeNumeric, mustBeNonempty}
        blocksize  (1,1) {mustBeInteger, mustBePositive}
        iterations (1,1) {mustBeInteger, mustBePositive}
        verbose    (1,1) logical = false
    end

    n_samples = numel(data);
    assert(blocksize < n_samples, 'blockbootstrp:blockTooLong', ...
           'The block size (%d) must be shorter than the data (%d samples).', ...
           blocksize, n_samples);

    sample_positions = (1:n_samples)';
    remainder        = rem(n_samples, blocksize);
    was_trimmed      = remainder ~= 0;

    if was_trimmed
        if verbose
            warning('blockbootstrp:notDivisible', ...
                    ['The length of the data (%d) is not divisible by the block size (%d). ' ...
                     '%d random samples (%.2f%%) are removed before blocking and the ' ...
                     'surrogates are linearly interpolated to replace them.'], ...
                    n_samples, blocksize, remainder, 100*remainder/n_samples);
        end
        discarded      = randperm(n_samples, remainder);
        kept_data      = data;              kept_data(discarded)      = [];
        kept_positions = sample_positions;  kept_positions(discarded) = [];
    else
        kept_data      = data;
        kept_positions = sample_positions;
    end

    numblocks = numel(kept_data) / blocksize; % integer by construction
    if verbose
        fprintf('Your number of blocks is %d\n', numblocks);
    end
    blocks = reshape(kept_data, [blocksize, numblocks])'; % one block per row

    resampled = bootstrp(iterations, @(x) x', blocks)';   % one surrogate per column

    if was_trimmed
        btstrppd = zeros(n_samples, iterations);
        for i_iteration = 1:iterations
            interpolant = griddedInterpolant(kept_positions, resampled(:, i_iteration));
            btstrppd(:, i_iteration) = interpolant(sample_positions);
        end
    else
        btstrppd = resampled;
    end
end


function scrambled = scrambling(x)
%SCRAMBLING  Circular shift of each column by a random amount.
%
%   scrambled = SCRAMBLING(x) rotates every column of x by an independent
%   random offset. Unlike a block bootstrap this preserves the
%   autocorrelation of the signal exactly, at the cost of introducing one
%   discontinuity per column at the wrap point. Used as an alternative
%   surrogate when cfg.surrogate_method is 'circular_shift'.
%
%   See also BLOCKBOOTSTRP.
    arguments
        x {mustBeNumeric}
    end

    scrambled = zeros(size(x));
    for i_column = 1:size(x, 2)
        offset = randi(size(x,1));
        scrambled(:, i_column) = [x(offset:end, i_column); x(1:offset-1, i_column)];
    end
end


function pair_results = analyse_pairs(traces, time_s, combinations, cfg)
%ANALYSE_PAIRS  Run SHIFTLINEARCORR on every pair, with progress reporting.
%
%   Surrogates are drawn afresh for each pair, as in the original analysis.
%   Precomputing one surrogate set per ROI and reusing it across pairs would
%   be much faster but would make the null distributions of different pairs
%   share draws, so it is deliberately not done here.
    n_pairs      = size(combinations, 1);
    pair_results = cell(n_pairs, 1);
    pair_timer   = tic;

    block_size = resolve_block_size(cfg, size(traces,1));

    for i_pair = 1:n_pairs
        pair_results{i_pair} = shiftlinearcorr( ...
            traces(:, combinations(i_pair,1)), ...
            traces(:, combinations(i_pair,2)), ...
            time_s, ...
            'Metric',          cfg.metric, ...
            'RunTest',         cfg.run_significance_test, ...
            'NullPercentiles', cfg.null_percentiles, ...
            'SurrogateMethod', string(cfg.surrogate_method), ...
            'BlockSize',       block_size, ...
            'Iterations',      cfg.n_iterations);

        if i_pair == 1 && cfg.run_significance_test
            fprintf('  about %.1f min for %d pairs\n', ...
                    toc(pair_timer) * n_pairs / 60, n_pairs);
        elseif mod(i_pair, 10) == 0
            fprintf('  pair %d of %d\n', i_pair, n_pairs);
        end
    end

    pair_results = vertcat(pair_results{:});
end


function block_size = resolve_block_size(cfg, n_samples)
%RESOLVE_BLOCK_SIZE  Block length in samples, with a warning on the block count.
%
%   An empty cfg.block_size_samples selects round(n^(1/5)). Note that the
%   earlier default was written round(n^1/5), which MATLAB evaluates as n/5:
%   that always yields exactly 5 blocks, whatever the length of the
%   recording. The analysis reported in the manuscript passed the block size
%   explicitly and was not affected.
    if isempty(cfg.block_size_samples)
        block_size = max(2, round(n_samples^(1/5)));
    else
        block_size = cfg.block_size_samples;
    end

    if block_size >= n_samples
        error('correlation:blockTooLong', ...
              'The block size (%d) must be shorter than the analysed window (%d samples).', ...
              block_size, n_samples);
    end

    n_blocks = floor(n_samples / block_size);
    fprintf('Block bootstrap: %d samples per block, %d blocks\n', block_size, n_blocks);
    if n_blocks < 25
        warning('correlation:tooFewBlocks', ...
                ['The window yields only %d blocks. The null distribution will be coarse; ' ...
                 'use a longer section or a shorter block.'], n_blocks);
    end
end


function [is_significant, q_values] = apply_multiple_comparison(pair_results, cfg)
%APPLY_MULTIPLE_COMPARISON  Significance verdict for the cross-correlation feature.
%
%   With cfg.multiple_comparison 'none' a pair is significant when its
%   coefficient falls outside the null bounds, which is the criterion used
%   in the original analysis. With 'fdr' the p-values are corrected by the
%   Benjamini-Hochberg procedure at cfg.fdr_alpha instead.
%
%   Correction matters here because the number of tests grows as the square
%   of the number of ROIs: 20 ROIs give 190 pairs, so at a nominal 5% about
%   10 pairs are expected to come out significant by chance alone.
    outside  = [pair_results.outside_max]';
    p_values = [pair_results.p_max]';
    q_values = nan(size(p_values));

    switch cfg.multiple_comparison
        case 'none'
            is_significant = outside ~= 0;

        case 'fdr'
            [is_significant, q_values] = benjamini_hochberg(p_values, cfg.fdr_alpha);
            fprintf('Benjamini-Hochberg at alpha = %.3f: %d of %d pairs significant\n', ...
                    cfg.fdr_alpha, sum(is_significant), numel(is_significant));
    end

    is_significant(isnan(outside)) = false; % the test did not run
end


function [is_significant, q_values] = benjamini_hochberg(p_values, alpha)
%BENJAMINI_HOCHBERG  Control of the false discovery rate.
%
%   Returns the significance verdict and the adjusted p-values. Note that a
%   bootstrap p-value cannot be smaller than 1/iterations, which limits how
%   small the adjusted values can become.
    n_tests = numel(p_values);
    [sorted_p, order] = sort(p_values(:));

    adjusted = sorted_p .* n_tests ./ (1:n_tests)';
    adjusted = min(1, flipud(cummin(flipud(adjusted)))); % enforce monotonicity

    q_values        = nan(n_tests, 1);
    q_values(order) = adjusted;
    is_significant  = q_values <= alpha;
end


%% ========================================================================
%  FIGURES
%  ========================================================================

function fig = plot_raw_traces(dataset)
%PLOT_RAW_TRACES  One tile per ROI, with the section boundaries marked.
    fig = figure('Name', 'Raw traces', 'WindowState', 'maximized', 'Color', 'w');
    layout = tiledlayout(fig, size(dataset.traces,2), 1, 'TileSpacing', 'tight');

    section_limits = [0; dataset.time_s(dataset.breakpoints); dataset.time_s(end)];

    for i_roi = 1:size(dataset.traces, 2)
        ax = nexttile(layout);
        plot(ax, dataset.time_s, dataset.traces(:, i_roi));
        set(ax, 'box', 'off')

        if ~isempty(dataset.breakpoints)
            xline(ax, dataset.time_s(dataset.breakpoints), 'r--', 'LineWidth', 1.5);
        end

        if i_roi == 1
            y_top = max(get(ax, 'YLim'));
            for i_section = 1:numel(section_limits)-1
                text(ax, mean(section_limits(i_section:i_section+1)), y_top, ...
                     dataset.section_names{i_section}, 'FontWeight', 'bold');
            end
        end
    end

    xlabel(layout, 'Time (s)')
    ylabel(layout, '\DeltaF/F')
end


function fig = plot_detrending(time_s, original, detrended, roi_names)
%PLOT_DETRENDING  Original trace, detrended trace and removed trend per ROI.
    fig = figure('Name', 'Detrending', 'WindowState', 'maximized', 'Color', 'w');
    layout = tiledlayout(fig, 'vertical', 'TileSpacing', 'tight');

    for i_roi = 1:size(original, 2)
        ax = nexttile(layout);
        plot(ax, time_s, original(:, i_roi));
        hold(ax, 'on')
        plot(ax, time_s, detrended(:, i_roi));
        plot(ax, time_s, original(:, i_roi) - detrended(:, i_roi), 'k--', 'LineWidth', 2);
        hold(ax, 'off')
        set(ax, 'box', 'off')
        ylabel(ax, roi_names(i_roi), 'Interpreter', 'none');
    end

    xlabel(layout, 'Time (s)')
    ylabel(layout, '\DeltaF/F')
    legend_handle = legend('Original', 'Detrended', 'Trend');
    legend_handle.Layout.Tile = 'east';
end


function fig = plot_correlogram_grid(lag_s, correlograms, labels, figure_title)
%PLOT_CORRELOGRAM_GRID  One tile per correlogram, shared lag axis.
    fig = figure('Name', figure_title, 'WindowState', 'maximized', 'Color', 'w');
    layout = tiledlayout(fig, 'vertical', 'TileSpacing', 'tight');

    for i_column = 1:size(correlograms, 2)
        ax = nexttile(layout);
        plot(ax, lag_s, correlograms(:, i_column));
        set(ax, 'box', 'off')
        axis(ax, 'tight')
        ylabel(ax, labels(i_column), 'Interpreter', 'none');
    end

    title(layout,  figure_title)
    xlabel(layout, 'Lag (s)')
    ylabel(layout, 'Correlation Coefficient')
end


function [fig, c_bar] = correlation_map(corr_max, combinations, is_significant, props, opts)
%CORRELATION_MAP  Connectivity map of the significant pairs.
%
%   fig = CORRELATION_MAP(corr_max, combinations, is_significant, props)
%   draws one node per ROI at its centroid and one edge per significant
%   pair, coloured by its correlation coefficient, optionally over the cell
%   map image.
%
%   Name-value arguments
%     'ScaleColorMap'  how the colours span the data:
%                        "AbsMax"    symmetric around 0, spanning the
%                                    largest absolute coefficient (default)
%                        "MinMax"    from the smallest to the largest
%                                    coefficient present
%                        "FullScale" the fixed range [-1 1]
%     'CMap'           256-by-3 colormap, blue to white to red by default
%     'CellMap'        true to draw the map image behind the network
%     'Path'           folder searched for that image
%     'Suffix'         pattern of the image file, default '*map.png'
%
%   See also RUN_ANALYSIS.
    arguments
        corr_max       (:,1) {mustBeNumeric}
        combinations   (:,2) {mustBeNumeric}
        is_significant (:,1) {mustBeNumericOrLogical}
        props          table
        opts.ScaleColorMap {mustBeMember(opts.ScaleColorMap, ["AbsMax", "MinMax", "FullScale"])} = "AbsMax"
        opts.CMap (256,3) {mustBeNumeric} = default_colormap()
        opts.CellMap (1,1) {mustBeNumericOrLogical} = false
        opts.Path {mustBeTextScalar} = ''
        opts.Suffix {mustBeTextScalar} = '*map.png'
    end

    is_significant = logical(is_significant);

    fig = figure('WindowState', 'maximized');
    if opts.CellMap
        draw_cellmap_background(opts.Path, opts.Suffix);
    end
    hold on
    set(gca, 'DataAspectRatio', [1 1 1])

    % nodes, in the colour assigned to the ROI by the segmentation
    for i_roi = 1:height(props)
        plot(props.CentroidX(i_roi), props.CentroidY(i_roi), 'o', ...
             'LineWidth', 2, 'MarkerSize', 10, ...
             'Color', [props.ColorR(i_roi), props.ColorG(i_roi), props.ColorB(i_roi)]/255);
    end

    ax = gca;
    ax.YDir = 'reverse'; % match image coordinates
    xlabel('X position (pixel)')
    ylabel('Y position (pixel)')
    title('Functional connectivity map')
    colormap(opts.CMap)

    significant_values = corr_max(is_significant);
    [colour_index, colour_limits] = map_to_colormap(significant_values, ...
                                                    opts.ScaleColorMap, size(opts.CMap,1));

    % edges, one per significant pair
    significant_pairs = find(is_significant);
    for i_edge = 1:numel(significant_pairs)
        i_pair = significant_pairs(i_edge);
        plot([props.CentroidX(combinations(i_pair,1)), props.CentroidX(combinations(i_pair,2))], ...
             [props.CentroidY(combinations(i_pair,1)), props.CentroidY(combinations(i_pair,2))], ...
             'Color', opts.CMap(colour_index(i_edge), :), 'LineWidth', 1.5);
    end
    hold off

    c_bar = colorbar(ax);
    c_bar.Label.String = 'Correlation Coefficient Peak';
    if ~isempty(colour_limits)
        c_bar.Ticks      = linspace(0, 1, 5);
        c_bar.TickLabels = round(linspace(colour_limits(1), colour_limits(2), 5), 3);
    end
end


function [colour_index, colour_limits] = map_to_colormap(values, scale_mode, n_colours)
%MAP_TO_COLORMAP  Index of each value in the colormap, and the labelled range.
    colour_index  = [];
    colour_limits = [];
    if isempty(values)
        return
    end

    switch scale_mode
        case "AbsMax"
            span = max(abs(values));
            if span == 0
                span = 1; % every coefficient is zero, avoid a degenerate range
            end
            colour_limits = [-span, span];

        case "MinMax"
            colour_limits = [min(values), max(values)];
            if diff(colour_limits) == 0
                colour_limits = colour_limits + [-1 1]*max(eps, abs(colour_limits(1))*0.01);
            end

        case "FullScale"
            colour_limits = [-1, 1];
    end

    colour_index = round(rescale(values, 1, n_colours, ...
                                 'InputMin', colour_limits(1), ...
                                 'InputMax', colour_limits(2)));
    colour_index = min(max(colour_index, 1), n_colours);
end


function cmap = default_colormap()
%DEFAULT_COLORMAP  Blue to white to red, white at the midpoint.
    cmap = [[linspace(0,1,128)'; ones(128,1)], ...
            [linspace(0,1,128)'; linspace(1,0,128)'], ...
            [ones(128,1); linspace(1,0,128)']];
end


function draw_cellmap_background(folder_path, suffix)
%DRAW_CELLMAP_BACKGROUND  Show the cell map image behind the network.
    candidates = dir(fullfile(folder_path, suffix));

    if isempty(candidates)
        warning('correlation:noCellMap', ...
                'No image matching "%s" was found in "%s".', suffix, folder_path);
        return
    elseif numel(candidates) > 1
        chosen = uigetfile(fullfile(folder_path, '*.png'), 'Select the cell map', ...
                           'MultiSelect', 'off');
        if isequal(chosen, 0)
            return
        end
        image_data = imread(fullfile(folder_path, chosen));
    else
        image_data = imread(fullfile(candidates.folder, candidates.name));
    end

    if size(image_data, 3) == 1
        image_data = repmat(image_data, [1 1 3]); % grayscale to RGB
    end
    imshow(image_data);
end


function save_significant_pair_figures(corr, traces, time_s, combinations, results, roi_names, output_dir, cfg)
%SAVE_SIGNIFICANT_PAIR_FIGURES  Correlogram and both traces, per significant pair.
    pair_dir = fullfile(output_dir, 'figures', 'correlograms');
    if ~isfolder(pair_dir)
        mkdir(pair_dir);
    end

    significant_pairs = find(results.significant_max_corr);
    for i_edge = 1:numel(significant_pairs)
        i_pair = significant_pairs(i_edge);
        roi_a  = combinations(i_pair, 1);
        roi_b  = combinations(i_pair, 2);

        fig    = figure('Visible', 'off', 'Color', 'w');
        layout = tiledlayout(fig, 3, 1, 'TileSpacing', 'tight');

        ax = nexttile(layout);
        plot(ax, corr.lag_s, corr.cross(:, i_pair));
        hold(ax, 'on')
        plot(ax, results.lag_s(i_pair), results.max_corr_coeff(i_pair), 'r*');
        if ~isnan(results.null_lower_max(i_pair))
            yline(ax, results.null_lower_max(i_pair), 'r--');
            yline(ax, results.null_upper_max(i_pair), 'r--', 'Null bounds');
        end
        hold(ax, 'off')
        ylabel(ax, 'Cross-correlation coefficient')
        xlabel(ax, 'Lag (s)')

        ax = nexttile(layout);
        plot(ax, time_s, traces(:, roi_a));
        ylabel(ax, '\DeltaF/F')
        title(ax, roi_names(roi_a), 'Interpreter', 'none')

        ax = nexttile(layout);
        plot(ax, time_s, traces(:, roi_b));
        ylabel(ax, '\DeltaF/F')
        title(ax, roi_names(roi_b), 'Interpreter', 'none')

        title(layout, sprintf('%s Vs. %s', roi_names(roi_a), roi_names(roi_b)), ...
              'Interpreter', 'none')
        xlabel(layout, 'Time (s)')

        file_stem = matlab.lang.makeValidName(sprintf('%s_Vs_%s', roi_names(roi_a), roi_names(roi_b)));
        save_figure(fig, fullfile(pair_dir, file_stem), cfg);
        close(fig);
    end
end


function save_figure(fig, base_path, cfg)
%SAVE_FIGURE  Export a figure as vector graphics in the configured format.
    if ~cfg.save_files || isempty(fig) || ~isgraphics(fig)
        return
    end
    exportgraphics(fig, [base_path '.' cfg.figure_format], 'ContentType', 'vector');
end


%% ========================================================================
%  OUTPUT
%  ========================================================================

function results = build_results_table(dataset, combinations, pair_results, is_significant, q_values, distances)
%BUILD_RESULTS_TABLE  One row per pair, with every measure and its verdict.
    results = table( ...
        dataset.roi_names(combinations(:,1))', ...
        dataset.roi_names(combinations(:,2))', ...
        [pair_results.max_corr_coeff]', ...
        [pair_results.lag_s]', ...
        [pair_results.linear_corr_coeff]', ...
        [pair_results.null_lower_max]', ...
        [pair_results.null_upper_max]', ...
        [pair_results.p_max]', ...
        q_values, ...
        is_significant, ...
        [pair_results.null_lower_linear]', ...
        [pair_results.null_upper_linear]', ...
        [pair_results.p_linear]', ...
        [pair_results.outside_linear]' ~= 0, ...
        distances, ...
        'VariableNames', {'roi_a', 'roi_b', 'max_corr_coeff', 'lag_s', ...
                          'linear_corr_coeff', 'null_lower_max', 'null_upper_max', ...
                          'p_max', 'q_max', 'significant_max_corr', ...
                          'null_lower_linear', 'null_upper_linear', 'p_linear', ...
                          'significant_linear_corr', 'distance_pixels'});
end


function write_outputs(output_dir, results, corr, combinations, dataset, cfg)
%WRITE_OUTPUTS  Workbook of results plus the correlograms as csv.
    workbook = fullfile(output_dir, 'correlation_results.xlsx');

    significant = results(results.significant_max_corr, :);
    summary = table(height(results), height(significant), ...
                    height(significant)/height(results), ...
                    numel(dataset.roi_names), dataset.time_s(end), ...
                    'VariableNames', {'n_pairs', 'n_significant', 'ratio_significant', ...
                                      'n_roi', 'recording_duration_s'});

    writetable(results,     workbook, 'Sheet', 'all_pairs',         'WriteMode', 'overwritesheet');
    writetable(significant, workbook, 'Sheet', 'significant_pairs', 'WriteMode', 'overwritesheet');
    writetable(summary,     workbook, 'Sheet', 'summary',           'WriteMode', 'overwritesheet');
    writetable(config_to_table(cfg), workbook, 'Sheet', 'parameters', 'WriteMode', 'overwritesheet');

    % correlograms are long and narrow, csv suits them better than a sheet
    cross_table = array2table([corr.lag_s, corr.cross], 'VariableNames', ...
                              ['lag_s', pair_labels(dataset.roi_names, combinations)]);
    writetable(cross_table, fullfile(output_dir, 'crosscorrelograms.csv'));

    auto_table = array2table([corr.lag_s, corr.auto], 'VariableNames', ...
                             ['lag_s', matlab.lang.makeValidName(cellstr(dataset.roi_names))]);
    writetable(auto_table, fullfile(output_dir, 'autocorrelograms.csv'));
end


function labels = pair_labels(roi_names, combinations)
%PAIR_LABELS  Valid column names of the form roiA_Vs_roiB.
    labels = matlab.lang.makeValidName( ...
        cellstr(roi_names(combinations(:,1))' + "_Vs_" + roi_names(combinations(:,2))'))';
end


function output_dir = prepare_output_folder(folder_path, cfg)
%PREPARE_OUTPUT_FOLDER  Create the results folder and its figures subfolder.
    output_dir = fullfile(folder_path, 'Correlation_results');
    if ~cfg.save_files
        return
    end

    figure_dir = fullfile(output_dir, 'figures');
    if ~isfolder(figure_dir)
        mkdir(figure_dir);
    end
end


function parameter_table = config_to_table(cfg)
%CONFIG_TO_TABLE  Flatten the configuration struct into a two-column table.
    cfg.matlab_version = version;
    cfg.run_date       = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    cfg.script_name    = mfilename;

    names  = fieldnames(cfg);
    values = cellfun(@(name) value_to_text(cfg.(name)), names, 'UniformOutput', false);

    parameter_table = table(names, values, 'VariableNames', {'parameter', 'value'});
end


function text_value = value_to_text(value)
%VALUE_TO_TEXT  Represent a configuration value as a single string.
    if ischar(value) || isstring(value)
        text_value = char(value);
    elseif islogical(value)
        text_value = mat2str(value);
    elseif isempty(value)
        text_value = '[]';
    elseif isnumeric(value)
        text_value = strtrim(sprintf('%g ', value));
    else
        text_value = class(value);
    end
end


%% ========================================================================
%  CONFIGURATION HELPERS
%  ========================================================================

function cfg = validate_config(cfg)
%VALIDATE_CONFIG  Check the settings and fail early with a clear message.
    cfg.time_section_mode = validatestring(cfg.time_section_mode, ...
        {'BlockAnalysis', 'DualCursor', 'SingleCursor', 'full'}, mfilename, 'cfg.time_section_mode');
    cfg.detrend_method = validatestring(cfg.detrend_method, ...
        {'polynomial', 'mean', 'none'}, mfilename, 'cfg.detrend_method');
    cfg.surrogate_method = validatestring(cfg.surrogate_method, ...
        {'block_bootstrap', 'circular_shift'}, mfilename, 'cfg.surrogate_method');
    cfg.multiple_comparison = validatestring(cfg.multiple_comparison, ...
        {'none', 'fdr'}, mfilename, 'cfg.multiple_comparison');
    cfg.figure_format = validatestring(cfg.figure_format, {'pdf', 'emf'}, ...
        mfilename, 'cfg.figure_format');

    cfg.metric       = validate_member(cfg.metric, ["AbsMax", "MaxClosestToZero"], 'cfg.metric');
    cfg.colour_scale = validate_member(cfg.colour_scale, ["AbsMax", "MinMax", "FullScale"], 'cfg.colour_scale');

    mustBeInteger(cfg.detrend_order);  mustBePositive(cfg.detrend_order);
    mustBeInteger(cfg.n_iterations);   mustBePositive(cfg.n_iterations);
    mustBePositive(cfg.time_window_s);

    if ~isempty(cfg.block_size_samples)
        mustBeInteger(cfg.block_size_samples);
        mustBePositive(cfg.block_size_samples);
    end

    if numel(cfg.null_percentiles) ~= 2 || any(cfg.null_percentiles < 0) || ...
            any(cfg.null_percentiles > 100)
        error('correlation:badConfig', ...
              'cfg.null_percentiles must be two values between 0 and 100.');
    end
    cfg.null_percentiles = sort(cfg.null_percentiles); % [lower upper]

    if cfg.fdr_alpha <= 0 || cfg.fdr_alpha >= 1
        error('correlation:badConfig', 'cfg.fdr_alpha must lie strictly between 0 and 1.');
    end

    if strcmp(cfg.multiple_comparison, 'fdr') && ~cfg.run_significance_test
        error('correlation:badConfig', ...
              'cfg.multiple_comparison is ''fdr'' but cfg.run_significance_test is false.');
    end

    if strcmp(cfg.figure_format, 'emf') && ~ispc
        warning('correlation:emfNotSupported', ...
                'EMF export is only available on Windows. Figures will be saved as PDF.');
        cfg.figure_format = 'pdf';
    end
end


function value = validate_member(value, allowed, name)
%VALIDATE_MEMBER  Check a string setting against the allowed values.
    value = string(value);
    if ~ismember(value, allowed)
        error('correlation:badConfig', '%s must be one of: %s.', name, strjoin(allowed, ', '));
    end
end


function check_dependencies(cfg)
%CHECK_DEPENDENCIES  Verify that the required toolbox functions are available.
    required = {'xcov', 'Signal Processing Toolbox'; ...
                'pdist2', 'Statistics and Machine Learning Toolbox'};

    if cfg.run_significance_test && strcmp(cfg.surrogate_method, 'block_bootstrap')
        required = [required; {'bootstrp', 'Statistics and Machine Learning Toolbox'}];
    end
    if cfg.run_significance_test
        required = [required; {'prctile', 'Statistics and Machine Learning Toolbox'}];
    end
    if cfg.show_cellmap
        required = [required; {'imread', 'Image Processing Toolbox'}];
    end
    if strcmp(cfg.metric, "MaxClosestToZero")
        required = [required; {'findpeaks', 'Signal Processing Toolbox'}];
    end

    for i_function = 1:size(required, 1)
        if ~exist(required{i_function,1}, 'file')
            error('correlation:missingToolbox', ...
                  '"%s" was not found. This script requires the %s.', ...
                  required{i_function,1}, required{i_function,2});
        end
    end
end


function set_random_seed(cfg)
%SET_RANDOM_SEED  Fix the random stream when a seed is configured.
    if ~cfg.run_significance_test
        return
    end

    if isempty(cfg.random_seed)
        fprintf(['Random seed not set: the surrogate draws differ between runs, so\n' ...
                 'p-values vary on the order of 1/n_iterations. Set cfg.random_seed\n' ...
                 'for an exactly repeatable run.\n']);
    else
        rng(cfg.random_seed);
        fprintf('Random seed set to %d.\n', cfg.random_seed);
    end
end
