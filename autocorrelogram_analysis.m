%% ========================================================================
%  AUTOCORRELOGRAM_ANALYSIS
%  Quantification of membrane-potential oscillations from whole-cell
%  current-clamp recordings stored as Axon Binary Files (.abf).
%
%  ------------------------------------------------------------------------
%  SUMMARY
%  ------------------------------------------------------------------------
%  A membrane-potential trace is mean-subtracted and its normalised
%  autocorrelation is computed. The most prominent local maximum of the
%  positive half of the autocorrelogram provides two measures: its height is
%  the autocorrelation coefficient, used as an index of oscillation
%  strength, and the reciprocal of its lag is the oscillation frequency. The
%  dominant frequency is also estimated independently from the amplitude
%  spectrum of the autocorrelogram.
%
%  Each coefficient is then tested against a null distribution built by
%  block bootstrap of the same analysed window: contiguous blocks of the
%  trace are resampled with replacement, which preserves the structure of
%  the signal within a block while destroying it across blocks. The
%  resulting distribution is a null for rhythmicity at timescales longer
%  than the block. The reported p-value is two-sided, the proportion of
%  surrogate coefficients whose absolute value is at least as large as the
%  empirical one.
%
%  ------------------------------------------------------------------------
%  PIPELINE
%  ------------------------------------------------------------------------
%  One pass over the recordings, one step per recording:
%
%    1. import a single channel of the .abf
%    2. optional lowpass filtering and downsampling, applied to the whole
%       recording with the same settings for every file (cfg.lowpass_hz,
%       cfg.downsample_factor)
%    3. selection of the analysed window (cfg.window_mode), by clicking on
%       the plotted trace or from typed start times
%    4. autocorrelogram, with optional manual override of the detected peak
%    5. amplitude spectrum of the autocorrelogram
%    6. optional block-bootstrap test (cfg.run_bootstrap)
%
%  Every measure is computed on the same analysed window, so the reported
%  coefficient and its p-value always refer to the same stretch of data.
%
%  ------------------------------------------------------------------------
%  QUICK START
%  ------------------------------------------------------------------------
%   1. Put abfload.m (see REQUIREMENTS) on the MATLAB path.
%   2. Edit the CONFIGURATION section below.
%   3. Run this file. A dialog asks for the folder holding the recordings.
%
%  Before a long run, set cfg.n_bootstrap to about 10 and check that the
%  windows, figures and tables are what you expect. The script prints a
%  runtime estimate for each recording after its first bootstrap iteration.
%
%  ------------------------------------------------------------------------
%  REQUIREMENTS
%  ------------------------------------------------------------------------
%   MATLAB R2020a or newer. The limiting functions are exportgraphics and
%   the WriteMode option of writetable, both introduced in R2020a;
%   argument-validation blocks require R2019b.
%   Signal Processing Toolbox   - lowpass, downsample, findpeaks, xcorr
%   Statistics and Machine Learning Toolbox - bootstrp, quantile
%   abfload.m by Harald Hentschke, third-party, not redistributed here:
%     https://www.mathworks.com/matlabcentral/fileexchange/6190
%   EMF figure export requires Windows; on other platforms the script warns
%   and saves PDF instead.
%
%  ------------------------------------------------------------------------
%  INPUT DATA
%  ------------------------------------------------------------------------
%   A folder of .abf files, searched recursively when cfg.search_subfolders
%   is true. Only one channel is analysed (cfg.channel). Traces are handled
%   as n-by-2 matrices: column 1 time in s starting at 0, column 2 membrane
%   potential in mV.
%
%  ------------------------------------------------------------------------
%  OUTPUT
%  ------------------------------------------------------------------------
%   A single folder, Autocorrelogram_results, containing
%     autocorrelogram_results.xlsx  sheet 'results', one row per recording,
%                                   and sheet 'parameters' with every
%                                   setting, the MATLAB version and the run
%                                   date, so the file documents its own
%                                   provenance
%     figures/                      per recording: preprocessing, analysed
%                                   window, autocorrelogram, amplitude
%                                   spectrum and, when the bootstrap runs,
%                                   the null distribution
%     Analysed_windows.pdf          overview of every analysed window
%   The results table is also returned to the workspace as "results".
%
%  ------------------------------------------------------------------------
%  REPRODUCIBILITY
%  ------------------------------------------------------------------------
%   Preprocessing is set once in the configuration and applied identically
%   to every recording, so it is fully described by the parameters sheet.
%   The analysed window depends on user clicks unless cfg.window_mode is
%   'fixed'; the limits actually used are written to the results table, so a
%   run can be repeated exactly by switching to 'fixed' and supplying those
%   limits. Set cfg.random_seed to an integer to make the bootstrap draws
%   repeatable; leaving it empty uses the default MATLAB stream.
%
%  ------------------------------------------------------------------------
%  AUTHORS, CITATION, LICENSE
%  ------------------------------------------------------------------------
%   Author:   Andrea Locarno, Stockholm University, ORCID: 0000-0003-0640-1510
%   Contact:  andrea.locarno@dbb.su.se; andrealocarno91@gmail.com
%   Citation: https://doi.org/10.64898/2026.07.16.737908
%   License:  CC-BY-4.0. abfload.m is distributed under its
%             own license and is not covered by it.
%  ========================================================================

clearvars
close all
clc


%% ========================== CONFIGURATION ==============================
% Every user-modifiable setting lives here. Nothing below this section
% needs to be edited to run the analysis.

cfg = struct();

% --- data --------------------------------------------------------------
cfg.channel           = 1;    % channel of the .abf to analyse
cfg.search_subfolders = true; % search the selected folder recursively
cfg.start_dir         = '';   % folder the selection dialog opens in

% --- preprocessing, applied identically to every recording -------------
cfg.lowpass_hz         = [];  % passband in Hz, [] to skip the filtering
cfg.downsample_factor  = 1; % keep 1 sample every n, 1 to skip

% --- analysed window ---------------------------------------------------
% 'click_start'  plot the trace, click once, take cfg.window_size_s from there
% 'click_range'  plot the trace, click the start and the end
% 'fixed'        no clicking, type the start time of each recording
% In the two click modes the user is asked first, per recording, and can
% decline; declining analyses the whole recording.
cfg.window_mode   = 'click_start';
cfg.window_size_s = 100; % window length for 'click_start' and 'fixed'

% --- autocorrelogram ---------------------------------------------------
cfg.allow_manual_peak = false; % offer manual override of the detected peak
cfg.spectrum_xlim_hz  = 10;    % x-axis limit of the spectrum figure only

% --- block bootstrap ---------------------------------------------------
cfg.run_bootstrap     = true;
cfg.n_bootstrap       = 1000;          % repetitions
cfg.ci_quantiles      = [0.025 0.975]; % lower and upper quantile, 95% CI
cfg.block_size_s      = 2;             % block length
cfg.min_blocks        = 25;            % warn below this many blocks
cfg.random_seed       = [];            % [] = do not touch the RNG
cfg.progress_interval = 100;           % print progress every n repetitions

% --- output ------------------------------------------------------------
cfg.save_files    = true;   % false to run without writing to disk
cfg.figure_format = 'pdf';  % 'emf' (Windows only) or 'pdf'


%% ============================== RUN ====================================

cfg = validate_config(cfg);
check_dependencies(cfg);

results = run_analysis(cfg);


%% ========================================================================
%  PIPELINE
%  ========================================================================

function results = run_analysis(cfg)
%RUN_ANALYSIS  Analyse every recording of a folder in a single pass.
%
%   results = RUN_ANALYSIS(cfg) asks for a folder of .abf files, processes
%   each of them in turn (import, preprocessing, window selection,
%   autocorrelogram, amplitude spectrum and optional block-bootstrap test)
%   and returns a table with one row per recording. Results, figures and the
%   settings that produced them are written to disk when cfg.save_files is
%   true.
%
%   Recordings are handled one at a time so that memory use does not grow
%   with the size of the folder.
%
%   See also MP_AUTOCORRELOGRAM, BOOTSTRAP_AUTOCORR.

    folder_path = uigetdir(cfg.start_dir, 'Select the folder with the .abf files');
    if isequal(folder_path, 0)
        error('autocorrelogram:noFolderSelected', 'No folder was selected.');
    end

    abf_list = list_recordings(folder_path, cfg.search_subfolders);
    n_files  = numel(abf_list);
    fprintf('Found %d recording(s) in "%s".\n', n_files, folder_path);

    % 'fixed' needs the start times before anything else, so that the whole
    % run can then proceed without further input
    fixed_start_s = nan(n_files, 1);
    if strcmp(cfg.window_mode, 'fixed')
        fixed_start_s = ask_fixed_start_times(abf_list);
    end

    output_dir = prepare_output_folder(folder_path, cfg);
    set_random_seed(cfg);

    records          = cell(n_files, 1);
    analysed_windows = cell(n_files, 1);

    for i_file = 1:n_files
        file_path = fullfile(abf_list(i_file).folder, abf_list(i_file).name);
        fprintf('\n[%d/%d] %s\n', i_file, n_files, abf_list(i_file).name);

        % --- import and preprocessing -----------------------------------
        raw_trace = load_recording(file_path, cfg.channel);
        prep      = preprocess_trace(raw_trace, cfg);

        % --- analysed window --------------------------------------------
        window = select_window(prep.trace, cfg, abf_list(i_file).name, fixed_start_s(i_file));
        analysed_windows{i_file} = window.trace;
        fprintf('  window %.2f - %.2f s (%.2f s, %d samples)\n', ...
                window.start_time_s, window.end_time_s, window.duration_s, size(window.trace,1));

        % --- autocorrelogram and spectrum --------------------------------
        autocorr = mp_autocorrelogram(window.trace, 'Interactive', cfg.allow_manual_peak);
        spectrum = autocorrelogram_freq_fft(autocorr.correlogram, ...
                                            'FrequencyLimitHz', cfg.spectrum_xlim_hz);

        % --- block bootstrap ---------------------------------------------
        if cfg.run_bootstrap
            boot = bootstrap_autocorr(window.trace, autocorr.coeff_automatic, cfg);
        else
            boot = empty_bootstrap_result();
        end

        records{i_file} = build_record(abf_list(i_file), prep, window, autocorr, spectrum, boot);

        % --- figures ------------------------------------------------------
        if cfg.save_files
            recording_name = erase(abf_list(i_file).name, '.abf');
            figure_dir     = fullfile(output_dir, 'figures');
            save_figure(prep.figure,     fullfile(figure_dir, [recording_name '_preprocessing']),  cfg);
            save_figure(window.figure,   fullfile(figure_dir, [recording_name '_analysed_window']), cfg);
            save_figure(autocorr.figure, fullfile(figure_dir, [recording_name '_correlogram']),     cfg);
            save_figure(spectrum.figure, fullfile(figure_dir, [recording_name '_spectrum']),        cfg);
            save_figure(boot.figure,     fullfile(figure_dir, [recording_name '_bootstrap_null']),  cfg);
        end
        close all
    end

    results = struct2table(vertcat(records{:}));
    warn_on_unequal_windows(results);

    if cfg.save_files
        overview_fig = plot_window_overview(analysed_windows, {abf_list.name});
        exportgraphics(overview_fig, fullfile(output_dir, 'Analysed_windows.pdf'), ...
                       'ContentType', 'vector');
        close(overview_fig);

        write_result_workbook(fullfile(output_dir, 'autocorrelogram_results.xlsx'), results, cfg);
        fprintf('\nResults written to "%s".\n', output_dir);
    end

    fprintf('Done.\n');
end


function record = build_record(file_info, prep, window, autocorr, spectrum, boot)
%BUILD_RECORD  Assemble one row of the results table.
%
%   Field order defines the column order of the output table. Measures that
%   were not computed are NaN rather than absent, so that every recording
%   yields the same columns.
    record = struct( ...
        'file_name',             string(file_info.name), ...
        'folder',                string(file_info.folder), ...
        'lowpass_hz',            prep.lowpass_hz, ...
        'downsample_factor',     prep.downsample_factor, ...
        'sampling_rate_hz',      prep.sampling_rate_hz, ...
        'window_start_s',        window.start_time_s, ...
        'window_end_s',          window.end_time_s, ...
        'window_duration_s',     window.duration_s, ...
        'coeff_automatic',       autocorr.coeff_automatic, ...
        'lag_automatic_s',       autocorr.lag_automatic_s, ...
        'freq_automatic_hz',     autocorr.freq_automatic_hz, ...
        'coeff_manual',          autocorr.coeff_manual, ...
        'lag_manual_s',          autocorr.lag_manual_s, ...
        'freq_manual_hz',        autocorr.freq_manual_hz, ...
        'freq_spectrum_peak_hz', spectrum.peak_frequency_hz, ...
        'block_size_s',          boot.block_size_s, ...
        'n_blocks',              boot.n_blocks, ...
        'ci_lower',              boot.ci_lower, ...
        'ci_upper',              boot.ci_upper, ...
        'p_value',               boot.p_value);
end


%% ========================================================================
%  DATA IMPORT
%  ========================================================================

function abf_list = list_recordings(folder_path, search_subfolders)
%LIST_RECORDINGS  List the .abf files of a folder, optionally recursively.
    if search_subfolders
        abf_list = dir(fullfile(folder_path, '**', '*.abf'));
    else
        abf_list = dir(fullfile(folder_path, '*.abf'));
    end

    if isempty(abf_list)
        error('autocorrelogram:noAbfFound', ...
              'No .abf file was found in "%s".', folder_path);
    end
end


function trace_data = load_recording(file_path, channel)
%LOAD_RECORDING  Import one channel of a single .abf file.
%
%   trace_data = LOAD_RECORDING(file_path, channel) returns an n-by-2
%   matrix: column 1 time in s starting at 0, column 2 the signal. The whole
%   recording is read; the analysed window is cut afterwards, so that
%   filtering is applied before windowing and no filter transient is
%   introduced at the window edges.
%
%   See also SELECT_WINDOW, ABFLOAD.
    arguments
        file_path (1,:) char
        channel   (1,1) {mustBeInteger, mustBePositive}
    end

    [signal, sampling_interval_us] = abfload(file_path);

    if channel > size(signal, 2)
        error('autocorrelogram:channelOutOfRange', ...
              'Channel %d was requested but "%s" contains %d channel(s).', ...
              channel, file_path, size(signal, 2));
    end
    signal = signal(:, channel);

    % abfload returns the sampling interval in microseconds
    time_s     = (0:size(signal,1)-1)' * sampling_interval_us * 1e-6;
    trace_data = [time_s, signal];
end


function start_times_s = ask_fixed_start_times(abf_list)
%ASK_FIXED_START_TIMES  One dialog collecting the start time of every recording.
    default_answer = repmat({'0'}, numel(abf_list), 1);
    answer = inputdlg({abf_list.name}, 'Write start time of each recording in s', ...
                      [1 40], default_answer);

    if isempty(answer)
        error('autocorrelogram:noStartTimes', 'No start times were provided.');
    end

    start_times_s = cellfun(@str2double, answer);
    if any(isnan(start_times_s)) || any(start_times_s < 0)
        error('autocorrelogram:badStartTime', ...
              'Start times must be non-negative numbers, in seconds.');
    end
end


%% ========================================================================
%  PREPROCESSING
%  ========================================================================

function prep = preprocess_trace(trace_data, cfg)
%PREPROCESS_TRACE  Lowpass filter and downsample a whole recording.
%
%   prep = PREPROCESS_TRACE(trace_data, cfg) applies, in this order, an IIR
%   lowpass at cfg.lowpass_hz and downsampling by cfg.downsample_factor.
%   Either step is skipped when its setting is [] or 1. The same settings
%   are used for every recording, so preprocessing is fully described by the
%   parameters sheet of the output workbook.
%
%   Lowpass filtering removes the action-potential band from the spectrum.
%   Downsampling speeds up the analysis, and the bootstrap in particular,
%   but introduces aliasing if pushed too far; it should be preceded by
%   filtering below the resulting Nyquist frequency.
%
%   prep is a struct with fields trace, lowpass_hz, downsample_factor,
%   sampling_rate_hz and figure.
%
%   See also SELECT_WINDOW, LOWPASS, DOWNSAMPLE.
    arguments
        trace_data (:,2) {mustBeNumeric, mustBeNonempty}
        cfg        (1,1) struct
    end

    original_rate_hz = 1 / (trace_data(2,1) - trace_data(1,1));
    processed        = trace_data;

    if ~isempty(cfg.lowpass_hz)
        if cfg.lowpass_hz >= original_rate_hz/2
            error('autocorrelogram:badPassband', ...
                  ['cfg.lowpass_hz (%g Hz) must be below the Nyquist frequency of the ' ...
                   'recording (%g Hz).'], cfg.lowpass_hz, original_rate_hz/2);
        end
        processed(:,2) = lowpass(processed(:,2), cfg.lowpass_hz, original_rate_hz, ...
                                 'ImpulseResponse', 'iir');
    end

    if cfg.downsample_factor > 1
        processed = downsample(processed, cfg.downsample_factor);

        new_rate_hz = 1 / (processed(2,1) - processed(1,1));
        if isempty(cfg.lowpass_hz) || cfg.lowpass_hz > new_rate_hz/2
            warning('autocorrelogram:aliasingRisk', ...
                    ['Downsampling to %.1f Hz without filtering below %.1f Hz may alias. ' ...
                     'Consider setting cfg.lowpass_hz.'], new_rate_hz, new_rate_hz/2);
        end
    end

    prep.trace             = processed;
    prep.lowpass_hz        = value_or_nan(cfg.lowpass_hz);
    prep.downsample_factor = cfg.downsample_factor;
    prep.sampling_rate_hz  = 1 / (processed(2,1) - processed(1,1));
    prep.figure            = gobjects(0);

    if isequal(processed, trace_data)
        return % nothing was done, no before/after figure to draw
    end
    prep.figure = plot_before_after(trace_data, processed, cfg);
end


%% ========================================================================
%  WINDOW SELECTION
%  ========================================================================

function window = select_window(trace_data, cfg, recording_name, fixed_start_s)
%SELECT_WINDOW  Define the stretch of trace on which every measure is computed.
%
%   window = SELECT_WINDOW(trace_data, cfg, recording_name, fixed_start_s)
%   returns the analysed window according to cfg.window_mode:
%
%     'click_start'  the trace is plotted and, if the user accepts, one
%                    click sets the start of a window of cfg.window_size_s.
%                    Keeping the length constant across recordings keeps the
%                    autocorrelation coefficients and the bootstrap null
%                    distributions comparable.
%     'click_range'  the trace is plotted and, if the user accepts, two
%                    clicks set the start and the end. Clicks outside the
%                    trace are clamped to its limits and clicks given in
%                    reverse order are swapped.
%     'fixed'        no plot and no clicking; the window starts at
%                    fixed_start_s and lasts cfg.window_size_s.
%
%   In the two click modes the user may decline, in which case the whole
%   recording is analysed.
%
%   window is a struct with fields trace (time restarted at 0),
%   start_time_s, end_time_s, duration_s and figure. start_time_s and
%   end_time_s are expressed in the time base of the recording, so they can
%   be reused later with cfg.window_mode = 'fixed'.
%
%   See also PREPROCESS_TRACE, MP_AUTOCORRELOGRAM.
    arguments
        trace_data     (:,2) {mustBeNumeric, mustBeNonempty}
        cfg            (1,1) struct
        recording_name (1,:) char
        fixed_start_s  (1,1) double
    end

    trace_start_s = trace_data(1,1);
    trace_end_s   = trace_data(end,1);

    switch cfg.window_mode
        case 'fixed'
            limits_s = [fixed_start_s, fixed_start_s + cfg.window_size_s];

        case {'click_start', 'click_range'}
            selection_fig = plot_for_selection(trace_data, recording_name);

            if ask_yes_no(sprintf('Select an analysis window in %s? (No analyses the whole recording)', ...
                                  recording_name), 'Window selection', 'Yes')
                limits_s = collect_clicks(trace_data, cfg.window_mode, cfg.window_size_s);
            else
                limits_s = [trace_start_s, trace_end_s];
            end

            pause(0.5)
            close(selection_fig);

        otherwise
            error('autocorrelogram:badConfig', 'Unknown cfg.window_mode "%s".', cfg.window_mode);
    end

    % clamp to the recording and warn if the requested window did not fit
    clamped_s = [max(limits_s(1), trace_start_s), min(limits_s(2), trace_end_s)];
    if clamped_s(2) - clamped_s(1) < (limits_s(2) - limits_s(1)) - eps(trace_end_s)
        warning('autocorrelogram:windowTruncated', ...
                ['The requested window (%.2f - %.2f s) extends beyond %s, which lasts ' ...
                 '%.2f s. It has been truncated to %.2f - %.2f s.'], ...
                limits_s(1), limits_s(2), recording_name, trace_end_s, clamped_s(1), clamped_s(2));
    end

    [~, index_start] = min(abs(trace_data(:,1) - clamped_s(1))); % snap to sampled points
    [~, index_end]   = min(abs(trace_data(:,1) - clamped_s(2)));
    selected = trace_data(index_start:index_end, :);

    if size(selected, 1) < 2
        error('autocorrelogram:emptyWindow', ...
              'The window selected in %s contains fewer than 2 samples.', recording_name);
    end

    window.start_time_s = selected(1,1);
    window.end_time_s   = selected(end,1);
    window.duration_s   = window.end_time_s - window.start_time_s;
    window.trace        = [selected(:,1) - selected(1,1), selected(:,2)]; % time restarted at 0

    window.figure = figure;
    plot(window.trace(:,1), window.trace(:,2));
    title(sprintf('Analysed window: %s', recording_name), 'Interpreter', 'none')
    xlabel('Time (s)')
    ylabel('Potential (mV)')
end


function limits_s = collect_clicks(trace_data, window_mode, window_size_s)
%COLLECT_CLICKS  Read the one or two clicks defining the window.
    waitfor(msgbox(click_instructions(window_mode, window_size_s), 'modal'));

    first_click = read_click(trace_data);
    xline(first_click, 'r--', 'LineWidth', 1.5);

    if strcmp(window_mode, 'click_start')
        limits_s = [first_click, first_click + window_size_s];
        xline(limits_s(2), 'r--', 'LineWidth', 1.5);
    else
        second_click = read_click(trace_data);
        xline(second_click, 'r--', 'LineWidth', 1.5);
        limits_s = sort([first_click, second_click]); % tolerate clicks in reverse order
    end

    waitfor(msgbox(sprintf('You have selected %.2f s', diff(limits_s))));
end


function instructions = click_instructions(window_mode, window_size_s)
%CLICK_INSTRUCTIONS  Text of the dialog explaining what to click.
    if strcmp(window_mode, 'click_start')
        instructions = sprintf(['Click once on the trace to set the start of the window. ' ...
                                'The window will last %g s from that point.'], window_size_s);
    else
        instructions = ['Click twice on the trace: 1st click for the start, 2nd for the ' ...
                        'end. If you click outside of the trace limits, its limit will be used.'];
    end
end


function time_point = read_click(trace_data)
%READ_CLICK  Read one click and clamp it to the limits of the trace.
    point = ginput(1);
    if isempty(point)
        error('autocorrelogram:noClick', 'No point was selected on the trace.');
    end

    time_point = min(max(point(1,1), trace_data(1,1)), trace_data(end,1));
end


%% ========================================================================
%  ANALYSIS
%  ========================================================================

function autocorr = mp_autocorrelogram(trace_data, opts)
%MP_AUTOCORRELOGRAM  Autocorrelation of a membrane-potential trace.
%
%   autocorr = MP_AUTOCORRELOGRAM(trace_data) computes the normalised
%   autocorrelation of the mean-subtracted trace, finds the local maxima of
%   the positive half of the autocorrelogram and keeps the most prominent
%   one. Its height is the autocorrelation coefficient, an index of
%   oscillation strength; the reciprocal of its lag is the oscillation
%   frequency. Mean subtraction removes the offset due to the holding
%   potential, and prominence rather than height is used so that the peak of
%   the oscillation is preferred over the shoulder of the central peak.
%
%   autocorr = MP_AUTOCORRELOGRAM(trace_data, 'Interactive', true) also
%   plots the autocorrelogram with the detected peak marked and offers the
%   user the option to mark a different one by clicking. The manual value is
%   reported alongside the automatic one and never replaces it: the
%   bootstrap compares the automatic coefficient of the recording with the
%   automatic coefficient of each surrogate, so that both sides of the test
%   use the same statistic.
%
%   trace_data is n-by-2: column 1 time in s starting at 0, column 2 signal.
%
%   autocorr is a struct with fields
%     coeff_automatic    height of the most prominent peak
%     lag_automatic_s    lag of that peak
%     freq_automatic_hz  1 / lag_automatic_s
%     coeff_manual       height at the manually marked lag, NaN if unused
%     lag_manual_s       manually marked lag, NaN if unused
%     freq_manual_hz     1 / lag_manual_s, NaN if unused
%     correlogram        2n-1 by 2, column 1 lag in s, column 2 correlation
%     figure             figure handle, empty when Interactive is false
%
%   See also AUTOCORRELOGRAM_FREQ_FFT, BOOTSTRAP_AUTOCORR, XCORR, FINDPEAKS.
    arguments
        trace_data       (:,2) {mustBeNumeric, mustBeNonempty}
        opts.Interactive (1,1) logical = false
    end

    % normalised autocorrelation of the mean-subtracted signal
    [correlation, lags] = xcorr(trace_data(:,2) - mean(trace_data(:,2)), 'coeff');

    % lags are returned in samples; convert to seconds using the duration of
    % the analysed window, the autocorrelogram being symmetric around zero
    lag_s = linspace(-trace_data(end,1), trace_data(end,1), numel(lags));

    % only the positive half is searched, the function being symmetric
    positive_half = round(numel(correlation)/2):numel(correlation);
    [peak_heights, peak_lags, ~, peak_prominence] = ...
        findpeaks(correlation(positive_half), lag_s(positive_half));

    if isempty(peak_heights)
        error('autocorrelogram:noPeak', ...
              ['No local maximum was found in the autocorrelogram. The window may be ' ...
               'too short or the signal may not be rhythmic.']);
    end

    [~, i_most_prominent] = max(peak_prominence);

    autocorr.coeff_automatic   = peak_heights(i_most_prominent);
    autocorr.lag_automatic_s   = peak_lags(i_most_prominent);
    autocorr.freq_automatic_hz = 1 / autocorr.lag_automatic_s;
    autocorr.coeff_manual      = NaN;
    autocorr.lag_manual_s      = NaN;
    autocorr.freq_manual_hz    = NaN;
    autocorr.correlogram       = [lag_s(:), correlation(:)];
    autocorr.figure            = gobjects(0);

    if ~opts.Interactive
        return
    end

    autocorr.figure = figure('WindowState', 'maximized');
    plot(lag_s, correlation);
    title('Autocorrelogram')
    xlabel('Time (s)')
    ylabel('Correlation Coefficient')
    hold on
    plot(autocorr.lag_automatic_s, autocorr.coeff_automatic + 0.025, ...
         'rv', 'MarkerFaceColor', 'r');
    legend('', 'Automatic CorrCoeff')

    keep_automatic = 'Keep the automatic detection';
    mark_manual    = 'Mark a different peak';
    answer = questdlg('Do you accept the automatically detected CorrCoeff, or do you want to mark a different peak?', ...
                      'AutomaticVsManual', keep_automatic, mark_manual, keep_automatic);

    if strcmp(answer, mark_manual)
        waitfor(msgbox('Enlarge the figure window, and click once on the peak', 'modal'));
        point = ginput(1);

        if isempty(point)
            warning('autocorrelogram:noManualPeak', ...
                    'No peak was marked, only the automatic value is reported.');
        else
            % ginput returns any point of the axes, so it is snapped to the
            % nearest sampled lag
            [~, i_nearest] = min(abs(lag_s - point(1,1)));

            autocorr.coeff_manual   = correlation(i_nearest);
            autocorr.lag_manual_s   = lag_s(i_nearest);
            autocorr.freq_manual_hz = 1 / autocorr.lag_manual_s;

            plot(autocorr.lag_manual_s, autocorr.coeff_manual + 0.025, ...
                 'gv', 'MarkerFaceColor', 'g');
            legend('', 'Automatic CorrCoeff', 'Manual CorrCoeff')
        end
    end

    hold off
    pause(0.5)
    set(autocorr.figure, 'Visible', 'off'); % kept for saving, not shown further
end


function spectrum = autocorrelogram_freq_fft(correlogram, opts)
%AUTOCORRELOGRAM_FREQ_FFT  Dominant frequency of an autocorrelogram by FFT.
%
%   spectrum = AUTOCORRELOGRAM_FREQ_FFT(correlogram) takes the positive half
%   of the autocorrelogram, which carries all of its information because the
%   function is symmetric, and computes its single-sided amplitude spectrum.
%   The frequency of the largest amplitude estimates the dominant
%   oscillation frequency, independently of the peak-detection estimate
%   returned by MP_AUTOCORRELOGRAM.
%
%   correlogram is m-by-2 as returned by MP_AUTOCORRELOGRAM: column 1 lag in
%   s, column 2 correlation.
%
%   spectrum is a struct with fields
%     peak_frequency_hz  frequency of the largest amplitude
%     peak_amplitude     amplitude at that frequency
%     spectrum           k-by-2, column 1 frequency in Hz, column 2 amplitude
%     figure             figure handle with the half correlogram and spectrum
%
%   Name-value argument
%     'FrequencyLimitHz'  x-axis limit of the plotted spectrum, default 10.
%                         It affects the figure only, not the returned peak.
%
%   See also MP_AUTOCORRELOGRAM, FFT.
    arguments
        correlogram           (:,2) {mustBeNumeric, mustBeNonempty}
        opts.FrequencyLimitHz (1,1) {mustBePositive} = 10
    end

    half_index  = round(size(correlogram,1)/2);
    lag_s       = correlogram(half_index:end, 1);
    correlation = correlogram(half_index:end, 2);

    sampling_rate_hz = 1 / (correlogram(2,1) - correlogram(1,1));
    n_samples        = numel(correlation);

    two_sided = abs(fft(correlation) / n_samples);

    % single-sided spectrum. floor() keeps the indexing valid for an odd
    % n_samples as well, and returns the same values as n_samples/2 + 1 when
    % n_samples is even
    n_single  = floor(n_samples/2) + 1;
    amplitude = two_sided(1:n_single);
    if mod(n_samples, 2) == 0
        amplitude(2:end-1) = 2 * amplitude(2:end-1); % DC and Nyquist are unique
    else
        amplitude(2:end)   = 2 * amplitude(2:end);   % odd length has no Nyquist bin
    end
    frequency_hz = sampling_rate_hz * (0:n_single-1)' / n_samples;

    [peak_amplitude, i_peak] = max(amplitude);

    spectrum.peak_frequency_hz = frequency_hz(i_peak);
    spectrum.peak_amplitude    = peak_amplitude;
    spectrum.spectrum          = [frequency_hz, amplitude];

    spectrum.figure = figure;
    subplot(2,1,1)
    plot(lag_s, correlation)
    title('Correlogram')
    xlabel('Time (s)')
    ylabel('Correlation Coefficient')

    subplot(2,1,2)
    plot(frequency_hz, amplitude)
    hold on
    plot(spectrum.peak_frequency_hz, peak_amplitude, 'r*')
    xlim([0 opts.FrequencyLimitHz])
    title('Single-Sided Amplitude Spectrum of Correlogram')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    hold off
end


function boot = bootstrap_autocorr(trace_data, empirical_coeff, cfg)
%BOOTSTRAP_AUTOCORR  Block-bootstrap test of an autocorrelation coefficient.
%
%   boot = BOOTSTRAP_AUTOCORR(trace_data, empirical_coeff, cfg) builds
%   cfg.n_bootstrap surrogates of the analysed window by resampling
%   contiguous blocks of cfg.block_size_s with replacement, and computes on
%   each of them the same statistic as on the data, that is the height of
%   the most prominent peak of its autocorrelogram. The two-sided p-value is
%   the proportion of surrogates whose absolute coefficient is at least as
%   large as the empirical one.
%
%   Because the peak is searched in every surrogate as it is in the data,
%   the null accounts for the search over lags and the test is conservative
%   in that respect.
%
%   The block must be long enough to contain the correlation imposed by any
%   preprocessing filter, and short enough for the window to be divided into
%   at least cfg.min_blocks blocks; a warning is issued otherwise.
%
%   boot is a struct with fields p_value, ci_lower, ci_upper, n_blocks,
%   block_size_s, null_coeff and figure.
%
%   See also BLOCKBOOTSTRP, MP_AUTOCORRELOGRAM.
    arguments
        trace_data      (:,2) {mustBeNumeric, mustBeNonempty}
        empirical_coeff (1,1) double
        cfg             (1,1) struct
    end

    sampling_interval_s = trace_data(2,1) - trace_data(1,1);
    block_size_samples  = round(cfg.block_size_s / sampling_interval_s);
    n_blocks            = floor(size(trace_data,1) / block_size_samples);

    if block_size_samples < 2
        error('autocorrelogram:blockTooShort', ...
              ['cfg.block_size_s (%g s) is shorter than 2 samples at the analysed ' ...
               'sampling rate.'], cfg.block_size_s);
    end
    if n_blocks < cfg.min_blocks
        warning('autocorrelogram:tooFewBlocks', ...
                ['The window yields only %d blocks of %g s (cfg.min_blocks = %d). The ' ...
                 'null distribution will be coarse; use a longer window or a shorter block.'], ...
                n_blocks, cfg.block_size_s, cfg.min_blocks);
    end
    if ~isempty(cfg.lowpass_hz) && cfg.block_size_s < 4/cfg.lowpass_hz
        warning('autocorrelogram:blockShorterThanFilter', ...
                ['The block (%g s) is short relative to the correlation imposed by the ' ...
                 '%g Hz lowpass (about %.2f s). The null may be too permissive; check ' ...
                 'the sensitivity of the p-values to a longer block.'], ...
                cfg.block_size_s, cfg.lowpass_hz, 1/cfg.lowpass_hz);
    end

    null_coeff  = nan(cfg.n_bootstrap, 1);
    time_column = trace_data(:,1);
    iteration_timer = tic;

    for i_rep = 1:cfg.n_bootstrap
        surrogate = blockbootstrp(trace_data(:,2), block_size_samples, 1);
        surrogate_autocorr = mp_autocorrelogram([time_column, surrogate]);
        null_coeff(i_rep)  = surrogate_autocorr.coeff_automatic;

        if i_rep == 1
            fprintf('  bootstrap: about %.1f min for %d iterations\n', ...
                    toc(iteration_timer) * cfg.n_bootstrap / 60, cfg.n_bootstrap);
        elseif mod(i_rep, cfg.progress_interval) == 0
            fprintf('  iteration %d of %d\n', i_rep, cfg.n_bootstrap);
        end
    end

    boundaries = quantile(null_coeff, cfg.ci_quantiles);

    boot.p_value      = sum(abs(null_coeff) >= abs(empirical_coeff)) / cfg.n_bootstrap;
    boot.ci_lower     = boundaries(1);
    boot.ci_upper     = boundaries(2);
    boot.n_blocks     = n_blocks;
    boot.block_size_s = cfg.block_size_s;
    boot.null_coeff   = null_coeff;
    boot.figure       = plot_bootstrap_null(null_coeff, empirical_coeff, boundaries, boot.p_value);

    fprintf('  coefficient %.4f, p = %.4f\n', empirical_coeff, boot.p_value);
end


function boot = empty_bootstrap_result()
%EMPTY_BOOTSTRAP_RESULT  Placeholder used when the bootstrap is switched off.
    boot = struct('p_value', NaN, 'ci_lower', NaN, 'ci_upper', NaN, ...
                  'n_blocks', NaN, 'block_size_s', NaN, ...
                  'null_coeff', [], 'figure', gobjects(0));
end


function [surrogates, n_blocks] = blockbootstrp(data, block_size, n_iterations, verbose)
%BLOCKBOOTSTRP  Block bootstrap of a time series.
%
%   surrogates = BLOCKBOOTSTRP(data, block_size, n_iterations) splits data
%   into contiguous non-overlapping blocks of block_size samples and draws
%   blocks with replacement to rebuild a series of the original length.
%   Resampling whole blocks preserves the structure of the signal within a
%   block while destroying it across blocks.
%
%   If the length of data is not an exact multiple of block_size, a number
%   of randomly chosen samples equal to the remainder is discarded before
%   blocking, and each surrogate is linearly interpolated back onto the
%   original sample positions so that all surrogates keep the input length.
%
%   surrogates is length(data)-by-n_iterations. n_blocks is the number of
%   blocks actually used.
%
%   See also BOOTSTRAP_AUTOCORR, BOOTSTRP.
    arguments
        data         (:,1) {mustBeNumeric, mustBeNonempty}
        block_size   (1,1) {mustBeInteger, mustBePositive}
        n_iterations (1,1) {mustBeInteger, mustBePositive}
        verbose      (1,1) logical = false
    end

    n_samples = numel(data);
    assert(block_size < n_samples, ...
           'blockbootstrp:blockTooLong', ...
           'The block size (%d) must be shorter than the data (%d samples).', ...
           block_size, n_samples);

    sample_positions = (1:n_samples)';
    remainder        = rem(n_samples, block_size);
    was_trimmed      = remainder ~= 0;

    if was_trimmed
        if verbose
            warning('blockbootstrp:notDivisible', ...
                    ['The length of the data (%d) is not divisible by the block size (%d). ' ...
                     '%d random samples (%.2f%%) are removed before blocking and the ' ...
                     'surrogates are linearly interpolated to replace them.'], ...
                    n_samples, block_size, remainder, 100*remainder/n_samples);
        end
        discarded      = randperm(n_samples, remainder);
        kept_data      = data;              kept_data(discarded)      = [];
        kept_positions = sample_positions;  kept_positions(discarded) = [];
    else
        kept_data      = data;
        kept_positions = sample_positions;
    end

    n_blocks = numel(kept_data) / block_size; % integer by construction
    if verbose
        fprintf('Number of blocks: %d\n', n_blocks);
    end
    blocks = reshape(kept_data, [block_size, n_blocks])'; % one block per row

    resampled = bootstrp(n_iterations, @(x) x', blocks)'; % one surrogate per column

    if was_trimmed
        surrogates = zeros(n_samples, n_iterations);
        for i_iteration = 1:n_iterations
            interpolant = griddedInterpolant(kept_positions, resampled(:, i_iteration));
            surrogates(:, i_iteration) = interpolant(sample_positions);
        end
    else
        surrogates = resampled;
    end
end


%% ========================================================================
%  FIGURES
%  ========================================================================

function fig = plot_for_selection(trace_data, recording_name)
%PLOT_FOR_SELECTION  Full-screen plot of a trace, used for clicking on it.
    fig = figure('WindowState', 'maximized');
    plot(trace_data(:,1), trace_data(:,2));
    title(sprintf('Select the desired time range: %s', recording_name), 'Interpreter', 'none');
    xlabel('Time (s)');
    ylabel('Potential (mV)');

    tick_values = xticks;
    xticks(min(tick_values):10:max(tick_values)); % denser ticks help the user aim
    hold on
end


function fig = plot_before_after(trace_before, trace_after, cfg)
%PLOT_BEFORE_AFTER  Two stacked panels showing the effect of preprocessing.
    fig = figure('Visible', 'on');

    subplot(2,1,1)
    plot(trace_before(:,1), trace_before(:,2));
    title('Original')
    xlabel('Time (s)')
    ylabel('Potential (mV)')
    annotation('textbox', [0.8 0.6 0.1 0.1], ...
               'String', sprintf('SR = %.2f kHz', sampling_rate_khz(trace_before)));

    subplot(2,1,2)
    plot(trace_after(:,1), trace_after(:,2));
    title('Preprocessed')
    xlabel('Time (s)')
    ylabel('Potential (mV)')
    annotation('textbox', [0.8 0.1 0.1 0.1], ...
               'String', sprintf('SR = %.2f kHz, lowpass = %s Hz, downsample = %d', ...
                                 sampling_rate_khz(trace_after), ...
                                 num2str(value_or_nan(cfg.lowpass_hz)), ...
                                 cfg.downsample_factor));
end


function fig = plot_bootstrap_null(null_coeff, empirical_coeff, boundaries, p_value)
%PLOT_BOOTSTRAP_NULL  Null distribution with the empirical value marked.
    fig = figure;
    histogram(null_coeff, 'Normalization', 'probability');
    hold on
    xline(empirical_coeff, 'r', 'LineWidth', 2, 'Label', 'recording');
    xline(boundaries(1), 'k--');
    xline(boundaries(2), 'k--');
    hold off

    title(sprintf('Bootstrap null distribution (p = %.4f)', p_value))
    xlabel('Autocorrelation Coefficient')
    ylabel('Probability')
    legend('Surrogates', 'Recording', 'Quantiles', 'Location', 'best')
end


function fig = plot_window_overview(windows, file_names)
%PLOT_WINDOW_OVERVIEW  One stacked panel per analysed window, for screening.
    fig    = figure('WindowState', 'maximized');
    layout = tiledlayout('vertical');

    for i_file = 1:numel(windows)
        nexttile
        plot(windows{i_file}(:,1), windows{i_file}(:,2));
        title(file_names{i_file}, 'Interpreter', 'none');
    end

    title(layout,  'Analysed Windows')
    xlabel(layout, 'Time (s)')
    ylabel(layout, 'Potential (mV)')
end


function save_figure(fig, base_path, cfg)
%SAVE_FIGURE  Export a figure as vector graphics in the configured format.
    if isempty(fig) || ~isgraphics(fig)
        return % the step that would have produced this figure did not run
    end
    exportgraphics(fig, [base_path '.' cfg.figure_format], ...
                   'BackgroundColor', 'none', 'ContentType', 'vector');
end


%% ========================================================================
%  OUTPUT
%  ========================================================================

function output_dir = prepare_output_folder(folder_path, cfg)
%PREPARE_OUTPUT_FOLDER  Create the results folder and its figures subfolder.
    output_dir = fullfile(folder_path, 'Autocorrelogram_results');
    if ~cfg.save_files
        return
    end

    figure_dir = fullfile(output_dir, 'figures');
    if ~isfolder(figure_dir)
        mkdir(figure_dir);
    end
end


function write_result_workbook(file_path, result_table, cfg)
%WRITE_RESULT_WORKBOOK  Write the results and the settings that produced them.
%
%   Sheet 'results'    one row per recording.
%   Sheet 'parameters' every field of cfg, the MATLAB version and the run
%                      date, so that the file documents its own provenance.
    writetable(result_table, file_path, 'Sheet', 'results', 'WriteMode', 'overwritesheet');
    writetable(config_to_table(cfg), file_path, 'Sheet', 'parameters', 'WriteMode', 'overwritesheet');
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


function warn_on_unequal_windows(results)
%WARN_ON_UNEQUAL_WINDOWS  Flag windows of different length across recordings.
%
%   The autocorrelation coefficient and the width of the bootstrap null
%   distribution both depend on the length of the analysed window, so
%   comparing recordings analysed over very different durations is not
%   straightforward.
    durations = results.window_duration_s;
    if numel(durations) < 2
        return
    end

    if (max(durations) - min(durations)) > 0.1 * median(durations)
        warning('autocorrelogram:unequalWindows', ...
                ['The analysed windows differ in length (%.2f to %.2f s). Coefficients ' ...
                 'and p-values are not directly comparable across recordings analysed ' ...
                 'over different durations.'], min(durations), max(durations));
    end
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


function value = value_or_nan(value)
%VALUE_OR_NAN  Replace an empty setting by NaN, for storage in a table.
    if isempty(value)
        value = NaN;
    end
end


%% ========================================================================
%  CONFIGURATION AND DIALOG HELPERS
%  ========================================================================

function cfg = validate_config(cfg)
%VALIDATE_CONFIG  Check the settings and fail early with a clear message.
    cfg.window_mode = validatestring(cfg.window_mode, ...
                                     {'click_start', 'click_range', 'fixed'}, ...
                                     mfilename, 'cfg.window_mode');
    cfg.figure_format = validatestring(cfg.figure_format, {'emf', 'pdf'}, ...
                                       mfilename, 'cfg.figure_format');

    mustBeInteger(cfg.channel);              mustBePositive(cfg.channel);
    mustBeInteger(cfg.downsample_factor);    mustBePositive(cfg.downsample_factor);
    mustBePositive(cfg.window_size_s);
    mustBePositive(cfg.spectrum_xlim_hz);
    mustBeInteger(cfg.n_bootstrap);          mustBePositive(cfg.n_bootstrap);
    mustBePositive(cfg.block_size_s);
    mustBeInteger(cfg.min_blocks);           mustBePositive(cfg.min_blocks);
    mustBeInteger(cfg.progress_interval);    mustBePositive(cfg.progress_interval);

    if ~isempty(cfg.lowpass_hz)
        mustBePositive(cfg.lowpass_hz);
    end

    if numel(cfg.ci_quantiles) ~= 2 || any(cfg.ci_quantiles < 0) || ...
            any(cfg.ci_quantiles > 1) || cfg.ci_quantiles(1) >= cfg.ci_quantiles(2)
        error('autocorrelogram:badConfig', ...
              'cfg.ci_quantiles must be [lower upper] with 0 <= lower < upper <= 1.');
    end

    if cfg.run_bootstrap && cfg.block_size_s >= cfg.window_size_s
        error('autocorrelogram:badConfig', ...
              'cfg.block_size_s (%g) must be shorter than cfg.window_size_s (%g).', ...
              cfg.block_size_s, cfg.window_size_s);
    end

    if ~isempty(cfg.random_seed)
        mustBeInteger(cfg.random_seed);
        mustBeNonnegative(cfg.random_seed);
    end

    if strcmp(cfg.figure_format, 'emf') && ~ispc
        warning('autocorrelogram:emfNotSupported', ...
                'EMF export is only available on Windows. Figures will be saved as PDF.');
        cfg.figure_format = 'pdf';
    end
end


function check_dependencies(cfg)
%CHECK_DEPENDENCIES  Verify that external functions and toolboxes are available.
    if ~exist('abfload', 'file')
        error('autocorrelogram:missingDependency', ...
              ['abfload.m was not found on the MATLAB path. It is third-party code ' ...
               'and must be downloaded separately, see the header of this file.']);
    end

    required = {'findpeaks',  'Signal Processing Toolbox'; ...
                'downsample', 'Signal Processing Toolbox'};
    if ~isempty(cfg.lowpass_hz)
        required = [required; {'lowpass', 'Signal Processing Toolbox'}];
    end
    if cfg.run_bootstrap
        required = [required; {'bootstrp', 'Statistics and Machine Learning Toolbox'}];
    end

    for i_function = 1:size(required, 1)
        if ~exist(required{i_function,1}, 'file')
            error('autocorrelogram:missingToolbox', ...
                  '"%s" was not found. This script requires the %s.', ...
                  required{i_function,1}, required{i_function,2});
        end
    end
end


function set_random_seed(cfg)
%SET_RANDOM_SEED  Fix the random stream when a seed is configured.
    if ~cfg.run_bootstrap
        return
    end

    if isempty(cfg.random_seed)
        fprintf(['Random seed not set: the bootstrap draws differ between runs, so\n' ...
                 'p-values vary on the order of 1/n_bootstrap. Set cfg.random_seed\n' ...
                 'for an exactly repeatable run.\n']);
    else
        rng(cfg.random_seed);
        fprintf('Random seed set to %d.\n', cfg.random_seed);
    end
end


function is_yes = ask_yes_no(question, dialog_title, default_answer)
%ASK_YES_NO  Yes/No dialog returning a logical. Closing it keeps the default.
    answer = questdlg(question, dialog_title, 'Yes', 'No', default_answer);
    if isempty(answer)
        answer = default_answer;
    end
    is_yes = strcmp(answer, 'Yes');
end


function rate_khz = sampling_rate_khz(trace_data)
%SAMPLING_RATE_KHZ  Sampling rate of a trace in kHz, from its first two samples.
    rate_khz = 1 / (trace_data(2,1) - trace_data(1,1)) / 1000;
end
