function lev0_readframes_fromTTL_v2(dirsel, overwrite, ttlProvider, opts)
% Build dataFrame from a TTL vector (e.g., 2 kHz DAQ) with robust pulse detection.
% Uses onset/offset detection (from your Detect_TTL_input), filters by width and refractory,
% and segments adaptively based on gaps relative to the median IFI.
%
% Inputs:
%   dirsel, overwrite    : as before (indices into global info, and overwrite flag)
%   ttlProvider          : function handle or struct
%                          - If function: S = ttlProvider(exptag) must return fields:
%                               S.ttl (vector), S.fs (Hz), optional S.name, S.segments (frame ranges), S.threshold_ON
%                          - If struct: fields above directly provided
%   opts (struct, optional fields):
%       .threshold_ON       : numeric threshold for TTL detection (overrides S.threshold_ON)
%       .minPulseWidth_ms   : default 1.0 ms (reject shorter)
%       .maxPulseWidth_ms   : default Inf   (no upper rejection)
%       .minIFI_ms          : default 15 ms (refractory between pulses)
%       .gapFactor          : default 5 (segment if gap > gapFactor * medianIFI)
%       .plotDebug          : default false (quick diagnostics)
%
% Output: saves dataFrame identical to old lev0_readframes, so lev1 can be reused.

global outputDirCardin
global info

analysis = 'lev0_readframes';
nDirs = length(info);
if nargin < 1 || isempty(dirsel), dirsel = 1:nDirs; end
if nargin < 2 || isempty(overwrite), overwrite = 0; end
if nargin < 4, opts = struct(); end

% Defaults
if ~isfield(opts,'minPulseWidth_ms'), opts.minPulseWidth_ms = 1.0; end
if ~isfield(opts,'maxPulseWidth_ms'), opts.maxPulseWidth_ms = Inf; end
if ~isfield(opts,'minIFI_ms'),        opts.minIFI_ms        = 15;  end
if ~isfield(opts,'gapFactor'),         opts.gapFactor        = 5;   end
if ~isfield(opts,'plotDebug'),         opts.plotDebug        = false; end

for iDir = dirsel
    exptag = info(iDir).dir;
    mouse  = exptag(1:6);
    outputDir = fullfile(outputDirCardin, analysis, mouse, exptag);
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    % ---- Get TTL and fs ----
    if isa(ttlProvider, 'function_handle')
        S = ttlProvider(exptag);
    else
        S = ttlProvider;
    end
    ttl = S.ttl(:);                 % column for convenience
    fs  = S.fs;
    if isfield(S,'name') && ~isempty(S.name), baseName = S.name; else, baseName = exptag; end

    outputFilename = fullfile(outputDir, baseName);
    if exist([outputFilename '.mat'], 'file') && overwrite == 0
        fprintf('skipping %s file %s \n', exptag, outputFilename);
        continue
    end

    % ---- Timebase ----
    N  = numel(ttl);
    t  = (0:N-1) / fs;

    % ---- Basic preconditioning ----
    % Replace NaNs with local median if any
    if any(isnan(ttl))
        nanmask = isnan(ttl);
        ttl(nanmask) = interp1(find(~nanmask), ttl(~nanmask), find(nanmask), 'linear','extrap');
    end
    % Optional: light detrend if there is slow drift around TTL baseline
    % ttl = detrend(ttl,1);

    % ---- Threshold selection ----
    if isfield(opts,'threshold_ON') && ~isempty(opts.threshold_ON)
        thrON = opts.threshold_ON;
    elseif isfield(S,'threshold_ON') && ~isempty(S.threshold_ON)
        thrON = S.threshold_ON;
    else
        % Auto: mid between modes (fallback to half-range)
        lo = prctile(ttl, 5); hi = prctile(ttl, 95);
        thrON = lo + 0.5*(hi - lo);
    end

    % ---- Detect onsets/offsets using your helper ----
    [Start_ON, Stop_ON] = local_detect_ttl_onoff(ttl, thrON);

    if isempty(Start_ON)
        warning('No TTL pulses detected for %s. Check threshold/polarity.', exptag);
        % still save an empty dataFrame for traceability
    end

    % ---- Pulse width filter ----
    pwSamp = Stop_ON - Start_ON + 1;
    minPW  = round((opts.minPulseWidth_ms/1000) * fs);
    maxPW  = isfinite(opts.maxPulseWidth_ms) * round((opts.maxPulseWidth_ms/1000) * fs) + ~isfinite(opts.maxPulseWidth_ms)*Inf;
    keepPW = (pwSamp >= max(1,minPW)) & (pwSamp <= maxPW);

    Start_ON = Start_ON(keepPW);
    Stop_ON  = Stop_ON(keepPW);
    pwSamp   = pwSamp(keepPW);

    % ---- Refractory (min IFI) ----
    minIFI = round((opts.minIFI_ms/1000) * fs);
    if numel(Start_ON) > 1
        isi = diff(Start_ON);
        keep = [true; isi >= minIFI];
        Start_ON = Start_ON(keep);
        Stop_ON  = Stop_ON(keep(1:end-1));  % best effort; last Stop kept implicitly
        % Recompute pwSamp for kept pulses
        pwSamp = Stop_ON - Start_ON(1:numel(Stop_ON)) + 1;
    end

    frameInd = Start_ON(:)';   % row vector of rising edges

    % ---- Segment detection (adaptive to your frame rate) ----
    segments_frame = [];
    if isfield(S,'segments') && ~isempty(S.segments)
        segments_frame = S.segments;
        segments_frame(:,1) = max(segments_frame(:,1), 1);
        segments_frame(:,2) = min(segments_frame(:,2), numel(frameInd));
    else
        if numel(frameInd) >= 3
            IFI = diff(frameInd) / fs;             % seconds
            medIFI = median(IFI);                   % robust estimate of frame period
            gapThresh_sec = opts.gapFactor * medIFI;
            gapThresh_samp = max(2, round(gapThresh_sec * fs));
            dFI = diff(frameInd);
            segStarts = [1, find(dFI > gapThresh_samp) + 1];
            segEnds   = [find(dFI > gapThresh_samp), numel(frameInd)];
            segments_frame = [segStarts(:), segEnds(:)];
        elseif ~isempty(frameInd)
            segments_frame = [1, numel(frameInd)];
        else
            segments_frame = zeros(0,2);
        end
    end

    % ---- Build imageInd (Nx4) & imageTime (Nx2) ----
    numMovies = size(segments_frame,1);
    imageInd  = zeros(numMovies, 4);
    imageTime = zeros(numMovies, 2);
    for k = 1:numMovies
        fStart = segments_frame(k,1);
        fEnd   = segments_frame(k,2);
        hsStart = frameInd(fStart);
        hsEnd   = frameInd(fEnd);
        imageInd(k,:)  = [hsStart, hsEnd, fStart, fEnd];
        imageTime(k,:) = [t(hsStart), t(hsEnd)];
    end

    % ---- Pack & save dataFrame ----
    dataFrame = struct();
    dataFrame.trial{1}           = ttl(:)'; % keep raw waveform (row)
    dataFrame.fsample            = fs;
    dataFrame.time{1}            = t;
    dataFrame.label              = {baseName};
    dataFrame.cfg                = [];
    dataFrame.time_timestamp{1}  = t;

    dataFrame.frameInd           = frameInd;
    dataFrame.imageInd           = imageInd;
    dataFrame.imageTime          = imageTime;
    dataFrame.imageIndLabel      = {'Start Index High Freq','End Index High Freq','Start Index Low Freq','End Index Low Freq'};
    dataFrame.imageTimeLabel     = {'Start Time','End Time'};
    dataFrame.numMovies          = numMovies;

    % Optional debug
    if opts.plotDebug
        local_plot_debug(ttl, t, frameInd, Start_ON, Stop_ON, fs, exptag);
    end

    save([outputFilename '.mat'], 'dataFrame', '-v7.3');
    fprintf('Saved dataFrame: %s.mat (frames=%d, movies=%d)\n', outputFilename, numel(frameInd), numMovies);
end
end

% --------- helpers ---------

function [Start_ON,Stop_ON] = local_detect_ttl_onoff(sg,thrON)
% Equivalent to your Detect_TTL_input (no figure by default)
lgON = sg(:) > thrON;
d = diff([false; lgON; false]);
Start_ON = find(d == 1);
Stop_ON  = find(d == -1) - 1;
end

function local_plot_debug(ttl, t, frameInd, Start_ON, Stop_ON, fs, titleTag)
if numel(ttl) > 2e5
    % don’t plot the whole session—plot first few seconds
    Tplot = min(5, t(end));
    idx = t <= Tplot;
    ttls = ttl(idx); ts = t(idx);
    so = Start_ON(Start_ON <= find(idx,1,'last'));
    eo = Stop_ON(Stop_ON <= find(idx,1,'last'));
else
    ttls = ttl; ts = t;
    so = Start_ON; eo = Stop_ON;
end
figure('Name',['TTL Debug - ' titleTag],'Color','w'); hold on
plot(ts, ttls, 'k-');
if ~isempty(so)
    stem(ts(so), max(ttls)*ones(size(so)), 'r','filled','Marker','none');
end
if ~isempty(eo)
    stem(ts(eo), max(ttls)*0.9*ones(size(eo)), 'g','filled','Marker','none');
end
xlabel('Time (s)'); ylabel('TTL (a.u.)');
legend({'TTL','Start','Stop'}); title(sprintf('fs=%.1f Hz, pulses=%d',fs,numel(frameInd)));
grid on
end