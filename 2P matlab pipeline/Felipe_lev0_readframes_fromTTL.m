    function Felipe_lev0_readframes_fromTTL(dirsel, overwrite, opts)
    % SIMPLE VERSION – reads TTL pulses (2 kHz) and produces dataFrame.
    %
    % The TTL file MUST be named:
    %           TTL_<exptag>.mat
    % Example:  info(1).dir = 'SD2_260218_01'
    %           TTL filename = 'TTL_SD2_260218_01.mat'
    %
    % The TTL .mat file MUST contain:
    %   ttl  - vector (1×N or N×1)
    %   fs   - sampling rate (2000)
    %
    
    %% ------------------------------------------------------------------------
    % Globals
    % -------------------------------------------------------------------------
    global outputDirCardin info
    
    %% Defaults for opts
    if nargin < 3 || isempty(opts)
        opts.threshold_ON      = [];
        opts.minPulseWidth_ms  = 1;
        opts.minIFI_ms         = 15;
        opts.gapFactor         = 5;
        opts.plotDebug         = false;
    end
    if ~isfield(opts,'threshold_ON'),      opts.threshold_ON = [];  end
    if ~isfield(opts,'minPulseWidth_ms'),  opts.minPulseWidth_ms = 1;  end
    if ~isfield(opts,'minIFI_ms'),         opts.minIFI_ms = 15;       end
    if ~isfield(opts,'gapFactor'),         opts.gapFactor = 5;        end
    if ~isfield(opts,'plotDebug'),         opts.plotDebug = false;    end
    
    %% ------------------------------------------------------------------------
    % TTL folder (CHANGE THIS TO YOUR LOCATION)
    % -------------------------------------------------------------------------
    
    % ---- TTL folder (respect opts) ----
    if isfield(opts,'ttlFolder') && ~isempty(opts.ttlFolder)
        ttlFolder = opts.ttlFolder;
    else
        ttlFolder = 'D:\TTL\';  % default fallback
    end

    
    %% ------------------------------------------------------------------------
    % Loop experiments
    % -------------------------------------------------------------------------
    for iDir = dirsel
    
        exptag = info(iDir).dir;
        mouse  = info(iDir).mouse;
        fsTTL  = info(iDir).fsample;
    
        fprintf('\n=============================================\n');
        fprintf('  Level 0 (TTL readframes) – %s\n', exptag);
        fprintf('=============================================\n');
    
        %---------- Build output folder -------------------------------
        outputDir = fullfile(outputDirCardin, 'lev0_readframes', mouse, exptag);
        if ~exist(outputDir, 'dir'), mkdir(outputDir); end
    
        outputFilename = fullfile(outputDir, exptag);
    
        if exist([outputFilename '.mat'], 'file') && overwrite == 0
            fprintf('Skipping (exists): %s\n', outputFilename);
            continue
        end
    
        %---------- Load TTL file --------------------------------------
        ttlFile = fullfile(ttlFolder, ['TTL_' exptag '.mat']);
    
        if ~exist(ttlFile, 'file')
            error('TTL file not found: %s', ttlFile);
        end
    
        tmp = load(ttlFile);   % must contain "ttl" and "fs"
    
        ttl = tmp.ttl(:);
        fs  = tmp.fs;
    
        if fs ~= fsTTL
            warning('TTL fs (%d Hz) != info fsample (%d Hz). Using TTL fs.', fs, fsTTL);
        end
    
        %---------- Timebase --------------------------------------------
        N = numel(ttl);
        t = (0:N-1) / fs;
    
        %---------- Thresholding -----------------------------------------
        if isempty(opts.threshold_ON)
            lo = prctile(ttl,5);
            hi = prctile(ttl,95);
            thr = lo + 0.5*(hi-lo);
        else
            thr = opts.threshold_ON;
        end
    
        % Detect onset/offset:
        lg = ttl > thr;
        d  = diff([false; lg; false]);
        Start_ON = find(d == 1);
        Stop_ON  = find(d == -1) - 1;
    
        %---------- Pulse width filter ------------------------------------
        pw = Stop_ON - Start_ON + 1;
        minPW = round((opts.minPulseWidth_ms/1000)*fs);
        keep = pw >= minPW;
        Start_ON = Start_ON(keep);
    
        %---------- Refractory filter -------------------------------------
        minIFI = round((opts.minIFI_ms/1000)*fs);
        if numel(Start_ON) > 1
            isi = diff(Start_ON);
            keep = [true; isi >= minIFI];
            Start_ON = Start_ON(keep);
        end
    
        frameInd = Start_ON(:)';
    
        %---------- Segment detection -------------------------------------
        if numel(frameInd) >= 3
            IFI = diff(frameInd) / fs;
            medIFI = median(IFI);
            gapThresh = opts.gapFactor * medIFI;
            gapSamp   = round(gapThresh * fs);
    
            dFI = diff(frameInd);
            segStarts = [1, find(dFI > gapSamp) + 1];
            segEnds   = [find(dFI > gapSamp), numel(frameInd)];
            segments = [segStarts(:), segEnds(:)];
        else
            segments = [1, numel(frameInd)];
        end
    
        %---------- Build imageInd and imageTime --------------------------
        numMovies = size(segments,1);
        imageInd  = zeros(numMovies,4);
        imageTime = zeros(numMovies,2);
    
        for k = 1:numMovies
            fStart = segments(k,1);
            fEnd   = segments(k,2);
    
            hsStart = frameInd(fStart);
            hsEnd   = frameInd(fEnd);
    
            imageInd(k,:)  = [hsStart hsEnd fStart fEnd];
            imageTime(k,:) = [t(hsStart) t(hsEnd)];
        end
    
        %---------- Pack dataFrame -----------------------------------------
        dataFrame = struct();
        dataFrame.trial{1}          = ttl';
        dataFrame.fsample           = fs;
        dataFrame.time{1}           = t;
        dataFrame.frameInd          = frameInd;
        dataFrame.imageInd          = imageInd;
        dataFrame.imageTime         = imageTime;
        dataFrame.numMovies         = numMovies;
    
        dataFrame.imageIndLabel  = {'Start_HS','End_HS','Start_Frame','End_Frame'};
        dataFrame.imageTimeLabel = {'Start_Time','End_Time'};
    
        %---------- Plot debug ------------------------------------------------
        if opts.plotDebug
            figure('Color','w'); hold on;
            plot(t, ttl);
            stem(t(frameInd), max(ttl)*ones(size(frameInd)), 'r');
            title(['TTL detection: ' exptag]);
            xlabel('Time (s)');
            ylabel('TTL (a.u.)');
        end
    
        %---------- Save -------------------------------------------------------
        save([outputFilename '.mat'], 'dataFrame','-v7.3');
        fprintf('Saved: %s.mat\n', outputFilename);
    
    end
    
    end