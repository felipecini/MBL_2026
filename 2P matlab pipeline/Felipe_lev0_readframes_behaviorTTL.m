function Felipe_lev0_readframes_behaviorTTL(dirsel, overwrite, opts)
% Reads ONLY behavior TTL pulses from TTL_<exptag> files and saves dataTTL_beh.
% Output saved at: fullfile(outputDirCardin, 'lev0_readframes_beh', mouse, exptag, [exptag '.mat'])

    %% Globals
    global outputDirCardin info

    % Defaults
    if nargin < 3 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'ttlFolder'),         opts.ttlFolder = 'D:\TTL'; end
    if ~isfield(opts,'threshold_ON'),      opts.threshold_ON = [];    end
    if ~isfield(opts,'minPulseWidth_ms'),  opts.minPulseWidth_ms = 1; end
    if ~isfield(opts,'minIFI_ms'),         opts.minIFI_ms = 15;       end
    if ~isfield(opts,'gapFactor'),         opts.gapFactor = 5;        end
    if ~isfield(opts,'plotDebug'),         opts.plotDebug = true;    end
    if ~isfield(opts,'behColumn'),         opts.behColumn = 2;        end % if TTL has [sample, ca, beh], beh is col 2

    analysis = 'lev0_readframes_beh';

    for iDir = dirsel
        exptag = info(iDir).dir;
        mouse  = info(iDir).mouse;
        fs_info= info(iDir).fsample;

        fprintf('\n=============================================\n');
        fprintf('  Level 0 (TTL behavior only) – %s\n', exptag);
        fprintf('=============================================\n');

        outDir = fullfile(outputDirCardin, analysis, mouse, exptag);
        if ~exist(outDir,'dir'), mkdir(outDir); end
        outFile = fullfile(outDir, exptag);

        if exist([outFile '.mat'],'file') && overwrite==0
            fprintf('Skipping (exists): %s\n', outFile);
            continue
        end

        % Locate TTL file
        ttlFileMat = fullfile(opts.ttlFolder, ['TTL_' exptag '.mat']);
        ttlFileCsv = fullfile(opts.ttlFolder, ['TTL_' exptag '.csv']);
        ttlFileTxt = fullfile(opts.ttlFolder, ['TTL_' exptag '.txt']);

  
        % ---------- Locate and load behavior TTL and fs ----------
        behTTL = []; fs = [];
        
        if exist(ttlFileMat,'file')
            S = load(ttlFileMat);
        
            % 1) Nested under Meta.*  (your case)
            if isfield(S,'meta') && isstruct(S.meta)
                % try common names for the behavior channel
                cand = {'behaviorTTL','behTTL','beh','behavior','ttl_beh'};
                for k = 1:numel(cand)
                    if isfield(S.meta, cand{k})
                        behTTL = double(S.meta.(cand{k})(:));
                        srcField = ['S.meta.' cand{k}];
                        break
                    end
                end
                % try sampling rate candidates in Meta
                fs_cand = {'behaviorFs','fs','TTLfs'};
                for k = 1:numel(fs_cand)
                    if isempty(fs) && isfield(S.meta, fs_cand{k}) && ~isempty(S.meta.(fs_cand{k}))
                        fs = double(S.meta.(fs_cand{k}));
                        break
                    end
                end
            end
        
            % 2) Nested under TTL.* (other rigs)
            if isempty(behTTL) && isfield(S,'TTL') && isstruct(S.TTL)
                cand = {'behaviorTTL','behTTL','beh','behavior','ttl_beh'};
                for k = 1:numel(cand)
                    if isfield(S.TTL, cand{k})
                        behTTL = double(S.TTL.(cand{k})(:));
                        srcField = ['S.TTL.' cand{k}];
                        break
                    end
                end
                fs_cand = {'fs','behaviorFs','TTLfs'};
                for k = 1:numel(fs_cand)
                    if isempty(fs) && isfield(S.TTL, fs_cand{k}) && ~isempty(S.TTL.(fs_cand{k}))
                        fs = double(S.TTL.(fs_cand{k}));
                        break
                    end
                end
            end
        
            % 3) Flat arrays / numeric matrices
            if isempty(behTTL) && isfield(S,'ttl_beh') && isfield(S,'fs')
                behTTL = double(S.ttl_beh(:)); fs = double(S.fs); srcField = 'S.ttl_beh';
            elseif isempty(behTTL) && isfield(S,'ttl') && isnumeric(S.ttl) && size(S.ttl,2) >= 2 && isfield(S,'fs')
                behTTL = double(S.ttl(:,2));   fs = double(S.fs); srcField = 'S.ttl(:,2)';
            elseif isempty(behTTL) && isfield(S,'TTL') && isnumeric(S.TTL) && size(S.TTL,2) >= 3 && isfield(S,'fs')
                behTTL = double(S.TTL(:,3));   fs = double(S.fs); srcField = 'S.TTL(:,3)';
            end
        
            % 4) Fallback fs if still empty
            if (isempty(fs) || fs <= 0)
                if isfield(S,'fs') && ~isempty(S.fs)
                    fs = double(S.fs);
                else
                    fs = fs_info; % <- from your info(iDir).fsample
                    warning('fs not found in MAT; using info(iDir).fsample = %g Hz', fs);
                end
            end
        
            if isempty(behTTL)
                error('Could not find behavior TTL in %s (checked Meta.*, TTL.*, ttl_beh, ttl/TLL matrices).', ttlFileMat);
            end
        
        elseif exist(ttlFileCsv,'file') || exist(ttlFileTxt,'file')
            if exist(ttlFileCsv,'file'), fpath = ttlFileCsv; else, fpath = ttlFileTxt; end
            M = readmatrix(fpath);
            if size(M,2) < opts.behColumn
                error('CSV/TXT TTL must have >= %d columns; got %d.', opts.behColumn, size(M,2));
            end
            behTTL = double(M(:, opts.behColumn));
            fs     = fs_info; % assume known from info
        
        else
            error('TTL file not found for %s in %s', exptag, opts.ttlFolder);
        end
        
        % Consistency notice if fs differs
        if fs ~= fs_info
            warning('TTL fs (%g) != info fsample (%g). Proceeding with TTL fs.', fs, fs_info);
        end
        
        % Helpful debug print
        fprintf('Loaded behavior TTL from %s | N=%d | fs=%.3f Hz | min=%.5f | med=%.5f | max=%.5f\n', ...
                exist('srcField','var')*srcField + ~exist('srcField','var'), ...
                numel(behTTL), fs, min(behTTL), median(behTTL), max(behTTL));
        
            
        % Detect behavior frame pulses
        [frameInd_beh, segments_beh, t_beh, thrUsed] = detect_frames_only(behTTL, fs, opts);

        dataTTL_beh = struct();
        dataTTL_beh.fs        = fs;
        dataTTL_beh.frameInd  = frameInd_beh(:)'; % sample indices where behavior frames occurred
        dataTTL_beh.segments  = segments_beh;
        dataTTL_beh.imageInd  = build_imageInd(frameInd_beh, segments_beh);
        dataTTL_beh.imageTime = [t_beh(dataTTL_beh.imageInd(:,1)) t_beh(dataTTL_beh.imageInd(:,2))];
        dataTTL_beh.numMovies = size(dataTTL_beh.imageInd,1);
        dataTTL_beh.threshold = thrUsed;

        if opts.plotDebug
            figure('Color','w','Name',['TTL behavior - ' exptag]); 
            ax1 = subplot(2,1,1); plot(t_beh, behTTL, 'k-'); hold on;
            yline(thrUsed, '--r'); title('Behavior TTL (raw)'); xlabel('Time (s)'); ylabel('TTL');
            ax2 = subplot(2,1,2); stem(t_beh(frameInd_beh), ones(size(frameInd_beh)), 'b.'); 
            title('Detected behavior frames'); xlabel('Time (s)'); ylabel('frames');
            linkaxes([ax1, ax2], 'x');
        end

        save([outFile '.mat'], 'dataTTL_beh', '-v7.3');
        fprintf('Saved behavior TTL: %s.mat\n', outFile);
    end
end

% --- helpers ---
function [frameInd, segments, t, thr] = detect_frames_only(ttl, fs, opts)
    N = numel(ttl);
    t = (0:N-1)./fs;

    if isempty(opts.threshold_ON)
        lo  = prctile(ttl,5);
        hi  = prctile(ttl,95);
        thr = lo + 0.5*(hi-lo);
    else
        thr = opts.threshold_ON;
    end

    lg = ttl > thr;
    d  = diff([false; lg; false]);
    Start_ON = find(d==1);
    Stop_ON  = find(d==-1) - 1;

    % pulse width filter
    minPW = round((opts.minPulseWidth_ms/1000)*fs);
    if ~isempty(Start_ON)
        pw = Stop_ON - Start_ON + 1;
        keep = pw >= minPW;
        Start_ON = Start_ON(keep);
    end

    % refractory IFI
    minIFI = round((opts.minIFI_ms/1000)*fs);
    if numel(Start_ON) > 1
        isi = diff(Start_ON);
        keep = [true; isi >= minIFI];
        Start_ON = Start_ON(keep);
    end

    frameInd = Start_ON(:)';

    % segment detection (gaps >> median IFI)
    if numel(frameInd) >= 3
        IFI = diff(frameInd)/fs;
        medIFI = median(IFI);
        gapThresh = opts.gapFactor * medIFI;
        gapSamp   = round(gapThresh * fs);
        dFI = diff(frameInd);
        segStarts = [1, find(dFI > gapSamp) + 1];
        segEnds   = [find(dFI > gapSamp), numel(frameInd)];
        segments  = [segStarts(:), segEnds(:)];
    else
        segments = [1, numel(frameInd)];
    end
end

function imageInd = build_imageInd(frameInd, segments)
    if isempty(frameInd)
        imageInd = zeros(0,4);
        return
    end
    imageInd = zeros(size(segments,1),4);
    for k=1:size(segments,1)
        fStart = segments(k,1);
        fEnd   = segments(k,2);
        hsStart= frameInd(fStart);
        hsEnd  = frameInd(fEnd);
        imageInd(k,:) = [hsStart hsEnd fStart fEnd];
    end
end