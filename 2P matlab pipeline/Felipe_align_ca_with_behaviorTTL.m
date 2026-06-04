function Felipe_align_ca_with_behaviorTTL(dirsel, overwrite)
% Align your EXISTING calcium TTL frames (dataFrame) to behavior TTL frames (dataTTL_beh).
% Saves align_beh with:
%   .fs
%   .ca.frameSamp
%   .beh.frameSamp
%   .ca2beh_frameIdx  (for each CA frame, the matched BEH frame index)
%   .ca.frameTime_s
%   .beh.frameTime_s

    global outputDirCardin info
    analysisCaTTL  = 'lev0_readframes';      % your existing calcium TTL output folder
    analysisBehTTL = 'lev0_readframes_beh';  % behavior TTL folder from function above

    for iDir = dirsel
        exptag = info(iDir).dir;
        parts  = split(exptag,'_'); mouse = parts{1};

        % Existing calcium TTL file (from your previous code)
        caTTL_file  = fullfile(outputDirCardin, analysisCaTTL,  mouse, exptag, [exptag '.mat']);
        % Behavior TTL file (from the behavior-only reader)
        behTTL_file = fullfile(outputDirCardin, analysisBehTTL, mouse, exptag, [exptag '.mat']);

        if ~exist(caTTL_file,'file')
            warning('Missing CALCIUM TTL file for %s. Run your calcium TTL extraction first. Skipping.', exptag);
            continue
        end
        if ~exist(behTTL_file,'file')
            warning('No behavior TTL for %s. Skipping alignment.', exptag);
            continue
        end

        % Load calcium TTL (expects dataFrame per your existing code)
        S_ca = load(caTTL_file);
        if ~isfield(S_ca,'dataFrame')
            warning('calcium TTL file for %s does not contain dataFrame. Skipping.', exptag);
            continue
        end
        dataFrame = S_ca.dataFrame;

        % Load behavior TTL (expects dataTTL_beh)
        S_beh = load(behTTL_file);
        if ~isfield(S_beh,'dataTTL_beh')
            warning('behavior TTL file for %s does not contain dataTTL_beh. Skipping.', exptag);
            continue
        end
        dataTTL_beh = S_beh.dataTTL_beh;

        % Extract frames
        if ~isfield(dataFrame,'frameInd') || isempty(dataFrame.frameInd)
            warning('Empty calcium frameInd for %s. Skipping.', exptag);
            continue
        end
        if ~isfield(dataTTL_beh,'frameInd') || isempty(dataTTL_beh.frameInd)
            warning('Empty behavior frameInd for %s. Skipping.', exptag);
            continue
        end

        fs_ca  = dataFrame.fsample;     % TTL fs used in your calcium extractor
        fs_beh = dataTTL_beh.fs;        % TTL fs used by behavior
        if fs_ca ~= fs_beh
            warning('TTL fs differ (CA=%g, BEH=%g) for %s. Proceeding with CA fs for timing.', fs_ca, fs_beh, exptag);
        end
        fs = fs_ca; % common sample rate for time conversion

        caSamp  = double(dataFrame.frameInd(:));
        behSamp = double(dataTTL_beh.frameInd(:));

        % nearest behavior frame index for each calcium frame
        behIdx = interp1(behSamp, (1:numel(behSamp))', caSamp, 'nearest', 'extrap');

        align_beh = struct();
        align_beh.fs               = fs;
        align_beh.ca.frameSamp     = caSamp;
        align_beh.beh.frameSamp    = behSamp;
        align_beh.ca2beh_frameIdx  = behIdx;
        align_beh.ca.frameTime_s   = caSamp ./ fs;
        align_beh.beh.frameTime_s  = behSamp ./ fs;

        % Monotonicity check (ideally non-decreasing)
        if any(diff(behIdx) < 0)
            warning('Non-monotonic CA->BEH mapping in %s. Check for dropped/extra pulses.', exptag);
        end

        % Save next to behavior TTL result
        outDir = fileparts(behTTL_file);
        save(fullfile(outDir, [exptag '_align_beh.mat']), 'align_beh', '-v7.3');
        fprintf('Saved CA↔BEH alignment: %s\n', fullfile(outDir, [exptag '_align_beh.mat']));
    end
end