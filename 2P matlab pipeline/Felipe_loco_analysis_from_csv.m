function Felipe_loco_analysis_from_csv(dirsel, overwrite, optsBout, optsPlot)
% ──────────────────────────────────────────────────────────────────────────────
%  Felipe_loco_analysis_from_csv
%
%  Performs:
%   1) Calcium ↔ Behavior mapping (size-safe: clamp & pad where needed)
%   2) Speed mapping
%   3) Locomotion bouts (your rules)
%   4) ΔF/F summaries:
%        - All stillness
%        - All locomotion
%        - Pre‑locomotion stillness (the still period immediately before each bout)
%   5) Peri‑locomotion onset ΔF/F plot (mean ± SEM/SD)
%   6) Extract exemplars of CA traces + speed for 5 locomotion bouts
%
%  Saves:
%     <exptag>_loco_summary.mat
%     <exptag>_box_still_vs_loco.png
%     <exptag>_peri_onset_mean_<sem|sd>.png
%     <exptag>_examples.mat
% ──────────────────────────────────────────────────────────────────────────────

    global info outputDirCardin

    analysisBeh = 'lev0_readframes_beh';
    analysisCA  = 'lev2_calculate_dFF';

    %% Defaults
    if nargin < 3 || isempty(optsBout), optsBout = struct(); end
    if nargin < 4 || isempty(optsPlot), optsPlot = struct(); end

    defBout = struct('thrOn_cm_s',1, 'peakMin_cm_s',10,'mergeGap_s',1,'minDur_s',5);
    defPlot = struct('peri_pre_s',5, 'peri_post_s',10, ...
                     'shadeType','sem', 'minOnsets',3, ...
                     'preStill_s', 2, ...           % NEW: stillness window immediately before onset (s)
                     'examples_N', 5, ...           % how many exemplars to save
                     'examples_strategy', 'longest' ... % 'longest' or 'first'
                     );

    optsBout = setDefaults(optsBout, defBout);
    optsPlot = setDefaults(optsPlot, defPlot);

    %% Loop experiments
    for iDir = dirsel

        exptag = info(iDir).dir;
        mouse  = info(iDir).mouse;

        behDir  = fullfile(outputDirCardin, analysisBeh, mouse, exptag);
        behTTLf = fullfile(behDir, [exptag '.mat']);
        behSpdf = fullfile(behDir, [exptag '_behSpeed.mat']);
        alignf  = fullfile(behDir, [exptag '_align_beh.mat']);
        dfff    = fullfile(outputDirCardin, analysisCA, mouse, exptag, [exptag '.mat']);

        if ~exist(behTTLf,'file') || ~exist(behSpdf,'file') || ~exist(alignf,'file') || ~exist(dfff,'file')
            warning('%s missing speed/TTL/alignment/dFF. Skipping.\n', exptag);
            continue
        end

        %% Load files
        S_ttl   = load(behTTLf,'dataTTL_beh'); behTTL    = S_ttl.dataTTL_beh;
        S_spd   = load(behSpdf,'behSpeed');    behSpeed  = S_spd.behSpeed;
        S_aln   = load(alignf,'align_beh');    align_beh = S_aln.align_beh;
        S_dff   = load(dfff,'datadFF');        datadFF   = S_dff.datadFF;

        %% Prepare calcium matrix (cells × frames), ensure orientation
        if ~isfield(datadFF,'dFF')
            error('datadFF.dFF not found in %s', dfff);
        end
        F = datadFF.dFF;
        % Orient to cells x frames
        if size(F,1) == numel(align_beh.ca.frameSamp) && size(F,1) > size(F,2)
            F = F.'; % transpose if frames x cells
        end
        good = true(size(F,1),1);
        if isfield(datadFF,'goodCells') && numel(datadFF.goodCells)==size(F,1)
            good = logical(datadFF.goodCells);
        end
        F_good = F(good,:); % only good cells

        %% Behavior domain
        nBehTTL = numel(behTTL.frameInd);
        behFs   = behSpeed.behFs; % must exist; created in speed builder

        %% Map calcium frames to behavior indices (CA→BEH)
        behIdx = align_beh.ca2beh_frameIdx(:);
        nCA    = numel(behIdx);

        % Clamp BEH indices
        behIdx(behIdx < 1)       = 1;
        behIdx(behIdx > nBehTTL) = nBehTTL;

        % Ensure behSpeed covers behavior indices needed by alignment
        if numel(behSpeed.speed_cm_s) < nBehTTL
            padN = nBehTTL - numel(behSpeed.speed_cm_s);
            behSpeed.speed_cm_s = [behSpeed.speed_cm_s(:); nan(padN,1)];
        end

        %% Recompute bouts on full TTL-length speed (handles NaN as non-loco)
        bouts = Felipe_find_locomotion_bouts(behSpeed.speed_cm_s, behFs, optsBout);

        %% Map speed and loco to calcium domain
        speed_ca = behSpeed.speed_cm_s(behIdx);  % N_CA x 1
        loco_ca  = bouts.locoMask(behIdx);       % N_CA x 1 (logical)

        %% Ensure calcium matrix has exactly N_CA columns
        if size(F_good,2) ~= nCA
            if size(F_good,2) > nCA
                F_good = F_good(:,1:nCA);
            else
                F_good = [F_good, nan(size(F_good,1), nCA - size(F_good,2))];
            end
        end

        %% Unconditional CA timebase (FIX for the 'caTime' error)
        caTime = align_beh.ca.frameTime_s(:);           % N_CA x 1
        caFs   = 1 / median(diff(caTime));              % Hz

        % Compute peri window frames ONCE (used by peri plot and exemplars)
        preF   = round(optsPlot.peri_pre_s  * caFs);
        postF  = round(optsPlot.peri_post_s * caFs);
        tperi  = (-preF:postF) / caFs;

        %% -------- Pre‑locomotion stillness window (behavior domain → calcium) --------
        preStillF  = round(optsPlot.preStill_s * behFs);
        preLoc_beh_mask = false(nBehTTL,1);
        for k = 1:numel(bouts.starts_idx)
            onset = bouts.starts_idx(k);
            st = max(1, onset - preStillF);
            if st <= onset-1
                preLoc_beh_mask(st:onset-1) = true;
            end
        end
        preLoc_ca_mask = preLoc_beh_mask(behIdx);
        preLoc_ca_mask = preLoc_ca_mask & ~loco_ca; % ensure pure stillness

        %% -------- Per-cell summaries --------
        nCells = size(F_good,1);
        r_speed   = nan(nCells,1);
        F_still   = nan(nCells,1);
        F_loco    = nan(nCells,1);
        F_preLoc  = nan(nCells,1);

        maskStill = ~loco_ca & isfinite(speed_ca);
        maskLoco  =  loco_ca & isfinite(speed_ca);
        maskPre   =  preLoc_ca_mask & isfinite(speed_ca);

        for c = 1:nCells
            f = F_good(c,:).';  % column

            m = isfinite(f) & isfinite(speed_ca);   % column logical
            if nnz(m) > 5
                r_speed(c) = corr(f(m), speed_ca(m), 'type','Pearson');
            end

            if any(maskStill), F_still(c) = mean(f(maskStill), 'omitnan'); end
            if any(maskLoco),  F_loco(c)  = mean(f(maskLoco),  'omitnan'); end
            if any(maskPre),   F_preLoc(c)= mean(f(maskPre),   'omitnan'); end
        end

        Delta_loco_minus_still = F_loco - F_still;
        Delta_loco_minus_pre   = F_loco - F_preLoc;

        summary = table(find(good), r_speed, F_still, F_preLoc, F_loco, ...
            Delta_loco_minus_still, Delta_loco_minus_pre, ...
            'VariableNames',{'cell','r_F_speed','F_still','F_preLoc','F_loco', ...
                             'DeltaF_loco_minus_still','DeltaF_loco_minus_preLoc'});

        %% -------- Box plot --------
        fig1 = figure('Color','w','Name',['Box ΔF/F still vs loco - ' exptag],'Position',[100 100 680 420]);
        boxplot([F_still(:), F_loco(:)], {'Still','Loco'});
        ylabel('Mean ΔF/F');
        title(sprintf('%s | nCells=%d', exptag, nCells));
        grid on
        saveas(fig1, fullfile(behDir,[exptag '_box_still_vs_loco.png']));

        %% -------- Peri-onset plot (only if enough onsets) --------
        onsets = bouts.starts_idx(:);
        if numel(onsets) >= optsPlot.minOnsets

            onset_t = onsets / behFs;
            onset_ca_idx = round(interp1(caTime, 1:nCA, onset_t, 'nearest','extrap'));

            periCells = nan(nCells, numel(tperi));
            for c = 1:nCells
                f = F_good(c,:).';
                M = nan(numel(onset_ca_idx), numel(tperi));
                for e = 1:numel(onset_ca_idx)
                    i0 = onset_ca_idx(e);
                    idx = (i0-preF):(i0+postF);
                    valid = idx>=1 & idx<=nCA;
                    row = nan(1,numel(tperi));
                    row(valid) = f(idx(valid));
                    M(e,:) = row;
                end
                periCells(c,:) = mean(M,1,'omitnan');
            end

            mu = mean(periCells,1,'omitnan');
            switch lower(optsPlot.shadeType)
                case 'sem'
                    sp = std(periCells,[],1,'omitnan') ./ sqrt(sum(~isnan(periCells),1));
                    lab = 'SEM';
                otherwise
                    sp = std(periCells,[],1,'omitnan');
                    lab = 'SD';
            end

            fig2 = figure('Color','w','Name',['Peri-movement ΔF/F - ' exptag],'Position',[120 120 720 460]);
            hold on
            fill([tperi fliplr(tperi)], [mu-sp fliplr(mu+sp)], [0.6 0.8 1.0], ...
                'FaceAlpha',0.35, 'EdgeColor','none');
            plot(tperi, mu, 'b-', 'LineWidth',2);
            xline(0,'k--');
            xlabel('Time from locomotion onset (s)'); ylabel('ΔF/F');
            title(sprintf('%s | mean ± %s (nCells=%d, nOnsets=%d)', exptag, lab, nCells, numel(onsets)));
            grid on; box on
            saveas(fig2, fullfile(behDir,[exptag '_peri_onset_mean_' lower(optsPlot.shadeType) '.png']));
        end

        %% -------- Exemplars: 5 bouts (speed + ΔF/F traces per cell) --------
        nExamples = min(optsPlot.examples_N, numel(onsets));
        ex = struct([]);
        if nExamples > 0
            % choose bout indices
            switch lower(optsPlot.examples_strategy)
                case 'longest'
                    [~, order] = sort(bouts.duration_s, 'descend');
                    pick = order(1:nExamples);
                otherwise % 'first'
                    pick = 1:nExamples;
            end

            ex = repmat(struct('tCA',[],'speed',[],'Fcells',[],'onset_idx_beh',[]), nExamples, 1);
            for ii = 1:nExamples
                k = pick(ii);
                onset_beh = bouts.starts_idx(k);
                onset_t   = onset_beh / behFs;
                i0 = round(interp1(caTime, 1:nCA, onset_t, 'nearest','extrap'));

                win   = (i0-preF):(i0+postF);
                valid = win>=1 & win<=nCA;

                ex(ii).tCA         = tperi;                 % peri axis in seconds
                ex(ii).onset_idx_beh = onset_beh;
                ex(ii).speed       = nan(1,numel(win));
                ex(ii).speed(valid)= speed_ca(win(valid));

                ex(ii).Fcells      = nan(nCells, numel(win));
                for c = 1:nCells
                    f = F_good(c,:).';
                    ex(ii).Fcells(c,valid) = f(win(valid));
                end
            end

            save(fullfile(behDir,[exptag '_examples.mat']), 'ex','-v7.3');
            % ----------- PLOT EXAMPLES (very minimal) -----------
for ii = 1:numel(ex)
    t   = ex(ii).tCA(:)';
    spd = ex(ii).speed(:)';
    F   = ex(ii).Fcells;      % nCells × T

    figure('Color','w','Name',sprintf('Example %d - %s', ii, exptag));

    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    % ==== SPEED ====
    nexttile;
    plot(t, spd, 'k', 'LineWidth',1.2); hold on;
    xline(0,'r','LineWidth',1.2);
    ylabel('Speed (cm/s)');
    title(sprintf('Bout %d (onset beh idx = %d)', ii, ex(ii).onset_idx_beh));
    grid on;

    % ==== ΔF/F traces stacked ====
    nexttile;
    nCells = size(F,1);
    colors = lines(nCells);

    % Vertical separation between cells
    spread = nanstd(F(:));
    if ~isfinite(spread) || spread==0, spread = 1; end
    dy = 1.5 * spread;

    offset = 0;
    for c = 1:nCells
        plot(t, F(c,:) + offset, 'Color', colors(c,:), 'LineWidth',1); hold on;
        offset = offset + dy;
    end

    xline(0,'r','LineWidth',1.2);
    xlabel('Time (s)');
    ylabel('\DeltaF/F (offset)');
    grid on;
end
        end

        %% -------- Save summary --------
        save(fullfile(behDir,[exptag '_loco_summary.mat']), ...
             'summary','bouts','behSpeed','align_beh','optsBout','optsPlot', '-v7.3');

        fprintf('\nSaved analysis, figures, and exemplars for %s\n', exptag);
    end
end

%% Utility: set defaults for missing fields
function S = setDefaults(S, D)
    f = fieldnames(D);
    for k = 1:numel(f)
        if ~isfield(S,f{k}) || isempty(S.(f{k}))
            S.(f{k}) = D.(f{k});
        end
    end
end