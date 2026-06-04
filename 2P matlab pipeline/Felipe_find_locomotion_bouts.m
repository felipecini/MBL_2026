function bouts = Felipe_find_locomotion_bouts(speed_cm_s, behFs, optsBout)
% Detect locomotion bouts:
%  1) Threshold at thrOn_cm_s (default 1 cm/s)
%  2) Remove bouts with peak < peakMin_cm_s (default 10 cm/s)
%  3) Merge bouts separated by <= mergeGap_s (default 1 s)
%  4) Remove bouts shorter than minDur_s (default 5 s)

    if nargin < 3 || isempty(optsBout), optsBout = struct(); end
    if ~isfield(optsBout,'thrOn_cm_s'),   optsBout.thrOn_cm_s   = 1;  end
    if ~isfield(optsBout,'peakMin_cm_s'), optsBout.peakMin_cm_s = 10; end
    if ~isfield(optsBout,'mergeGap_s'),   optsBout.mergeGap_s   = 1;  end
    if ~isfield(optsBout,'minDur_s'),     optsBout.minDur_s     = 5;  end

    v = speed_cm_s(:);
    N = numel(v);
    if N==0 || ~isfinite(behFs) || behFs<=0
        bouts = struct('starts_idx',[],'ends_idx',[],'duration_s',[], ...
                       'locoMask',false(0,1),'behFs',behFs,'starts_t_s',[],'ends_t_s',[]);
        return
    end

    isLoco = v >= optsBout.thrOn_cm_s;
    isLoco(~isfinite(v)) = false;

    d = diff([false; isLoco; false]);
    starts = find(d==1);
    ends   = find(d==-1) - 1;

    % peak filter
    keep = true(size(starts));
    for i=1:numel(starts)
        seg = v(starts(i):ends(i));
        if ~any(isfinite(seg)) || max(seg,[], 'omitnan') < optsBout.peakMin_cm_s
            keep(i) = false;
        end
    end
    starts = starts(keep);
    ends   = ends(keep);

    % merge close bouts
    if numel(starts)>1
        gapMax = round(optsBout.mergeGap_s * behFs);
        newS = starts(1);
        newE = ends(1);
        for i=2:numel(starts)
            if (starts(i) - newE(end) - 1) <= gapMax
                newE(end) = ends(i);
            else
                newS(end+1) = starts(i); %#ok<AGROW>
                newE(end+1) = ends(i);   %#ok<AGROW>
            end
        end
        starts = newS(:);
        ends   = newE(:);
    end

    % min duration
    minLen = round(optsBout.minDur_s * behFs);
    durF = (ends - starts + 1);
    keep = durF >= minLen;
    starts = starts(keep);
    ends   = ends(keep);
    durF   = durF(keep);

    % final mask
    locoMask = false(N,1);
    for i=1:numel(starts), locoMask(starts(i):ends(i)) = true; end

    bouts = struct();
    bouts.starts_idx = starts;
    bouts.ends_idx   = ends;
    bouts.duration_s = durF / behFs;
    bouts.locoMask   = locoMask;
    bouts.behFs      = behFs;
    bouts.starts_t_s = starts / behFs;
    bouts.ends_t_s   = ends   / behFs;
end