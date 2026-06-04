function bout = detect_bouts_simple(t, v, varargin)
% DETECT_BOUTS_SIMPLE
% Behavior-only bout detection with your rules. Returns time-intervals.
% - sleep: v < 0.35
% - low  : 0.35 <= v < 2
% - high : v >= 2, then keep only bouts with peak >= PeakMinHigh (default 4)
% - Merge gaps <= MaxGapSec (default 1 s), drop bouts < MinDurSec (default 5 s)
% - Optional manual sleep union.
%
% INPUTS
%   t : Nx1 time (s), numeric (duplicates allowed)
%   v : Nx1 speed (units/s), numeric
%
% NAME-VALUE (optional)
%   'MaxGapSec'    (default 1)
%   'MinDurSec'    (default 5)
%   'PeakMinHigh'  (default 4)
%   'ManualSleep'  (default [], Kx2 [on off] in seconds)
%
% OUTPUT
%   bout.sleep / bout.low / bout.high  --> [on_s off_s] (after merge & filters)
%   bout.params                         --> struct of parameters used

% ---- params
p = inputParser;
p.addRequired('t', @(x)isnumeric(x)&&isvector(x));
p.addRequired('v', @(x)isnumeric(x)&&isvector(x));
p.addParameter('MaxGapSec',   1,  @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinDurSec',   5,  @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('PeakMinHigh', 4,  @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('ManualSleep', [], @(x)isnumeric(x) && (isempty(x) || size(x,2)==2));
p.parse(t, v, varargin{:});
prm = p.Results;

% ---- clean inputs (no asserts)
t = t(:); v = v(:);
ok = isfinite(t) & isfinite(v);
t = t(ok); v = v(ok);

% enforce strictly increasing time by keeping first of duplicates
[ t, keepIdx ] = unique(t, 'stable');
v = v(keepIdx);

if numel(t) < 2
    bout.sleep=[]; bout.low=[]; bout.high=[]; bout.params=prm;
    return
end

% ---- masks (your thresholds)
sleep_mask = (v < 0.3);
low_mask   = (v >= 0.3) & (v < 2.0);
high_mask  = (v >= 2.0);

% ---- raw bouts
sleep_b = mask_to_bouts(t, sleep_mask);
low_b   = mask_to_bouts(t, low_mask);
high_b  = mask_to_bouts(t, high_mask);

% ---- merge & min-duration
sleep_b = merge_and_minlen(sleep_b, prm.MaxGapSec, prm.MinDurSec);
low_b   = merge_and_minlen(low_b,   prm.MaxGapSec, prm.MinDurSec);
high_b  = merge_and_minlen(high_b,  prm.MaxGapSec, prm.MinDurSec);

% ---- keep high bouts with peak >= PeakMinHigh
if ~isempty(high_b) && prm.PeakMinHigh > 0
    keep = false(size(high_b,1),1);
    for i = 1:size(high_b,1)
        seg = (t >= high_b(i,1)) & (t <= high_b(i,2));
        if any(seg) && max(v(seg)) >= prm.PeakMinHigh
            keep(i) = true;
        end
    end
    high_b = high_b(keep,:);
end

% ---- manual sleep union (optional)
if ~isempty(prm.ManualSleep)
    sleep_b = union_bouts(sleep_b, prm.ManualSleep);
    sleep_b = merge_and_minlen(sleep_b, 0, prm.MinDurSec);
end

% ---- output
bout.sleep  = sleep_b;
bout.low    = low_b;
bout.high   = high_b;
bout.params = prm;

end % main


% ---------- helpers ----------
function bouts = mask_to_bouts(t, mask)
    mask = logical(mask(:));
    on  = find(diff([false; mask])==1);
    off = find(diff([mask; false])==-1);
    bouts = [t(on) t(off)];
end

function bouts2 = merge_and_minlen(bouts, max_gap_s, min_dur_s)
    if isempty(bouts), bouts2=bouts; return; end
    bouts = sortrows(bouts,1);
    cur = bouts(1,:); out = zeros(0,2);
    for i = 2:size(bouts,1)
        if (bouts(i,1) - cur(2)) <= max_gap_s
            cur(2) = max(cur(2), bouts(i,2));  % extend
        else
            out = [out; cur]; %#ok<AGROW>
            cur = bouts(i,:);
        end
    end
    out = [out; cur];
    dur = out(:,2) - out(:,1);
    bouts2 = out(dur >= min_dur_s, :);
end

function U = union_bouts(A, B)
    if isempty(A), U=B; return; end
    if isempty(B), U=A; return; end
    C = sortrows([A; B],1);
    cur=C(1,:); out=zeros(0,2);
    for i=2:size(C,1)
        if C(i,1) <= cur(2)     % touch/overlap
            cur(2)=max(cur(2), C(i,2));
        else
            out=[out; cur]; %#ok<AGROW>
            cur=C(i,:);
        end
    end
    U=[out; cur];
end