function out = detect_and_plot_bouts(t, v, varargin)
% DETECT_AND_PLOT_BOUTS
% Detects sleep, low, and high movement bouts from time (t) and speed (v),
% applies merging and minimum duration rules, optional peak filter for high,
% optional manual sleep additions, and (optionally) plots shaded bouts.
%
% INPUTS
%   t : Nx1 time vector (s), strictly increasing (duplicates allowed; handled)
%   v : Nx1 speed vector
%
% NAME-VALUE (optional)
%   'MaxGapSec'    (default 1)   - merge bouts whose gap <= MaxGapSec
%   'MinDurSec'    (default 5)   - drop bouts with duration < MinDurSec
%   'PeakMinHigh'  (default 4)   - keep "high" bouts only if peak >= this
%   'ManualSleep'  (default [])  - Kx2 [on off] seconds to union with sleep
%   'DoPlot'       (default true)
%   'Colors'       (default: sleep=[0.45 0.10 0.65], low=[0.95 0.60 0.10], high=[0.10 0.40 0.90])
%   'Alpha'        (default 0.85)
%   'LineWidth'    (default 1.4)
%
% OUTPUT (struct)
%   out.sleep, out.low, out.high : [on_s off_s] in seconds (after merging/filters)
%   out.params : params used
%   out.axes   : handle to behavior (speed) axes if plotted, [] otherwise
%
% -------------------------------------------------------------------------

% ---- Parse inputs
p = inputParser;
p.addRequired('t', @(x)isvector(x) && isnumeric(x));
p.addRequired('v', @(x)isvector(x) && isnumeric(x));
p.addParameter('MaxGapSec',   1,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinDurSec',   5,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('PeakMinHigh', 4,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('ManualSleep', [],   @(x)isnumeric(x) && (isempty(x) || size(x,2)==2));
p.addParameter('DoPlot',      true, @(x)islogical(x)&&isscalar(x));
defaultColors.sleep = [0.45 0.10 0.65]; % purple
defaultColors.low   = [0.95 0.60 0.10]; % orange
defaultColors.high  = [0.10 0.40 0.90]; % blue
p.addParameter('Colors', defaultColors, @(s)isstruct(s) && all(isfield(s,{'sleep','low','high'})));
p.addParameter('Alpha', 0.85, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('LineWidth', 1.4, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(t, v, varargin{:});
prm = p.Results;

% ---- Clean inputs
t = t(:); v = v(:);
valid = isfinite(t) & isfinite(v);
t = t(valid); v = v(valid);
[t, ia] = unique(t, 'stable');  % enforce monotonic time; drop dup times
v = v(ia);

if numel(t) < 2
    warning('Input time/speed too short after cleaning.');
    out = make_out([],[],[],prm,[]);
    return
end

% ---- Threshold masks (your logic)
sleep_mask = (v < 0.35);
low_mask   = (v >= 0.35) & (v < 2.0);
high_mask  = (v >= 2.0);

% ---- Raw bouts from masks (in time)
sleep_bouts_raw = mask_to_bouts(t, sleep_mask);
low_bouts_raw   = mask_to_bouts(t, low_mask);
high_bouts_raw  = mask_to_bouts(t, high_mask);

% ---- Merge gaps and drop short bouts
sleep_bouts = merge_and_minlen(sleep_bouts_raw, prm.MaxGapSec, prm.MinDurSec);
low_bouts   = merge_and_minlen(low_bouts_raw,   prm.MaxGapSec, prm.MinDurSec);
high_bouts  = merge_and_minlen(high_bouts_raw,  prm.MaxGapSec, prm.MinDurSec);

% ---- Peak criterion for "high"
high_bouts  = keep_if_peak_at_least(high_bouts, t, v, prm.PeakMinHigh);

% ---- Manually union sleep (optional)
if ~isempty(prm.ManualSleep)
    sleep_bouts = union_bouts(sleep_bouts, prm.ManualSleep);
    sleep_bouts = merge_and_minlen(sleep_bouts, 0, prm.MinDurSec); % enforce min dur after union
end

% ---- Optional speed plot
ax = [];
if prm.DoPlot
    figure; set(gcf,'Renderer','opengl');
    ax = axes; hold(ax,'on');
    h = plot(ax, t, v, 'k', 'LineWidth', prm.LineWidth);
    yL = ylim(ax);
    draw_bouts(ax, sleep_bouts, prm.Colors.sleep, prm.Alpha, yL);
    draw_bouts(ax, low_bouts,   prm.Colors.low,   prm.Alpha, yL);
    draw_bouts(ax, high_bouts,  prm.Colors.high,  prm.Alpha, yL);
    uistack(h, 'top');
    xlabel(ax,'Time (s)'); ylabel(ax,'Speed (units/s)');
    title(ax,'Sleep (<0.35), Low (0.35–2), High (≥2 & peak≥4) movement bouts');
    box(ax,'off');
end

% ---- Output
out = make_out(sleep_bouts, low_bouts, high_bouts, prm, ax);

end % ====== main ======


% -------------------- helpers --------------------
function bouts = mask_to_bouts(t, mask)
    mask = logical(mask(:));
    on_idx  = find(diff([false; mask]) == 1);
    off_idx = find(diff([mask; false]) == -1);
    bouts = [t(on_idx) t(off_idx)];
end

function bouts_out = merge_and_minlen(bouts, max_gap_s, min_dur_s)
    if isempty(bouts), bouts_out = bouts; return; end
    bouts = sortrows(bouts,1);
    cur = bouts(1,:); out = zeros(0,2);
    for k = 2:size(bouts,1)
        gap = bouts(k,1) - cur(2);
        if gap <= max_gap_s
            cur(2) = max(cur(2), bouts(k,2));
        else
            out = [out; cur]; %#ok<AGROW>
            cur = bouts(k,:);
        end
    end
    out = [out; cur];
    dur = out(:,2) - out(:,1);
    bouts_out = out(dur >= min_dur_s, :);
end

function bouts_out = keep_if_peak_at_least(bouts, t, v, peak_min)
    if isempty(bouts), bouts_out = bouts; return; end
    keep = false(size(bouts,1),1);
    for i = 1:size(bouts,1)
        seg = (t >= bouts(i,1)) & (t <= bouts(i,2));
        if any(seg) && max(v(seg)) >= peak_min
            keep(i) = true;
        end
    end
    bouts_out = bouts(keep,:);
end

function bouts_u = union_bouts(A, B)
    if isempty(A), bouts_u = B; return; end
    if isempty(B), bouts_u = A; return; end
    C = sortrows([A; B], 1);
    cur = C(1,:); out = zeros(0,2);
    for i = 2:size(C,1)
        if C(i,1) <= cur(2)    % overlap/touch
            cur(2) = max(cur(2), C(i,2));
        else
            out = [out; cur]; %#ok<AGROW>
            cur = C(i,:);
        end
    end
    bouts_u = [out; cur];
end

function draw_bouts(ax, bouts, colorRGB, alphaVal, yL)
    if isempty(bouts), return; end
    for i = 1:size(bouts,1)
        x = [bouts(i,1) bouts(i,2) bouts(i,2) bouts(i,1)];
        y = [yL(1) yL(1) yL(2) yL(2)];
        patch('Parent', ax, 'XData', x, 'YData', y, ...
              'FaceColor', colorRGB, 'FaceAlpha', alphaVal, 'EdgeColor','none');
    end
end

function out = make_out(sleep_b, low_b, high_b, prm, ax)
    out.sleep  = sleep_b;
    out.low    = low_b;
    out.high   = high_b;
    out.params = prm;
    out.axes   = ax;
end