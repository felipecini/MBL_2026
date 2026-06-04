function out = detect_movement_bouts(t, v, varargin)
% DETECT_MOVEMENT_BOUTS
% Build sleep/low/high movement bouts from time (t) and speed (v), merge gaps,
% enforce minimum duration and (optionally) a peak criterion for "high".
% Returns bouts, transitions (Q↔H), a flat events table, and (optionally) plots
% the shaded labels over the speed trace. Can also save MAT/CSV outputs.
%
% INPUTS
%   t : Nx1 time (s)              -- strictly increasing (duplicates ok; handled)
%   v : Nx1 speed (units/s)       -- numeric vector, same length as t
%
% NAME-VALUE PARAMS (all optional, flexible)
%   'SleepMax'     (default 0.5)     : speed < SleepMax          → sleep mask
%   'LowRange'     (default [0.5 1]) : LowRange(1) ≤ v < LowRange(2) → low mask
%   'HighMin'      (default 1.0)     : v ≥ HighMin               → high mask
%   'PeakMinHigh'  (default 10)      : keep high bouts only if peak ≥ PeakMinHigh
%   'MaxGapSec'    (default 1)       : merge bouts if gap ≤ MaxGapSec
%   'MinDurSec'    (default 5)       : drop bouts shorter than this after merging
%   'ManualSleep'  (default [])      : K×2 [on off] (s) to union with sleep
%   'DoPlot'       (default false)   : plot shaded bouts over v(t)
%   'Alpha'        (default 0.90)    : patch opacity (0–1)
%   'Colors'       (struct with fields sleep/low/high; dark defaults provided)
%   'LineWidth'    (default 1.4)     : speed line width
%   'SaveBase'     (default '')      : if non-empty, save MAT + CSV sidecars
%
% OUTPUT (struct)
%   out.sleep/low/high : table with columns:
%        on_s, off_s, duration_s, on_idx, off_idx, peak_value, peak_time_s
%   out.transitions    : struct with fields:
%        Q2H_onsets_s, H2Q_onsets_s, Q2H_idx, H2Q_idx
%   out.events_table   : table with rows of atomic events:
%        t_s, idx, label  (labels like 'sleep_on','high_off','Q2H_onset', etc.)
%   out.params         : the final parameter values
%   out.axes           : axes handle (if DoPlot=true), else []
%
% -------------------------------------------------------------------------

% ---------- Parse inputs ----------
p = inputParser;
p.addRequired('t', @(x)isvector(x) && isnumeric(x));
p.addRequired('v', @(x)isvector(x) && isnumeric(x));

p.addParameter('SleepMax',    0.5,    @(x)isnumeric(x)&&isscalar(x));
p.addParameter('LowRange',    [0.5 1],@(x)isnumeric(x)&&numel(x)==2&&x(1)<x(2));
p.addParameter('HighMin',     1.0,    @(x)isnumeric(x)&&isscalar(x));
p.addParameter('PeakMinHigh', 10,     @(x)isnumeric(x)&&isscalar(x)&&x>=0);

p.addParameter('MaxGapSec',   1,      @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinDurSec',   5,      @(x)isnumeric(x)&&isscalar(x)&&x>=0);

p.addParameter('ManualSleep', [],     @(x)isnumeric(x) && (isempty(x) || size(x,2)==2));

p.addParameter('DoPlot',      false,  @(x)islogical(x)&&isscalar(x));
p.addParameter('Alpha',       0.90,   @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
C.sleep = [0.45 0.10 0.65];  C.low = [0.95 0.60 0.10]; C.high = [0.10 0.40 0.90];
p.addParameter('Colors',      C,      @(s)isstruct(s) && all(isfield(s,{'sleep','low','high'})));
p.addParameter('LineWidth',   1.4,    @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.addParameter('SaveBase',    '',     @(x)ischar(x) || isstring(x));

p.parse(t, v, varargin{:});
prm = p.Results;

% ---------- Sanitize / normalize ----------
t = t(:); v = v(:);
valid = isfinite(t) & isfinite(v);
t = t(valid); v = v(valid);

% strictly increasing time
[t, ia] = unique(t, 'stable');
v = v(ia);

if numel(t) < 2
    warning('Time vector too short after cleaning.');
    out = pack_out(empty_bout_table(), empty_bout_table(), empty_bout_table(), ...
                   struct('Q2H_onsets_s',[],'H2Q_onsets_s',[],'Q2H_idx',[],'H2Q_idx',[]), ...
                   table(), prm, []);
    return;
end

% ---------- Build masks ----------
sleep_mask = (v < prm.SleepMax);
low_mask   = (v >= prm.LowRange(1)) & (v < prm.LowRange(2));
high_mask  = (v >= prm.HighMin);

% ---------- Extract raw bouts from masks (by index) ----------
sleep_idx = mask_to_idx_blocks(sleep_mask);
low_idx   = mask_to_idx_blocks(low_mask);
high_idx  = mask_to_idx_blocks(high_mask);

% ---------- Merge gaps & drop short ----------
sleep_idx = merge_and_minlen_by_idx(t, sleep_idx, prm.MaxGapSec, prm.MinDurSec);
low_idx   = merge_and_minlen_by_idx(t, low_idx,   prm.MaxGapSec, prm.MinDurSec);
high_idx  = merge_and_minlen_by_idx(t, high_idx,  prm.MaxGapSec, prm.MinDurSec);

% ---------- Peak filter for "high" bouts ----------
if ~isempty(high_idx) && prm.PeakMinHigh > -inf
    keep = false(size(high_idx,1),1);
    for i = 1:size(high_idx,1)
        seg = high_idx(i,1):high_idx(i,2);
        [pk, pk_rel] = max(v(seg));
        if pk >= prm.PeakMinHigh, keep(i) = true; end %#ok<*AGROW>
    end
    high_idx = high_idx(keep,:);
end

% ---------- Manual sleep union (optional) ----------
if ~isempty(prm.ManualSleep)
    sleep_idx = union_idx_blocks(sleep_idx, time_to_idx_blocks(t, prm.ManualSleep));
    % re-apply min duration after union (no extra gap merge)
    sleep_idx = drop_short_by_idx(t, sleep_idx, prm.MinDurSec);
end

% ---------- Convert idx blocks to tables with times & peaks ----------
sleep_tbl = idx_blocks_to_table('sleep', t, v, sleep_idx);
low_tbl   = idx_blocks_to_table('low',   t, v, low_idx);
high_tbl  = idx_blocks_to_table('high',  t, v, high_idx);

% ---------- Transitions (quiescence ↔ high) by edges ----------
% define quiescence as "sleep"; you can change to (sleep|low) if desired
Q_mask = false(size(t)); Q_mask(sleep_idx_to_linear(sleep_idx)) = true;
H_mask = false(size(t)); H_mask(high_idx_to_linear(high_idx))   = true;

Q2H_idx = find(Q_mask(1:end-1)==true & H_mask(2:end)==true) + 1;
H2Q_idx = find(H_mask(1:end-1)==true & Q_mask(2:end)==true) + 1;

transitions.Q2H_onsets_s = t(Q2H_idx);
transitions.H2Q_onsets_s = t(H2Q_idx);
transitions.Q2H_idx      = Q2H_idx;
transitions.H2Q_idx      = H2Q_idx;

% ---------- Flat events table (sorted) ----------
events_table = make_events_table(t, sleep_tbl, low_tbl, high_tbl, transitions);

% ---------- Optional plot ----------
ax = [];
if prm.DoPlot
    figure; set(gcf,'Renderer','opengl');
    ax = axes; hold(ax,'on');
    h = plot(ax, t, v, 'k', 'LineWidth', prm.LineWidth);
    yL = ylim(ax);
    draw_patches_from_idx(ax, t, sleep_idx, prm.Colors.sleep, prm.Alpha, yL);
    draw_patches_from_idx(ax, t, low_idx,   prm.Colors.low,   prm.Alpha, yL);
    draw_patches_from_idx(ax, t, high_idx,  prm.Colors.high,  prm.Alpha, yL);
    uistack(h,'top');
    xlabel(ax,'Time (s)'); ylabel(ax,'Speed');
    title(ax, sprintf('Sleep < %.3g, Low in [%.3g, %.3g), High ≥ %.3g (peak ≥ %.3g)', ...
        prm.SleepMax, prm.LowRange(1), prm.LowRange(2), prm.HighMin, prm.PeakMinHigh));
    box(ax,'off');
end

% ---------- Package & (optionally) save ----------
out = pack_out(sleep_tbl, low_tbl, high_tbl, transitions, events_table, prm, ax);

if ~isempty(prm.SaveBase)
    base = char(prm.SaveBase);
    save([base '_bouts.mat'], 'out', '-v7.3');

    % Bouts CSV
    writetable(add_bout_type(sleep_tbl,'sleep'), [base '_bouts_sleep.csv']);
    writetable(add_bout_type(low_tbl,  'low'),   [base '_bouts_low.csv']);
    writetable(add_bout_type(high_tbl, 'high'),  [base '_bouts_high.csv']);

    % Transitions CSV
    Ttr = table( ...
        transitions.Q2H_onsets_s(:), transitions.H2Q_onsets_s(:), ...
        'VariableNames', {'Q2H_onset_s','H2Q_onset_s'});
    writetable(Ttr, [base '_transitions.csv']);

    % All events CSV
    writetable(events_table, [base '_events.csv']);
end

end % ===== END MAIN FUNCTION =====


% --------------------- Local helpers ---------------------

function tbl = empty_bout_table()
    tbl = table( ...
        zeros(0,1), zeros(0,1), zeros(0,1), ...
        zeros(0,1,'uint32'), zeros(0,1,'uint32'), ...
        zeros(0,1), zeros(0,1), ...
        'VariableNames', {'on_s','off_s','duration_s','on_idx','off_idx','peak_value','peak_time_s'});
end

function out = pack_out(sleep_tbl, low_tbl, high_tbl, transitions, events_table, prm, ax)
    out.sleep        = sleep_tbl;
    out.low          = low_tbl;
    out.high         = high_tbl;
    out.transitions  = transitions;
    out.events_table = events_table;
    out.params       = prm;
    out.axes         = ax;
end

function A = add_bout_type(tbl, typ)
    if isempty(tbl)
        A = tbl;
    else
        A = addvars(tbl, repmat(string(typ),height(tbl),1), 'Before',1, 'NewVariableNames','bout_type');
    end
end

function blocks = mask_to_idx_blocks(mask)
    mask = logical(mask(:));
    on  = find(diff([false; mask])==1);
    off = find(diff([mask; false])==-1);
    blocks = [on off];
end

function blocks = time_to_idx_blocks(t, bouts_time)
    if isempty(bouts_time), blocks = zeros(0,2); return; end
    on_idx  = arrayfun(@(x)find_first_ge(t,x), bouts_time(:,1));
    off_idx = arrayfun(@(x)find_last_le(t,x),  bouts_time(:,2));
    blocks  = [on_idx off_idx];
    blocks  = blocks(blocks(:,1) <= blocks(:,2), :); % keep valid
end

function idx = find_first_ge(t, x)
    k = find(t >= x, 1, 'first'); if isempty(k), k = numel(t); end; idx = k;
end
function idx = find_last_le(t, x)
    k = find(t <= x, 1, 'last');  if isempty(k), k = 1;        end; idx = k;
end

function blocks = merge_and_minlen_by_idx(t, blocks, max_gap_s, min_dur_s)
    if isempty(blocks), return; end
    % sort by start
    blocks = sortrows(blocks,1);
    out = [];
    cur = blocks(1,:);
    for i = 2:size(blocks,1)
        gap = t(blocks(i,1)) - t(cur(2));
        if gap <= max_gap_s
            cur(2) = max(cur(2), blocks(i,2)); % extend
        else
            out = [out; cur]; %#ok<AGROW>
            cur = blocks(i,:);
        end
    end
    out = [out; cur];

    % min duration filter
    dur = t(out(:,2)) - t(out(:,1));
    blocks = out(dur >= min_dur_s, :);
end

function blocks = union_idx_blocks(A, B)
    if isempty(A), blocks = B; return; end
    if isempty(B), blocks = A; return; end
    C = [A; B];
    C = sortrows(C,1);
    out = C(1,:);
    for i = 2:size(C,1)
        if C(i,1) <= out(end,2)+1 % touching/overlapping in index space
            out(end,2) = max(out(end,2), C(i,2));
        else
            out = [out; C(i,:)]; %#ok<AGROW>
        end
    end
    blocks = out;
end

function blocks = drop_short_by_idx(t, blocks, min_dur_s)
    if isempty(blocks), return; end
    dur = t(blocks(:,2)) - t(blocks(:,1));
    blocks = blocks(dur >= min_dur_s, :);
end

function tbl = idx_blocks_to_table(kind, t, v, blocks)
    if isempty(blocks), tbl = empty_bout_table(); return; end
    n = size(blocks,1);
    on_idx  = blocks(:,1);
    off_idx = blocks(:,2);
    on_s    = t(on_idx);
    off_s   = t(off_idx);
    duration_s = off_s - on_s;
    peak_value  = zeros(n,1);
    peak_time_s = zeros(n,1);
    for i = 1:n
        seg = on_idx(i):off_idx(i);
        [pk, r] = max(v(seg));
        peak_value(i)  = pk;
        peak_time_s(i) = t(seg(r));
    end
    tbl = table(on_s, off_s, duration_s, uint32(on_idx), uint32(off_idx), peak_value, peak_time_s, ...
        'VariableNames', {'on_s','off_s','duration_s','on_idx','off_idx','peak_value','peak_time_s'});
end

function lin = sleep_idx_to_linear(blocks)
    lin = [];
    for i = 1:size(blocks,1)
        lin = [lin, blocks(i,1):blocks(i,2)]; %#ok<AGROW>
    end
    lin = unique(lin);
end
function lin = high_idx_to_linear(blocks)
    lin = sleep_idx_to_linear(blocks); % same logic
end

function draw_patches_from_idx(ax, t, blocks, col, alpha, yL)
    if isempty(blocks), return; end
    for i = 1:size(blocks,1)
        x1 = t(blocks(i,1)); x2 = t(blocks(i,2));
        patch('Parent',ax, ...
              'XData',[x1 x2 x2 x1], 'YData',[yL(1) yL(1) yL(2) yL(2)], ...
              'FaceColor',col, 'FaceAlpha',alpha, 'EdgeColor','none');
    end
end

function T = make_events_table(t, sleep_tbl, low_tbl, high_tbl, transitions)
    % Build a simple chronological list of events with labels
    ev_t = []; ev_idx = []; ev_lab = strings(0,1);

    append_bouts = @(tbl, label_on, label_off) ...
        [tbl.on_s, tbl.on_idx, repmat(string(label_on),height(tbl),1); ...
         tbl.off_s, tbl.off_idx, repmat(string(label_off),height(tbl),1)];

    S = append_bouts(sleep_tbl,'sleep_on','sleep_off');
    L = append_bouts(low_tbl,  'low_on','low_off');
    H = append_bouts(high_tbl, 'high_on','high_off');

    ev = [S; L; H; ...
          transitions.Q2H_onsets_s(:), transitions.Q2H_idx(:), repmat("Q2H_onset", numel(transitions.Q2H_idx),1); ...
          transitions.H2Q_onsets_s(:), transitions.H2Q_idx(:), repmat("H2Q_onset", numel(transitions.H2Q_idx),1)];

    if isempty(ev)
        T = table([], [], strings(0,1), 'VariableNames', {'t_s','idx','label'});
    else
        [~, order] = sort(ev(:,1));
        ev_sorted = ev(order,:);
        T = table(ev_sorted(:,1), uint32(ev_sorted(:,2)), string(ev_sorted(:,3)), ...
            'VariableNames', {'t_s','idx','label'});
    end
end
``