function bouts = detect_locomotion_bouts(t, v)
% DETECT_LOCOMOTION_BOUTS  Apply bout rules to behavior speed.
% Inputs:
%   t : behavior time vector (s), length N
%   v : speed vector (cm/s),    length N
% Output:
%   bouts : table with columns
%       on_s, off_s, duration_s, peak_cm_s, peak_time_s, on_idx, off_idx
%
% Rules:
%   - Binary locomotion = (v >= 1 cm/s); rising/falling edges -> on/off
%   - Remove bouts with peak < 10 cm/s
%   - Merge bouts whose gap <= 1 s
%   - Remove merged bouts with duration < 5 s

    % --- sanitize inputs ---
    t = t(:); v = v(:);
    valid = isfinite(t) & isfinite(v);
    t = t(valid); v = v(valid);

    % enforce strictly increasing time (drop duplicates if any)
    [t, ia] = unique(t, 'stable'); 
    v = v(ia);

    if numel(t) < 2
        bouts = empty_table(); 
        return
    end

    % --- 1) threshold at 1 cm/s ---
    locomask = (v >= 2);

    % rising (on) and falling (off) edges
    on_idx  = find(diff([false; locomask]) == 1);
    off_idx = find(diff([locomask; false]) == -1);

    if isempty(on_idx)
        bouts = empty_table(); 
        return
    end

    on_t = t(on_idx);
    off_t = t(off_idx);

    % --- 2) remove bouts with peak < 10 cm/s ---
    keep = true(numel(on_idx),1);
    peak_val  = zeros(numel(on_idx),1);
    peak_time = zeros(numel(on_idx),1);

    for k = 1:numel(on_idx)
        seg = on_idx(k):off_idx(k);
        [peak_val(k), rel] = max(v(seg));
        peak_time(k) = t(seg(rel));
        if peak_val(k) < 4
            keep(k) = false;
        end
    end

    on_idx  = on_idx(keep);
    off_idx = off_idx(keep);
    on_t    = on_t(keep);
    off_t   = off_t(keep);
    peak_val  = peak_val(keep);
    peak_time = peak_time(keep);

    if isempty(on_idx)
        bouts = empty_table(); 
        return
    end

    % --- 3) merge bouts if gap <= 1 s ---
    % We'll sweep and merge consecutive bouts whose start - previous_end <= 1
    m_on_idx  = on_idx(1);
    m_off_idx = off_idx(1);
    m_on_t    = on_t(1);
    m_off_t   = off_t(1);

    % compute peak over the current merged window
    seg = m_on_idx:m_off_idx;
    [m_peak, rel] = max(v(seg));
    m_peak_time = t(seg(rel));

    out_on_idx = []; out_off_idx = [];
    out_on_t   = []; out_off_t   = [];
    out_peak   = []; out_peak_t  = [];

    for k = 2:numel(on_idx)
        gap = on_t(k) - m_off_t;  % time gap between previous end and this start
        if gap <= 1
            % extend current merged bout
            m_off_idx = off_idx(k);
            m_off_t   = off_t(k);
            seg = m_on_idx:m_off_idx;
            [m_peak, rel] = max(v(seg));
            m_peak_time = t(seg(rel));
        else
            % commit previous merged bout
            out_on_idx(end+1,1) = m_on_idx;
            out_off_idx(end+1,1)= m_off_idx;
            out_on_t(end+1,1)   = m_on_t;
            out_off_t(end+1,1)  = m_off_t;
            out_peak(end+1,1)   = m_peak;
            out_peak_t(end+1,1) = m_peak_time;

            % start a new merged bout
            m_on_idx  = on_idx(k);
            m_off_idx = off_idx(k);
            m_on_t    = on_t(k);
            m_off_t   = off_t(k);
            seg = m_on_idx:m_off_idx;
            [m_peak, rel] = max(v(seg));
            m_peak_time = t(seg(rel));
        end
    end
    % commit the last merged bout
    out_on_idx(end+1,1) = m_on_idx;
    out_off_idx(end+1,1)= m_off_idx;
    out_on_t(end+1,1)   = m_on_t;
    out_off_t(end+1,1)  = m_off_t;
    out_peak(end+1,1)   = m_peak;
    out_peak_t(end+1,1) = m_peak_time;

    % --- 4) remove merged bouts with duration < 5 s ---
    dur = out_off_t - out_on_t;
    keep = dur >= 5;

    bouts = table( ...
        out_on_t(keep), out_off_t(keep), dur(keep), ...
        out_peak(keep), out_peak_t(keep), ...
        out_on_idx(keep), out_off_idx(keep), ...
        'VariableNames', {'on_s','off_s','duration_s','peak_cm_s','peak_time_s','on_idx','off_idx'} ...
    );

    % nested helper for empty output
    function T = empty_table()
        T = table( ...
            zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1,'uint32'), zeros(0,1,'uint32'), ...
            'VariableNames', {'on_s','off_s','duration_s','peak_cm_s','peak_time_s','on_idx','off_idx'} ...
        );
    end
end