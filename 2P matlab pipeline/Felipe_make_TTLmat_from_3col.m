function outFile = Felipe_make_TTLmat_from_3col_EX(exptag, srcFile, outFolder, varargin)
% Create a simple TTL MAT file for Level-0 readframes from a 3-column source:
%   col1 = sample index (2 kHz)
%   col2 = calcium TTL (analog or 0/1)
%   col3 = behavior TTL (stored in meta only)
%
% INPUTS
%   exptag   : string for logging (e.g., 'SD2_260218_01') — not used in naming unless you want
%   srcFile  : FULL path to the 3-column source file (.mat/.csv/.txt)
%   outFolder: folder to save the preprocessed TTL .mat
%
% NAME-VALUE PAIRS (optional):
%   'Fs'          : sampling rate Hz (default x)
%   'ColSample'   : which column is sample index (default 1)
%   'ColCa'       : which column is Ca TTL (default 2)
%   'ColBeh'      : which column is behavior TTL (default 3)
%   'Threshold'   : numeric threshold for Ca TTL (default = auto)
%   'OutBaseName' : custom output base name (default = basename of srcFile)
%   'PlotDebug'   : true/false (default false)
%
% OUTPUT
%   outFile : full path to the saved MAT file containing:
%             ttl (1xN), fs (scalar), meta (struct)
%
% EXAMPLE
%   out = Felipe_make_TTLmat_from_3col_EX('SD2_260218_01', ...
%         'D:\PX1171b\SD2\my_3col_ttl.csv', ...
%         'D:\PX1171b\SD2', ...
%         'Fs',2000,'OutBaseName','TTL_SD2_260218_01','PlotDebug',true);

% -------------------- Parse inputs --------------------
p = inputParser;
p.addRequired('exptag',   @(x)ischar(x) || isstring(x));
p.addRequired('srcFile',  @(x)ischar(x) || isstring(x));
p.addRequired('outFolder',@(x)ischar(x) || isstring(x));
p.addParameter('Fs',           2000,     @(x)isnumeric(x) && isscalar(x) && x>0);
p.addParameter('ColSample',    1,        @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('ColCa',        2,        @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('ColBeh',       3,        @(x)isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('Threshold',    [],       @(x)isnumeric(x) && (isempty(x) || isscalar(x)));
p.addParameter('OutBaseName',  '',       @(x)ischar(x) || isstring(x));
p.addParameter('PlotDebug',    false,    @(x)islogical(x) || isnumeric(x));
p.parse(exptag, srcFile, outFolder, varargin{:});
S = p.Results;

if ~exist(S.srcFile,'file')
    error('Source TTL file not found: %s', S.srcFile);
end
if ~exist(S.outFolder,'dir')
    mkdir(S.outFolder);
end

% -------------------- Load the 3-col matrix --------------------
[~, srcBase, ext] = fileparts(S.srcFile);
switch lower(ext)
    case '.mat'
        tmp = load(S.srcFile);
        vn  = fieldnames(tmp);
        M   = [];
        for vi = 1:numel(vn)
            V = tmp.(vn{vi});
            if isnumeric(V) && ismatrix(V) && size(V,2) >= 3
                M = V; break
            end
        end
        if isempty(M)
            error('No N x >=3 numeric matrix found in MAT file: %s', S.srcFile);
        end
    otherwise
        M = readmatrix(S.srcFile);
        if size(M,2) < 3
            error('TTL file must have >= 3 columns: %s', S.srcFile);
        end
end

% -------------------- Extract columns --------------------
colSample = S.ColSample; colCa = S.ColCa; colBeh = S.ColBeh;

sampleCol = double(M(:, colSample));
ttlCaRaw  = double(M(:, colCa));
ttlBhRaw  = double(M(:, colBeh));

% Sanity
if any(diff(sampleCol) <= 0)
    warning('Sample column not strictly increasing. Proceed with caution.');
end

fs = S.Fs;
t  = sampleCol(:)' / fs;  %#ok<NASGU> % used only for optional plotting

% -------------------- Threshold Ca TTL to 0/1 --------------------
if isempty(S.Threshold)
    lo  = prctile(ttlCaRaw,5);
    hi  = prctile(ttlCaRaw,95);
    thr = lo + 0.8*(hi - lo);
else
    thr = S.Threshold;
end
ttlBinary = ttlCaRaw(:)' > thr;  % row vector logical

% (Optional) quick QC plot
if S.PlotDebug
    figure('Color','w'); 
    subplot(2,1,1); plot(sampleCol/fs, ttlCaRaw, 'k-'); hold on; yline(thr,'r--');
    title(sprintf('Ca TTL (raw) with threshold — %s', char(S.exptag)));
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;
    subplot(2,1,2); plot(sampleCol/fs, double(ttlBinary), 'b-');
    title('Ca TTL (binary)'); xlabel('Time (s)'); ylabel('0/1'); grid on;
    drawnow;
end

% -------------------- Prepare outputs --------------------
ttl  = double(ttlBinary);  % 1xN double 0/1
fs   = fs;                 %#ok<NASGU>

meta = struct();
meta.exptag        = char(S.exptag);
meta.srcFile       = char(S.srcFile);
meta.columns       = struct('sample', colSample, 'ca', colCa, 'beh', colBeh);
meta.fs            = S.Fs;
meta.sampleIndex   = sampleCol;    % keep for traceability
meta.behaviorTTL   = ttlBhRaw;     % stored for reference (not used by Level-0 reader)

% -------------------- Output name & save --------------------
if isempty(S.OutBaseName)
    % default: TTL_<basename_of_src>
    outBase = ['TTL_' srcBase];
else
    outBase = char(S.OutBaseName);
    % ensure it starts with TTL_ for compatibility with the reader (optional)
    if ~startsWith(outBase, 'TTL_')
        outBase = ['TTL_' outBase];
    end
end

outFile = fullfile(S.outFolder, [outBase '.mat']);
save(outFile, 'ttl', 'fs', 'meta', '-v7.3');
fprintf('Saved preprocessed TTL MAT: %s\n', outFile);
end