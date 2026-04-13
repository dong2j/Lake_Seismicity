clc; clear; close all

%% --------- SETTINGS ----------
rootFolder = './Raw_Data/GLWS v1.1/Global database of lake water storage GLWS time series v1.1/All';   % change to your folder path, e.g. './Poly' or './ConstantArea'
pattern    = '*.csv';
% If your CSVs have extra non-data lines at the top (metadata),
% set this to the number of lines to skip. Otherwise keep 0.
numHeaderLinesToSkip = 0;

%% --------- FIND FILES ----------
files = dir(fullfile(rootFolder, pattern));
if isempty(files)
    error('No CSV files found in: %s', fullfile(rootFolder, pattern));
end

%% --------- OUTPUT CONTAINER ----------
outDir = './Processed_Data/EachLakeTimeSeries/';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% --------- LOOP ----------
p = gcp('nocreate');
if isempty(p)
    parpool;
end

nFiles = numel(files);

% Preallocate output struct array with fixed indices
LakeSeries(nFiles,1) = struct( ...
    'ID', [], ...
    'Name', "", ...
    'File', "", ...
    'tDec', [], ...
    'level', [], ...
    'storage', [], ...
    'level_extreme', [], ...
    'storage_extreme', [], ...
    'storage_extreme_time', [], ...
    'storage_beta', [], ...
    'storage_p', [] );

parfor k = 1:nFiles
    fprintf('%d of %d processed\n', k, nFiles);

    fpath = fullfile(files(k).folder, files(k).name);
    fname = files(k).name;

    % Create a local result for this iteration
    S = struct( ...
        'ID', [], ...
        'Name', "", ...
        'File', "", ...
        'tDec', [], ...
        'level', [], ...
        'storage', [], ...
        'level_extreme', [], ...
        'storage_extreme', [], ...
        'storage_extreme_time', [], ...
        'storage_beta', [], ...
        'storage_p', [] );

    tok = regexp(fname, '^ID(\d+)([A-Za-z]+)', 'tokens', 'once');
    if isempty(tok)
        warning('Could not parse ID/Name from filename: %s (skipping)', fname);
        LakeSeries(k) = S;
        continue
    end

    lakeID   = str2double(tok{1});
    lakeName = string(tok{2});

    % ---- Read CSV ----
    opts = detectImportOptions(fpath, 'NumHeaderLines', numHeaderLinesToSkip);
    T = readtable(fpath, opts);

    timeRaw    = T.Time;
    lvlRaw     = T.Mean_Levels;
    storageRaw = T.rws;

    % ---- Convert to numeric ----
    timeNum    = toNumeric(timeRaw);
    lvlNum     = toNumeric(lvlRaw);
    storageNum = toNumeric(storageRaw);

    good = isfinite(timeNum) & isfinite(lvlNum) & isfinite(storageNum);
    timeNum    = timeNum(good);
    lvlNum     = lvlNum(good);
    storageNum = storageNum(good);

    if isempty(timeNum)
        warning('No valid data in file: %s (skipping)', fname);
        LakeSeries(k) = S;
        continue
    end

    % ---- Parse YYYYDDD ----
    YYYY = floor(timeNum / 1000);
    DDD  = round(timeNum - YYYY * 1000);

    % ---- Convert to decimal year ----
    tDec = nan(size(timeNum));
    for i = 1:numel(timeNum)
        y = YYYY(i);
        d = DDD(i);
        daysInYear = 365 + isLeapYear(y);
        tDec(i) = y + (d - 1) / daysInYear;
    end

    % ---- Metrics ----
    [max_sto, idx_max] = max(storageNum);
    [min_sto, idx_min] = min(storageNum);

    [betaS, pS] = trend_sen_nullp_from_decimal(tDec, YYYY, DDD, storageNum);

    % ---- Fill local struct ----
    S.ID    = lakeID;
    S.Name  = lakeName;
    S.File  = string(fname);
    S.tDec  = tDec;
    S.level = lvlNum - lvlNum(1);
    S.storage = storageNum;
    S.level_extreme = max(lvlNum) - min(lvlNum);
    S.storage_extreme = max_sto - min_sto;
    S.storage_extreme_time = abs(tDec(idx_max) - tDec(idx_min));
    S.storage_beta = betaS;
    S.storage_p    = pS;

    % ---- Plot ----
    fig = figure('Visible', 'off', 'Color', 'w');

    plot(S.tDec, S.storage, 'k-', 'LineWidth', 1.1);
    hold on;
    plot(S.tDec, (S.tDec - S.tDec(1)) * betaS, 'r-', 'LineWidth', 1.1);

    title(sprintf('Lake %d | %s | p = %.3g', lakeID, lakeName, pS), ...
        'Interpreter', 'none');
    xlabel('Time (year)');
    ylabel('Storage (Gt)');
    set(gca, 'FontSize', 10);

    outFile = fullfile(outDir, sprintf('%d.png', lakeID));
    print(fig, outFile, '-dpng', '-r72');
    close(fig);

    % ---- Assign once to sliced output ----
    LakeSeries(k) = S;
end

% Remove skipped/empty entries if desired
keep = ~cellfun(@isempty, {LakeSeries.ID});
LakeSeries = LakeSeries(keep);

nAll = numel(LakeSeries);
%% --------- SAVE OUTPUT ----------
save('./Processed_Data/LakeSeries_all.mat', 'LakeSeries');
disp('Saved: LakeSeries_all.mat');


inLakeTxt  = './Processed_Data/Processed_Lake.txt';
outLakeTxt = './Processed_Data/Processed_Lake_correct_trend.txt';

Lake_out   = readtable(inLakeTxt);
ID_Lake      = Lake_out.LakeID;
TrendP_Lake  = Lake_out.TrendPval;
Trend_Lake   = Lake_out.Trend;
% Replace Trend and TrendPval using LakeSeries
for i = 1:nAll
    jkf = find(ID_Lake == LakeSeries(i).ID, 1);
    if isempty(jkf)
        continue
    end
    Trend_Lake(jkf)  = LakeSeries(i).storage_beta;
    TrendP_Lake(jkf)= LakeSeries(i).storage_p;
end

% Write updated columns back into table
Lake_out.Trend     = Trend_Lake;
Lake_out.TrendPval = TrendP_Lake;

% (Optional) reorder columns for consistency
Lake_out = Lake_out(:, {'LakeID','LakeName','LakeArea','TypeName', ...
                        'Trend','TrendPval','X','Y'});

% Save as txt
writetable(Lake_out, outLakeTxt, 'Delimiter','\t');

%% =========================
% Helper: plot all series normalized to [-1,1] + median
%% =========================

function x = toNumeric(v)
    if isnumeric(v)
        x = double(v);
    elseif iscell(v)
        % cell array of strings/numbers
        x = nan(size(v));
        for i = 1:numel(v)
            if isnumeric(v{i})
                x(i) = double(v{i});
            else
                x(i) = str2double(string(v{i}));
            end
        end
    else
        % string/categorical/char
        x = str2double(string(v));
    end
end

function tf = isLeapYear(y)
    % Gregorian leap year rule
    tf = (mod(y,4)==0 & mod(y,100)~=0) | (mod(y,400)==0);
end

function [beta_obs, pval] = trend_sen_nullp_from_decimal(tDec, YYYY, DDD, x)
% Trend beta = Sen's slope on DESEASONALIZED monthly medians
% p-value from "null randomization more extreme than observed beta":
%    p = mean( |beta_null| >= |beta_obs| )  (two-sided)
%
% Null randomization: permute YEAR blocks (keeps within-year month pattern, breaks long-term ordering)
% Optional CI: block bootstrap resampling YEAR blocks (with replacement)

    % ---- 1) aggregate to monthly medians (faster + robust) ----
    dt = datetime(YYYY,1,1) + days(DDD-1);
    mo = month(dt);
    yr = YYYY;

    % group by (yr,mo)
    [G, ~] = findgroups(yr(:), mo(:));
    xMon = splitapply(@median, x(:),    G);
    tMon = splitapply(@median, tDec(:), G);
    yrMon = splitapply(@(v)v(1), yr(:), G);
    moMon = splitapply(@(v)v(1), mo(:), G);

    % sort by time
    [tMon, isrt] = sort(tMon);
    xMon  = xMon(isrt);
    yrMon = yrMon(isrt);
    moMon = moMon(isrt);

    % ---- 2) deseasonalize by month (monthly median climatology) ----
    S = nan(12,1);
    for m = 1:12
        ii = (moMon == m);
        S(m) = median(xMon(ii), 'omitnan');
    end
    A = xMon - S(moMon);

    % ---- 3) observed Sen slope ----
    beta_obs = sen_slope(tMon, A);

    % ---- 4) null distribution by permuting YEAR blocks ----
    Bnull = 10000;
    yearsSorted = unique(yrMon, 'sorted');
    nY = numel(yearsSorted);
    beta_null = nan(Bnull,1);

    % Pre-index for speed
    idxByYear = cell(nY,1);
    for k = 1:nY
        idxByYear{k} = find(yrMon == yearsSorted(k));
    end

    for b = 1:Bnull
        permYears = yearsSorted(randperm(nY));  % a permutation (no replacement)
        Aperm = A;

        % map: destination year block gets values from source permYears(k)
        for k = 1:nY
            dstIdx = idxByYear{k};
            srcIdx = find(yrMon == permYears(k));
            % lengths should match if each year has same set of months;
            % but if missing months differ by year, we handle by matching by month:
            if numel(dstIdx) == numel(srcIdx)
                Aperm(dstIdx) = A(srcIdx);
            else
                % robust matching by month within year
                dstMonths = moMon(dstIdx);
                srcMonths = moMon(srcIdx);
                tmp = Aperm(dstIdx);
                for mm = unique(dstMonths(:))'
                    di = dstIdx(dstMonths==mm);
                    si = srcIdx(srcMonths==mm);
                    if ~isempty(di) && ~isempty(si)
                        % if multiple points for same month (rare after monthly median), align by order
                        tmp(dstMonths==mm) = A(si(1:min(numel(si),sum(dstMonths==mm))));
                    end
                end
                Aperm(dstIdx) = tmp;
            end
        end

        beta_null(b) = sen_slope(tMon, Aperm);
    end

    % "more extreme than observed"
    if beta_obs>0
        pval = (sum(beta_null >= beta_obs)) / Bnull;
    else
        pval = (sum(beta_null < beta_obs)) / Bnull;
    end

end

function beta = sen_slope(t, a)
% Median of all pairwise slopes (a_j-a_i)/(t_j-t_i), j>i
% Robust, returns NaN if insufficient data.

    t = t(:); 
    a = a(:);

    n = numel(t);  
    
    DT = t.' - t;          % n x n
    DA = a.' - a;          % n x n

    mask = triu(true(n),1) & (DT ~= 0) & isfinite(DT) & isfinite(DA);
    S = DA(mask) ./ DT(mask);

    beta = median(S, 'omitnan');
end



