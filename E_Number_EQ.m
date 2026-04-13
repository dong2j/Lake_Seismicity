clc; clear; close all

%% =======================
%  INPUTS
%% =======================
inCatMat    = './Processed_Data/Declustered_Events.txt';
inLakeTxt   = './Processed_Data/Processed_Lake_correct_trend.txt';
p_thre      = 0.1;   % (optional) lake significance threshold
maxDepth_km = 30;

Dis_thres = [10 20 30 40 50 100 150];   % km
t_end     = 2020;    % force curves to extend to 2020

titlestr  = 'Increase + Decrease + Stable'; %'Increase'; %'Decrease';

%% =======================
%  LOAD EARTHQUAKE CATALOG
%% =======================
Cat_out = readtable(inCatMat);
Cat_out = Cat_out(Cat_out.Depth < maxDepth_km,:);

Time_EQ = Cat_out.Time;
Lat_EQ  = Cat_out.Latitude;
Lon_EQ  = Cat_out.Longitude;
Mag_EQ  = Cat_out.Magnitude;
Dep_EQ  = Cat_out.Depth;

% robust total time span for rate calculation:
tYear0   = floor(min(Time_EQ));
tYear1   = floor(max(Time_EQ));
total_yr = tYear1 - tYear0 + 1;   % inclusive years

%% =======================
%  LOAD LAKES
%% =======================
Lake_out    = readtable(inLakeTxt);

Type_Lake   = string(Lake_out.TypeName);
ID_Lake     = Lake_out.LakeID;
TrendP_Lake = Lake_out.TrendPval;
Trend_Lake  = Lake_out.Trend;
X_Lake      = Lake_out.X;   % lon
Y_Lake      = Lake_out.Y;   % lat

% ---- optional: keep only significant lakes ----
%{
jkf = TrendP_Lake < p_thre & Trend_Lake > 0;
Type_Lake = Type_Lake(jkf);
ID_Lake   = ID_Lake(jkf);
X_Lake    = X_Lake(jkf);
Y_Lake    = Y_Lake(jkf);
%}

%% =======================
%  EXPORT 1: Closest lake for EVERY earthquake (no cutoff)
%% =======================
bigR = Inf;
[~, closestDist_km_all, closestLakeID_all, ~] = assignEQtoLakesWithin( ...
    Lat_EQ, Lon_EQ, ID_Lake, Y_Lake, X_Lake, bigR);

Tclosest = table(Time_EQ, Mag_EQ, Lat_EQ, Lon_EQ, Dep_EQ, ...
    double(closestLakeID_all), double(closestDist_km_all), ...
    'VariableNames', {'Time','Magnitude','Latitude','Longitude','Depth', ...
                      'ClosestLakeID','ClosestLakeDist_km'});

writetable(Tclosest, ['./Processed_Data/Declustered_Events_withClosestLake_depth',num2str(maxDepth_km),'.txt'], ...
    'Delimiter','\t');

%% =======================
%  GROUPS: All / Reservoir / Natural
%% =======================
idx_all = true(size(ID_Lake));
idx_res = contains(lower(Type_Lake), "reserv");   % matches Reservoir/reservoir/...
idx_nat = contains(lower(Type_Lake), "natural");  % matches Natural/natural/...

groups(1).name = "All";
groups(1).idx  = idx_all;

groups(2).name = "Reservoir";
groups(2).idx  = idx_res;

groups(3).name = "Natural";
groups(3).idx  = idx_nat;

%% =======================
%  FIGURE 1 + EXPORT 2 (merged): cumulative curves + wide tables
%% =======================
figure('Units','normalized','Position',[0.2 0.2 0.6 0.3],'Color','w');
colo = jet(length(Dis_thres));

for gi = 1:3
    subplot(1,3,gi); hold on

    IDg = ID_Lake(groups(gi).idx);
    Xg  = X_Lake(groups(gi).idx);
    Yg  = Y_Lake(groups(gi).idx);

    outFile = "./Processed_Data/Cumulative_EQ_vs_Time_" + groups(gi).name + "_depth" + num2str(maxDepth_km) + ".txt";

    % One function does BOTH: (i) plot, (ii) export wide table
    plotAndExportCumulativeWide( ...
        Time_EQ, Lat_EQ, Lon_EQ, ...
        IDg, Yg, Xg, ...
        Dis_thres, t_end, ...
        colo, true, outFile);

    title(groups(gi).name);
    xlabel('Time');
    ylabel('Cumulative # EQ');
    grid on;
    set(gca,'FontSize',16,'YScale','log');
    ylim([1 1e3]);
end
sgtitle(titlestr);

%% =======================
%  FIGURE 2: EQ RATE DENSITY vs DISTANCE (cumulative within R)
%  density(R) = N(d_closest <= R) / ( total_yr * (pi*R^2*NL) )
%% =======================
figure('Units','normalized','Position',[0.2 0.55 0.6 0.3],'Color','w');

for gi = 1:3
    subplot(1,3,gi); hold on

    IDg = ID_Lake(groups(gi).idx);
    Xg  = X_Lake(groups(gi).idx);
    Yg  = Y_Lake(groups(gi).idx);
    NLg = numel(IDg);

    % Closest distance from each EQ to the group lakes (km)
    bigR = Inf;
    [~, closestDist_km, ~, ~] = assignEQtoLakesWithin( ...
        Lat_EQ, Lon_EQ, IDg, Yg, Xg, bigR);

    % cumulative counts within each outer radius
    Rvec   = 0:10:150;
    counts = arrayfun(@(r) sum(closestDist_km <= r), Rvec);

    % cumulative "search area": pi R^2 for each lake, summed over lakes
    area_km2 = pi * (Rvec.^2) * NLg;

    % rate density: events / (yr * km^2)
    density = counts ./ area_km2 ./ total_yr;

    bar(Rvec, density);

    title(groups(gi).name);
    xlabel('Distance to closest lake (km)');
    ylabel('EQ rate density (yr^{-1} km^{-2})');
    grid on;
    set(gca,'FontSize',16);
    ylim([0 1e-6]);
end
sgtitle(titlestr);

%% =======================
%  EXPORT: EQs within 20 km of closest lake (ALL lakes)
%% =======================
if maxDepth_km == 70 
    d_export = 20; % km

    IDg = ID_Lake;
    Xg  = X_Lake;
    Yg  = Y_Lake;

    % lake type for each lake (aligned with IDg, Xg, Yg)
    LakeType_g = strings(numel(IDg),1);
    LakeType_g(idx_res) = "Reservoir";
    LakeType_g(idx_nat) = "Natural";

    bigR = Inf;
    [~, closestDist_km, closestLakeID, closestLakeIdx] = assignEQtoLakesWithin( ...
        Lat_EQ, Lon_EQ, IDg, Yg, Xg, bigR);

    in = (closestDist_km <= d_export);

    LakeName_g = string(Lake_out.LakeName);

    Texp = table(Time_EQ(in), Mag_EQ(in), Lat_EQ(in), Lon_EQ(in), Dep_EQ(in), ...
        double(closestDist_km(in)), double(closestLakeID(in)), ...
        LakeName_g(closestLakeIdx(in)), ...
        LakeType_g(closestLakeIdx(in)), ...
        Yg(closestLakeIdx(in)), Xg(closestLakeIdx(in)), ...
        'VariableNames', {'Time','Magnitude','Latitude','Longitude','Depth_km', ...
                        'ClosestDist_km','ClosestLakeID','LakeName','LakeType','LakeLat','LakeLon'});

    Texp.Group = repmat("All", height(Texp), 1);

    fname = "./Processed_Data/EQ_within_" + num2str(d_export) + "km_All.txt";
    writetable(Texp, fname, 'Delimiter','\t');

    %% =================================================
    % Export lake time series for lakes referenced by Texp
    %% =================================================
    load('./Processed_Data/LakeSeries_all.mat');
    % expects struct array LakeSeries_all with fields: ID, tDec, storage

    outDir = './Processed_Data/LakeSeries_export';
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    allIDs   = double(Texp.ClosestLakeID);
    allNames = string(Texp.LakeName);

    [uIDs, ia] = unique(allIDs, 'stable');
    uNames = allNames(ia);

    % --- build a fast lookup of LakeSeries_all IDs ---
    seriesIDs = [LakeSeries.ID];

    fprintf('Exporting %d lake time series to: %s\n', numel(uIDs), outDir);

    for k = 1:numel(uIDs)
        thisID   = uIDs(k);
        thisName = uNames(k);

        j = find(seriesIDs == thisID, 1, 'first');
        if isempty(j)
            warning('Lake ID %g not found in LakeSeries_all. Skipping.', thisID);
            continue
        end

        tDec = LakeSeries(j).tDec(:);
        stor = LakeSeries(j).storage(:);

        safeName = regexprep(thisName, '[^\w\-]+', '_');
        safeName = regexprep(safeName, '_+', '_');
        safeName = regexprep(safeName, '^_|_$', '');

        if strlength(safeName) == 0
            safeName = "UnknownName";
        end

        fname = fullfile(outDir, sprintf('ID%d_%s.txt', thisID, safeName));

        Tout = table(tDec, stor, 'VariableNames', {'tDec','storage'});
        writetable(Tout, fname, 'Delimiter', '\t');
    end
end

%% =======================
%  FUNCTIONS
%% =======================

function [within, closestDist_km, closestLakeID, closestLakeIdx] = ...
    assignEQtoLakesWithin(Lat_EQ, Lon_EQ, LakeID, LakeLat, LakeLon, Dis_thre_km)

    R = 6371.0; % km
    Lat_EQ  = Lat_EQ(:);
    Lon_EQ  = Lon_EQ(:);
    LakeLat = LakeLat(:);
    LakeLon = LakeLon(:);
    LakeID  = LakeID(:);

    NE = numel(Lat_EQ);

    phiL = deg2rad(LakeLat);
    lamL = deg2rad(LakeLon);
    phiE = deg2rad(Lat_EQ);
    lamE = deg2rad(Lon_EQ);

    closestDist_km = nan(NE,1);
    closestLakeID  = nan(NE,1);
    closestLakeIdx = nan(NE,1);

    chunkSize = 5000;
    nChunks = ceil(NE / chunkSize);

    for c = 1:nChunks
        i1 = (c-1)*chunkSize + 1;
        i2 = min(c*chunkSize, NE);
        idx = i1:i2;

        dphi = phiE(idx) - phiL.';   % (nChunk x NL)
        dlam = lamE(idx) - lamL.';   % (nChunk x NL)

        a = sin(dphi/2).^2 + cos(phiE(idx)).*cos(phiL.').*sin(dlam/2).^2;
        d = 2*R*asin(min(1, sqrt(a)));  % km

        [dmin, jmin] = min(d, [], 2);

        closestDist_km(idx) = dmin;
        closestLakeIdx(idx) = jmin;
        closestLakeID(idx)  = LakeID(jmin);
    end

    within = closestDist_km <= Dis_thre_km;

    if ~isinf(Dis_thre_km)
        closestDist_km(~within) = nan;
        closestLakeID(~within)  = nan;
        closestLakeIdx(~within) = nan;
    end
end

function Tout = plotAndExportCumulativeWide( ...
    Time_EQ, Lat_EQ, Lon_EQ, ...
    ID_Lake, Y_Lake, X_Lake, ...
    Dis_thres, t_end, ...
    colo, doPlot, outFile)

    nD  = numel(Dis_thres);
    leg = cell(nD,1);

    % store curves (only one main loop for compute+plot)
    TimeCols = cell(nD,1);
    CumCols  = cell(nD,1);
    maxLen   = 0;

    for ii = 1:nD
        Dis_thre = Dis_thres(ii);

        % ---- compute step curve for this distance (inline; no 2nd function) ----
        [within, ~, ~, ~] = assignEQtoLakesWithin( ...
            Lat_EQ, Lon_EQ, ID_Lake, Y_Lake, X_Lake, Dis_thre);

        tSort = sort(Time_EQ(within));

        if isempty(tSort)
            t_step = [];
            n_step = [];
        else
            nCum = (1:numel(tSort))';

            t_step = nan(2*numel(tSort)-1,1);
            n_step = nan(2*numel(tSort)-1,1);

            t_step(1) = tSort(1);
            n_step(1) = 1;

            for k = 2:numel(tSort)
                t_step(2*k-2) = tSort(k);
                n_step(2*k-2) = nCum(k-1);

                t_step(2*k-1) = tSort(k);
                n_step(2*k-1) = nCum(k);
            end

            t_step = [t_step; t_end];
            n_step = [n_step; n_step(end)];
        end

        TimeCols{ii} = t_step(:);
        CumCols{ii}  = n_step(:);
        maxLen = max(maxLen, numel(t_step));

        if doPlot && ~isempty(t_step)
            semilogy(t_step, n_step, 'Color', colo(ii,:), 'LineWidth', 2);
        end

        leg{ii} = ['Within ' num2str(Dis_thre) ' km'];
    end

    if doPlot
        legend(leg,'Location','northwest');
    end

    % Build wide table (pad with NaN so columns have equal length)
    Tout = table();
    for ii = 1:nD
        t = TimeCols{ii};
        n = CumCols{ii};

        t = [t; nan(maxLen - numel(t), 1)];
        n = [n; nan(maxLen - numel(n), 1)];

        Tout.("Time_"  + Dis_thres(ii) + "km") = t;
        Tout.("CumEQ_" + Dis_thres(ii) + "km") = n;
    end

    writetable(Tout, outFile, 'Delimiter','\t');
end


