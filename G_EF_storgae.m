clc; clear; close all

%% =======================
%  INPUTS
%% =======================
inCatMat    = './Processed_Data/Declustered_Events.txt';
inLakeTxt   = './Processed_Data/Processed_Lake_correct_trend.txt';
inLakeSeriesMat = './Processed_Data/LakeSeries_all.mat';  % expects LakeSeries (struct) or LakeSeries_all

%% =======================
%  SETTINGS
%% =======================
Dis_thres = 20:5:150;      % dmax (km), cumulative [0, dmax]
dtbin     = 1/365.25;      % years (~daily)

% EQ filter
maxDepth_km = 30;

% Jackknife over lakes: omit 10% lakes each realization
nJack    = 10000;
omitFrac = 0.10;

% Null model: shuffle EQ times (locations fixed)
nNull = 10000;

% Parallel toggle
useParallel = true;

%% =======================
%  LOAD EARTHQUAKE CATALOG
%% =======================
EQ = readtable(inCatMat);
EQ = EQ(EQ.Depth < maxDepth_km,:);

Time_EQ = EQ.Time;        % decimal year assumed
Mag_EQ  = EQ.Magnitude;   %#ok<NASGU>
Lat_EQ  = EQ.Latitude;
Lon_EQ  = EQ.Longitude;
Dep_EQ  = EQ.Depth;       %#ok<NASGU>

NE = numel(Lat_EQ);

%% =======================
%  LOAD LAKES (METADATA)
%% =======================
Lake_out = readtable(inLakeTxt);

Type_Lake = Lake_out.TypeName;
isRe = strcmp(Type_Lake,'Reservoir');
isNa = strcmp(Type_Lake,'Natural lake');

Group(1).name = 'All';       Group(1).mask = true(height(Lake_out),1);
Group(2).name = 'Reservoir'; Group(2).mask = isRe;
Group(3).name = 'Natural';   Group(3).mask = isNa;

%% =======================
%  LOAD LAKE TIME SERIES (GLWS)
%% =======================
S = load(inLakeSeriesMat);

% accept either LakeSeries or LakeSeries_all
if isfield(S,'LakeSeries')
    LakeSeries = S.LakeSeries;
elseif isfield(S,'LakeSeries_all')
    LakeSeries = S.LakeSeries_all;
else
    error('Cannot find LakeSeries or LakeSeries_all in %s', inLakeSeriesMat);
end

seriesIDs = [LakeSeries.ID];
[tf, loc] = ismember(Lake_out.LakeID, seriesIDs);
LakeRow_to_SeriesIndex = nan(height(Lake_out),1);
LakeRow_to_SeriesIndex(tf) = loc(tf);

%% =======================
%  PRECOMPUTE PER-LAKE THRESHOLD + HIGH/LOW DURATION + OBS WINDOW
%  ALSO STORE A GRIDDED SERIES FOR FAST INTERPOLATION LATER
%% =======================
nL = height(Lake_out);
Sthr_lake = nan(nL,1);
TH_lake   = nan(nL,1);    % years
TL_lake   = nan(nL,1);    % years
tmin_lake = nan(nL,1);
tmax_lake = nan(nL,1);

tgrid_cell = cell(nL,1);
sgrid_cell = cell(nL,1);

for j = 1:nL
    sidx = LakeRow_to_SeriesIndex(j);
    if ~isfinite(sidx), continue; end

    time    = LakeSeries(sidx).tDec(:);
    storage = LakeSeries(sidx).storage(:);

    if numel(time) < 2, continue; end
    good0 = isfinite(time) & isfinite(storage);
    time = time(good0); storage = storage(good0);
    if numel(time) < 2, continue; end

    % ensure monotonic time for interp
    [time, isort] = sort(time);
    storage = storage(isort);

    tmin_lake(j) = time(1);
    tmax_lake(j) = time(end);

    tgrid = (time(1):dtbin:time(end))';
    sgrid = interp1(time, storage, tgrid, 'linear', NaN);

    goodg = isfinite(sgrid);
    if sum(goodg) < 2, continue; end

    Sthr = median(sgrid(goodg));

    Sthr_lake(j) = Sthr;
    TH_lake(j)   = sum(sgrid(goodg) >  Sthr) * dtbin;
    TL_lake(j)   = sum(sgrid(goodg) <= Sthr) * dtbin;

    tgrid_cell{j} = tgrid;
    sgrid_cell{j} = sgrid;
end

%% =======================
%  PRECOMPUTE DISTANCE MATRICES (ONCE) PER GROUP
%% =======================
fprintf('\nPrecomputing EQ-to-lake distance matrices (main cost, done once)...\n');

D_all = haversineDistMatrix_single(Lat_EQ, Lon_EQ, Lake_out.Y,                Lake_out.X);
D_re  = haversineDistMatrix_single(Lat_EQ, Lon_EQ, Lake_out.Y(Group(2).mask), Lake_out.X(Group(2).mask));
D_na  = haversineDistMatrix_single(Lat_EQ, Lon_EQ, Lake_out.Y(Group(3).mask), Lake_out.X(Group(3).mask));

% Group column -> global row index in Lake_out
idxAll_global = (1:height(Lake_out))';
idxRe_global  = find(Group(2).mask);
idxNa_global  = find(Group(3).mask);

%% =======================
%  OBSERVED EF(dmax): EF_storage = lambda_high / lambda_low
%  IMPORTANT: NO "CLOSEST LAKE" ASSIGNMENT.
%  Each EQ contributes to EVERY lake within dmax.
%% =======================
nD = numel(Dis_thres);

EF_obs_all = nan(nD,1);
EF_obs_re  = nan(nD,1);
EF_obs_na  = nan(nD,1);

fprintf('Computing observed EF(dmax)...\n');
for ii = 1:nD
    dmax = Dis_thres(ii);

    EF_obs_all(ii) = EF_storage_allLakesWithin( ...
        dmax, D_all, Time_EQ, idxAll_global, ...
        Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
        tgrid_cell, sgrid_cell);

    EF_obs_re(ii)  = EF_storage_allLakesWithin( ...
        dmax, D_re, Time_EQ, idxRe_global, ...
        Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
        tgrid_cell, sgrid_cell);

    EF_obs_na(ii)  = EF_storage_allLakesWithin( ...
        dmax, D_na, Time_EQ, idxNa_global, ...
        Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
        tgrid_cell, sgrid_cell);
end

%% =======================
%  START PARPOOL (optional)
%% =======================
if useParallel
    p = gcp('nocreate');
    if isempty(p), parpool; end
end

%% =======================
%  JACKKNIFE: omit 10% lakes each realization
%% =======================
fprintf('Jackknife %d realizations (omit %.0f%% lakes)...\n', nJack, 100*omitFrac);

EF_all_j = nan(nJack, nD);
EF_re_j  = nan(nJack, nD);
EF_na_j  = nan(nJack, nD);

nAll = size(D_all,2);
nRe  = size(D_re,2);
nNa  = size(D_na,2);

if useParallel
    parfor r = 1:nJack
        % --- sample kept lakes (indices are within each group-matrix) ---
        nOmit_all = round(omitFrac * nAll);
        idxKeep_all = setdiff(1:nAll, randperm(nAll, nOmit_all));

        nOmit_re = round(omitFrac * nRe);
        idxKeep_re = setdiff(1:nRe, randperm(nRe, nOmit_re));

        nOmit_na = round(omitFrac * nNa);
        idxKeep_na = setdiff(1:nNa, randperm(nNa, nOmit_na));

        % --- kept global lake indices and kept distance matrices ---
        keepAll_global = idxAll_global(idxKeep_all);
        keepRe_global  = idxRe_global(idxKeep_re);
        keepNa_global  = idxNa_global(idxKeep_na);

        D_all_keep = D_all(:, idxKeep_all);
        D_re_keep  = D_re(:,  idxKeep_re);
        D_na_keep  = D_na(:,  idxKeep_na);

        tmp_all = nan(1,nD);
        tmp_re  = nan(1,nD);
        tmp_na  = nan(1,nD);

        for ii = 1:nD
            dmax = Dis_thres(ii);

            tmp_all(ii) = EF_storage_allLakesWithin( ...
                dmax, D_all_keep, Time_EQ, keepAll_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);

            tmp_re(ii)  = EF_storage_allLakesWithin( ...
                dmax, D_re_keep, Time_EQ, keepRe_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);

            tmp_na(ii)  = EF_storage_allLakesWithin( ...
                dmax, D_na_keep, Time_EQ, keepNa_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);
        end

        EF_all_j(r,:) = tmp_all;
        EF_re_j(r,:)  = tmp_re;
        EF_na_j(r,:)  = tmp_na;
    end
else
    for r = 1:nJack
        nOmit_all = round(omitFrac * nAll);
        idxKeep_all = setdiff(1:nAll, randperm(nAll, nOmit_all));

        nOmit_re = round(omitFrac * nRe);
        idxKeep_re = setdiff(1:nRe, randperm(nRe, nOmit_re));

        nOmit_na = round(omitFrac * nNa);
        idxKeep_na = setdiff(1:nNa, randperm(nNa, nOmit_na));

        keepAll_global = idxAll_global(idxKeep_all);
        keepRe_global  = idxRe_global(idxKeep_re);
        keepNa_global  = idxNa_global(idxKeep_na);

        D_all_keep = D_all(:, idxKeep_all);
        D_re_keep  = D_re(:,  idxKeep_re);
        D_na_keep  = D_na(:,  idxKeep_na);

        for ii = 1:nD
            dmax = Dis_thres(ii);

            EF_all_j(r,ii) = EF_storage_allLakesWithin( ...
                dmax, D_all_keep, Time_EQ, keepAll_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);

            EF_re_j(r,ii)  = EF_storage_allLakesWithin( ...
                dmax, D_re_keep, Time_EQ, keepRe_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);

            EF_na_j(r,ii)  = EF_storage_allLakesWithin( ...
                dmax, D_na_keep, Time_EQ, keepNa_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);
        end
    end
end

CI_ratio_all = prctile(EF_all_j, [5 95], 1).';
CI_ratio_re  = prctile(EF_re_j,  [5 95], 1).';
CI_ratio_na  = prctile(EF_na_j,  [5 95], 1).';

%% =======================
%  NULL MODEL: shuffle EQ times (locations fixed)
%% =======================
fprintf('Null model %d realizations (shuffle EQ times)...\n', nNull);

EF_all_null = nan(nNull, nD);
EF_re_null  = nan(nNull, nD);
EF_na_null  = nan(nNull, nD);

if useParallel
    parfor r = 1:nNull
        Time_shuf = Time_EQ(randperm(NE));

        tmp_all = nan(1,nD);
        tmp_re  = nan(1,nD);
        tmp_na  = nan(1,nD);

        for ii = 1:nD
            dmax = Dis_thres(ii);

            tmp_all(ii) = EF_storage_allLakesWithin( ...
                dmax, D_all, Time_shuf, idxAll_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);

            tmp_re(ii)  = EF_storage_allLakesWithin( ...
                dmax, D_re, Time_shuf, idxRe_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);

            tmp_na(ii)  = EF_storage_allLakesWithin( ...
                dmax, D_na, Time_shuf, idxNa_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);
        end

        EF_all_null(r,:) = tmp_all;
        EF_re_null(r,:)  = tmp_re;
        EF_na_null(r,:)  = tmp_na;
    end
else
    for r = 1:nNull
        Time_shuf = Time_EQ(randperm(NE));

        for ii = 1:nD
            dmax = Dis_thres(ii);

            EF_all_null(r,ii) = EF_storage_allLakesWithin( ...
                dmax, D_all, Time_shuf, idxAll_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);

            EF_re_null(r,ii)  = EF_storage_allLakesWithin( ...
                dmax, D_re, Time_shuf, idxRe_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);

            EF_na_null(r,ii)  = EF_storage_allLakesWithin( ...
                dmax, D_na, Time_shuf, idxNa_global, ...
                Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
                tgrid_cell, sgrid_cell);
        end
    end
end

CI_null_all = prctile(EF_all_null, [5 95], 1).';
CI_null_re  = prctile(EF_re_null,  [5 95], 1).';
CI_null_na  = prctile(EF_na_null,  [5 95], 1).';

%% =======================
%  PLOT (null band + jackknife CI + observed)
%% =======================
col_all  = [150  75   0] / 255;   % brown
col_re   = [120  60 160] / 255;   % purple
col_na   = [ 60 150  80] / 255;   % green
col_null = [0.7 0.7 0.7];

alpha_CI   = 0.35;
alpha_null = 0.45;

figure('Units','normalized','Position',[0.18 0.2 0.72 0.32]);

subplot(1,3,1); hold on;
shadedCI(Dis_thres, CI_null_all(:,1), CI_null_all(:,2), col_null, alpha_null);
shadedCI(Dis_thres, CI_ratio_all(:,1), CI_ratio_all(:,2), col_all, alpha_CI);
plot(Dis_thres, EF_obs_all, '-', 'Color', col_all, 'LineWidth', 1.8);
yline(1,'--','Color',[0.4 0.4 0.4]); grid on;
xlabel('d_{max} (km)'); ylabel('EF_{storage}=\lambda_{high}/\lambda_{low}');
title('All'); ylim([0 4]); set(gca,'FontSize',14);

subplot(1,3,2); hold on;
shadedCI(Dis_thres, CI_null_re(:,1), CI_null_re(:,2), col_null, alpha_null);
shadedCI(Dis_thres, CI_ratio_re(:,1), CI_ratio_re(:,2), col_re, alpha_CI);
plot(Dis_thres, EF_obs_re, '-', 'Color', col_re, 'LineWidth', 1.8);
yline(1,'--','Color',[0.4 0.4 0.4]); grid on;
xlabel('d_{max} (km)'); ylabel('EF_{storage}');
title('Reservoir'); ylim([0 4]); set(gca,'FontSize',14);

subplot(1,3,3); hold on;
shadedCI(Dis_thres, CI_null_na(:,1), CI_null_na(:,2), col_null, alpha_null);
shadedCI(Dis_thres, CI_ratio_na(:,1), CI_ratio_na(:,2), col_na, alpha_CI);
plot(Dis_thres, EF_obs_na, '-', 'Color', col_na, 'LineWidth', 1.8);
yline(1,'--','Color',[0.4 0.4 0.4]); grid on;
xlabel('d_{max} (km)'); ylabel('EF_{storage}');
title('Natural lake'); ylim([0 4]); set(gca,'FontSize',14);

%% =======================
%  SAVE RESULTS TABLE
%% =======================
Distance_km = Dis_thres(:);

T = table( ...
    Distance_km, ...
    EF_obs_all(:), CI_ratio_all(:,1), CI_ratio_all(:,2), ...
    EF_obs_re(:),  CI_ratio_re(:,1),  CI_ratio_re(:,2), ...
    EF_obs_na(:),  CI_ratio_na(:,1),  CI_ratio_na(:,2), ...
    CI_null_all(:,1), CI_null_all(:,2), ...
    CI_null_re(:,1),  CI_null_re(:,2), ...
    CI_null_na(:,1),  CI_null_na(:,2), ...
    'VariableNames', { ...
        'dmax_km', ...
        'EF_All', 'Jack90_All_Low','Jack90_All_High', ...
        'EF_Reservoir','Jack90_Re_Low','Jack90_Re_High', ...
        'EF_Natural',  'Jack90_Na_Low','Jack90_Na_High', ...
        'Null90_All_Low','Null90_All_High', ...
        'Null90_Re_Low','Null90_Re_High', ...
        'Null90_Na_Low','Null90_Na_High' ...
    });

fname = ['./Processed_Data/EQ_EF_Storage_Jackknife_Null_depth',num2str(maxDepth_km),'km.txt'];
writetable(T, fname, 'Delimiter', '\t');
fprintf('Saved table to %s\n', fname);



%% =======================
%  EXPORT ALL EQ-RESERVOIR PAIRS WITHIN 20 KM
%% =======================
if maxDepth_km == 30
    d_export = 20;   % km

    % Reservoir-only metadata
    Lake_out_Re = Lake_out(Group(2).mask, :);

    % Find ALL earthquake-reservoir pairs within threshold
    [iEQ, iLake] = find(D_re <= d_export);

    % Extract pairwise information
    pairTime   = Time_EQ(iEQ);
    pairMag    = Mag_EQ(iEQ);
    pairLat    = Lat_EQ(iEQ);
    pairLon    = Lon_EQ(iEQ);
    pairDep    = Dep_EQ(iEQ);
    pairDist   = D_re(sub2ind(size(D_re), iEQ, iLake));

    pairLakeID   = Lake_out_Re.LakeID(iLake);
    pairLakeName = Lake_out_Re.LakeName(iLake);
    pairLakeTrend = Lake_out_Re.Trend(iLake);
    pairLakeTrendPval = Lake_out_Re.TrendPval(iLake);
    pairLakeLat  = Lake_out_Re.Y(iLake);
    pairLakeLon  = Lake_out_Re.X(iLake);

    % Build output table
    T_re = table(pairTime, pairMag, pairLat, pairLon, pairDep, ...
        pairDist, pairLakeID, pairLakeName, pairLakeTrend, pairLakeTrendPval, pairLakeLat, pairLakeLon, ...
        'VariableNames', {'Time','Magnitude','Latitude','Longitude','Depth_km', ...
                          'Dist_km','LakeID','LakeName','Trend','TrendPval','LakeLat','LakeLon'});

    % Optional: sort by earthquake time, then distance
    T_re = sortrows(T_re, {'Time','Dist_km'});

    % Optional: write out
    writetable(T_re, './Processed_Data/EQ_Reservoir_Pairs_within20km.txt', ...
        'Delimiter', '\t');

end


%% ================== LOCAL FUNCTIONS ==================

function EF = EF_storage_allLakesWithin( ...
    dmax, D, Time_EQ, lakeGlobalIdx, ...
    Sthr_lake, TH_lake, TL_lake, tmin_lake, tmax_lake, ...
    tgrid_cell, sgrid_cell)

% - For each lake i, consider EQs within dmax of that lake AND within obs window.
% - Classify each EQ by storage(tEQ) > median threshold or <=.
% - Sum N_high,i, N_low,i, T_high,i, T_low,i over lakes.
% - EF = (sum N_high / sum T_high) / (sum N_low / sum T_low).

    EF = nan;

    nL = numel(lakeGlobalIdx);

    NHsum = 0; NLsum = 0;
    THsum = 0; TLsum = 0;

    for k = 1:nL
        j = lakeGlobalIdx(k); % global row in Lake_out

        inDist = (D(:,k) <= dmax);
        if ~any(inDist), continue; end

        t = Time_EQ(inDist);
        inWin = (t >= tmin_lake(j)) & (t <= tmax_lake(j));
        if ~any(inWin), continue; end
        t = t(inWin);

        tg = tgrid_cell{j};
        sg = sgrid_cell{j};

        sEQ = interp1(tg, sg, t, 'linear', NaN);
        good = isfinite(sEQ);
        if ~any(good), continue; end
        sEQ = sEQ(good);

        thr = Sthr_lake(j);
        NH = sum(sEQ >  thr);
        NL = sum(sEQ <= thr);

        NHsum = NHsum + NH;
        NLsum = NLsum + NL;
        THsum = THsum + TH_lake(j);
        TLsum = TLsum + TL_lake(j);
    end

    if THsum<=0 || TLsum<=0 || NLsum==0
        return
    end

    lambdaH = NHsum / THsum;
    lambdaL = NLsum / TLsum;
    EF = lambdaH / max(eps, lambdaL);
end

function D = haversineDistMatrix_single(latE, lonE, latL, lonL)
% Return NE x NL matrix of great-circle distances (km), stored as single.
% Uses chunking over EQs to keep memory stable.
    R = 6371.0;

    latE = latE(:); lonE = lonE(:);
    latL = latL(:); lonL = lonL(:);

    NE = numel(latE);
    NL = numel(latL);

    if NL == 0
        D = zeros(NE, 0, 'single');
        return
    end

    phiE = deg2rad(latE);
    lamE = deg2rad(lonE);
    phiL = deg2rad(latL);
    lamL = deg2rad(lonL);

    D = zeros(NE, NL, 'single');

    chunkSize = 3000;
    nChunks = ceil(NE/chunkSize);

    for c = 1:nChunks
        i1 = (c-1)*chunkSize + 1;
        i2 = min(c*chunkSize, NE);
        idx = i1:i2;

        dphi = phiE(idx) - phiL.';   % (nChunk x NL)
        dlam = lamE(idx) - lamL.';   % (nChunk x NL)

        a = sin(dphi/2).^2 + cos(phiE(idx)).*cos(phiL.').*sin(dlam/2).^2;
        d = 2*R*asin(min(1, sqrt(a)));

        D(idx,:) = single(d);
    end
end

function shadedCI(x, ylow, yhigh, faceColor, faceAlpha)
    x = x(:); ylow = ylow(:); yhigh = yhigh(:);
    fill([x; flipud(x)], [ylow; flipud(yhigh)], faceColor, ...
        'EdgeColor','none', 'FaceAlpha',faceAlpha);hold on;
    plot(x,ylow,'--','color',faceColor,'linewidth',2.5);
    plot(x,yhigh,'--','color',faceColor,'linewidth',2.5);
end

