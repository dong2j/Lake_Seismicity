clc; clear; close all

%% =======================
%  INPUTS
%% =======================
inCatMat   = './Processed_Data/Declustered_Events.txt';
inLakeTxt  = './Processed_Data/Processed_Lake_correct_trend.txt';

% EQ filter
maxDepth_km = 30;

%% =======================
%  LOAD EARTHQUAKE CATALOG
%% =======================
Cat_out = readtable(inCatMat);
Cat_out = Cat_out(Cat_out.Depth<maxDepth_km,:);
Time_EQ = Cat_out.Time;        % decimal year assumed
Mag_EQ  = Cat_out.Magnitude;
Lat_EQ  = Cat_out.Latitude;
Lon_EQ  = Cat_out.Longitude;
Dep_EQ  = Cat_out.Depth;

NE = numel(Lat_EQ);

%% =======================
%  LOAD LAKES
%% =======================
Lake_out = readtable(inLakeTxt);

Type_Lake   = Lake_out.TypeName;
TrendP_Lake = Lake_out.TrendPval;
Trend_Lake  = Lake_out.Trend;

isSig = TrendP_Lake < 0.1;
isRe  = strcmp(Type_Lake,'Reservoir');
isNa  = strcmp(Type_Lake,'Natural lake');

keepPos_Re = isSig & isRe & (Trend_Lake > 0);
keepNeg_Re = isSig & isRe & (Trend_Lake < 0);
keepPos_Na = isSig & isNa & (Trend_Lake > 0);
keepNeg_Na = isSig & isNa & (Trend_Lake < 0);

Lake_outPos_Re = Lake_out(keepPos_Re,:);
Lake_outNeg_Re = Lake_out(keepNeg_Re,:);
Lake_outPos_Na = Lake_out(keepPos_Na,:);
Lake_outNeg_Na = Lake_out(keepNeg_Na,:);

Lake_outPos_All = unique([Lake_outPos_Re; Lake_outPos_Na], 'rows');
Lake_outNeg_All = unique([Lake_outNeg_Re; Lake_outNeg_Na], 'rows');

fprintf('Lakes kept: Pos_Re=%d, Neg_Re=%d, Pos_Na=%d, Neg_Na=%d, Pos_All=%d, Neg_All=%d\n', ...
    height(Lake_outPos_Re), height(Lake_outNeg_Re), height(Lake_outPos_Na), height(Lake_outNeg_Na), ...
    height(Lake_outPos_All), height(Lake_outNeg_All));

%% ================= USER SETTINGS =================
Dis_thres = 20:5:150;   % dmax (km), cumulative: [0, dmax]
t1 = 1992; t2 = 2006; t3 = 2020;

% Jackknife over lakes: omit 10% lakes each realization
nJack    = 10000;
omitFrac = 0.10;

% Null model: shuffle EQ times (locations fixed)
nNull = 10000;

% Export EQ list for one chosen dmax
d_export = 20;  % km
%% =================================================

nD = numel(Dis_thres);

isEarlier_obs = (Time_EQ>=t1 & Time_EQ<t2);
isLater_obs   = (Time_EQ>=t2 & Time_EQ<t3);

%% =================================================
% 1) PRECOMPUTE EQ-to-LAKE DISTANCE MATRICES (ONCE)
%    (for speed in jackknife + null)
%% =================================================
fprintf('\nPrecomputing distance matrices (this is the main cost, done once)...\n');

D_Pos_Re  = haversineDistMatrix_single(Lat_EQ, Lon_EQ, Lake_outPos_Re.Y,  Lake_outPos_Re.X);
D_Neg_Re  = haversineDistMatrix_single(Lat_EQ, Lon_EQ, Lake_outNeg_Re.Y,  Lake_outNeg_Re.X);
D_Pos_Na  = haversineDistMatrix_single(Lat_EQ, Lon_EQ, Lake_outPos_Na.Y,  Lake_outPos_Na.X);
D_Neg_Na  = haversineDistMatrix_single(Lat_EQ, Lon_EQ, Lake_outNeg_Na.Y,  Lake_outNeg_Na.X);

% NEW: ALL (unique union)
D_Pos_All = haversineDistMatrix_single(Lat_EQ, Lon_EQ, Lake_outPos_All.Y, Lake_outPos_All.X);
D_Neg_All = haversineDistMatrix_single(Lat_EQ, Lon_EQ, Lake_outNeg_All.Y, Lake_outNeg_All.X);

%% =================================================
% 2) OBSERVED EF(dmax): using closest distance within each group
% =================================================
fprintf('Computing observed EF(dmax)...\n');

% Helper: safe min over columns (if group empty -> Inf distances)
dmin_Pos_Re  = safeMinDist(D_Pos_Re,  NE);
dmin_Neg_Re  = safeMinDist(D_Neg_Re,  NE);
dmin_Pos_Na  = safeMinDist(D_Pos_Na,  NE);
dmin_Neg_Na  = safeMinDist(D_Neg_Na,  NE);
dmin_Pos_All = safeMinDist(D_Pos_All, NE);
dmin_Neg_All = safeMinDist(D_Neg_All, NE);

EQ_enhance_ratio    = nan(nD,1);
EQ_enhance_ratio_Re = nan(nD,1);
EQ_enhance_ratio_Na = nan(nD,1);

for ii = 1:nD
    dmax = Dis_thres(ii);

    % --- ALL ---
    inPosAll = (dmin_Pos_All <= dmax);
    inNegAll = (dmin_Neg_All <= dmax);

    bPosAll = sum(inPosAll & isEarlier_obs);  aPosAll = sum(inPosAll & isLater_obs);
    bNegAll = sum(inNegAll & isEarlier_obs);  aNegAll = sum(inNegAll & isLater_obs);

    EQ_enhance_ratio(ii) = (bNegAll + aPosAll) / max(eps, (aNegAll + bPosAll));

    % --- Reservoir-only ---
    inPosRe = (dmin_Pos_Re <= dmax);
    inNegRe = (dmin_Neg_Re <= dmax);

    bPosRe = sum(inPosRe & isEarlier_obs);  aPosRe = sum(inPosRe & isLater_obs);
    bNegRe = sum(inNegRe & isEarlier_obs);  aNegRe = sum(inNegRe & isLater_obs);

    EQ_enhance_ratio_Re(ii) = (bNegRe + aPosRe) / max(eps, (aNegRe + bPosRe));

    % --- Natural-only ---
    inPosNa = (dmin_Pos_Na <= dmax);
    inNegNa = (dmin_Neg_Na <= dmax);

    bPosNa = sum(inPosNa & isEarlier_obs);  aPosNa = sum(inPosNa & isLater_obs);
    bNegNa = sum(inNegNa & isEarlier_obs);  aNegNa = sum(inNegNa & isLater_obs);

    EQ_enhance_ratio_Na(ii) = (bNegNa + aPosNa) / max(eps, (aNegNa + bPosNa));
end

%% =================================================
% 3) JACKKNIFE (FAST): omit lakes, recompute min over kept columns
%    No haversine recomputation, only min() + counts.
%% =================================================
fprintf('Jackknife %d realizations (omit %.0f%% lakes)...\n', nJack, 100*omitFrac);

EF_all_j = nan(nJack, nD);
EF_re_j  = nan(nJack, nD);
EF_na_j  = nan(nJack, nD);

nPosRe  = size(D_Pos_Re,2);
nNegRe  = size(D_Neg_Re,2);
nPosNa  = size(D_Pos_Na,2);
nNegNa  = size(D_Neg_Na,2);
nPosAll = size(D_Pos_All,2);
nNegAll = size(D_Neg_All,2);

for r = 1:nJack

    % number to omit
    nOmit_PosRe  = max(1, round(omitFrac * nPosRe));
    nOmit_NegRe  = max(1, round(omitFrac * nNegRe));
    nOmit_PosNa  = max(1, round(omitFrac * nPosNa));
    nOmit_NegNa  = max(1, round(omitFrac * nNegNa));
    nOmit_PosAll = max(1, round(omitFrac * nPosAll));
    nOmit_NegAll = max(1, round(omitFrac * nNegAll));

    % indices to KEEP
    idxKeep_PosRe  = setdiff(1:nPosRe,  randperm(nPosRe,  nOmit_PosRe));
    idxKeep_NegRe  = setdiff(1:nNegRe,  randperm(nNegRe,  nOmit_NegRe));
    idxKeep_PosNa  = setdiff(1:nPosNa,  randperm(nPosNa,  nOmit_PosNa));
    idxKeep_NegNa  = setdiff(1:nNegNa,  randperm(nNegNa,  nOmit_NegNa));
    idxKeep_PosAll = setdiff(1:nPosAll, randperm(nPosAll, nOmit_PosAll));
    idxKeep_NegAll = setdiff(1:nNegAll, randperm(nNegAll, nOmit_NegAll));
    
    % --- closest distance to kept lakes in each group ---
    dPosRe  = min(D_Pos_Re(:,  idxKeep_PosRe),  [], 2);
    dNegRe  = min(D_Neg_Re(:,  idxKeep_NegRe),  [], 2);
    dPosNa  = min(D_Pos_Na(:,  idxKeep_PosNa),  [], 2);
    dNegNa  = min(D_Neg_Na(:,  idxKeep_NegNa),  [], 2);
    dPosAll = min(D_Pos_All(:, idxKeep_PosAll), [], 2);
    dNegAll = min(D_Neg_All(:, idxKeep_NegAll), [], 2);


    for ii = 1:nD
        dmax = Dis_thres(ii);

        % --- ALL  ---
        inPosAll = (dPosAll <= dmax);
        inNegAll = (dNegAll <= dmax);

        bPosAll = sum(inPosAll & isEarlier_obs);  aPosAll = sum(inPosAll & isLater_obs);
        bNegAll = sum(inNegAll & isEarlier_obs);  aNegAll = sum(inNegAll & isLater_obs);

        EF_all_j(r,ii) = (bNegAll + aPosAll) / max(eps, (aNegAll + bPosAll));

        % --- Reservoir ---
        inPosRe = (dPosRe <= dmax);
        inNegRe = (dNegRe <= dmax);

        bPosRe = sum(inPosRe & isEarlier_obs);  aPosRe = sum(inPosRe & isLater_obs);
        bNegRe = sum(inNegRe & isEarlier_obs);  aNegRe = sum(inNegRe & isLater_obs);

        EF_re_j(r,ii) = (bNegRe + aPosRe) / max(eps, (aNegRe + bPosRe));

        % --- Natural ---
        inPosNa = (dPosNa <= dmax);
        inNegNa = (dNegNa <= dmax);

        bPosNa = sum(inPosNa & isEarlier_obs);  aPosNa = sum(inPosNa & isLater_obs);
        bNegNa = sum(inNegNa & isEarlier_obs);  aNegNa = sum(inNegNa & isLater_obs);

        EF_na_j(r,ii) = (bNegNa + aPosNa) / max(eps, (aNegNa + bPosNa));
    end
end

CI_ratio    = prctile(EF_all_j, [5 95], 1).';
CI_ratio_Re = prctile(EF_re_j,  [5 95], 1).';
CI_ratio_Na = prctile(EF_na_j,  [5 95], 1).';

%% =================================================
% 4) NULL MODEL (no influence): shuffle EQ times
%    Distances fixed, only time masks change.
%% =================================================
fprintf('Null model %d realizations (shuffle EQ times)...\n', nNull);

EF_all_null = nan(nNull, nD);
EF_re_null  = nan(nNull, nD);
EF_na_null  = nan(nNull, nD);

for r = 1:nNull
    
%     tmin = 1992;
%     tmax = 2020;
%     Time_shuf = tmin + (tmax - tmin) * rand(NE,1);

    Time_shuf = Time_EQ(randperm(NE));
    isEarlier = (Time_shuf>=t1 & Time_shuf<t2);
    isLater   = (Time_shuf>=t2 & Time_shuf<t3);

    for ii = 1:nD
        dmax = Dis_thres(ii);

        % --- ALL (correct, no double counting) ---
        inPosAll = (dmin_Pos_All <= dmax);
        inNegAll = (dmin_Neg_All <= dmax);

        bPosAll = sum(inPosAll & isEarlier);  aPosAll = sum(inPosAll & isLater);
        bNegAll = sum(inNegAll & isEarlier);  aNegAll = sum(inNegAll & isLater);

        EF_all_null(r,ii) = (bNegAll + aPosAll) / max(eps, (aNegAll + bPosAll));

        % --- Reservoir ---
        inPosRe = (dmin_Pos_Re <= dmax);
        inNegRe = (dmin_Neg_Re <= dmax);

        bPosRe = sum(inPosRe & isEarlier);  aPosRe = sum(inPosRe & isLater);
        bNegRe = sum(inNegRe & isEarlier);  aNegRe = sum(inNegRe & isLater);

        EF_re_null(r,ii) = (bNegRe + aPosRe) / max(eps, (aNegRe + bPosRe));

        % --- Natural ---
        inPosNa = (dmin_Pos_Na <= dmax);
        inNegNa = (dmin_Neg_Na <= dmax);

        bPosNa = sum(inPosNa & isEarlier);  aPosNa = sum(inPosNa & isLater);
        bNegNa = sum(inNegNa & isEarlier);  aNegNa = sum(inNegNa & isLater);

        EF_na_null(r,ii) = (bNegNa + aPosNa) / max(eps, (aNegNa + bPosNa));
    end
end

CI_null_all = prctile(EF_all_null, [5 95], 1).';
CI_null_re  = prctile(EF_re_null,  [5 95], 1).';
CI_null_na  = prctile(EF_na_null,  [5 95], 1).';

%% =================================================
% 5) PLOT (null band + jackknife CI + observed)
%% =================================================

% ---------- COLORS ----------
col_all   = [150  75   0] / 255;   % brown
col_re    = [120  60 160] / 255;   % purple
col_na    = [ 60 150  80] / 255;   % green
col_null  = [0.7 0.7 0.7];         % grey

alpha_CI   = 0.35;   % jackknife CI transparency
alpha_null = 0.45;   % null band transparency

figure('Units','normalized','Position',[0.18 0.2 0.72 0.32]);

% ===================== ALL =====================
subplot(1,3,1); hold on;

% null (grey)
shadedCI(Dis_thres, CI_null_all(:,1), CI_null_all(:,2), col_null, alpha_null);

% jackknife CI (light brown)
shadedCI(Dis_thres, CI_ratio(:,1), CI_ratio(:,2), col_all, alpha_CI);

% observed
plot(Dis_thres, EQ_enhance_ratio, '-', ...
    'Color', col_all, 'LineWidth', 1.8);

yline(1,'--','Color',[0.4 0.4 0.4]);
grid on;
xlabel('d_{max} (km)');
ylabel('EF');
title('All');
ylim([0 6]);
set(gca,'FontSize',14);

% ===================== RESERVOIR =====================
subplot(1,3,2); hold on;

shadedCI(Dis_thres, CI_null_re(:,1), CI_null_re(:,2), col_null, alpha_null);
shadedCI(Dis_thres, CI_ratio_Re(:,1), CI_ratio_Re(:,2), col_re, alpha_CI);

plot(Dis_thres, EQ_enhance_ratio_Re, '-', ...
    'Color', col_re, 'LineWidth', 1.8);

yline(1,'--','Color',[0.4 0.4 0.4]);
grid on;
xlabel('d_{max} (km)');
ylabel('EF');
title('Reservoir');
ylim([0 6]);
set(gca,'FontSize',14);

% ===================== NATURAL LAKE =====================
subplot(1,3,3); hold on;

shadedCI(Dis_thres, CI_null_na(:,1), CI_null_na(:,2), col_null, alpha_null);
shadedCI(Dis_thres, CI_ratio_Na(:,1), CI_ratio_Na(:,2), col_na, alpha_CI);

plot(Dis_thres, EQ_enhance_ratio_Na, '-', ...
    'Color', col_na, 'LineWidth', 1.8);

yline(1,'--','Color',[0.4 0.4 0.4]);
grid on;
xlabel('d_{max} (km)');
ylabel('EF');
title('Natural lake');
ylim([0 6]);
set(gca,'FontSize',14);

%% =================================================
% 6) SAVE RESULTS TABLE
%% =================================================
Distance_km = Dis_thres(:);

T = table( ...
    Distance_km, ...
    EQ_enhance_ratio(:), CI_ratio(:,1), CI_ratio(:,2), ...
    EQ_enhance_ratio_Re(:), CI_ratio_Re(:,1), CI_ratio_Re(:,2), ...
    EQ_enhance_ratio_Na(:), CI_ratio_Na(:,1), CI_ratio_Na(:,2), ...
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

fname = './Processed_Data/EQ_EF_LongTerm_Jackknife_Null_depth30km.txt';
writetable(T, fname, 'Delimiter', '\t');
fprintf('Saved table to %s\n', fname);
%}
%% =================================================
% 7) EXPORT EQ LISTS WITHIN d_export (FIXED)
%    (count-once within each group = closest lake in that group)
%% =================================================
if maxDepth_km == 70
    fprintf('Exporting EQ lists within %g km...\n', d_export);

    % --- POS RE ---
    [dmin, idxMin] = min(D_Pos_Re, [], 2);
    in = (dmin <= d_export);
    T_pos_re = table(Time_EQ(in), Mag_EQ(in), Lat_EQ(in), Lon_EQ(in), Dep_EQ(in), ...
        double(dmin(in)), Lake_outPos_Re.LakeID(idxMin(in)), Lake_outPos_Re.LakeName(idxMin(in)), Lake_outPos_Re.Y(idxMin(in)),  Lake_outPos_Re.X(idxMin(in)), ...
        'VariableNames', {'Time','Magnitude','Latitude','Longitude','Depth_km', ...
                          'ClosestDist_km','ClosestLakeID','LakeName','LakeLat','LakeLon'});
    T_pos_re.Group = repmat("Pos_Re", height(T_pos_re), 1);

    % --- POS NA ---
    [dmin, idxMin] = min(D_Pos_Na, [], 2);
    in = (dmin <= d_export);
    T_pos_na = table(Time_EQ(in), Mag_EQ(in), Lat_EQ(in), Lon_EQ(in), Dep_EQ(in), ...
        double(dmin(in)), Lake_outPos_Na.LakeID(idxMin(in)), Lake_outPos_Na.LakeName(idxMin(in)), Lake_outPos_Na.Y(idxMin(in)),  Lake_outPos_Na.X(idxMin(in)), ...
        'VariableNames', {'Time','Magnitude','Latitude','Longitude','Depth_km', ...
                          'ClosestDist_km','ClosestLakeID','LakeName','LakeLat','LakeLon'});
    T_pos_na.Group = repmat("Pos_Na", height(T_pos_na), 1);

    T_pos = [T_pos_re; T_pos_na];
    fname_pos = ['./Processed_Data/EQ_within_' num2str(d_export) 'km_IncreaseLakes.txt'];
    writetable(T_pos, fname_pos, 'Delimiter','\t');


    % --- NEG RE ---
    [dmin, idxMin] = min(D_Neg_Re, [], 2);
    in = (dmin <= d_export);
    T_neg_re = table(Time_EQ(in), Mag_EQ(in), Lat_EQ(in), Lon_EQ(in), Dep_EQ(in), ...
        double(dmin(in)), Lake_outNeg_Re.LakeID(idxMin(in)), Lake_outNeg_Re.LakeName(idxMin(in)), Lake_outNeg_Re.Y(idxMin(in)),  Lake_outNeg_Re.X(idxMin(in)), ...
        'VariableNames', {'Time','Magnitude','Latitude','Longitude','Depth_km', ...
                          'ClosestDist_km','ClosestLakeID','LakeName','LakeLat','LakeLon'});
    T_neg_re.Group = repmat("Neg_Re", height(T_neg_re), 1);

    % --- NEG NA ---
    [dmin, idxMin] = min(D_Neg_Na, [], 2);
    in = (dmin <= d_export);
    T_neg_na = table(Time_EQ(in), Mag_EQ(in), Lat_EQ(in), Lon_EQ(in), Dep_EQ(in), ...
        double(dmin(in)), Lake_outNeg_Na.LakeID(idxMin(in)),Lake_outNeg_Na.LakeName(idxMin(in)), Lake_outNeg_Na.Y(idxMin(in)),  Lake_outNeg_Na.X(idxMin(in)), ...
        'VariableNames', {'Time','Magnitude','Latitude','Longitude','Depth_km', ...
                          'ClosestDist_km','ClosestLakeID','LakeName','LakeLat','LakeLon'});
    T_neg_na.Group = repmat("Neg_Na", height(T_neg_na), 1);

    T_neg = [T_neg_re; T_neg_na];
    fname_neg = ['./Processed_Data/EQ_within_' num2str(d_export) 'km_DecreaseLakes.txt'];
    writetable(T_neg, fname_neg, 'Delimiter','\t');


    %% =================================================
    % Export lake time series for lakes referenced by T_pos / T_neg

    load('./Processed_Data/LakeSeries_all.mat');   % expects struct array LakeSeries_all with fields: ID, tDec, storage (or storgae)

    outDir = './Processed_Data/LakeSeries_export';
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    % --- collect referenced lakes from T_pos / T_neg (if they exist) ---
    allIDs   = [T_pos.ClosestLakeID;T_neg.ClosestLakeID];
    allNames = [string(T_pos.LakeName);string(T_neg.LakeName)];


    % --- unique by lake ID ---
    [uIDs, ia] = unique(allIDs, 'stable');
    uNames = allNames(ia);   % representative name for each ID

    % --- build a fast lookup of LakeSeries_all IDs ---
    seriesIDs = [LakeSeries.ID];

    %
    fprintf('Exporting %d lake time series to: %s\n', numel(uIDs), outDir);

    for k = 1:numel(uIDs)
        thisID = uIDs(k);
        thisName = uNames(k);

        % find matching series entry
        j = find(seriesIDs == thisID, 1, 'first');
        if isempty(j)
            warning('Lake ID %g not found in LakeSeries_all. Skipping.', thisID);
            continue
        end

        tDec = LakeSeries(j).tDec;
        stor = LakeSeries(j).storage;
        tDec = tDec(:);
        stor = stor(:);

        % build safe filename: ID_LakeName.txt
        safeName = regexprep(thisName, '[^\w\-]+', '_');  % keep letters/numbers/_/-
        safeName = regexprep(safeName, '_+', '_');        % collapse ____
        safeName = regexprep(safeName, '^_|_$', '');      % trim leading/trailing _

        if strlength(safeName) == 0
            safeName = "UnknownName";
        end

        fname = fullfile(outDir, sprintf('ID%d_%s.txt', thisID, safeName));

        % write as 2-column table
        Tout = table(tDec, stor, 'VariableNames', {'tDec','storage'});
        writetable(Tout, fname, 'Delimiter', '\t');

    end
end

%}
%% ================== LOCAL FUNCTIONS ==================

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

    chunkSize = 3000; % adjust if needed
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

function dmin = safeMinDist(D, NE)
% min over columns, but if D has 0 columns -> Inf
    if isempty(D) || size(D,2)==0
        dmin = inf(NE,1);
    else
        dmin = min(D, [], 2);
    end
end

function dmin = safeMinDistCols(D, keepMask, NE)
% min over kept columns, but safe for empty
    if isempty(D) || size(D,2)==0 || ~any(keepMask)
        dmin = inf(NE,1);
    else
        dmin = min(D(:, keepMask), [], 2);
    end
end

function shadedCI(x, ylow, yhigh, faceColor, faceAlpha)
    x = x(:); ylow = ylow(:); yhigh = yhigh(:);
    fill([x; flipud(x)], [ylow; flipud(yhigh)], faceColor, ...
        'EdgeColor','none', 'FaceAlpha',faceAlpha);hold on;
    plot(x,ylow,'--','color',faceColor,'linewidth',2.5);
    plot(x,yhigh,'--','color',faceColor,'linewidth',2.5);
end

