clc,clear
close all

S=load('./Processed_Data/LakeSeries_all.mat');
LakeSeries=S.LakeSeries;
nAll = numel(LakeSeries);

inLakeTxt_orig  = './Processed_Data/Processed_Lake.txt';
Lake_out_orig   = readtable(inLakeTxt_orig);
ID_Lake_orig      = Lake_out_orig.LakeID;
TrendP_Lake_orig  = Lake_out_orig.TrendPval;
Trend_Lake_orig   = Lake_out_orig.Trend;


inLakeTxt  = './Processed_Data/Processed_Lake_correct_trend.txt';
Lake_out   = readtable(inLakeTxt);
Type_Lake = Lake_out.TypeName;
ID_Lake      = Lake_out.LakeID;
TrendP_Lake  = Lake_out.TrendPval;
Trend_Lake   = Lake_out.Trend;


%{
jkf=find(Trend_Lake<0&TrendP_Lake<0.1);
Type_Lake=Type_Lake(jkf);
ID_Lake=ID_Lake(jkf);
TrendP_Lake=TrendP_Lake(jkf);
Trend_Lake=Trend_Lake(jkf);
%}

NL = length(ID_Lake);
Re_flag=zeros(1,NL);
Na_flag=zeros(1,NL);
for i = 1:NL
    Re_flag(i) = strcmp(Type_Lake(i), 'Reservoir');
    Na_flag(i) = strcmp(Type_Lake(i), 'Natural lake');
end

jkf1=find(Na_flag==0);
Trend_Lake_Re  = Trend_Lake(jkf1);
Lake_out_Re = Lake_out(jkf1,:);
Trend_Lake_Re_orig  = Trend_Lake_orig(jkf1);
Lake_out_Re_orig = Lake_out_orig(jkf1,:);


jkf2=find(Na_flag==1);
Trend_Lake_Na  = Trend_Lake(jkf2);
Lake_out_Na = Lake_out(jkf2,:);
Trend_Lake_Na_orig  = Trend_Lake_orig(jkf2);
Lake_out_Na_orig = Lake_out_orig(jkf2,:);


[num_inc_Re,num_sta_Re,num_dec_Re]=get_each_cat_num(Lake_out_Re)
[num_inc_Re_orig,num_sta_Re_orig,num_dec_Re_orig]=get_each_cat_num(Lake_out_Re_orig)

[num_inc_Na,num_sta_Na,num_dec_Na]=get_each_cat_num(Lake_out_Na)
[num_inc_Na_orig,num_sta_Na_orig,num_dec_Na_orig]=get_each_cat_num(Lake_out_Na_orig)


figure('Units','normalized','Position',[0.2 0.2 0.4 0.2])

edges = -2:0.1:2;

% =========================
% (1) Relative trend change
% =========================
subplot(1,2,1)
histogram( (Trend_Lake(jkf1) - Trend_Lake_orig(jkf1)) ./ abs(Trend_Lake_orig(jkf1)), ...
           edges, 'Normalization','probability'); 
hold on;
histogram( (Trend_Lake(jkf2) - Trend_Lake_orig(jkf2)) ./ abs(Trend_Lake_orig(jkf2)), ...
           edges, 'Normalization','probability');

legend('Reservoir','Natural','Location','best');
xlabel('(Recalculated trend − original trend) / abs(original trend)');
ylabel('Probability');
grid on;

set(gca,'fontsize',16);

% =========================
% (2) Sign consistency
% =========================
% --- Signs (original -> recalculated) ---
s0_1 = sign(Trend_Lake_orig(jkf1));   % original (Reservoir)
s1_1 = sign(Trend_Lake(jkf1));        % recalculated

s0_2 = sign(Trend_Lake_orig(jkf2));   % original (Natural)
s1_2 = sign(Trend_Lake(jkf2));

% --- Define category order as STRINGS ---
order = ["+ → +", "+ → -", "- → +", "- → -"];

% --- Reservoir categories ---
cat1 = strings(numel(s0_1),1);
cat1(s0_1>0 & s1_1>0) = "+ → +";
cat1(s0_1>0 & s1_1<0) = "+ → -";
cat1(s0_1<0 & s1_1>0) = "- → +";
cat1(s0_1<0 & s1_1<0) = "- → -";

% --- Natural categories ---
cat2 = strings(numel(s0_2),1);
cat2(s0_2>0 & s1_2>0) = "+ → +";
cat2(s0_2>0 & s1_2<0) = "+ → -";
cat2(s0_2<0 & s1_2>0) = "- → +";
cat2(s0_2<0 & s1_2<0) = "- → -";

% --- Convert to categorical with specified order ---
cat1 = categorical(cat1, order, 'Ordinal', true);
cat2 = categorical(cat2, order, 'Ordinal', true);

% --- Plot (probability) ---
subplot(1,2,2)
histogram(cat1, 'Normalization','probability'); hold on;
histogram(cat2, 'Normalization','probability');

legend('Reservoir','Natural','Location','best');
xlabel('Original trend → Recalculated trend');
ylabel('Probability');
grid on;

set(gca,'fontsize',16);




titlestr=['Increase + Decrease + Stable'];%['Decrease'];%['Increase'];%

%{
[~,idx]=max(Trend_Lake_Re);
Lake_out_Re(idx,:)
targetID = Lake_out_Re.LakeID(idx);          % extract the numeric ID
jkf = find([LakeSeries.ID] == targetID); % compare against struct array IDs
figure;
plot(LakeSeries(jkf).tDec,LakeSeries(jkf).level);
figure;
plot(LakeSeries(jkf).tDec,LakeSeries(jkf).storage);
%}


% -------- Build masks for positive and negative significant trends --------
keepPos_Na = false(nAll,1);
keepNeg_Na = false(nAll,1);
keepEq_Na = false(nAll,1);
keepEq_Re = false(nAll,1);
keepPos_Re = false(nAll,1);
keepNeg_Re = false(nAll,1);

for i = 1:nAll
    jkf = find(ID_Lake == LakeSeries(i).ID, 1);
    if isempty(jkf), continue; end

    sig = TrendP_Lake(jkf) < 0.1;
    sig_n = TrendP_Lake(jkf) >= 0.1;
    Re_flag = strcmp(Type_Lake(jkf), 'Reservoir');
    Na_flag = strcmp(Type_Lake(jkf), 'Natural lake');

    if sig && Re_flag && Trend_Lake(jkf) > 0
        keepPos_Re(i) = true;
    elseif sig && Re_flag&& Trend_Lake(jkf) < 0
        keepNeg_Re(i) = true;
    elseif sig && Na_flag && Trend_Lake(jkf) > 0
        keepPos_Na(i) = true;
    elseif sig && Na_flag&& Trend_Lake(jkf) < 0
        keepNeg_Na(i) = true;
    elseif sig_n && Na_flag
        keepEq_Na(i) = true;
    elseif sig_n && Re_flag
        keepEq_Re(i) = true;
    end
end

LakeSeriesPos_Re = LakeSeries(keepPos_Re);
LakeSeriesNeg_Re = LakeSeries(keepNeg_Re);
LakeSeriesPos_Na = LakeSeries(keepPos_Na);
LakeSeriesNeg_Na = LakeSeries(keepNeg_Na);
LakeSeriesEq_Na = LakeSeries(keepEq_Na);
LakeSeriesEq_Re = LakeSeries(keepEq_Re);

N_Pos_Re  = numel(LakeSeriesPos_Re);
N_Neg_Re  = numel(LakeSeriesNeg_Re);
N_Pos_Na  = numel(LakeSeriesPos_Na);
N_Neg_Na  = numel(LakeSeriesNeg_Na);
N_Eq_Re = numel(LakeSeriesEq_Re);
N_Eq_Na = numel(LakeSeriesEq_Na);

N_Pos_All = N_Pos_Re + N_Pos_Na;
N_Neg_All = N_Neg_Re + N_Neg_Na;


LakeSeries_Re=[LakeSeriesPos_Re;LakeSeriesNeg_Re;LakeSeriesEq_Re];
LakeSeries_Na=[LakeSeriesPos_Na;LakeSeriesNeg_Na;LakeSeriesEq_Na];




%% ---- Convert to strings once ----
sN_Pos_Re  = num2str(N_Pos_Re);
sN_Eq_Re = num2str(N_Eq_Re);
sN_Neg_Re  = num2str(N_Neg_Re);
sN_Pos_Na  = num2str(N_Pos_Na);
sN_Eq_Na = num2str(N_Eq_Na);
sN_Neg_Na  = num2str(N_Neg_Na);
sN_Pos_All = num2str(N_Pos_All);
sN_Neg_All = num2str(N_Neg_All);


% -------- Plot both groups --------
% ================= Positive / Negative – Reservoir =================
if N_Pos_Re > 0
    [tCommon_Pos_Re, medianCurve_level_Pos_Re, medianCurve_storage_Pos_Re, ...
        area_Pos_Re, COV_area_Pos_Re] = GroupWithMedian( ...
        LakeSeriesPos_Re, ...
        ['Significant positive trend, Reservoir, N = ', sN_Pos_Re]);
else
    tCommon_Pos_Re = [];
    medianCurve_level_Pos_Re = [];
    medianCurve_storage_Pos_Re = [];
    area_Pos_Re = [];
    COV_area_Pos_Re = [];
end

if N_Neg_Re > 0
    [tCommon_Neg_Re, medianCurve_level_Neg_Re, medianCurve_storage_Neg_Re, ...
        area_Neg_Re, COV_area_Neg_Re] = GroupWithMedian( ...
        LakeSeriesNeg_Re, ...
        ['Significant negative trend, Reservoir, N = ', sN_Neg_Re]);
else
    tCommon_Neg_Re = [];
    medianCurve_level_Neg_Re = [];
    medianCurve_storage_Neg_Re = [];
    area_Neg_Re = [];
    COV_area_Neg_Re = [];
end


% ================= Positive / Negative – Natural Lakes =================
if N_Pos_Na > 0
    [tCommon_Pos_Na, medianCurve_level_Pos_Na, medianCurve_storage_Pos_Na, ...
        area_Pos_Na, COV_area_Pos_Na] = GroupWithMedian( ...
        LakeSeriesPos_Na, ...
        ['Significant positive trend, Natural Lake, N = ', sN_Pos_Na]);
else
    tCommon_Pos_Na = [];
    medianCurve_level_Pos_Na = [];
    medianCurve_storage_Pos_Na = [];
    area_Pos_Na = [];
    COV_area_Pos_Na = [];
end

if N_Neg_Na > 0
    [tCommon_Neg_Na, medianCurve_level_Neg_Na, medianCurve_storage_Neg_Na, ...
        area_Neg_Na, COV_area_Neg_Na] = GroupWithMedian( ...
        LakeSeriesNeg_Na, ...
        ['Significant negative trend, Natural Lake, N = ', sN_Neg_Na]);
else
    tCommon_Neg_Na = [];
    medianCurve_level_Neg_Na = [];
    medianCurve_storage_Neg_Na = [];
    area_Neg_Na = [];
    COV_area_Neg_Na = [];
end


% ================= Positive / Negative – All Lakes =================
if N_Pos_All > 0
    [tCommon_Pos_all, medianCurve_level_Pos_all, medianCurve_storage_Pos_all, ...
        ~, COV_area_Pos_all] = GroupWithMedian( ...
        [LakeSeriesPos_Na; LakeSeriesPos_Re], ...
        ['Significant positive trend, All, N = ', sN_Pos_All]);
else
    tCommon_Pos_all = [];
    medianCurve_level_Pos_all = [];
    medianCurve_storage_Pos_all = [];
    COV_area_Pos_all = [];
end

if N_Neg_All > 0
    [tCommon_Neg_all, medianCurve_level_Neg_all, medianCurve_storage_Neg_all, ...
        ~, COV_area_Neg_all] = GroupWithMedian( ...
        [LakeSeriesNeg_Na; LakeSeriesNeg_Re], ...
        ['Significant negative trend, All, N = ', sN_Neg_All]);
else
    tCommon_Neg_all = [];
    medianCurve_level_Neg_all = [];
    medianCurve_storage_Neg_all = [];
    COV_area_Neg_all = [];
end


[~, ~, ~, area_Re, ~] = GroupWithMedian( LakeSeries_Re,' ');
[~, ~, ~, area_Na, ~] = GroupWithMedian( LakeSeries_Na,' ');

        
    
radius_Re=sqrt(area_Pos_Re*1e9/pi)/1e3;
radius_Na=sqrt(area_Na*1e9/pi)/1e3;


%% Plot histogram

figure('Units','normalized','Position',[0.2 0.2 0.8 0.2])
subplot(1,4,1)
edges = 0:1:20;
plt_histogram_RevsNa(radius_Re,radius_Na,edges,titlestr)
xlabel('Radius (km)')

subplot(1,4,2)
edges = 0:0.005:0.4;
plt_histogram_RevsNa(abs(Trend_Lake_Re),abs(Trend_Lake_Na),edges,titlestr)
xlabel('Absolute Trend of storage (Gt/yr)');

subplot(1,4,3)
edges = 0:1:20;
plt_histogram_RevsNa([LakeSeries_Re.storage_extreme],[LakeSeries_Na.storage_extreme],edges,titlestr)
xlabel('Peak-to-peak difference of storage (Gt)');

subplot(1,4,4)
edges = 0:2:30;
plt_histogram_RevsNa([LakeSeries_Re.storage_extreme_time],[LakeSeries_Na.storage_extreme_time],edges,titlestr)
xlabel('Time between peak-to-peak difference of storage (yr)');



%%
% Concatenate data
ID      = [LakeSeries_Re.ID,      LakeSeries_Na.ID];
Extreme = [LakeSeries_Re.storage_extreme, LakeSeries_Na.storage_extreme];
Area = [area_Re,area_Na];

% Build lake-type labels (same length as ID)
Type = [ ...
    repmat("Reservoir",    numel(LakeSeries_Re), 1); ...
    repmat("NaturalLake",  numel(LakeSeries_Na), 1)  ...
];

% Build table
T = table(ID', Type, Extreme',Area', ...
    'VariableNames', {'LakeID','LakeType','MaxStorageChange_Gt','Area_km'});

% Write to txt (tab-delimited, easy to read in any editor)
fname = './Processed_Data/Lake_maximum_storgae_change.txt';
writetable(T, fname, 'Delimiter', '\t');

fprintf('Saved table to %s\n', fname);


 %% ======================= CDF of COV_area_* =======================
figure('Units','normalized','Position',[0.2 0.2 0.8 0.3])
subplot(1,3,1)
hold on;
% Plot only if data exist
if N_Pos_Re > 0
    [f_Pos_Re, x_Pos_Re] = ecdf(COV_area_Pos_Re);
    plot(x_Pos_Re, f_Pos_Re, 'LineWidth', 2);
end

if N_Neg_Re > 0
    [f_Neg_Re, x_Neg_Re] = ecdf(COV_area_Neg_Re);
    plot(x_Neg_Re, f_Neg_Re, 'LineWidth', 2);
end

if N_Pos_Na > 0
    [f_Pos_Na, x_Pos_Na] = ecdf(COV_area_Pos_Na);
    plot(x_Pos_Na, f_Pos_Na, 'LineWidth', 2);
end

if N_Neg_Na > 0
    [f_Neg_Na, x_Neg_Na] = ecdf(COV_area_Neg_Na);
    plot(x_Neg_Na, f_Neg_Na, 'LineWidth', 2);
end

xlim([0 1.5])
ylim([0 1])
xlabel('std of area / mean of area');
ylabel('CDF')

% Legend: include only the groups that were plotted
leg = {};
if N_Pos_Re > 0, leg{end+1} = ['Reservoir (Increase), N = ', sN_Pos_Re]; end
if N_Neg_Re > 0, leg{end+1} = ['Reservoir (Decrease), N = ', sN_Neg_Re]; end
if N_Pos_Na > 0, leg{end+1} = ['Natural (Increase), N = ', sN_Pos_Na]; end
if N_Neg_Na > 0, leg{end+1} = ['Natural (Decrease), N = ', sN_Neg_Na]; end
if ~isempty(leg)
    legend(leg, 'Location','best');
end

set(gca,'FontSize',14)


%% ======================= Median curves: LEVEL =======================
subplot(1,3,2)
hold on;
if N_Pos_Re > 0
    plot(tCommon_Pos_Re, medianCurve_level_Pos_Re,'-*','color',1/255*[213,118,107],'linewidth',2.5);
end
if N_Neg_Re > 0
    plot(tCommon_Neg_Re, medianCurve_level_Neg_Re,'-*','color',1/255*[168,200,228],'linewidth',2.5);
end
if N_Pos_Na > 0
    plot(tCommon_Pos_Na, medianCurve_level_Pos_Na,'--','color',1/255*[213,118,107],'linewidth',2.5);
end
if N_Neg_Na > 0
    plot(tCommon_Neg_Na, medianCurve_level_Neg_Na,'--','color',1/255*[168,200,228],'linewidth',2.5);
end
if N_Pos_All > 0
    plot(tCommon_Pos_all, medianCurve_level_Pos_all,'-','color',1/255*[213,118,107],'linewidth',2.5);
end
if N_Neg_All > 0
    plot(tCommon_Neg_all, medianCurve_level_Neg_all,'-','color',1/255*[168,200,228],'linewidth',2.5);
end

leg = {};
if N_Pos_Re > 0,  leg{end+1} = ['Reservoir, Increase, N = ', sN_Pos_Re]; end
if N_Neg_Re > 0,  leg{end+1} = ['Reservoir, Decrease, N = ', sN_Neg_Re]; end
if N_Pos_Na > 0,  leg{end+1} = ['Natural Lake, Increase, N = ', sN_Pos_Na]; end
if N_Neg_Na > 0,  leg{end+1} = ['Natural Lake, Decrease, N = ', sN_Neg_Na]; end
if N_Pos_All > 0, leg{end+1} = ['All, Increase, N = ', sN_Pos_All]; end
if N_Neg_All > 0, leg{end+1} = ['All, Decrease, N = ', sN_Neg_All]; end
if ~isempty(leg)
    legend(leg, 'location','northwest');
end

xlabel('Year');
ylabel('Normalized Water Level');
set(gca,'fontsize',16);
ylim([-1 1]);
xlim([1992 2020]);


%% ======================= Median curves: STORAGE =======================
subplot(1,3,3)
hold on;

if N_Pos_Re > 0
    plot(tCommon_Pos_Re, medianCurve_storage_Pos_Re,'-*','color',1/255*[213,118,107],'linewidth',2.5);
end
if N_Neg_Re > 0
    plot(tCommon_Neg_Re, medianCurve_storage_Neg_Re,'-*','color',1/255*[168,200,228],'linewidth',2.5);
end
if N_Pos_Na > 0
    plot(tCommon_Pos_Na, medianCurve_storage_Pos_Na,'--','color',1/255*[213,118,107],'linewidth',2.5);
end
if N_Neg_Na > 0
    plot(tCommon_Neg_Na, medianCurve_storage_Neg_Na,'--','color',1/255*[168,200,228],'linewidth',2.5);
end
if N_Pos_All > 0
    plot(tCommon_Pos_all, medianCurve_storage_Pos_all,'-','color',1/255*[213,118,107],'linewidth',2.5);
end
if N_Neg_All > 0
    plot(tCommon_Neg_all, medianCurve_storage_Neg_all,'-','color',1/255*[168,200,228],'linewidth',2.5);
end

leg = {};
if N_Pos_Re > 0,  leg{end+1} = ['Reservoir, Increase, N = ', sN_Pos_Re]; end
if N_Neg_Re > 0,  leg{end+1} = ['Reservoir, Decrease, N = ', sN_Neg_Re]; end
if N_Pos_Na > 0,  leg{end+1} = ['Natural Lake, Increase, N = ', sN_Pos_Na]; end
if N_Neg_Na > 0,  leg{end+1} = ['Natural Lake, Decrease, N = ', sN_Neg_Na]; end
if N_Pos_All > 0, leg{end+1} = ['All, Increase, N = ', sN_Pos_All]; end
if N_Neg_All > 0, leg{end+1} = ['All, Decrease, N = ', sN_Neg_All]; end
if ~isempty(leg)
    legend(leg, 'location','northwest');
end

xlabel('Year');
ylabel('Normalized Water Storage');
set(gca,'fontsize',16);
ylim([-1 1]);
xlim([1992 2020]);


%% ======================= Time vector assert + Year =======================
% Only assert for groups that actually exist (N>0)
tList = {};
if N_Pos_Re > 0,  tList{end+1} = tCommon_Pos_Re;  end
if N_Neg_Re > 0,  tList{end+1} = tCommon_Neg_Re;  end
if N_Pos_Na > 0,  tList{end+1} = tCommon_Pos_Na;  end
if N_Neg_Na > 0,  tList{end+1} = tCommon_Neg_Na;  end
if N_Pos_All > 0, tList{end+1} = tCommon_Pos_all; end
if N_Neg_All > 0, tList{end+1} = tCommon_Neg_all; end

if numel(tList) >= 2
    t0 = tList{1};
    for ii = 2:numel(tList)
        assert(isequal(t0, tList{ii}), 'Time vectors are not identical');
    end
end

% Pick a valid Year vector from the first available group
if N_Pos_Re > 0
    Year = tCommon_Pos_Re(:);
elseif N_Neg_Re > 0
    Year = tCommon_Neg_Re(:);
elseif N_Pos_Na > 0
    Year = tCommon_Pos_Na(:);
elseif N_Neg_Na > 0
    Year = tCommon_Neg_Na(:);
elseif N_Pos_All > 0
    Year = tCommon_Pos_all(:);
elseif N_Neg_All > 0
    Year = tCommon_Neg_all(:);
else
    Year = [];
end


%% save file
T_median = table( ...
    Year, ...
    medianCurve_storage_Pos_Re(:), ...
    medianCurve_storage_Neg_Re(:), ...
    medianCurve_storage_Pos_Na(:), ...
    medianCurve_storage_Neg_Na(:), ...
    medianCurve_storage_Pos_all(:), ...
    medianCurve_storage_Neg_all(:), ...
    'VariableNames', { ...
        'Year', ...
        'Reservoir_Increase', ...
        'Reservoir_Decrease', ...
        'Natural_Increase', ...
        'Natural_Decrease', ...
        'All_Increase', ...
        'All_Decrease' ...
    });


fname = './Processed_Data/WaterStorageTrend.txt';

% ---- write metadata header ----
fid = fopen(fname, 'w');

fprintf(fid, '# Median water storage curves\n');
fprintf(fid, '# Sample sizes\n');
fprintf(fid, '# N_Reservoir_Increase\t%d\n', N_Pos_Re);
fprintf(fid, '# N_Reservoir_Notrend\t%d\n', N_Eq_Re);
fprintf(fid, '# N_Reservoir_Decrease\t%d\n', N_Neg_Re);
fprintf(fid, '# N_Natural_Increase\t%d\n',   N_Pos_Na);
fprintf(fid, '# N_Natural_Notrend\t%d\n', N_Eq_Na);
fprintf(fid, '# N_Natural_Decrease\t%d\n',   N_Neg_Na);
fprintf(fid, '# N_All_Increase\t%d\n',       N_Pos_All);
fprintf(fid, '# N_All_Decrease\t%d\n',       N_Neg_All);
fprintf(fid, '\n');   % blank line before table

fclose(fid);
writetable(T_median, fname, ...
           'Delimiter', '\t', ...
           'WriteMode', 'append', ...
           'WriteVariableNames', true);
       
 %}      
       
       
function []=plt_histogram_RevsNa(Index_Re,Index_Na,edges,titlestr)

    % Define colors
    col_Re = [0.5 0 0.5];   % purple
    col_Na = [0 0.6 0];     % green

    % ---- Histogram (PDF / probability) ----
    yyaxis left
    ax = gca;
    ax.YColor = 'k';

    h1 = histogram(Index_Re, edges, ...
        'Normalization','probability', ...
        'DisplayStyle','stairs', ...
        'LineWidth',1.5, ...
        'EdgeColor',col_Re);
    hold on
    h2 = histogram(Index_Na, edges, ...
        'Normalization','probability', ...
        'DisplayStyle','stairs', ...
        'LineWidth',1.5, ...
        'EdgeColor',col_Na);
    ylabel('Probability');

    % ---- CDF ----
    yyaxis right
    ax = gca;
    ax.YColor = 'k';

    [f_Re, x_Re] = ecdf(Index_Re);
    [f_Na, x_Na] = ecdf(Index_Na);

    jkf1 = find(f_Na <= 0.05);
    jkf2 = find(f_Re <= 0.05);
    xmin = min([x_Na(jkf1(end)), x_Re(jkf2(end))]);

    jkf1 = find(f_Na >= 0.95);
    jkf2 = find(f_Re >= 0.95);
    xmax = max([x_Na(jkf1(1)), x_Re(jkf2(1))]);

    p1 = plot(x_Re, f_Re, '-','LineWidth',2, 'Color',col_Re);
    hold on
    p2 = plot(x_Na, f_Na,'-', 'LineWidth',2, 'Color',col_Na);
    ylabel('CDF');

    xlim([min(edges) max(edges)]);

    legend([h1 h2 p1 p2], ...
           {'Reservoir PDF','Lake PDF','Reservoir CDF','Lake CDF'}, ...
           'fontsize',10,'Location','best');

    grid on;
    set(gca,'fontsize',16);
    title(titlestr);
end
%}
function [tCommon, medianCurve1, medianCurve2, area_mean, COV_area]=GroupWithMedian(LakeSeriesGroup, plotTitleStr)

nLake = numel(LakeSeriesGroup);
hold on

if nLake == 0
    title([plotTitleStr ' (no lakes)']);
    xlabel('Year'); ylabel('Normalized level (–1 to 1)');
    set(gca,'fontsize',14); box on
    return
end

% Find global time range across this group (FIXED)
tMin = inf; tMax = -inf;
for i = 1:nLake
    tMin = min(tMin, min(LakeSeriesGroup(i).tDec));
    tMax = max(tMax, max(LakeSeriesGroup(i).tDec));
end

% Common monthly axis
tCommon = (floor(tMin*12)/12) : (1/12) : (ceil(tMax*12)/12);
Y = nan(numel(tCommon), nLake);

for i = 1:nLake
    t = LakeSeriesGroup(i).tDec(:);
    y1 = LakeSeriesGroup(i).level(:);
    y2 = LakeSeriesGroup(i).storage(:);

    good = isfinite(t) & isfinite(y1) & isfinite(y2);
    t = t(good);
    y1 = y1(good);
    y2 = y2(good);
    
    % area
    area=y2./y1;
    area = area(isfinite(area));
    area_mean(i) = mean(area);
    COV_area(i) = std(area)/mean(area);


    % Normalize to [-1, 1]
    y1min = min(y1);
    y1max = max(y1);
    
    y2min = min(y2);
    y2max = max(y2);
  
    y1_norm = 2*(y1 - y1min)/(y1max - y1min) - 1;
    y2_norm = 2*(y2 - y2min)/(y2max - y2min) - 1;

    % Interpolate to common time
    Y1(:,i) = interp1(t, y1_norm, tCommon, 'linear', NaN);
    Y2(:,i) = interp1(t, y2_norm, tCommon, 'linear', NaN);
    
    
    if max(y1_norm)>100
        LakeSeriesGroup(i).Name
    end
end

% Median curve (black)
medianCurve1 = median(Y1, 2, 'omitnan');
medianCurve2 = median(Y2, 2, 'omitnan');
end


function [num_inc,num_sta,num_dec]=get_each_cat_num(Lake_out)
    TrendP_Lake  = Lake_out.TrendPval;
    Trend_Lake   = Lake_out.Trend;
    num_inc = length(find(TrendP_Lake<0.1&Trend_Lake>0));
    num_dec = length(find(TrendP_Lake<0.1&Trend_Lake<0));
    num_sta = length(TrendP_Lake)-num_inc-num_dec;
end
