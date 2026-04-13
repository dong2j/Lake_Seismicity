
clc; clear; close all;
%{
load('events_plus_M6.mat');
Lat = T_filtered.lat;
Lon = T_filtered.lon;
Mag=T_filtered.mag;
Time0  = T_filtered.TimeDT;
t0 = datetime(1976,1,1,0,0,0,'TimeZone', Time0.TimeZone);
days_since_1976 = days(Time0 - t0);
Time=1976+days_since_1976/365.25;


jkf=find(Time>=1980&Mag>=(6-1e-12));
%}
D=load('./Raw_Data/Global.txt');
Time=D(:,1);
Mag=D(:,5);
Lat=D(:,2);
Lon=D(:,3);
Depth=D(:,4);
jkf=find(Time>=1987&Mag>=(5.5-1e-12)&Depth<=70);

% Reorder everything else using the same indices
Mag_sel = Mag(jkf);
Lat_sel = Lat(jkf);
Lon_sel = Lon(jkf);
Depth_sel = Depth(jkf);
Time_sel = Time(jkf);


[parent,lgRmin,lgTmin,lgETAmin,gmm,lgETAmin_threshold,cluster_ratio,median_cluster_lgTmin,median_cluster_lgRmin]=NND(Mag_sel,Time_sel,2,Lat_sel,Lon_sel,1,1.6);
color_curve0=1/255*[122 146 158;211 120 120;240 209 130; 60 110 185 ];

%%
figure('units','normalized','position',[0.1,0.1,0.4,0.3])
subplot(1,2,1)
plt_density_scatter(lgTmin,lgRmin)
xx=floor(min(lgTmin)):0.2:ceil(max(lgTmin));
yy=lgETAmin_threshold-xx;
plot(xx,yy,'--k','linewidth',2.5);

subplot(1,2,2)
histogram(lgETAmin, 'Normalization', 'pdf', 'FaceColor', color_curve0(1,:), 'EdgeColor', 'none');
hold on;
mu1 = gmm.mu(1);
mu2 = gmm.mu(2);
sigma1 = sqrt(gmm.Sigma(1));
sigma2 = sqrt(gmm.Sigma(2));
xx = linspace(min(lgETAmin), max(lgETAmin), 1000);
pdf1 = pdf('Normal', xx, mu1, sigma1);
pdf2 = pdf('Normal', xx, mu2, sigma2);
plot(xx, pdf1 * gmm.ComponentProportion(1),'color',color_curve0(2,:), 'LineWidth', 2);hold on;
plot(xx, pdf2 * gmm.ComponentProportion(2), 'color',color_curve0(3,:), 'LineWidth', 2);hold on;
plot(xx, pdf1 * gmm.ComponentProportion(1)+pdf2 * gmm.ComponentProportion(2),'--','color',color_curve0(4,:), 'LineWidth', 2);hold on;
plot([lgETAmin_threshold, lgETAmin_threshold], [0 0.4], 'k--', 'LineWidth', 2);
xlabel('$\log_{10} \eta$', 'Interpreter', 'latex')
ylabel('Probability Density');
set(gca,'fontsize',16);grid on;box on;grid minor;

%%
figure('units','normalized','position',[0.1,0.1,0.4,0.8])
back_idx = find(lgETAmin>=lgETAmin_threshold);
back_eta=lgETAmin(back_idx);
back_lgRmin=lgRmin(back_idx);
back_lgTmin=lgTmin(back_idx);

subplot(2,2,1)
plt_density_scatter(back_lgTmin,back_lgRmin)


Time_sel_new = Time_sel(2:end);
Mag_sel_new = Mag_sel(2:end);
Lat_sel_new = Lat_sel(2:end);
Lon_sel_new = Lon_sel(2:end);
Depth_sel_new = Depth_sel(2:end);

subplot(2,2,2)
Time_min=1987;
Time_max=2020;
tt=linspace(Time_min,Time_max,100);
for jj=1:length(tt)
    jkf=find(Time_sel_new<=tt(jj));
    cul_EQ0(jj)=length(jkf);
end
plot(tt,cul_EQ0,'-r','linewidth',2.5);hold on;

Time_back=Time_sel_new(back_idx);
Mag_back=Mag_sel_new(back_idx);
Lat_back=Lat_sel_new(back_idx);
Lon_back=Lon_sel_new(back_idx);
Depth_back = Depth_sel_new(back_idx);

for jj=1:length(tt)
    jkf=find(Time_back<=tt(jj));
    cul_EQ(jj)=length(jkf);
end
plot(tt,cul_EQ,'-','color',[0.75 0.75 0.75],'linewidth',2.5);hold on;

subplot(2,2,3)
[~,lgRmin,lgTmin,lgETAmin,~,~,~,~,~]=NND(Mag_back,Time_back,2,Lat_back,Lon_back,1,1.6);
plt_density_scatter(lgTmin,lgRmin);

%% Build catalog
Cat_out = [Time_back, Mag_back, Lat_back, Lon_back, Depth_back];

% Time window
jkf = Cat_out(:,1) >= 1992 & Cat_out(:,1) < 2020;
Cat_out = Cat_out(jkf,:);

% Convert to table with column names
T = table( ...
    Cat_out(:,1), ... % Time
    Cat_out(:,2), ... % Magnitude
    Cat_out(:,3), ... % Latitude
    Cat_out(:,4), ... % Longitude
    Cat_out(:,5), ... % Depth
    'VariableNames', {'Time','Magnitude','Latitude','Longitude','Depth'} );

% Output path
outTxt = './Processed_Data/Declustered_Events.txt';

% Write as tab-delimited text
writetable(T, outTxt, 'Delimiter','\t');

disp("Wrote: " + outTxt);

%%
all_back_eta = [];
all_back_lgRmin = [];
all_back_lgTmin = [];

num_rand = 100;    
for i=1:length(num_rand)

    n=length(back_idx);
    pMag = randperm(n);
    pLoc = randperm(n);

    % Apply shuffles
    Mag_tmp = Mag_back(pMag);     % shuffle magnitudes
    Lat_tmp = Lat_back(pLoc);     % shuffle locations
    Lon_tmp = Lon_back(pLoc);     % keep Lat/Lon paired


    [~,lgRmin,lgTmin,lgETAmin,~,~,~,~,~]=NND(Mag_tmp,Time_back,2,Lat_tmp,Lon_tmp,1,1.6);
    
   
    all_back_eta=[all_back_eta;lgETAmin];
    all_back_lgRmin=[all_back_lgRmin;lgRmin];
    all_back_lgTmin=[all_back_lgTmin;lgTmin];
end

subplot(2,2,4)
plt_density_scatter(all_back_lgTmin,all_back_lgRmin)

%{

%% stochasticly get background events
figure('units','normalized','position',[0.1,0.1,0.4,0.8])

PDF_back=pdf1 * gmm.ComponentProportion(1);
PDF_all=pdf1 * gmm.ComponentProportion(1)+pdf2 * gmm.ComponentProportion(2);
P_back=PDF_back./PDF_all;

n = numel(lgETAmin);
num_rand = 100;
back_idxs = cell(num_rand,1);

all_back_eta = [];
all_back_lgRmin = [];
all_back_lgTmin = [];
for jj = 1:num_rand
    P_retain = interp1(xx, P_back, lgETAmin, 'nearest', 'extrap');  % nx1

    retain_idx = rand(n,1) < P_retain;
    back_idxs{jj} = find(retain_idx);
    all_back_eta=[all_back_eta;lgETAmin(back_idxs{jj})];
    all_back_lgRmin=[all_back_lgRmin;lgRmin(back_idxs{jj})];
    all_back_lgTmin=[all_back_lgTmin;lgTmin(back_idxs{jj})];

end
subplot(2,2,1)
plt_density_scatter(all_back_lgTmin,all_back_lgRmin)


%%
Time_sel_new = Time_sel(2:end);
Mag_sel_new = Mag_sel(2:end);
Lat_sel_new = Lat_sel(2:end);
Lon_sel_new = Lon_sel(2:end);
subplot(3,2,4)
Time_min=1976;
Time_max=2020;
tt=linspace(Time_min,Time_max,100);
for jj=1:length(tt)
    jkf=find(Time_sel_new<=tt(jj));
    cul_EQ0(jj)=length(jkf);
end
plot(tt,cul_EQ0,'-r','linewidth',2.5);hold on;


for i=1:length(num_rand)
    back_idx=back_idxs{i};
    Time_tmp=Time_sel_new(back_idx);
    for jj=1:length(tt)
        jkf=find(Time_tmp<=tt(jj));
        cul_EQ(i,jj)=length(jkf);
    end
    plot(tt,cul_EQ(i,:),'-','color',[0.75 0.75 0.75],'linewidth',2.5);hold on;
end

%%
    

all_back_eta = [];
all_back_lgRmin = [];
all_back_lgTmin = [];

for i=1:length(num_rand)
    back_idx=back_idxs{i};
    Mag_tmp=Mag_sel_new(back_idx);
    Time_tmp=Time_sel_new(back_idx);
    Lat_tmp=Lat_sel_new(back_idx);
    Lon_tmp=Lon_sel_new(back_idx);

    [~,lgRmin,lgTmin,lgETAmin,~,~,~,~,~]=NND(Mag_tmp,Time_tmp,2,Lat_tmp,Lon_tmp,1,1.6);
    
   
    all_back_eta=[all_back_eta;lgETAmin];
    all_back_lgRmin=[all_back_lgRmin;lgRmin];
    all_back_lgTmin=[all_back_lgTmin;lgTmin];
end 

subplot(3,2,5)
plt_density_scatter(all_back_lgTmin,all_back_lgRmin)


%%
all_back_eta = [];
all_back_lgRmin = [];
all_back_lgTmin = [];

    
for i=1:length(num_rand)
    back_idx=back_idxs{i};
    Mag_tmp=Mag_sel_new(back_idx);
    Time_tmp=Time_sel_new(back_idx);
    Lat_tmp=Lat_sel_new(back_idx);
    Lon_tmp=Lon_sel_new(back_idx);

    n=length(back_idx);
    pMag = randperm(n);
    pLoc = randperm(n);

    % Apply shuffles
    Mag_tmp = Mag_tmp(pMag);     % shuffle magnitudes
    Lat_tmp = Lat_tmp(pLoc);     % shuffle locations
    Lon_tmp = Lon_tmp(pLoc);     % keep Lat/Lon paired


    [~,lgRmin,lgTmin,lgETAmin,~,~,~,~,~]=NND(Mag_tmp,Time_tmp,2,Lat_tmp,Lon_tmp,1,1.6);
    
   
    all_back_eta=[all_back_eta;lgETAmin];
    all_back_lgRmin=[all_back_lgRmin;lgRmin];
    all_back_lgTmin=[all_back_lgTmin;lgTmin];
end

subplot(3,2,6)
plt_density_scatter(all_back_lgTmin,all_back_lgRmin)

%% output

idx_back=back_idxs{1}+1;

Mag_sel_new=Mag_sel_new(back_idx);
Time_sel_new=Time_sel_new(back_idx);
Lat_sel_new=Lat_sel_new(back_idx);
Lon_sel_new=Lon_sel_new(back_idx);
Cat_out=[Time_sel_new,Mag_sel_new,Lat_sel_new,Lon_sel_new];
jkf=find(Cat_out(:,1)>=1992&Cat_out(:,1)<2021);
Cat_out = Cat_out(jkf,:);
outMat = './Processed_Data/Declustered_Events.mat';
save(outMat, 'Cat_out');
disp("Wrote: " + outMat);
%{
T_filtered=T_filtered(idx_back,:);
jkf=find(year(T_filtered.TimeDT)>1992);
T_filtered=T_filtered(jkf,:);
outMat = 'events_1992plus_M6plus_background.mat';
save(outMat, 'T_filtered');
disp("Wrote: " + outMat);
%}

%}
%}
%% helper
function [parent,lgRmin,lgTmin,lgETAmin,gmm,lgETAmin_threshold,cluster_ratio,median_cluster_lgTmin,median_cluster_lgRmin] = ...
    NND(Dm, Dt, dayoryear, Dlat, Dlon, b, df)

    % Ensure column vectors
    Dm   = Dm(:);
    Dt   = Dt(:);
    Dlat = Dlat(:);
    Dlon = Dlon(:);

    n = length(Dm);

    % ---- Time unit handling ----
    if dayoryear == 1
        dayfac = 1/365.25; % days -> years
    else
        dayfac = 1;        % already in years (or keep as-is)
    end

    % ---- Fix zero/duplicate time differences ----
    difft = diff(Dt);
    if any(difft > 0)
        min_nonzero = min(difft(difft > 0));
        difft(difft == 0) = min_nonzero;
        dt_floor = min(difft)/2;
    else
        dt_floor = eps; % fallback if all times are identical
    end

    % ---- Precompute radians ----
    latRad = deg2rad(Dlat);
    lonRad = deg2rad(Dlon);

    % ---- Outputs ----
    parent   = zeros(n-1, 1);
    lgRmin   = zeros(n-1, 1);
    lgTmin   = zeros(n-1, 1);
    lgETAmin = zeros(n-1, 1);

    R_earth_km = 6371; % km

    % ---- Main loop: compute distance vector only (1:j-1) -> j ----
    for j = 2:n
        idx = 1:(j-1);

        % Haversine distance from all previous points to point j
        dlat = latRad(j) - latRad(idx);
        dlon = lonRad(j) - lonRad(idx);

        % Wrap to [-pi, pi] so dateline works
        dlon = mod(dlon + pi, 2*pi) - pi;

        a = sin(dlat/2).^2 + cos(latRad(j)).*cos(latRad(idx)).*sin(dlon/2).^2;
        c = 2 * atan2(sqrt(a), sqrt(1 - a));
        rij = R_earth_km * c; % (j-1)x1, km

        % Avoid exact zeros (shouldn't happen except duplicates)
        rij(rij == 0) = min(rij(rij > 0))/2;

        % R term
        R = rij(:).^df .* 10.^(-b .* Dm(idx) / 2);

        % T term
        tij = Dt(j) - Dt(idx);
        tij = max(dt_floor, tij);
        tij = dayfac * tij;
        T = tij(:) .* 10.^(-b .* Dm(idx) / 2);

        % ETA
        ETA = R .* T;

        % parent = argmin ETA
        [~, k] = min(ETA);
        parent(j-1)   = k;
        lgRmin(j-1)   = log10(R(k));
        lgTmin(j-1)   = log10(T(k));
        lgETAmin(j-1) = log10(ETA(k));
    end

    % ---- GMM thresholding / clustering ----
    options = statset('MaxIter', 1000);
    gmm = [];
    lgETAmin_threshold = NaN;
    cluster_ratio = NaN;
    median_cluster_lgTmin = NaN;
    median_cluster_lgRmin = NaN;

    try
        gmm = fitgmdist(lgETAmin, 2, 'Options', options);

        mu1 = gmm.mu(1);
        mu2 = gmm.mu(2);
        sigma1 = sqrt(gmm.Sigma(1));
        sigma2 = sqrt(gmm.Sigma(2));

        f = @(x) pdf('Normal', x, mu1, sigma1) * gmm.ComponentProportion(1) - ...
                 pdf('Normal', x, mu2, sigma2) * gmm.ComponentProportion(2);

        if sign(f(mu1)) ~= sign(f(mu2))
            lgETAmin_threshold = fzero(f, [min(mu1,mu2), max(mu1,mu2)]);

            jkf = find(lgETAmin < lgETAmin_threshold);
            cluster_ratio = length(jkf) / length(lgETAmin);

            if ~isempty(jkf)
                median_cluster_lgTmin = median(lgTmin(jkf));
                median_cluster_lgRmin = median(lgRmin(jkf));
            end
        end
    catch
        % keep NaNs
    end
end



function []=plt_density_scatter(lgTmin,lgRmin)
        lgTmin_min=-9;
        lgTmin_max=-1;
        lgRmin_min=-4;
        lgRmin_max=4;
        density_2D=density2C(lgTmin,lgRmin,lgTmin_min:0.2:lgTmin_max,lgRmin_min:0.2:lgRmin_max);
        scatter(lgTmin,lgRmin,20,density_2D,'filled');hold on;
        colormap(gca, cool);
        xlabel('$\log_{10} T$', 'Interpreter', 'latex')
        ylabel('$\log_{10} R$', 'Interpreter', 'latex')
        set(gca,'fontsize',16);
        grid on;box on;grid minor;
        xlim([lgTmin_min,lgTmin_max]);
        ylim([lgRmin_min,lgRmin_max]);
end

function [h]=density2C(X,Y,XList,YList)
    [XMesh,YMesh]=meshgrid(XList,YList);
    XYi=[XMesh(:) YMesh(:)];
    F=ksdensity([X,Y],XYi);
    ZMesh=zeros(size(XMesh));
    ZMesh(1:length(F))=F;
    h=interp2(XMesh,YMesh,ZMesh,X,Y);
end