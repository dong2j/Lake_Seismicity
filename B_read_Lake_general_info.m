clc;
clear;
close all;

% Read shapefile and convert to table
S = shaperead('./Raw_Data/GLWS v1.1/Global database of lake water storage GLWS shapefile v1.1/Global database of lake water storage GLWS shapefile v1.1.shp');
T = struct2table(S);

% Read centroid data
Location_centroid = readtable('./Raw_Data/GLWS_v1.1_Lakes_Centroid.xlsx');

centroid_data_ID = Location_centroid.LakeID;
centroid_data_X  = Location_centroid.Lon;
centroid_data_Y  = Location_centroid.Lat;

% Lake IDs from shapefile
shp_ID = T.LakeID;

% Preallocate X and Y if they do not already exist
if ~ismember('X', T.Properties.VariableNames)
    T.X = nan(height(T), 1);
end
if ~ismember('Y', T.Properties.VariableNames)
    T.Y = nan(height(T), 1);
end

% Match centroid coordinates to shapefile table
for i = 1:length(shp_ID)
    idx = find(centroid_data_ID == shp_ID(i));

    if length(idx) == 1
        dloc(i) = sqrt((centroid_data_X(idx)-T.X(i)).^2+(centroid_data_Y(idx)-T.Y(i)).^2);
        T.X(i) = centroid_data_X(idx);
        T.Y(i) = centroid_data_Y(idx);
    elseif isempty(idx)
        fprintf('No match for LakeID %d\n', shp_ID(i));
    else
        fprintf('More than one match for LakeID %d\n', shp_ID(i));
    end
end

% Select output columns
Lake_out = T(:, {'LakeID','LakeName','LakeArea','TypeName', ...
                 'Trend','TrendPval','X','Y'});

% Write output
writetable(Lake_out, './Processed_Data/Processed_Lake.txt');



% --- marker size (clip) ---
sz = max(abs(Lake_out.Trend), 0.01);
sz = min(sz, 0.3);
ms = 200 * sz;

% --- color by sign: red = positive, blue = negative ---
isPos = Lake_out.Trend > 0;
isNeg = Lake_out.Trend < 0;

figure; hold on

% --- global basemap (Mapping Toolbox) ---
axesm('robinson','Frame','on','Grid','on');   % or 'eqdcylin'
geoshow('landareas.shp','FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.6 0.6 0.6]);

% --- scatter on top (X=lon, Y=lat) ---
scatterm(Lake_out.Y(isPos), Lake_out.X(isPos), ms(isPos), 'r', 'filled');
scatterm(Lake_out.Y(isNeg), Lake_out.X(isNeg), ms(isNeg), 'b', 'filled');

figure;
histogram(dloc);
