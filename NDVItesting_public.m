clear;clc;
% load HDF4 file from MODIS
info = hdfinfo("MOD13Q1.A2023353.h10v06.061.2024005132830.hdf");

% Access the 250m 16 days NDVI data from the first SDS in the first Vgroup
dataset_name = info.Vgroup.Vgroup(1).SDS(1).Name;  % '250m 16 days NDVI' (name of the dataset)

% Read the NDVI data using hdfread
NDVI = hdfread(info.Filename, dataset_name);
NDVI = double(NDVI) / 10000;
values = struct2cell(info.Attributes); % view attributes if needed
disp(values)

% Extract information from file name
[filepath,name,ext] = fileparts(info.Filename);
C = strsplit(name,'.');

% Get starting date (year, day of year) of composite image
Yr = strsplit(string(C(2)),'A');
Yr = char(Yr(2));
DOY = Yr(5:end);
Yr = Yr(1:4); 

Yr = str2double(Yr);
DOY = str2double(DOY);

% Combine year and DOY to get datetime
dt = datetime('today');
dt = dateshift(dt,'start','year');
dt.Year = Yr;
dt = dt+caldays(DOY)-1;

% Get horizontal and vertical tile coordinates
HV = strsplit(string(C(3)),{'h','v'});
h = str2double(HV(2));
v = str2double(HV(3));
%% Save variables from StructMetadata.0
structure = info.Attributes;
val = structure(2).Value;

% Split the string into lines
lines = strsplit(val, newline);

% Initialize an empty structure to store the variables
data = struct();

% Loop through each line and extract the variable name and value
for i = 1:length(lines)
    % Skip empty lines
    if isempty(strtrim(lines{i}))
        continue;
    end
    
    % Find the index of the equal sign to split the name and value
    eqIndex = strfind(lines{i}, '=');
    
    if ~isempty(eqIndex)
        % Extract the variable name (before the equal sign)
        varName = strtrim(lines{i}(1:eqIndex-1));
        
        % Extract the value (after the equal sign)
        valueStr = strtrim(lines{i}(eqIndex+1:end));
        
        % Attempt to convert the value to numeric if possible
        try
            % Handle cases with parentheses (e.g., tuples)
            if valueStr(1) == '(' && valueStr(end) == ')'
                valueStr = valueStr(2:end-1); % Remove the parentheses
                values = str2double(strsplit(valueStr, ','));
                data.(varName) = values;
            else
                % Convert to numeric value if it's not a string
                data.(varName) = str2double(valueStr);
                if isnan(data.(varName)) % if conversion fails, store as string
                    data.(varName) = valueStr;
                end
            end
        catch
            % If unable to convert to numeric, store the value as a string
            data.(varName) = valueStr;
        end
    end
end

% Display the structure to see the saved data
disp(data);

%% Get coordinate values for each datapoint
% Extract the UpperLeftPointMtrs and LowerRightMtrs from the metadata
UpperLeftPointMtrs = data.UpperLeftPointMtrs;
LowerRightMtrs = data.LowerRightMtrs;

% Grid dimensions
XDim = 4800; % Number of columns
YDim = 4800; % Number of rows

%% Reproject from Sinusoidal to lat and lon
radius = data.ProjParams(1);
T = 1111950; % height and width of each MODIS tile in the projection plane (meters);
xmin = -20015109; % western limit of projection plane (m)
ymax = 10007555; % northern limit of projection plane (m)
w = 231.65635826; % actual size of 250-m MODIS sinusoidal grid cell (m)

% Compute the position of the center of the sinusoidal grid cell

m = XDim;%rows
n = YDim;%columns

x = ones(XDim);
y = ones(XDim);

% Use formula from Appendix B of MODIS User Guide
% https://modis-land.gsfc.nasa.gov/pdf/MODIS_C61_BA_User_Guide_1.0.pdf
for i = 1:m
    for j = 1:n
        x(i,j) = (j + 0.5)*w + (h*T) + xmin;
        y(i,j) = ymax - (i + 0.5)*w - (v*T);
    end
end

lat = (y./radius)*180/pi;
lon = 180/pi*(x./(radius*cos(y./radius)));

%% Cropping 
% Load the KML file with the kml2struct function from https://www.mathworks.com/matlabcentral/fileexchange/35642-kml2struct
kmlFile = "MikeIslandPolygon.kml";  % Replace with your actual KML file path
kmlData = kml2struct(kmlFile);

% Extract polygon coordinates from the KML data
lon_kml = kmlData.Lon(~isnan(kmlData.Lon));
lat_kml = kmlData.Lat(~isnan(kmlData.Lat));

% Convert KML polygon to a mask for cropping the grid
inPolygon = inpolygon(lon, lat, lon_kml, lat_kml);

% Apply the mask to the gridded data to crop it
cropped_data = NDVI(inPolygon);
cropped_lat = lat(inPolygon);
cropped_lon = lon(inPolygon);