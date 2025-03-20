% ptp_analysis.m
% Purpose: Generate point-to-point (PTP) diagnostic comparison plots of outputs
% from IMPACTS WRF to EXRAD and NEXRAD observations
% 
% Author(s): Jon Seibert
% Last updated: 26 January 2023
% 
% Inputs: Various
% Outputs: *.png [ ]
%
% Process: (For each comparable timestep)
%   1) Read in obs, background, and analysis steps for a given datetime
%   2) Calculate basic ensemble mean, use mean for comparisons
%   3) Regrid WRF to each obs for direct comparison
%   4) Compare background-analysis, background-obs and analysis-obs and dif between
%   5) Calculate RMSE and Bias (mean of dif) for each
%   6) Plot all
%
% TODO:
%   - ADD BACKGROUND INPUT TO COMPARE IT TO ANALYSIS AND/OR OBS (previous
%   timestep) Format: wrfinput_d02_2022-01-29_13:00:00_001, in 1200 folder,
%   for 1300
%
% NOTES:
%  - TIME_COUNT is not robust when crossing month boundaries- replace with
%    manual time-tick count if using it to do so.
%  - NOT ROBUST TO USE IN THE SOUTHERN HEMISPHERE LATITUDES!

script_name = 'p3_plotting.m';

local_wd = "C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code";
server_wd = "/storage/home/jjs5895/projects/IMPACTS/code";

% Check for correct execution mode
working_dir = strrep(pwd,"\","/");
if(working_dir == local_wd)
    local = 1;
elseif(working_dir == server_wd)
    local = 0;
else
    error('Unknown environment detected.')
end
   

%% 0. General Settings

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/input';
    input_path_gis = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    input_path_nas = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/input/exrad_data/';
    intermediate_path = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/output/';
    path_to_code = "C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code";
    path_to_extra_code = './downloaded_code';
    input_path_ecmwf = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Data/ECMWF';
    input_path_st = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Data/IMPACTS_field_Catalog';
    % Experiment paths
    exp_path_1 = 'AIR'; % TESTING
    exp_path_2 = 'CONV/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/';
    input_path_gis = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/GIS';
    intermediate_path = '/storage/home/jjs5895/projects/IMPACTS/intermediate';
    output_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/ptp';
    path_to_code = "/storage/home/jjs5895/projects/IMPACTS/code";
    path_to_extra_code = '/storage/home/jjs5895/projects/IMPACTS/code/downloaded_code'; % Specify path to borders.m
    input_path_ecmwf = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/ECMWF';
    input_path_st = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/IMPACTS_field_catalog';
    % Experiment paths
    exp_path_1 = 'AIRCFT/fc';
    exp_path_2 = 'CONV/fc/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
end

addpath(path_to_extra_code);

% Date & time
% Range: 2022-01-29-0300 : 2022-01-30-0300
cur_year = 2020;
cur_month = 2;
cur_day = 7;
%ex_hour_list = ["1412","1443","1501","1516"];
%gis_hour_list = ["1410","1445","1500","1515"];
%gis_hour_list = ["1300","1400","1500","1600","1700","1800"];
%ex_hour_list = gis_hour_list;

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 900; %800
fig_height = 900; %750
fig_width_small = 350;
fig_height_small = 350;

limit_borders = true;

% Spatial limits of analysis and plot (degrees lat/lon)
% These are overridden by latlon_wrf_trim.mat
w_lim = -79;
e_lim = -69.75;
s_lim = 36;
n_lim = 46;
w_lim_d01 = -97;
e_lim_d01 = -67;
s_lim_d01 = 29.5;
n_lim_d01 = 48;
limit_borders = true; % Whether to apply spatial limits to plot
trim_e = false; % Whether to trim the eastern edge

% Figure font sizes
title_font_size = 20;
label_font_size = 20;
axes_font_size = 20;

data_type = "refl"; % Variable being analyzed
units = "dBZ";

cmap = 'reflmap'; % Choice of colorbar colormap
clim_lower = -30; % Colorbar limits
clim_upper = 75;

dbz_nodata = -30;
dbz_nodata_gis = -66;

plot_boxes = true;

%% 1A. Settings- P3 Data

input_path_p3 = sprintf('%s/%s',input_path_base,'p3_data/');

data_src_p3 = "p3";

input_filename = '02072020_temp.csv';

%lon_name_p3 = "longitude";
%lat_name_p3 = "latitude";
%data_name_p3 = "Ku_band_reflectivity";

%use_files_ex = 4:7; % Specify list of file indices to analyze
%use_files_ex = [4 4 4 4 4 4]; % Specify list of file indices to analyze
%file_list = string(ls(input_path_p3));
%filename_list_ex = file_list(use_files_ex);

% NOTES: Lat/lon coverage different for each file/timestamp, comes in 3D,
% needs adjustment from 0:360 to -180:180

%% 1B. Settings- GIS Data (NEXRAD, but in images)

input_path_gis = strcat(input_path_base,'/gis_data/2020');

data_src_gis = "nex";

file_format_gis = '%s/n0q_%04.f%02.f%02.f%02.f%02.f.png';

lat_gis_raw = 49.9975:-0.005:23.0025;
lon_gis_raw = -125.9975:0.005:-65.0025;
[lon_gis,lat_gis] = meshgrid(lon_gis_raw,lat_gis_raw);

% NOTES: Comes in image format, does not contain own latlon data
% Conversion formula: data_gis = double(imread())*0.5 - 32.5; 
        
%% 2. Preliminary Processing & Setup

% Lay out title and filename formats
datetime_format_file = '%04.f%02.f%02.f%04.f'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%04.f'; % year, month, day, hour [as numbers]

%time_count = length(gis_hour_list);

load('latlon_gis_trim.mat'); % lon-gis_trim, lat_gis_trim, border idx


%% Read in data

input_file_full = sprintf('%s/%s',input_path_p3,input_filename);

p_data = readtable(input_file_full);

obs_type = p_data.obs_type_degC_;
obs_lon = p_data.obs_lon_deg_;
obs_lat = p_data.obs_lat_deg_;
obs_pressure = p_data.obs_pressure_level_mb_;
obs_value = p_data.obs_value;
obs_error = zeros(size(obs_type));
obs_error(:) = 2.6; % Error = 2.6 K
obs_time = p_data.obs_time;

obs_year = year(obs_time);
obs_month = month(obs_time);
obs_day = day(obs_time);
obs_hour = hour(obs_time);
obs_min = minute(obs_time);
obs_sec = second(obs_time);

obs_timestamp = obs_hour+(obs_min/60)+(obs_sec/3600);
obs_da_step = round(obs_timestamp);

%lon_gis_trim = 

%%
% Cutting out one of the greys makes the colors line up to the numbers properly on caxis([14 21])
steps_cmap =     [0 255 255;	  % 0 <= x < 5 dbz (Cyan)
        30 144 255;   % 10 <= x < 15 dbz (Blue)
        0 139 40; 	  % 30 <= x < 35 dbz (Dark Green)
        0 255 0; 	  % 20 <= x < 25 dbz (Lime)
        255 255 0; 	  % 35 <= x < 40 dbz (Yellow)
        255 127 0; 	  % 45 <= x < 50 dbz (Orange)
        170 0 0]; 	  % 55 <= x < 60 dbz (Red)

steps_cmap = steps_cmap/255;

%%

obs_timestamp_a = obs_timestamp(obs_timestamp < 18.5);
obs_timestamp_b = obs_timestamp(obs_timestamp >= 18.5);
obs_lon_a = obs_lon(obs_timestamp < 18.5);
obs_lon_b = obs_lon(obs_timestamp >= 18.5);
obs_lat_a = obs_lat(obs_timestamp < 18.5);
obs_lat_b = obs_lat(obs_timestamp >= 18.5);

cmap = jet;
%cmap = steps_cmap;

plot_color_a = obs_timestamp_a;
plot_color_b = obs_timestamp_b;
%plot_color = obs_da_step;

% Plot LARGE
f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
h = scatter(obs_lon_a,obs_lat_a,[],plot_color_a,'filled');
hold on;
h = scatter(obs_lon_b,obs_lat_b,[],plot_color_b,'x');
c = colorbar('FontSize',axes_font_size); % Make colorbar
colormap(cmap); % Set colors
caxis([14 20]);

% Plot state borders
borders('continental us','black','linewidth',1); 

% Focus on desired area, remove whitespace
if(limit_borders)
    xlim([w_lim e_lim]);
    ylim([s_lim n_lim]);
end

%plot(bounding_sequence_x,bounding_sequence_y,'r','linewidth',2.2);
hold off;

% Apply labels
title(sprintf('P-3 Flight Path (Hour UTC) | 2020-02-07'),'FontSize',title_font_size);
%title(sprintf('P3 Flight Path and Time Step of Assimilation (hour Z) | 2020-02-07'),'FontSize',title_font_size);
xlabel('Longitude (deg)','FontSize',label_font_size);
ylabel('Latitude (deg)','FontSize',label_font_size);
set(gca,'Fontsize',axes_font_size);
saveas(h,sprintf('%s/p3/p3_path.png',output_path_base)); % Save as .png
%saveas(h,sprintf('%s/p3/p3_da_steps.png',output_path_base)); % Save as .png

% Plot SMALL
        f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
        %title(sprintf(title_format_small,exp_name_a,exp_name_b,plot_type,data_type,timestamp_file),'FontSize',title_font_size); % Replace title
        saveas(h,sprintf('%s/p3/p3_path_small.png',output_path_base)); % Save as .png=

close('all');

%% Pressure level over time

cmap = jet;
%cmap = steps_cmap;

plot_color = obs_timestamp;
%plot_color = obs_da_step;

% Plot LARGE
f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
h = scatter(obs_timestamp,obs_pressure,'filled');
set(gca, 'YDir','reverse')


% Apply labels
title(sprintf('P3 Flight Pressure Levels vs Time | 2020-02-07'),'FontSize',title_font_size);
%title(sprintf('P3 Flight Path and Time Step of Assimilation (hour Z) | 2020-02-07'),'FontSize',title_font_size);
xlabel('Time (UTC)','FontSize',label_font_size);
ylabel('Pressure (mb)','FontSize',label_font_size);
set(gca,'Fontsize',axes_font_size);
saveas(h,sprintf('%s/p3/p3_pressure.png',output_path_base)); % Save as .png
%saveas(h,sprintf('%s/p3/p3_da_steps.png',output_path_base)); % Save as .png

% Plot SMALL
        %f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
        %title(sprintf(title_format_small,exp_name_a,exp_name_b,plot_type,data_type,timestamp_file),'FontSize',title_font_size); % Replace title
        %saveas(h,sprintf('%s/p3/p3_pressure_small.png',output_path_base)); % Save as .png=

close('all');


%% Plot temperature fields

cmap = 'reflmap';



        
    % Plot LARGE
    f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
    h = pcolor(lon_gis_trim,lat_gis_trim,data_to_plot); % Plot the data
    set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
    shading interp; % Smooth out plot from grid-boxes
    c = colorbar('FontSize',axes_font_size); % Make colorbar
    colormap(cmap); % Set colors
    caxis([clim_lower clim_upper]);

    % Plot state borders
    hold on;
    borders('continental us','black','linewidth',1); 
    hold off;

    % Focus on desired area, remove whitespace
    if(limit_borders)
        xlim([w_lim e_lim]);
        ylim([s_lim n_lim]);
    end

    % Apply labels
    title(sprintf(title_format_raw,upper('wrf'),domain,member_string,'Raw',data_type,units,timestamp_title),'FontSize',title_font_size); 
    xlabel('Longitude (deg)','FontSize',label_font_size);
    ylabel('Latitude (deg)','FontSize',label_font_size);
    set(gca,'Fontsize',axes_font_size);
    saveas(h,sprintf(output_format_raw,output_path_large,'wrf',domain,member_string,'raw',data_type,timestamp_file)); % Save as .png

