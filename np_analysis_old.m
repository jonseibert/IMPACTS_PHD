% np_analysis.m
%
% Purpose: Generate neighborhood-probability (NP) and neighborhood ensemble 
% probability (NEP) diagnostic plots of outputs from an IMPACTS WRF Case, 
% for a specified variable
% 
% Author(s): Jon Seibert
% Last updated: 9 June 2023
% 
% Inputs: *.nc [Various NetCDF files]
% Outputs: *.png [Set of NP plots with varied ROI/cutoff]
%
% Dependencies: neighborhood_prob_2.m, haversine_distance.m, borders.m, reflmap.m
%
% Process: 
% 1) Convert to binary
% 2) ROI checking:
% For each point no more than x gridpoints away:
%   a) Retrieve lat/lon coords. Use Haversine distance to measure distance D between centerpoint and target.
%   b) If D <= ROI, add binary value to list
% 3) Then compute mean of all values within radius
%
% Usage: Edit sections "0. General Settings" and  "1. Settings" to specify the inputs, run entire script
%
% TODO:
%   - Finish adding ability to NP the observations
%
% NOTES:
%   - CURRENTLY SET UP TO IGNORE THE ENSEMBLE MEAN AND MAKE ITS OWN 
%     (implemented for Reflectivity analysis)
%   - WILL BREAK IF CROSSING MONTH BOUNDARIES WITH DATA (find variable
%     'time_count' to adjust)
%   - Assumes lat/lon consistent across all times and members

%% 0. General Settings

% Local or server?
local = false; % true = running on local device, false = running on server
%local = true;

% Plot font sizes
title_font_size = 16;
label_font_size = 16;
axes_font_size = 14;

max_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 30, 31, 31]; % Number of days in each month

% Lay out title and filename formats
% Input filenames, plot titles, output filenames, datetime strings
% Compressed = for small plots
% Raw = for direct WRF member plotting
filename_format = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #
filename_format_alt = 'wrfinput_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00_%03.f'; % domain, year, month, day, hour, member #
title_format = '[%s-d%02.f-%s] Neigh. Prob. of %s > %d %s, ROI %d km (%s)'; % Data source (caps), domain, member #, data type, cutoff value, units, ROI/1000, datetime
title_format_compressed = 'NP>%d|%s|%s'; % cutoff, member #, datetime
title_format_raw =  '[%s|d0%d|%s] %s: %s(%s)|%s'; % data source, domain, plot type, member #/mean, data type, units, timestamp
output_format = '%s/np_%s_d%02.f_%s_%s_%d%s_%dkm_%s.png'; % output path, data source, domain, member #, data type, cutoff value, units, ROI/1000, datetime
output_format_raw = '%s/%s_d0%d_%s_%s_%s_%s.png'; % output path, domain, mem/mean, fit/dif, data type, timestamp
datetime_format_file = '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]
nep_save_format = 'nep_storage/nep_%s_%s_%d%s_%dkm'; % data_type, timestamp_file, cutoff, units, roi/1000
input_format_base_refl = '%s/%s_%s_d%02.f_%s_%s_%s.mat'; % output path, data source, plot type, domain, member #, data type, datetime

%% 1. Settings- WRF Data (2020 Case)

% Reference filename: wrfinput_d01_2020-02-07_14:00:00_001

% Denotes output folder name
run_name = 'np_air_2';
exp_name = 'AIR';
use_obs = false;
alt_input_format = false;
use_base_ref = true;

% Basic details
domain = 2;
num_members = 40;
data_src = lower(exp_name); % Dataset label
data_type = "refl"; % Variable being analyzed
units = "dBZ";

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020';
    %input_path_gis = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    intermediate_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/testing';
    path_to_borders = './borders';
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/AIRCFT/fc';
    input_path_br = '/storage/home/jjs5895/projects/IMPACTS/data/2020/';
    %input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/AIRCFT-F-18Z';
    %input_path_gis = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/GIS';
    intermediate_path = '/storage/home/jjs5895/projects/IMPACTS/intermediate';
    output_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/np';
    path_to_borders = '/storage/home/jjs5895/projects/IMPACTS/code/borders'; % Specify path to borders.m
end

% NetCDF retrieval
lon_name = "XLONG"; 
lat_name = "XLAT"; 
data_name = "REFL_10CM"; % Variable being analyzed

% Note whether transposition or compositing is necessary to properly align
% the data with desired grid
transpose_data = true;
flip_data = true;
composite_data = false;

% Specific binary conversion cutoff values
% Units = dBZ
cut_min = 25;
cut_step = 5;
cut_max = 30;

% Specific radius of influence values
% Units = m
roi_min = 30*1000;
roi_step = 10*1000;
roi_max = 30*1000;

% Date & time
% Range: 2020-02-07-1300 : 2020-02-07-1800
% But 1300 is the same as 1400
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 7;
start_hour = 14;
end_hour = 18;

% Hour increment size
hour_step = 1; 

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 925;
fig_height = 900;
fig_width_small = 350;
fig_height_small = 350;

% Spatial limits of analysis and plot (degrees lat/lon)
w_lim = -79;
e_lim = -69.75;
s_lim = 36;
n_lim = 46;
limit_borders = true; % Whether to apply spatial limits to plot
trim_e = false; % Whether to trim the eastern edge

cmap = 'reflmap'; % Choice of colorbar colormap
clim_lower = -30; % Colorbar limits
clim_upper = 75;

jet_blue_percent = 0.2; % Cut off coldest (Darkest blues) _% of jet colorbar for ease of analysis

%% OBS

data_src_gis = "nex";

input_path_gis = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/GIS';
file_format_gis = 'n0q_%04.f%02.f%02.f%02.f%02.f.png';

lat_gis_raw = 49.9975:-0.005:23.0025;
lon_gis_raw = -125.9975:0.005:-65.0025;
[lon_gis,lat_gis] = meshgrid(lon_gis_raw,lat_gis_raw);

%% 2. Main

tic; % Start script timer

if(use_base_ref)
    data_src_print = strcat(data_src,'-base');
else
    data_src_print = data_src;
end

% Directory structure:
% output/np/(run_name)/(variable_name)/(plot_size)
output_path_large = sprintf('%s/%s/%s/%s',output_path_base,run_name,data_type,'large');
output_path_small = sprintf('%s/%s/%s/%s',output_path_base,run_name,data_type,'small');

addpath(path_to_borders);

% Time counters [NOTE: NOT ROBUST TO CROSSING MONTH BOUNDARIES]
time_count = (end_day - start_day + 1)*24 - (start_hour) - (24 - end_hour) + 1;

year = start_year;
month = start_month;
day = start_day;
hour = start_hour;
member = 1;

% Define sample file to read from
timestamp_file = sprintf(datetime_format_file,year,month,day,hour+1);
filename = sprintf(filename_format,domain,member);

% CHANGE THIS IF TRYING TO OPERATE ON OBS

if(~use_obs)
    load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'lon_wrf_trim','lat_wrf_trim'); % Load in preestablished latlon values
    lon = lon_wrf_trim;
    lat = lat_wrf_trim;
else
    t
end

dims = size(lon); % Determine grid dimensions

% Define value lists for cutoffs and rois
cutoffs = cut_min:cut_step:cut_max;
rois = roi_min:roi_step:roi_max;
num_cutoffs = numel(cutoffs);
num_rois = numel(rois);
num_perms = num_cutoffs*num_rois;

% If output folders do not exist, create them
if(~isfolder(output_path_large) || ~isfolder(output_path_small))
    mkdir(output_path_large);
    mkdir(output_path_small);
end

% Chop down custom version of "jet" colormap
jet_dims = size(jet);
jet_lims = [round(jet_dims(1)*jet_blue_percent) jet_dims(1)];
jet_modded = jet;
jet_modded = jet_modded(jet_lims(1):jet_lims(2),:);

% MAIN LOOP
% For each timestamp being analyzed:
for time_idx = 1:time_count
   
% Define current loop's time strings
timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
timestamp_file = sprintf(datetime_format_file,year,month,day,hour);

% Array to hold the running sum sum of member NPs for NEP
sum_np = zeros(dims(1),dims(2),num_perms);

    % For each member:
    for member = 1:num_members
        
        member_string = sprintf('%03.f',member);
        
        if(use_base_ref) % UNDER CONSTRUCTION
            %input_format_base = '%s/%s_%s_d%02.f_%s_%s_%s.mat'; % output path, data source, plot type, domain, member #, data type, datetime
            input_path_base_refl = sprintf('%s/BASE_REF/%s',input_path_br,exp_name);
            load(sprintf(input_format_base_refl,input_path_base_refl,data_src,'base',domain,member_string,data_type,timestamp_file),'wrf_base_refl');
            data = wrf_base_refl;
        elseif(alt_input_format)
            filename = sprinf(filename_format_alt,domain,year,month,day,hour);
            data = ncread(sprintf('%s/%s/%s',input_path_base,member_string,filename),data_name); % Retrieve data
        else
            filename = sprintf(filename_format_wrf,domain,member_idx);       
            data = ncread(sprintf('%s/%s/%s',input_path_base,timestamp_file,filename),data_name); % Retrieve data
        end

        %data = ncread(sprintf('%s/%s/%s',input_path_base,timestamp_file,filename),data_name); % Retrieve data

        if(~use_base_ref)
            % If data needs to be composited down to 2D, do so
            if(composite_wrf)
                data = max(data,[],3);
            end

            % If data needs to be transposed to align with lat/lon grid, do so
            if(transpose_wrf) 
                data = data';
            end

            if(flip_wrf)
                data = flip(data);
            end
        end
     
        cutoff_idx = 1;
        
        

        % Run analysis and plot each variation
        for cutoff = cut_min:cut_step:cut_max % For each cutoff value:

            bin = (data >= cutoff); % Convert data to binary: does it meet the cutoff value or not
            roi_idx = 1;

            for roi = roi_min:roi_step:roi_max % For each ROI value: 
                
                np = neighborhood_prob_2(bin,lon,lat,roi); % Run NP analysis on binary data
                
                % Store for mean (storage only lasts across the current
                % time index)
                np_idx = num_rois*(cutoff_idx-1) + roi_idx;
                sum_np(:,:,np_idx) = sum_np(:,:,np_idx) + np;

                % Plot LARGE version of figure
                f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                h = pcolor(lon,lat,np); % Plot the data
                set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                shading interp; % Smooth out plot from grid-boxes
                c = colorbar('FontSize',axes_font_size); % Make colorbar
                colormap(jet_modded); % Set colors
                caxis([0 1]) % Make sure probability map goes all the way from 0 to 1 for comparison
                
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
                title(sprintf(title_format,upper(data_src),domain,member_string,data_type,cutoff,units,roi/1000,timestamp_title),'FontSize',title_font_size); 
                xlabel('Longitude','FontSize',label_font_size);
                ylabel('Latitude','FontSize',label_font_size);
                set(gca,'Fontsize',axes_font_size);
                saveas(h,sprintf(output_format,output_path_large,data_src_print,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png
                
                % Plot SMALL version of figure
                f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
                title(sprintf(title_format_compressed,cutoff,member_string,timestamp_title),'FontSize',title_font_size);
                saveas(h,sprintf(output_format,output_path_small,data_src_print,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png
                close('all');
                
                roi_idx = roi_idx + 1;
            end
            cutoff_idx = cutoff_idx + 1;
        end
        
        % Plot WRF on WRF grid for reference
        
        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        h = pcolor(lon,lat,data); % Plot the data
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
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);
        saveas(h,sprintf(output_format_raw,output_path_large,'wrf',domain,member_string,'raw',data_type,timestamp_file)); % Save as .png

    end

    %% 3. Make mean plots

    % MEAN OF THE NP (NEP)
    member_string = 'nep';

    cutoff_idx = 1;

    % Run analysis and plot each variation
    for cutoff = cut_min:cut_step:cut_max % For each cutoff value:

        roi_idx = 1;

        for roi = roi_min:roi_step:roi_max % For each ROI value: 

            % Extract correct index for the running sum
            np_idx = num_rois*(cutoff_idx-1) + roi_idx;
            
            % Compute mean of member NPs
            nep = sum_np(:,:,np_idx)/num_members;

            % Plot LARGE
            f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
            h = pcolor(lon,lat,nep); % Plot the data
            set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
            shading interp; % Smooth out plot from grid-boxes
            c = colorbar('FontSize',axes_font_size); % Make colorbar
            colormap(jet_modded); % Set colors
            caxis([0 1]) % Make sure probability map goes all the way from 0 to 1 for comparison

            % Plot state borders
            hold on;
            borders('continental us','black','linewidth',1); 
            hold off;

            % Focus on desired area, remove whitespace
            xlim([w_lim e_lim]);
            ylim([s_lim n_lim]);

            % Apply labels
            title(sprintf(title_format,upper(data_src),domain,upper(member_string),data_type,cutoff,units,roi/1000,timestamp_title),'FontSize',title_font_size); 
            xlabel('Longitude','FontSize',label_font_size);
            ylabel('Latitude','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);
            saveas(h,sprintf(output_format,output_path_large,data_src_print,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png

            % Plot SMALL
            f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
            title(sprintf(title_format_compressed,cutoff,member_string,timestamp_title),'FontSize',title_font_size);
            saveas(h,sprintf(output_format,output_path_small,data_src_print,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png
            close('all');

            roi_idx = roi_idx + 1;
        end

        cutoff_idx = cutoff_idx + 1;
    end
    
    % Time increment
    hour = hour + hour_step;

    if(hour > 23)
        hour = 0;
        day = day + 1;
    end

    if(day > max_day(month))
        day = 1;
        month = month + 1;
    end

    if (month > 12)
        month = 1;
        year = year + 1;
    end

end


% air_err = [22.482 20.470 19.213 18.357]
% conv_err = [22.465 20.605 19.419 18.399]

% Get total runtime and print to stdout
runtime = toc;
hours = floor(toc/3600);
mins = floor((toc/60) - (hours*60));
secs = toc - (hours*3600) - (mins*60);
fprintf('Done. Total script runtime = %02.f:%02.f:%02.f\n',hours,mins,secs)
