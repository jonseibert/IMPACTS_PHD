% wrf_3d_ref_to_base.m
%
% Purpose: Interpolate WRF output 3D reflectivity to NEXRAD 0.5 degree base reflectivity
% 
% Author(s): Jon Seibert
% Last updated: 26 Sept 2024
% 
% Inputs: [Ens Mem WRF Analysis Step].nc
% Outputs: NEXRAD_radar_heights_wrf_grid_trim.mat, [Ens Mem Base Refl].mat
%
% Usage: Edit sections "0. Script Controls" and "1. General Settings" to specify the inputs, run entire script
%
% TODO: Update documentation
%  
% Dependencies: nexrad_stations.mat, latlon_wrf_trim.mat,
%   linear_interpolate.m, borders.m, reflmap.m
%
% Local functions: align_data()
%
% NOTES:
%   - Make sure all units match! Was a temporary issue with elevations in ft
script_name = 'wrf_3d_ref_to_base.m';

%% 0. Script Controls

% Local or server? [Now handled automatically]
%local = 0; % 1 = running on local device, 0 = running on server

% Denotes output folder name
run_name = 'air_base_v3';

% Which experiment to use
exp_choice = 1; 

% Date & time
% Ranges: 2020-02-07-1400 : 2020-02-07-1800 : 2020-02-08-0000
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 7;
start_hour = 14;
end_hour = 18;

num_members = 40;
%num_members = 2; %TESTING

remake_beam_heights = 0;
remake_mat_files = 0;

show_plots = 0;
plot_heights = 0;

%% 1A. General Settings + WRF
% Shouldn't have to change the below for most runs

% Reference filenames: wrf_enkf_output_d02_001_2020-02-07_14-00-00.nc
%                      data/2020/AIRCFT-F-18Z/003/wrfout_d01_2020-02-07_19:00:00
%                      data/2020/AIRCFT/fc/202002071400/wrf_enkf_output_d02_016

% Detect working environment
local_wd = "C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code";
server_wd = "/storage/home/jjs5895/projects/IMPACTS/code";
working_dir = strrep(pwd,"\","/");
if(working_dir == local_wd)
    local = 1;
elseif(working_dir == server_wd)
    local = 0;
else
    error('Unknown environment detected.')
end

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020';
    input_path_gis = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    intermediate_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/testing';
    path_to_code = "C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code";
    path_to_extra_code = './downloaded_code';
    input_path_ecmwf = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Data/ECMWF';
    storage_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/output/base_ref/mats';
    % Experiment paths
    exp_path_1 = 'AIR';
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
    storage_path = '/storage/home/jjs5895/projects/IMPACTS/output/base_ref/mats/';
    % Experiment paths
    exp_path_1 = 'AIRCFT/fc';
    exp_path_2 = 'CONV/fc/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
end

addpath(path_to_extra_code); % Make sure Matlab can see the addon code

% Experiment names
exp_name_1 = 'AIR';
exp_name_2 = 'CONV';
exp_name_3 = 'AIR-F';
exp_name_4 = 'CONV-F';
exp_name_5 = 'NODA';

alt_input_exps = [3 4 5]; % Input is set up differently for these experiments

% WRF Details
domain = 2;
data_type = "REFL"; % Variable being analyzed
units = "dBZ";

% NetCDF retrieval
lon_name = "XLONG"; 
lat_name = "XLAT"; 
data_name = "REFL_10CM"; % Variable being analyzed
data_name_gp = "PHB"; % Base Geopotential (m^2/s^2)
data_name_gp_prime = "PH"; % Perturbation Geopotential (m^2/s^2)
data_name_terrain_height = "HGT"; % Terrain height above sea level (m)

% Note whether transposition or compositing is necessary to properly align
% the data with desired grid
transpose_data = true;
flip_data = true;
z_composite_data = false;

% Hour increment size
hour_step = 1; 

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 925;
fig_height = 900;
fig_width_small = 525;
fig_height_small = 500;

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

%grid_res = 3000; % 3 km

if(local)
    % Figure font sizes
    title_font_size = 18;
    label_font_size = 18;
    axes_font_size = 16;
    contour_font_size = 14;
else
    % Figure font sizes
    title_font_size = 24;
    label_font_size = 24;
    axes_font_size = 22;
    contour_font_size = 20;
end

% Lay out title and filename formats
% Input filenames, plot titles, output filenames, datetime strings
% Compressed = for small plots
filename_format = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #
filename_format_alt = 'wrfout_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00'; % domain, year, month, day, hour
output_format = '%s/%s_%s_d%02.f_%s_%s_%s.png'; % output path, data source, plot type, domain, member #, data type, datetime
output_format_data = '%s/%s_%s_d%02.f_%s_%s_%s.mat'; % output path, data source, plot type, domain, member #, data type, datetime
datetime_format_file = '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

%% 1B. Settings: NEXRAD Base Refl

% N0R product is angle 0.5, 16 levels/ 230 km (short range base refl.) N0Q
% is long range: 460 km, still 0.5. we're using N0Q. 8 bit Base Reflectivity 0.5dbz
beam_angle = 0.5;
max_dist = 230000; % 230 km NOW USING SHORTER MAX TO AVOID RADAR FOLDING

g = 9.8; % m/s^2;
Re = 6371000; % radius of Earth in m

radar_heights_filename = sprintf('%s/NEXRAD_radar_heights_wrf_grid_trim.mat',intermediate_path);

%% 2. Setup

if(~show_plots)
    set(groot,'DefaultFigureVisible','off') % Turn off figure popups for local
end

% Check for correct execution mode
working_dir = strrep(pwd,"\","/");
if(working_dir ~= path_to_code)
    error("Operating mode does not match execution directory.")
end

% Define run names and flags
switch exp_choice
    case 1
        input_path = sprintf('%s/%s',input_path_base,exp_path_1);
        exp_name = exp_name_1;
    case 2
        input_path = sprintf('%s/%s',input_path_base,exp_path_2);
        exp_name = exp_name_2;
    case 3
        input_path = sprintf('%s/%s',input_path_base,exp_path_3);
        exp_name = exp_name_3;
    case 4
        input_path = sprintf('%s/%s',input_path_base,exp_path_4);
        exp_name = exp_name_4;
    case 5
        input_path = sprintf('%s/%s',input_path_base,exp_path_5);
        exp_name = exp_name_5;
    otherwise
        error('Unrecognized experiment choice');
end

data_src = lower(exp_name);

% Announce operating mode
fprintf('Starting %s in %s mode (%s)\n',upper(script_name),mode,exp_name);
fprintf('Run designation: (%s)\n',run_name);
tic;

if(ismember(exp_choice,alt_input_exps))
    alt_input_format =  true;
else
    alt_input_format = false;
end

% Directory structure:
% output/np/(run_name)/(variable_name)/(plot_size)
output_path_large = sprintf('%s/%s/%s/%s',output_path_base,run_name,upper(data_type),'large');
output_path_small = sprintf('%s/%s/%s/%s',output_path_base,run_name,upper(data_type),'small');

% If output folders do not exist, create them
if(~isfolder(output_path_large) || ~isfolder(output_path_small))
    mkdir(output_path_large);
    mkdir(output_path_small);
end

% Load NEXRAD station and WRF grid data
load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat')); % lon_wrf_trim, lat_wrf_trim, n_idx,s_idx,w_idx,e_idx
load(sprintf('%s/%s',intermediate_path,'nexrad_stations.mat')); % wban, station_ids,_station_names,lon_deg,lat_deg,elevations,tower_heights

% Time counters [NOTE: NOT ROBUST TO CROSSING MONTH BOUNDARIES]
max_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 30, 31, 31]; % Number of days in each month
time_count = (end_day - start_day + 1)*24 - (start_hour) - (24 - end_hour) + 1;
loop_count = 1;

dims = size(lon_wrf_trim);
ylen_wrf = dims(1);
xlen_wrf = dims(2);
zlen_wrf = 50;
num_radars = length(station_ids);

%% 3. NEXRAD stations: compute radar grids

if(~remake_beam_heights && (exist(radar_heights_filename,'file') == 2))
    load(radar_heights_filename);
else
    radar_height_grid = zeros(num_radars,ylen_wrf,xlen_wrf); % Radar idx, lon, lat = beam height value
    radar_height_grid(:,:,:) = NaN; % Any that are set to non-NaN are within the radius

    for station_idx = 1:num_radars
        % Retrieve data
        rad_lon = lon_deg(station_idx);
        rad_lat = lat_deg(station_idx);
        rad_id = station_ids(station_idx);
        rad_name = station_names(station_idx);
        rad_elev = elevations(station_idx);
        rad_th = tower_heights(station_idx);

        % Compute wrf grid indices and beam height
        for y = 1:ylen_wrf
            for x = 1:xlen_wrf
                wrf_lon_val = lon_wrf_trim(y,x);
                wrf_lat_val = lat_wrf_trim(y,x);
                D = haversine_distance(rad_lon,rad_lat,wrf_lon_val,wrf_lat_val);
                if(D <= max_dist)
                    % Compute beam height H from slant distance R and ground distance D
                    R = D/(cosd(beam_angle));
                    %H = R*sind(beam_angle) + (R^2)/(2*(4/3)*Re) + rad_elev + rad_th; % Jon version
                    H = sqrt(R^2 + ((4/3)*Re)^2 + 2*R*(4/3)*Re*sind(beam_angle)) - (4/3)*Re + rad_elev + rad_th; % Keenan version adjusted for m
                    radar_height_grid(station_idx,y,x) = H;
                end
            end
        end
    end

    save(sprintf('%s/%s',intermediate_path,radar_heights_filename),'radar_height_grid');
end

%% Set up holder array for ensemble members
ens_mem_base_refl = zeros(ylen_wrf,xlen_wrf,num_members);

storage_path_full = sprintf('%s/%s',storage_path,lower(exp_name));

year = start_year;
month = start_month;
day = start_day;
hour = start_hour;

first_loop = true;

% For each timestamp being analyzed:
for time_idx = 1:time_count

    timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
    timestamp_file = sprintf(datetime_format_file,year,month,day,hour);

    % For each ensemble member:
    for member_idx = 1:num_members+1
    
        %% 4A. Read in data 
        plot_type = 'base';
        
        % Define current loop's time strings
        if(member_idx == num_members+1)
            member_string = "mean";
            wrf_base_refl = mean(ens_mem_base_refl,3);
            fprintf('Processing: Ensemble Mean, %s\n',timestamp_title);
        else
        
            member_string = sprintf('%03.f',member_idx);
    
            if(~alt_input_format)
                filename = sprintf(filename_format,domain,member_idx);
                input_path_full = sprintf('%s/%s/%s',input_path,timestamp_file,filename);
                storage_path_full = sprintf('%s/%s',storage_path_full,storage_filename); % WIP LOCATION
            else
                filename = sprintf(filename_format_alt,domain,year,month,day,hour);
                input_path_full = sprintf('%s/%03.f/%s',input_path,member_idx,filename);
            end
    
            % Retrieve data
            refl = ncread(input_path_full,data_name);
            gp_base = ncread(input_path_full,data_name_gp); % Base geopotential (m^2/s^2)
            gp_pr = ncread(input_path_full,data_name_gp_prime); % Perturbation geopotential
            terrain_height = ncread(input_path_full,data_name_terrain_height); % 
    
            fprintf('Processing: Member %d, %s\n',member_idx,timestamp_title);
    
            % Process input data to have correct limits and orientation
            refl = align_data(refl,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            gp_base = align_data(gp_base,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            gp_pr = align_data(gp_pr,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            terrain_height = align_data(terrain_height,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            
            refl_composite = max(refl,[],3);
    
            % Convert to model height
            model_height = ((gp_base+gp_pr)./g) - terrain_height; % model height above the *ground* (m)
            model_height_sl = model_height + terrain_height; % model height above sea level (m)

            %% 4B. Compute base reflectivity

            
    
            % For all latlon:
            %   For all radar grids:
            %       If there is a non-NaN beam height here, interpolate two vertically
            %       closest WRF refl points to get the value here, add to list
            %   Take max value of list, set into grid
    
            wrf_base_refl = zeros(ylen_wrf,xlen_wrf);
            wrf_base_refl(:) = -30;
            wrf_height_map = zeros(ylen_wrf,xlen_wrf);
            
            save_row_idx = 277;
            save_col_idx = 112;
    
            for y = 1:ylen_wrf
                for x = 1:xlen_wrf
                    interp_dbz = zeros(num_radars,1);
                    interp_dbz(:) = NaN;
                    for station_idx = 1:num_radars
                        beam_height = radar_height_grid(station_idx,y,x);
                        if(~isnan(beam_height))
                            for z = 1:zlen_wrf
                                mh = model_height(y,x,z);
                                if(mh == beam_height)
                                    interp_dbz(station_idx) = refl(y,x,z);
                                    break;
                                elseif(mh > beam_height)
                                    if(z == 1)
                                        interp_dbz(station_idx) = refl(y,x,z);
                                    else
                                        interp_dbz(station_idx) = linear_interpolate(model_height(y,x,z-1),mh,beam_height,refl(y,x,z-1),refl(y,x,z));
                                    end
                                    break;
                                end
                            end
                        end
                    end
                    [max_val,max_station_idx] = max(interp_dbz);
                    wrf_base_refl(y,x) = max_val;
                    wrf_height_map(y,x) = radar_height_grid(max_station_idx,y,x);
                    
                    if(y == save_row_idx && x == save_col_idx)
                        interp_dbz_save = interp_dbz;
                    end
                end
            end
            ens_mem_base_refl(:,:,member_idx) = wrf_base_refl;
            
            save(sprintf(output_format_data,output_path_large,data_src,plot_type,domain,member_string,upper(data_type),timestamp_file),'wrf_base_refl');
        end

        

        %% 5A. Plot
        
        % Plot new base reflectivity
        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        h = pcolor(lon_wrf_trim,lat_wrf_trim,wrf_base_refl); % Plot the data
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
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);
        title_format_raw = '[%s|d0%d|%s] Interpolated Base Reflectivity: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
        title(sprintf(title_format_raw,exp_name,domain,member_string,upper(data_type),units,timestamp_title),'FontSize',title_font_size);

        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,upper(data_type),timestamp_file)); % Save as .png
        
        % Plot SMALL
        f.Position = [fig_x fig_y 525 500]; % Shrink figure
        title_format_compressed = '[%s] IBR %s(%s)'; % exp name, domain, member #/mean, data type, units, timestamp
        title(sprintf(title_format_compressed,exp_name,upper(data_type),units),'FontSize',title_font_size);
        saveas(h,sprintf(output_format,output_path_small,data_src,plot_type,domain,member_string,upper(data_type),timestamp_file)); % Save as .png
        
        close('all');
        
        %% 5B. Plot vertical composite for reference

        plot_type = 'raw';

        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        h = pcolor(lon_wrf_trim,lat_wrf_trim,refl_composite); % Plot the data
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
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);
        title_format_raw = '[%s|d0%d|%s] Vertical Composite: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
        title(sprintf(title_format_raw,exp_name,domain,member_string,upper(data_type),units,timestamp_title),'FontSize',title_font_size);

        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,upper(data_type),timestamp_file)); % Save as .png
        
        close('all');
        
        %% 5C. Test plot: 2D "used" beam height map
        
        if(plot_heights)
        
            plot_type = 'hgt';

            wrf_height_min = 0;
            wrf_height_max = max(model_height,[],'all');

            f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
            h = pcolor(lon_wrf_trim,lat_wrf_trim,wrf_height_map); % Plot the data
            set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
            shading interp; % Smooth out plot from grid-boxes
            c = colorbar('FontSize',axes_font_size); % Make colorbar
            colormap(jet); % Set colors
            caxis([wrf_height_min wrf_height_max]);

            % Plot state borders
            hold on;
            borders('continental us','black','linewidth',1); 

            % Plot radar locations
            scatter(lon_deg,lat_deg,100,'pentagram','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]);
            hold off;

            % Focus on desired area, remove whitespace
            if(limit_borders)
                xlim([w_lim e_lim]);
                ylim([s_lim n_lim]);
            end

            % Apply labels
            xlabel('Longitude','FontSize',label_font_size);
            ylabel('Latitude','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);

            title_format_height = '[%s|d0%d|%s] Utilized beam height (%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
            title(sprintf(title_format_height,exp_name,domain,member_string,'m',timestamp_title),'FontSize',title_font_size);

            saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,'hgt',timestamp_file)); % Save as .png
            
            %% 5D. Vertical Profile at chosen location
            
            refl_to_plot = squeeze(refl(save_row_idx,save_col_idx,:));
            height_staggered = squeeze(model_height(save_row_idx,save_col_idx,:));
            height_to_plot = zeros(size(refl_to_plot));
            for idx = 1:length(height_to_plot)
                height_to_plot(idx) = (height_staggered(idx) + height_staggered(idx+1))/(2*1000);
            end

            f = figure('Position',[fig_x fig_y 800 900]); % Create initial blank figure
            h = plot(refl_to_plot,height_to_plot);
            hold on;
            h = scatter(interp_dbz_save,(radar_height_grid(:,save_row_idx,save_col_idx)/1000),'filled');
            %h = pcolor(lon_wrf_trim,lat_wrf_trim,wrf_height_map); % Plot the data
            %set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
            %shading interp; % Smooth out plot from grid-boxes
            %c = colorbar('FontSize',axes_font_size); % Make colorbar
            %colormap(cmap); % Set colors
            %caxis([clim_lower clim_upper]);
            ylim([0 20]);

            % Apply labels
            xlabel('Reflectivity(dBZ)','FontSize',label_font_size);
            ylabel('Elevation (km)','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);
            legend('WRF dBZ','Interpolated dBZ');

            plot_type = 'vpf';

            title_format_height = '[%s|d0%d|%s] VPF Big Stone Beach DE: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
            title(sprintf(title_format_height,exp_name,domain,member_string,upper(data_type),units,timestamp_title),'FontSize',title_font_size);

            saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,upper(data_type),timestamp_file)); % Save as .png
            hold off;
            
            close('all');
        end
    end

    %% 6. Time increment
    
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

set(groot,'DefaultFigureVisible','on');
fprintf('%s run complete.\n',upper(script_name));
% END

%% Local functions

% Process input data
function out_data = align_data(in_data,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx)

    temp_data = in_data;

    % If data needs to be composited down to 2D, do so
    if(z_composite_data)
        temp_data = max(temp_data,[],3);
    end

    % If data needs to be transposed and/or flipped to align with lat/lon grid, do so
    if(transpose_data) 
        temp_data = permute(temp_data,[2 1 3]);
    end

    if(flip_data)
        temp_data = flip(temp_data);
    end

    % Trim data to desired spatial limits
    out_data = temp_data(n_idx:s_idx,w_idx:e_idx,:);
end
