% np_analysis.m
% Purpose: Generate neighborhood-probability (NP) diagnostic plots of outputs
% from an IMPACTS WRF Case
% 
% Author(s): Jon Seibert
% Last updated: 8 September 2022
% 
% Inputs: *.nc [Various NetCDF files]
% Outputs: *.png [Set of NP plots with varied ROI/cutoff]
%
% Process: 
% 1) Convert to binary
% 2) ROI checking:
% For each point no more than x gridpoints away:
%   a) Retrieve lat/lon coords. Use Haversine distance to measure distance D between centerpoint and target.
%   b) If D <= ROI, add binary value to list
% 3) Then compute mean of all values within radius
%
% Usage: Edit 1. Settings section to specify the inputs, run entire script
%
% TODO:
%   - Adding T and P contour / colorfield plots
%   - Add functionality for NP analysis of T/P zones below certain
%   threshold for map tracking
%
% NOTES:
%   - CURRENTLY SET UP TO IGNORE THE ENSEMBLE MEAN AND MAKE ITS OWN 
%     (implemented for Reflectivity analysis)
%   - WILL BREAK IF CROSSING MONTH BOUNDARIES WITH DATA (find variable
%     'time_count' to adjust)
%   - Assumes lat/lon consistent across all times and members

%% 0. General Settings

%addpath('/storage/home/jjs5895/IMPACTS/code/borders');

% Plot font sizes
title_font_size = 16;
label_font_size = 14;
axes_font_size = 12;

% Whether to use a shorter title (for gridding smaller versions of plots)
compress_title = false;

max_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 30, 31, 31]; % Number of days in each month

g = 9.80665; % m/s^2 standard gravity at sea level
R = 287.0600676; % J/kg/K, specific gas constant

%% 1. Settings- WRF Data (2020 Case)

% Basic details
domain = 2;
num_members = 1;
data_src = "wrf"; % Dataset label
data_type = "refl"; % Variable being analyzed
units = "dBZ";

% Filepaths
input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2022/1200_test';
%input_path_base = '/storage/home/jjs5895/IMPACTS/data/2020/AIRCFT/fc';
%input_path_base = '/storage/home/jjs5895/IMPACTS/data/2020/AIRCFT900/fc';
output_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/np/Testing';

% NetCDF retrieval
lon_name = "XLONG"; 
lat_name = "XLAT"; 
data_name_refl = "REFL_10CM"; % Variable being analyzed
data_name_temp = "T"; % PERTURBATION POTENTIAL TEMPERATURE
data_name_pres = "PB"; % Base Pressure
data_name_pres_2 = "P"; % Perturbation pressure
data_name_psfc = "PSFC"; % Surface pressure?
data_name_t2m = "T2"; % 2m Temperature
data_name_q = "QVAPOR"; % Water vapor mixing ratio, kg/kg
data_name_gp = "PHB"; % base Geopotential (m^2/s^2)
data_name_gp_2 = "PH"; % Perturbation geopotential
time_name = "Times";

% Note whether transposition or compositing is necessary
transpose_data = true;
flip_data = true;
composite_data = true;

% Choose subset of members
%member_list = [1, 10, 20, 30, 40];

% Specific binary conversion cutoff values
% Units = dBZ
% WARNING: ERROR DISCOVERED IN CODE. NEP NOT RELIABLE IF MULTIPLE CUTOFFS
% ARE USED IN ONE RUN.
cut_min = 40;
cut_step = 5;
cut_max = 40;

% Specific radius of influence values
% Units = m
% WARNING: ERROR DISCOVERED IN CODE. NEP NOT RELIABLE IF MULTIPLE ROIS
% ARE USED IN ONE RUN.
roi_min = 30*1000;
roi_step = 10*1000;
roi_max = 30*1000;

% Date & time
% Range: 2022-01-29-0400 : 2022-01-30-0300
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 7;
start_hour = 4;
end_hour = 4; %18

hour_step = 1; % Hour increment size

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 900; %800
fig_height = 900; %750
fig_width_small = 350;
fig_height_small = 350;
w_lim = -77; % Degrees lat/lon
e_lim = -69.75;
s_lim = 36;
n_lim = 46;

%time_count = 3; %TEMPORARY

cmap = 'reflmap'; % Choice of colorbar colormap
clim_lower = -30; % Colorbar limits
clim_upper = 75;
limit_borders = true;

jet_blue_percent = 0.2; % Cut off coldest (Darkest blues) _% of jet colorbar for ease of analysis

run_name = 'testing_pt';
%run_name = '2020_conv'; % 0 - 1800
%run_name = '2020_aircft'; % %1300 -1800
%run_name = '2020_aircft900';

%wrfinput_d01_2020-02-07_14:00:00_001

%% 2. Main

% Lay out title and filename formats
title_format = '[%s-d%02.f-%s] Neigh. Prob. of %s > %d %s, ROI %d km (%s)'; % Data source (caps), domain, member #, data type, cutoff value, units, ROI/1000, datetime
title_format_compressed = 'NP>%d|%s|%s'; % cutoff, member #, datetime
output_format = '%s/np_%s_d%02.f_%s_%s_%d%s_%dkm_%s.png'; % output path, data source, domain, member #, data type, cutoff value, units, ROI/1000, datetime
datetime_format_file = '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]
nep_save_format = 'nep_storage/nep_%s_%s_%d%s_%dkm'; % data_type, timestamp_file, cutoff, units, roi/1000
filename_format = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #
filename_format_alt = 'wrfinput_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00_%03.f'; % domain, year, month, day, hour, member #

% Directory structure:
% output/np/(run_name)/(variable_name)/(plot_size)
%output_path_large = sprintf('%s/%s/%s/%s',output_path_base,run_name,data_type,'large');
%output_path_small = sprintf('%s/%s/%s/%s',output_path_base,run_name,data_type,'small');
output_path_large = output_path_base;
output_path_small = output_path_base;

% Time counters [NOTE: NOT ROBUST TO CROSSING MONTH BOUNDARIES]
time_count = (end_day - start_day + 1)*24 - (start_hour) - (24 - end_hour) + 1;

year = start_year;
month = start_month;
day = start_day;
hour = start_hour;
member = 1;

% Define sample file to read from
timestamp_file = sprintf(datetime_format_file,year,month,day,hour+1);
%filename = sprintf(filename_format,domain,member);
filename = 'wrf_enkf_output_d02_020';

% Retrieve lat/lon data
%raw_lon = ncread(sprintf('%s/%s/%s', input_path_base,timestamp_file,filename),lon_name);
%raw_lat = ncread(sprintf('%s/%s/%s', input_path_base,timestamp_file,filename),lat_name);
raw_lon = ncread(sprintf('%s/%s', input_path_base,filename),lon_name);
raw_lat = ncread(sprintf('%s/%s', input_path_base,filename),lat_name);

dims = size(raw_lon); % Determine grid dimensions

% If lat/lon data is 1D, make a 2D grid from them (NP code requires 2D)
if(dims(1) == 1 || dims(2) == 1)
    [lon,lat] = meshgrid(raw_lon,raw_lat);
else
    lon = raw_lon;
    lat = raw_lat;
end

% If lat/lon needs to be transposed and/or flipped (both = rotated) to 
% align with proper compass orientation, do so
if(transpose_data) 
    lon = lon';
    lat = lat';
end

if(flip_data)
    lon = flip(lon);
    lat = flip(lat);
end

cutoffs = cut_min:cut_step:cut_max;
rois = roi_min:roi_step:roi_max;
num_cutoffs = numel(cutoffs);
num_rois = numel(rois);
num_perms = num_cutoffs*num_rois*time_count; % 1 for each permutation: NOT YET PROPERLY HANDLED

% Compute sum of member NPs for NEP
sum_np = zeros(dims(1),dims(2),num_perms);

% If output folders do not exist, create them
%if(~isfolder(output_path_large) || ~isfolder(output_path_small))
%    mkdir(output_path_large);
%    mkdir(output_path_small);
%end

jet_dims = size(jet);
jet_lims = [round(jet_dims(1)*jet_blue_percent) jet_dims(1)];
jet_modded = jet;
jet_modded = jet_modded(jet_lims(1):jet_lims(2),:);
%%
% For each timestamp being analyzed:
for time_idx = 1:time_count
    
    %time_idx = 1;
   
timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
timestamp_file = sprintf(datetime_format_file,year,month,day,hour);

    % For each member:
    for member = 1:num_members
        
        %member = 1;
        
        %if(time_idx == 1)
        %    filename = sprintf(filename_format_alt,domain,year,month,day,hour+1,member);
        %else
        %    filename = sprintf(filename_format,domain,member);
        %end
        filename = 'wrf_enkf_output_d02_020';
        member_string = sprintf('%03.f',member);

        %data = ncread(sprintf('%s/%s/%s',input_path_base,timestamp_file,filename),data_name); % Retrieve data
        data = ncread(sprintf('%s/%s',input_path_base,filename),data_name_refl); % Reflectivity
        data_p_base = ncread(sprintf('%s/%s',input_path_base,filename),data_name_pres); % Base pressure
        data_p_pert = ncread(sprintf('%s/%s',input_path_base,filename),data_name_pres_2); % Perturbation pressure
        data_pt_pert = ncread(sprintf('%s/%s',input_path_base,filename),data_name_temp); % Perturbation potential temperature
        data_psfc = ncread(sprintf('%s/%s',input_path_base,filename),data_name_psfc); % PSFC = Surface Level Pressure?? MSLP??
        data_t2m = ncread(sprintf('%s/%s',input_path_base,filename),data_name_t2m); % 2m Temperature
        data_q = ncread(sprintf('%s/%s',input_path_base,filename),data_name_q); % 
        data_gp_base = ncread(sprintf('%s/%s',input_path_base,filename),data_name_gp); % 
        data_gp_pert = ncread(sprintf('%s/%s',input_path_base,filename),data_name_gp_2); % 
        
        data_p = data_p_base + data_p_pert; % 3D Pressure field
        data_gp = data_gp_base + data_gp_pert; % 3D Geopotential field
        data_gph = data_gp./g;
        
        data_pt = data_pt_pert + 300; % Potential temperature
        data_t = data_pt.*((data_p./100000).^(0.287)); % 3D Temperature field
        
        % If data needs to be composited down to 2D, do so
        if(composite_data)
            data = max(data,[],3);
        end

        % If data needs to be transposed to align with lat/lon grid, do so
        if(transpose_data) 
            data = data';
            data_p = permute(data_p,[2 1 3]);
            data_t = permute(data_t,[2 1 3]);
            data_t2m = permute(data_t2m,[2 1 3]);
            data_gp = permute(data_gp,[2 1 3]);
            data_gph = permute(data_gph,[2 1 3]);
            data_q = permute(data_q,[2 1 3]);
        end
        
        if(flip_data)
            data = flip(data);
            data_p = flip(data_p); % Flip DOES behave as intended when applied to a 3D array
            data_t = flip(data_t);
            data_t2m = flip(data_t2m);
            data_gp = flip(data_gp);
            data_gph = flip(data_gph);
            data_q = flip(data_q);
        end
        
        % Compute MSLP
        % equation: P1 = P2*e^(((g/(R*Tv))*(z2-z1))
        % (mean Tv)
        
        Re = 6.3781e6;
        z1 = 0; % m
        P2 = data_p(:,:,1); % Lowest level of P
        Tv = data_t2m.*(1 + 0.608*data_q(:,:,1)); % THIS SHOULD BE A MEAN TV- HOW GET?
        
        %z2 = data_gph(:,:,1)/g; % Geopotential height
        z2 = (data_gp.*Re)/((g*Re) - data_gp); % Height above mean sea level
        
        mslp = P2.*exp((g./(R.*Tv)).*(z2-z1)); % MSLP
        
        % Convert T from K to C
        data_t = data_t - 273.15;

        % Convert P from Pascals (Pa) to Millibars (mb)
        data_p = data_p./100;
        
        data_t_save = data_t;
        cutoff_idx = 1;
%%
        % Run analysis and plot each variation
        for cutoff = cut_min:cut_step:cut_max % For each cutoff value:

            bin = (data >= cutoff); % Convert data to binary
            roi_idx = 1;

            for roi = roi_min:roi_step:roi_max % For each ROI value: 
                
                % Determine if the NEP for this permutation exists.
                %savename = sprintf(nep_save_format,data_type,timestamp_file,cutoff,units,roi/1000);
                %nep_exists = isfile(sprintf('%s.mat',savename));
                
                % If the NEP is already calculated AND this member isn't in the include list, don't process it. (Still needed to slice it for mean calc)
                %if(nep_exists && (~any(ismember(member_list,member)))) 
                %    continue;
                %end

                np = neighborhood_prob_2(bin,lon,lat,roi); % Run NP analysis
                
                % Store for mean
                %np_idx = num_rois*(cutoff_idx-1) + roi_idx;
                np_idx = time_idx; %TEMPORARY
                sum_np(:,:,np_idx) = sum_np(:,:,np_idx) + np;

                % Plot LARGE
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
                xlim([w_lim e_lim]);
                ylim([s_lim n_lim]);

                % Apply labels
                title(sprintf(title_format,upper(data_src),domain,member_string,data_type,cutoff,units,roi/1000,timestamp_title),'FontSize',title_font_size); 
                xlabel('Longitude','FontSize',label_font_size);
                ylabel('Latitude','FontSize',label_font_size);
                saveas(h,sprintf(output_format,output_path_large,data_src,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png
                
                % Plot SMALL
                f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
                title(sprintf(title_format_compressed,cutoff,member_string,timestamp_title),'FontSize',title_font_size);
                saveas(h,sprintf(output_format,output_path_small,data_src,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png
                %close('all');
                
                roi_idx = roi_idx + 1;
            end
            cutoff_idx = cutoff_idx + 1;
        end
        
        % Plot WRF on WRF grid for reference
        
        title_format_raw =  '[%s|d0%d|%s] %s: %s(%s)|%s'; % data source, domain, plot type, member #/mean, data type, units, timestamp
        output_format_raw = '%s/%s_d0%d_%s_%s_%s_%s.png'; % output path, domain, mem/mean, fit/dif, data type, timestamp
    
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
        saveas(h,sprintf(output_format_raw,output_path_large,'wrf',domain,member_string,'raw',data_type,timestamp_file)); % Save as .png
        
        %%
        
        data_t = data_t_save;
        
        %%
        
        data_t = imgaussfilt(data_t,3);
        
        %% Plot T and P contours on WRF grid

        height_index = 5;
        
        
        % Formats
        title_format_raw =  '[%s|d0%d|%s] %s:%s(%s)|z%d|%s'; % data source, domain, plot type, member #/mean, data type, units, timestamp
        output_format_raw = '%s/%s_d0%d_%s_%s_%s_z%d_%s.png'; % output path, domain, mem/mean, fit/dif, data type, timestamp
    
        % Plot LARGE - T
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        [M,c] = contour(lon,lat,data_t(:,:,height_index),[-20 -16 -12 -8 -4 0 4 8 12 16 20]); % Plot the data
        %clabel(M,c);
        %h = contour(lon,lat,data_t(:,:,height_index),'ShowText','on'); % Plot the data
        cb = colorbar('FontSize',axes_font_size); % Make colorbar
        caxis([-20 20]);
        colorbar('Ticks',[-20 -16 -12 -8 -4 0 4 8 12 16 20],...
         'TickLabels',{'-20', '-16', '-12', '-8', '-4', '0', '4', '8', '12', '16', '20'})
        colormap(tmap); % Set colors
        c.LineWidth = 2;
        
        % Plot state borders
        hold on;
        %borders('continental us','white','linewidth',1); 
        borders('continental us','Color',[0.5 0.5 0.5],'linewidth',1); 
        hold off;

        % Focus on desired area, remove whitespace
        if(limit_borders)
            xlim([w_lim e_lim]);
            ylim([s_lim n_lim]);
        end

        % Apply labels
        title(sprintf(title_format_raw,upper('wrf'),domain,member_string,'Contours','T','C',height_index,timestamp_title),'FontSize',title_font_size); 
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        %set(gcf, 'InvertHardCopy', 'off');
        %set(gca,'color','Black');
        saveas(f,sprintf(output_format_raw,output_path_large,'wrf',domain,member_string,'raw','T',height_index,timestamp_file)); % Save as .png
        %% plot LARGE - P
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        [M,c] = contour(lon,lat,data_p(:,:,height_index)); % Plot the data
        %h = contour(lon,lat,data_p(:,:,height_index),'ShowText','on'); % Plot the data
        cb = colorbar('FontSize',axes_font_size); % Make colorbar
        %colormap(cmap); % Set colors
        c.LineWidth = 1;

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
        title(sprintf(title_format_raw,upper('wrf'),domain,member_string,'Contours','P','mb',height_index,timestamp_title),'FontSize',title_font_size); 
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        saveas(f,sprintf(output_format_raw,output_path_large,'wrf',domain,member_string,'P','mb',height_index,timestamp_file)); % Save as .png
        

    end

    %% 3. Make mean plots

    % MEAN OF THE NP (NEP)
    member_string = 'nep';

    cutoff_idx = 1;

    % Run analysis and plot each variation
    for cutoff = cut_min:cut_step:cut_max % For each cutoff value:

        roi_idx = 1;

        for roi = roi_min:roi_step:roi_max % For each ROI value: 

            %np_idx = num_rois*(cutoff_idx-1) + roi_idx;
            np_idx = time_idx; %TEMPORARY
            % Determine if the NEP for this permutation exists. If so, load it.
            % Otherise, compute from sum and save it as a file. 
            % (Sum should be 0 if the saved file exists.)
            savename = sprintf(nep_save_format,data_type,timestamp_file,cutoff,units,roi/1000);
            %nep_exists = isfile(sprintf('%s.mat',savename));
            %if(nep_exists)
            %    nep = load(savename,'nep');
            %else
                nep = sum_np(:,:,np_idx)/num_members;
                save(savename,'nep');
            %end

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
            title(sprintf(title_format,upper(data_src),domain,upper(member_string),data_type,cutoff,units,roi/1000,timestamp_title)); 
            xlabel('Longitude','FontSize',label_font_size);
            ylabel('Latitude','FontSize',label_font_size);
            saveas(h,sprintf(output_format,output_path_large,data_src,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png

            % Plot SMALL
            f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
            title(sprintf(title_format_compressed,cutoff,member_string,timestamp_title));
            saveas(h,sprintf(output_format,output_path_small,data_src,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png
            %close('all');

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
