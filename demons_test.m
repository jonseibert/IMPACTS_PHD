%% 0. General Settings

% Local or server?
%local = false; % true = running on local device, false = running on server
local = true;

% Print out processing notes?
print_warnings = false;
progress_bar = false;

script_name = 'snowband_tracking.m';

% Plot font sizes
title_font_size = 16;
label_font_size = 14;
axes_font_size = 14;

% Lay out title and filename formats
% Input filenames, plot titles, output filenames, datetime strings
% Compressed = for small plots
filename_format = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #
filename_format_alt = 'wrfinput_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00_%03.f'; % domain, year, month, day, hour, member #
%title_format = '[%s-d%02.f-%s] Neigh. Prob. of %s > %d %s, ROI %d km (%s)'; % Data source (caps), domain, member #, data type, cutoff value, units, ROI/1000, datetime
%title_format_compressed = 'NP>%d|%s|%s'; % cutoff, member #, datetime
%output_format = '%s/np_%s_d%02.f_%s_%s_%d%s_%dkm_%s.png'; % output path, data source, domain, member #, data type, cutoff value, units, ROI/1000, datetime
datetime_format_file = '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

%% 1. Settings- Data (WRF 2020 Case)

% Reference filenames: wrf_enkf_output_d02_001_2020-02-07_14-00-00.nc

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020';
    output_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/tracking';
    path_to_borders = './borders';
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/AIRCFT/fc';
    output_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/tracking';
    path_to_borders = '/storage/home/jjs5895/projects/IMPACTS/code/borders'; % Specify path to borders.m
end

% Denotes output folder name
run_name = '2020_aircft';
exp_name = 'AIR';

% Basic details
domain = 2;
%num_members = 40;
num_members = 1; %TESTING
data_src = "wrf"; % Dataset label
data_type = "refl"; % Variable being analyzed
units = "dBZ";

% Date & time
% Range: 2022-01-29-0400 : 2022-01-30-0300
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 7;
start_hour = 14;
end_hour = 15;

% Hour increment size
hour_step = 1; 

% Figure specs
fig_x = 100;
fig_y = 100;
%fig_width = 740; % Use with -77 limit
fig_width = 925;
fig_height = 900;
fig_width_small = 350;
fig_height_small = 350;

% Spatial limits of analysis and plot (degrees lat/lon)
%w_lim = -77;
w_lim = -79;
e_lim = -69.75;
s_lim = 36;
n_lim = 46;
limit_borders = true; % Whether to apply spatial limits to plot
trim_e = false; % Whether to trim the eastern edge

cmap = 'reflmap'; % Choice of colorbar colormap
clim_lower = -30; % Colorbar limits
clim_upper = 75;




%%

output_path_large = sprintf('%s/demons',output_path_base);
datetime_format_file = '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]

year = 2020;
month = 2;
day = 7;
hour = 14;

member_string = '001';
timestamp_file = sprintf(datetime_format_file,year,month,day,hour);


load(sprintf('intermediate/clusters_%s_%s',member_string,timestamp_file));

lon_a = lon_wrf_trim;
lat_a = lat_wrf_trim;
clusters_a = pt_clusters;
data_a = data_thresh;

hour = 15;
timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
load(sprintf('intermediate/clusters_%s_%s',member_string,timestamp_file));

lon_b = lon_wrf_trim;
lat_b = lat_wrf_trim;
clusters_b = pt_clusters;
data_b = data_thresh;

a_sum = sum(data_a,"all");
b_sum = sum(data_b,"all");

data_a_norm = data_a./a_sum;
data_b_norm = data_b./b_sum;

tformEstimate = imregcorr(data_a,data_b,"translation");
Rfixed = imref2d(size(data_b));
a_corr = imwarp(data_a,tformEstimate,"OutputView",Rfixed);

%imshow(a_corr)

[D,a_trans] = imregdemons(a_corr,data_b);
%[D,a_trans] = imregdemons(data_a_norm,data_b_norm);
%[D,a_trans] = imregdemons(data_a,data_b,'AccumulatedFieldSmoothing',1.5);
%[D,a_trans] = imregdemons(data_a,data_b,'PyramidLevels',2);
%[D,a_trans] = imregdemons(data_a,data_b,'AccumulatedFieldSmoothing',1.5,'PyramidLevels',2);

%a_trans = a_trans*a_sum;

data_a(data_a == 0) = -30;
data_b(data_b == 0) = -30;
a_trans(a_trans == 0) = -30;

%%

 for plot_idx = 1:3
        if(plot_idx == 1)
            plot_type = 'base14';  % Specify fit, dif, trim
            data_to_plot = data_a;
            hour = 14;
            timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
            timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
        elseif(plot_idx == 2)
            plot_type = '14-demonshift-15';  % Specify fit, dif, trim
            data_to_plot = a_trans;
            hour = 15;
            timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
            timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
        elseif(plot_idx == 3)
            plot_type = 'base15';
            data_to_plot = data_b;
            hour = 15;
            timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
            timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
        end
    
        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        h = pcolor(lon_wrf_trim,lat_wrf_trim,data_to_plot); % Plot the data
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

        title_format_raw =   '[%s|d0%d|%s] Raw: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
        output_format = '%s/%s_%s_d0%d_%s_%s_%s_%s_%s.png'; % output path, exp name 1, exp name 2, domain, mem/mean, bao_short, plot_type, data type, timestamp
        title(sprintf(title_format_raw,exp_name,domain,member_string,data_type,units,timestamp_title),'FontSize',title_font_size);
        output_name_b = 'raw';

        saveas(h,sprintf(output_format,output_path_large,lower(exp_name),output_name_b,domain,member_string,'a',plot_type,data_type,timestamp_file)); % Save as .png

        close("all");
end
