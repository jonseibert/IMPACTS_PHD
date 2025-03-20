% % plot_error_tables.m
%
% Purpose: 
% 
% Author(s): Jon Seibert
% Last updated: 26 March 2023
% 
% Inputs: 
% Outputs: 
% Dependencies: 
%
% NOTES:
%  - TIME_COUNT is not robust when crossing month boundaries- replace with
%    manual time-tick count if using it to do so.
%
% TODO:
%  - Add cross-experiment mean comparison plot- include stdevs?

script_name = 'check_error_tables.m';

%% 0A. Script Controls

first_pass = 1;

bias_max = 0;
bias_min = 0;
bs_max = 0;
bs_min = 0;
csi_max = 0;
csi_min = 0;
ets_max = 0;
ets_min = 0;
rmse_max = 0;
rmse_min = 0;
sd_max = 0;
sd_min = 0;
sde_max = 0;
sde_min = 0;

%%
for subdomain = 0:2
%subdomain = 0;
batch_name = 'value_check';
%batch_name = 'alpha_sbd0';
%batch_name = 'beta_sbd0';
exclude_noda = 0;

% Local or server execution mode
%local = false; % true = running on local device, false = running on server
local = true;

% Date & time
% Range: 2020-02-07-1400 : 2020-02-07-1800; 14-00, 18-00
% 1400 is same across all experiments
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 8;
start_hour = 14;
end_hour = 0;

show_plots = false;
show_progress = true; % If true, prints out progress bar/debug messages

% Plots

var_name_list = ["RMSE","Bias","ETS20","ETS35","CSI20","CSI35","SD","SDE","BS20","BS35"];

%run_name_list = ["air_obs_bulk","conv_obs_bulk","air-f_obs_bulk","conv-f_obs_bulk","noda_obs_bulk"];
%exp_name_list = ["AIR","CONV","AIR-F","CONV-F","NODA"];
%run_name_list = ["air_obs_bulk","conv_obs_bulk","noda_obs_bulk"];
%exp_name_list = ["AIR","CONV","NODA"];
%run_name_list = ["air-f_obs_bulk","conv-f_obs_bulk","noda_obs_bulk"];
%exp_name_list = ["AIR-F","CONV-F","NODA"];

if(subdomain == 0)
    %run_name_list = ["air_obs_sbd0","conv_obs_sbd0","air-f_obs_sbd0","conv-f_obs_sbd0","noda_obs_sbd0"];
    run_name_list = ["air_bulk_2_sbd0","conv_bulk_2_sbd0","air-f_bulk_2_sbd0","conv-f_bulk_2_sbd0","noda_bulk_2_sbd0"];
elseif(subdomain == 1)
    run_name_list = ["air_bulk_2_sbd1","conv_bulk_2_sbd1","air-f_bulk_2_sbd1","conv-f_bulk_2_sbd1","noda_bulk_2_sbd1"];
elseif(subdomain == 2)
    run_name_list = ["air_bulk_2_sbd2","conv_bulk_2_sbd2","air-f_bulk_2_sbd2","conv-f_bulk_2_sbd2","noda_bulk_2_sbd2"];
end
exp_name_list = ["AIR","CONV","AIR-F","CONV-F","NODA"];

exclude_comp_vars = ["BS20","BS35","SD","SDE"];
zero_line_vars = ["ETS20","ETS35","Bias"];

%run_name = "air_obs_bulk";
%exp_name = "AIR";

plot_value_mean = true;
plot_ens_mean = true;

combine_series_plumes = true;
squish = true;

hard_axes_caps = 1;
exclude_rmse_cap = 1;

run_all_exps = true;

override_f18 = 1; % Replace 18z error values on AIRF snd CONVF with NaN, sets equal to prior for exp comp. 
extend_noda = 1; % correct missing NoDA 14z values to shared equal field

sd_max = 0.7;

plume_linewidth = 1.2;
addon_linewidth = 2;
standard_linewidth = 3.5;
expc_linewidth = 4;
standard_ebsize = 12;





%% 0B. General Settings

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/ev_tables';
    intermediate_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = sprintf('C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/ev_plots/%s',batch_name);
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/ptp/ev_tables';
    intermediate_path = '/storage/home/jjs5895/projects/IMPACTS/intermediate';
    output_path_base = sprintf('/storage/home/jjs5895/projects/IMPACTS/output/ev_plots/%s',batch_name);
end

hour_step = 1; % Hour increment size

domain = 2;
num_members = 40;

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 925;
fig_height = 900;
fig_width_small = 800;
fig_height_small = 400;
fig_width_small_plume = 800;
fig_height_small_plume = 400;
fig_width_tile = 650;
fig_height_tile = 294;



% Figure font sizes
title_font_size = 18;
label_font_size = 18;
axes_font_size = 16;
title_font_size_tile = 18;
label_font_size_tile = 18;
axes_font_size_tile = 16;

data_type = "refl"; % Variable being analyzed
units = "dBZ";

err_rowNames = 1:(num_members+6);
err_rowNames = string(err_rowNames);
err_rowNames(num_members+1) = "EMean";
err_rowNames(num_members+2) = "Median";
err_rowNames(num_members+3) = "MeanEV";
err_rowNames(num_members+4) = "Min";
err_rowNames(num_members+5) = "Max";
err_rowNames(num_members+6) = "StdEV";

ens_mean_idx = 41;
ens_med_idx = 42;
mean_ev_idx = 43;
min_idx = 44;
max_idx = 45;
std_idx = 46;

num_rows = 46;
num_vars = 10;

comp_line_colors = ["#0072BD","#D95319","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];
error_dec = 3;

%% Formats

% Lay out title and filename formats
title_format_fit =   '[%s-%s|d0%d|%s] %s: %s(%s)|%s'; % exp name 1, exp name 2, domain, member #/mean, background/analysis/obs (bao), bao, plot_type, data type, units, timestamp
title_format_dif =   '[%s-%s|d0%d|%s] %s: %s(%s)|%s\nRMSE = %0.3f, Bias = %0.3f'; % exp name 1, exp name 2, domain, member #/mean, background/analysis/obs (bao), bao, plot_type, data type, units, timestamp, err, bias
title_format_trim =  '[%s] %s: %s(%s)|%s'; % data source, plot_type, data type, units, timestamp
title_format_raw =   '[%s|d0%d|%s] Raw: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
title_format_sd =   '[%s|d0%d] %s: %s(%s)|%s'; % exp name, domain, sd_type, data type, units, timestamp
title_format_small = '[%s-%s] %s-%s-%s'; % exp name 1, exp name 2,plot_type,data type,timestamp
output_format =         '%s/%s_%s_d0%d_%s_%s_%s_%s_%s.png'; % output path, exp name 1, exp name 2, domain, mem/mean, bao_short, plot_type, data type, timestamp
output_format_small =   '%s/%s_%s_d0%d_%s_%s_%s_%s_%s_small.png'; % output path, exp name 1, exp name 2, domain, mem/mean, bao_short, plot_type, data type, timestamp
input_format_base_refl = '%s/%s_%s_d%02.f_%s_%s_%s.mat'; % output path, data source, plot type, domain, member #, data type, datetime
datetime_format_file =  '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

%% Initial Setup

if(~show_plots)
    set(groot,'DefaultFigureVisible','off') % Turn off figure popups for local
end

% Announce operating mode
fprintf('Starting %s in %s mode.\n',script_name,mode);

max_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 30, 31, 31]; % Number of days in each month

% Time counters [NOTE: NOT ROBUST TO CROSSING MONTH BOUNDARIES! Hard-code if needed]
time_count = (end_day - start_day + 1)*24 - (start_hour) - (24 - end_hour) + 1;

load(sprintf('%s/%s',intermediate_path,'latlon_gis_trim.mat'),'lon_gis_trim','lat_gis_trim'); % Load in preestablished latlon values
load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'lon_wrf_trim','lat_wrf_trim'); % Load in preestablished latlon values

dimensions_wrf = size(lon_wrf_trim);
ylen_wrf = dimensions_wrf(1);
xlen_wrf = dimensions_wrf(2);

dimensions_gis = size(lon_gis_trim);
ylen_gis = dimensions_gis(1);
xlen_gis = dimensions_gis(2);

mean_tag_string = sprintf('em%d_vm%d',plot_ens_mean,plot_value_mean);

if(~run_all_exps)
    run_name_list = [run_name];
    exp_name_list = [exp_name];
end

output_path_large = sprintf('%s/%s',output_path_base,'large');
output_path_small = sprintf('%s/%s',output_path_base,'small');
output_path_tile = sprintf('%s/%s',output_path_base,'tile');

% If output folders do not exist, create them
if(~isfolder(output_path_base))
    mkdir(output_path_large);
    mkdir(output_path_small);
    mkdir(output_path_tile);
end

force_continuity = false;

%% Main Loop

num_exps = length(exp_name_list);

mem_storage = zeros([num_exps num_vars time_count num_members]);
vm_storage = zeros([num_exps num_vars time_count]);
em_storage = zeros([num_exps num_vars time_count]);
sd_storage = zeros([num_exps num_vars time_count]);

for exp_idx = 1:num_exps
    
    exp_name = exp_name_list(exp_idx);
    run_name = run_name_list(exp_idx);
    
    time_axis = strings(time_count,1);
    
    for var_idx = 1:num_vars

        year = start_year;
        month = start_month;
        day = start_day;
        hour = start_hour;

        plume_series = zeros([time_count,num_rows]);

        for time_idx = 1:time_count

            timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
            timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
            time_axis(time_idx) = sprintf('%02d00',hour);

            filename = sprintf('ev_table_%s_%s.txt',run_name,timestamp_file);
            filepath_full = sprintf('%s/%s',input_path_base,filename);
            if(isfile(filepath_full))
                
                ev_table = tdfread(filepath_full);

                switch var_idx
                    case 1
                        var_name = 'RMSE';
                        plume_series(time_idx,:) = ev_table.RMSE;
                    case 2
                        var_name = 'Bias';
                        plume_series(time_idx,:) = ev_table.Bias;
                    case 3
                        var_name = 'ETS20';
                        plume_series(time_idx,:) = ev_table.ETS20;
                    case 4
                        var_name = 'ETS35';
                        plume_series(time_idx,:) = ev_table.ETS35;
                    case 5
                        var_name = 'CSI20';
                        plume_series(time_idx,:) = ev_table.CSI20;
                    case 6
                        var_name = 'CSI35';
                        plume_series(time_idx,:) = ev_table.CSI35;
                    case 7
                        var_name = 'SD';
                        plume_series(time_idx,:) = ev_table.SD;
                    case 8
                        var_name = 'SDE';
                        plume_series(time_idx,:) = ev_table.SDE;
                        %continue;
                    case 9
                        var_name = 'BS20';
                        plume_series(time_idx,:) = ev_table.BS20;
                    case 10
                        var_name = 'BS35';
                        plume_series(time_idx,:) = ev_table.BS35;
                end  
            else % If the file does not exist, use all NaN
                %fprintf('\n');
                var_name = var_name_list(var_idx);
                plume_series(time_idx,:) = NaN;
            end
            
            if(ismember(exp_name,["AIR-F","CONV-F"]) && timestamp_file == "202002071800")
                if(override_f18)
                    plume_series(time_idx,:) = NaN;
                    force_continuity = true;
                end
            end
            
            if(exp_name == "NODA" && timestamp_file == "202002071400")
                plume_series(time_idx,:) = NaN;
            end

            %% Time increment
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

        %% Plot the current variable time series

        ens_mem = plume_series(:,1:num_members);
        ens_mean = plume_series(:,ens_mean_idx);
        mean_ev = plume_series(:,mean_ev_idx);
        ens_std = plume_series(:,std_idx);
        % min, max, med also exist
        
        % Store for comp plots
        mem_storage(exp_idx,var_idx,:,:) = ens_mem;
        em_storage(exp_idx,var_idx,:) = ens_mean;
        vm_storage(exp_idx,var_idx,:) = mean_ev;
        sd_storage(exp_idx,var_idx,:) = ens_std;
        
        plot_18z_line = false;
        
        if(combine_series_plumes)
            if(ismember(exp_idx,[1 2]))
                continue;
            elseif(ismember(exp_idx,[3 4])) % Combine
                plot_18z_line = true;
                ens_mem(1:5,:) = mem_storage(exp_idx-2,var_idx,1:5,:);
                ens_mean(1:5) = em_storage(exp_idx-2,var_idx,1:5);
                mean_ev(1:5) = vm_storage(exp_idx-2,var_idx,1:5);
                ens_std(1:5) = sd_storage(exp_idx-2,var_idx,1:5);
            end
        end
        
        min_value = min([min(ens_mem,[],'all'),min(mean_ev,[],'all'),min(ens_mean,[],'all')]);
        max_value = max([max(ens_mem,[],'all'),max(mean_ev,[],'all'),max(ens_mean,[],'all')]);

    end
end

%% Plot experiment comparisons

if(extend_noda && start_hour == 14)
    em_storage(5,:,1) = mean(em_storage(1:2,:,1));
    vm_storage(5,:,1) = mean(vm_storage(1:2,:,1));
    sd_storage(5,:,1) = mean(sd_storage(1:2,:,1));
end

for var_idx = 1:num_vars
    
    var_name = var_name_list(var_idx);
    if(ismember(var_name,["BS20","BS35","SD","SDE"]))
        plot_error_bars = false;
    else
        plot_error_bars = true;
    end
    
    for plot_idx = 1:3
        
        
        is_sd = false;
        if(plot_idx == 1)
            plot_type = 'em';  % Specify fit, dif, trim, raw
            data_to_plot = em_storage;
        elseif(plot_idx == 2)
            plot_type = 'vm';
            data_to_plot = vm_storage;
        else
            plot_type = 'sd';
            continue;
        end
        
        if(force_continuity && start_hour == 14) % Force lines to connect
            data_to_plot(3:4,var_idx,5) = data_to_plot(1:2,var_idx,5);
        end
        
        if(exclude_noda)
            num_exps_to_plot = num_exps - 1;
            data_to_plot = data_to_plot(1:num_exps_to_plot,:,:);
        else
            num_exps_to_plot = num_exps;
        end
        
        if(exclude_noda)
            max_value = max(data_to_plot(:,var_idx,:) + sd_storage(1:num_exps_to_plot,var_idx,:),[],'all');
            min_value = min(data_to_plot(:,var_idx,:) - sd_storage(1:num_exps_to_plot,var_idx,:),[],'all');
        else
            max_value = max(data_to_plot(:,var_idx,:) + sd_storage(:,var_idx,:),[],'all');
            min_value = min(data_to_plot(:,var_idx,:) - sd_storage(:,var_idx,:),[],'all');
        end

        if(first_pass && plot_idx == 1)
            if(var_name == "RMSE")
                rmse_max = max_value;
                rmse_min = min_value;
            elseif(var_name == "Bias")
                bias_max = max_value;
                bias_min = min_value;
            elseif(var_name == "SD")
                sd_max = max_value;
                sd_min = min_value;
            elseif(var_name == "SDE")
                sde_max = max_value;
                sde_min = min_value;
            elseif(var_name == "CSI20" || var_name == "CSI35")
                csi_max = max_value;
                csi_min = min_value;
            elseif(var_name == "ETS20" || var_name == "ETS35")
                ets_max = max_value;
                ets_min = min_value;
            elseif(var_name == "BS20" || var_name == "BS35")
                bs_max = max_value;
                bs_min = min_value;
            end
        else
            if(var_name == "RMSE")
                if(max_value > rmse_max)
                    rmse_max = max_value;
                end
                if(min_value < rmse_min)
                    rmse_min = min_value;
                end
            elseif(var_name == "Bias")
                if(max_value > bias_max)
                    bias_max = max_value;
                end
                if(min_value < bias_min)
                    bias_min = min_value;
                end
            elseif(var_name == "SD")
                if(max_value > sd_max)
                    sd_max = max_value;
                end
                if(min_value < sd_min)
                    sd_min = min_value;
                end
            elseif(var_name == "SDE")
                if(max_value > sde_max)
                    sde_max = max_value;
                end
                if(min_value < sde_min)
                    sde_min = min_value;
                end
            elseif(var_name == "CSI20" || var_name == "CSI35")
                if(max_value > csi_max)
                    csi_max = max_value;
                end
                if(min_value < csi_min)
                    csi_min = min_value;
                end
            elseif(var_name == "ETS20" || var_name == "ETS35")
                if(max_value > ets_max)
                    ets_max = max_value;
                end
                if(min_value < ets_min)
                    ets_min = min_value;
                end
            elseif(var_name == "BS20" || var_name == "BS35")
                if(max_value > bs_max)
                    bs_max = max_value;
                end
                if(min_value < bs_min)
                    bs_min = min_value;
                end
            end  
        end
    end
end

first_pass = 0;

%save(sprintf('intermediate/ev_minmax_sbd%d.mat',subdomain),'rmse_min','rmse_max','bias_min','bias_max',"bs_min","bs_max","ets_min","ets_max","csi_min","csi_max","sd_min","sd_max","sde_max","sde_min");
  

end

save('intermediate/ev_minmax_allsbd.mat','rmse_min','rmse_max','bias_min','bias_max',"bs_min","bs_max","ets_min","ets_max","csi_min","csi_max","sd_min","sd_max","sde_max","sde_min");
    

fprintf('Done.\n');
