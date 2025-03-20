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

script_name = 'significance_test.m';

%% 0A. Script Controls

subdomain = 0;
batch_name = 't_testing';
%batch_name = 'alpha_sbd0';
%batch_name = 'beta_sbd0';
exclude_noda = 0;

combine_domain_limits = 0;
single_domain_limits = 0;

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

show_plots = 0;
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

if(exclude_noda)
    bias_max = 11.5;
    bias_min = -1.2;
    bs_max = 0.3;
    bs_min = 0;
    csi_max = 0.7;
    csi_min = 0;
    ets_max = 0.3;
    ets_min = -0.2;
    rmse_max_vm = 28.5;
    rmse_min_vm = 18;
    rmse_max_em = 26;
    rmse_min_em = 16;
    %sd_max = 0.2;
    sd_min = 0;
else
    bias_max = 11.5;
    bias_min = -1.2;
    bs_max = 0.3;
    bs_min = 0;
    csi_max = 0.7;
    csi_min = 0;
    ets_max = 0.3;
    ets_min = -0.2;
    rmse_max_vm = 28.5;
    rmse_min_vm = 18;
    rmse_max_em = 26;
    rmse_min_em = 16;
    %sd_max = 0.2;
    sd_min = 0.1;
end

%% 0B. General Settings

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/output/ev_tables';
    intermediate_path = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = sprintf('C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/output/sig_tables/%s',batch_name);
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
plot_units = "dB";

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

output_path_large = sprintf('%s',output_path_base);


% If output folders do not exist, create them
if(~isfolder(output_path_base))
    mkdir(output_path_large);
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
                fprintf('File [%s] not found.\n',filepath_full);
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
        
    end
end

%% Plot experiment comparisons

ktest_results = zeros(2,num_vars,time_count,2).*NaN; % Slot 1 is AIR-CONV (1) or AIRF-CONVF (2); Slot 4 is Y/N (1) or P-value (2)
pos_dif_results = zeros(num_vars,time_count);

%mem_storage = zeros([num_exps num_vars time_count num_members]);
for var_idx = 1:num_vars
    
    var_name = var_name_list(var_idx);
    
    data_to_plot = mem_storage;
    if(force_continuity && start_hour == 14) % Force lines to connect
        data_to_plot(3:4,var_idx,5,:) = data_to_plot(1:2,var_idx,5,:);
    end

    for time_idx = 1:time_count
        if(time_idx < 6)
            data_a = squeeze(data_to_plot(1,var_idx,time_idx,:));
            data_b = squeeze(data_to_plot(2,var_idx,time_idx,:));
    
            % Conduct significance t-test (Kolmogorov-Smirnov 2-sample test)
            %p_value = 0.05; % Defaults to this already
            [h,p] = kstest2(data_a,data_b);
            ktest_results(1,var_idx,time_idx,:) = [h,p];
            
        elseif(time_idx >= 6)
            data_a = squeeze(data_to_plot(3,var_idx,time_idx,:));
            data_b = squeeze(data_to_plot(4,var_idx,time_idx,:));
            [h,p] = kstest2(data_a,data_b);
            ktest_results(2,var_idx,time_idx,:) = [h,p];
            
        end
        % Record which ensemble was closer to observations
        if(ismember(var_idx,[3,4,5,6])) % ETS,CSI higher better
            pos_dif_results(var_idx,time_idx) = mean(data_a,"omitnan") > mean(data_b,"omitnan");
        elseif(ismember(var_idx,[1,8,9,10])) % RMSE,SDE,BS lower better
            pos_dif_results(var_idx,time_idx) = mean(data_a,"omitnan") < mean(data_b,"omitnan");
        elseif(var_idx == 2) % Bias closer to 0 better
            pos_dif_results(var_idx,time_idx) = abs(mean(data_a,"omitnan")) < abs(mean(data_b,"omitnan"));
        end
    end

end

comparison_strings = ["A-C","AF-CF"];
var_strings = var_name_list;
time_strings = [14:23,0];


T = combinations(comparison_strings,var_strings,time_strings);
num_rows = size(T,1);
T.P_VALUE = zeros(num_rows,1).*NaN;
T.SIG = zeros(num_rows,1).*NaN;
T.POS_DIF = zeros(num_rows,1).*NaN;
T.SIG_SIGN = zeros(num_rows,1);

row_idx = 1;
for comp_idx = 1:2
    for var_idx = 1:num_vars
        var_name = var_name_list(var_idx);
        for time_idx = 1:time_count
            T{row_idx,"P_VALUE"} = ktest_results(comp_idx,var_idx,time_idx,2);
            T{row_idx,"SIG"} = max(ktest_results(comp_idx,var_idx,time_idx,1),0);
            T{row_idx,"POS_DIF"} = pos_dif_results(var_idx,time_idx);
            if(T{row_idx,"SIG"} & T{row_idx,"POS_DIF"})
                T{row_idx,"SIG_SIGN"} = 1;
            elseif(T{row_idx,"SIG"} & ~T{row_idx,"POS_DIF"})
                T{row_idx,"SIG_SIGN"} = -1;
            end
            row_idx = row_idx+1;
        end
    end
end

T.Properties.VariableNames = ["EXPS","VAR","TIME","P_VALUE","SIG","POS_DIF","SIG_SIGN"];

to_delete = (T.EXPS == "A-C" & ((T.TIME > 18) | T.TIME == 0)) | (T.EXPS == "AF-CF" & (T.TIME <= 18 & T.TIME ~= 0)); 
T(to_delete,:) = [];

% Save table of test results
%sig_table = array2table(ktest_results,'VariableNames',{'A','B','C','D'});
writetable(T,sprintf('%s/sig_test_table_%s.txt',output_path_base,batch_name),'Delimiter','\t','WriteRowNames',true);
    

fprintf('Done.\n');
