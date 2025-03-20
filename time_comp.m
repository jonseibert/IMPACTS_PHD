filename = 'wrf_enkf_output_d02_001';
input_path_wrf = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020/202002071800/';
filename_2 = 'wrfout_d02_18z_airf';

air_wind = ncread(sprintf('%s/%s',input_path_wrf,filename),'U');
airf_wind = ncread(sprintf('%s/%s',input_path_wrf,filename),'U');
air_temp = ncread(sprintf('%s/%s',input_path_wrf,filename),'T');
airf_temp = ncread(sprintf('%s/%s',input_path_wrf,filename),'T');

wind_dif = air_wind-airf_wind;
temp_dif = air_temp-airf_temp;

any((wind_dif ~= 0),'all')
any((temp_dif ~= 0),'all')