input_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/testing';
output_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/testing';
f1 = 'wrf_enkf_output_d02_001_18z';
f2 = 'wrfout_d02_18z_airf_test2';
%D1 = ncread(sprintf('%s/%s',input_path,f1),'REFL_10CM');
%D2 = ncread(sprintf('%s/%s',input_path,f2),'REFL_10CM');
%D1 = ncread(sprintf('%s/%s',input_path,f1),'P');
%D2 = ncread(sprintf('%s/%s',input_path,f2),'P');
%D1_Z = max(D1,[],3);
%D2_Z = max(D2,[],3);
%D1_Z = D1(:,:,35);
%D2_Z = D2(:,:,35);

g = 9.8;





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

%%
input_path_base = input_path;
%D1
filename = f1;

D1_refl = ncread(sprintf('%s/%s',input_path_base,filename),data_name_refl); % Reflectivity
D1_p_base = ncread(sprintf('%s/%s',input_path_base,filename),data_name_pres); % Base pressure
D1_p_pert = ncread(sprintf('%s/%s',input_path_base,filename),data_name_pres_2); % Perturbation pressure
D1_pt_pert = ncread(sprintf('%s/%s',input_path_base,filename),data_name_temp); % Perturbation potential temperature
D1_psfc = ncread(sprintf('%s/%s',input_path_base,filename),data_name_psfc); % PSFC = Surface Level Pressure?? MSLP??
D1_t2m = ncread(sprintf('%s/%s',input_path_base,filename),data_name_t2m); % 2m Temperature
D1_q = ncread(sprintf('%s/%s',input_path_base,filename),data_name_q); % 2m Temperature
D1_gp_base = ncread(sprintf('%s/%s',input_path_base,filename),data_name_gp); % 2m Temperature
D1_gp_pert = ncread(sprintf('%s/%s',input_path_base,filename),data_name_gp_2); % 2m Temperature
D1_U = ncread(sprintf('%s/%s',input_path_base,filename),'U'); % U
D1_V = ncread(sprintf('%s/%s',input_path_base,filename),'V'); % U
D1_W = ncread(sprintf('%s/%s',input_path_base,filename),'W'); % U
D1_THM = ncread(sprintf('%s/%s',input_path_base,filename),'THM');
D1_QV = ncread(sprintf('%s/%s',input_path_base,filename),'QVAPOR');

D1_p = D1_p_base + D1_p_pert; % 3D Pressure field
D1_gp = D1_gp_base + D1_gp_pert; % 3D Geopotential field
D1_gph = D1_gp./g;

D1_pt = D1_pt_pert + 300; % Potential temperature
D1_t = D1_pt.*((D1_p./100000).^(0.287)); % 3D Temperature field

%D2
filename = f2;

D2_refl = ncread(sprintf('%s/%s',input_path_base,filename),data_name_refl); % Reflectivity
D2_p_base = ncread(sprintf('%s/%s',input_path_base,filename),data_name_pres); % Base pressure
D2_p_pert = ncread(sprintf('%s/%s',input_path_base,filename),data_name_pres_2); % Perturbation pressure
D2_pt_pert = ncread(sprintf('%s/%s',input_path_base,filename),data_name_temp); % Perturbation potential temperature
D2_psfc = ncread(sprintf('%s/%s',input_path_base,filename),data_name_psfc); % PSFC = Surface Level Pressure?? MSLP??
D2_t2m = ncread(sprintf('%s/%s',input_path_base,filename),data_name_t2m); % 2m Temperature
D2_q = ncread(sprintf('%s/%s',input_path_base,filename),data_name_q); % 2m Temperature
D2_gp_base = ncread(sprintf('%s/%s',input_path_base,filename),data_name_gp); % 2m Temperature
D2_gp_pert = ncread(sprintf('%s/%s',input_path_base,filename),data_name_gp_2); % 2m Temperature
D2_U = ncread(sprintf('%s/%s',input_path_base,filename),'U'); % U
D2_V = ncread(sprintf('%s/%s',input_path_base,filename),'V'); % U
D2_W = ncread(sprintf('%s/%s',input_path_base,filename),'W'); % U
D2_THM = ncread(sprintf('%s/%s',input_path_base,filename),'THM');
D2_QV = ncread(sprintf('%s/%s',input_path_base,filename),'QVAPOR');

D2_p = D2_p_base + D2_p_pert; % 3D Pressure field
D2_gp = D2_gp_base + D2_gp_pert; % 3D Geopotential field
D2_gph = D2_gp./g;

D2_pt = D2_pt_pert + 300; % Potential temperature
D2_t = D2_pt.*((D2_p./100000).^(0.287)); % 3D Temperature field

%%

D1_Z = D1_refl;
D2_Z = D2_refl;

plot_filename = 'test_ZB.png';
units = 'dBZ';
var_name = 'Refl';

composite = 1;
level = 35;

plot_base_ref = 0;
%cmap = 'redblue';
cmap = 'reflmap';
use_clim = 1;
clim_lower = -30;
clim_upper = 75;
%%

load('intermediate/latlon_wrf_trim.mat')
%load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'lon_wrf_trim','lat_wrf_trim'); % Load in preestablished latlon values
lon_wrf = lon_wrf_trim;
lat_wrf = lat_wrf_trim;
clearvars lon_wrf_trim lat_wrf_trim;

plot_type = 'dif';


if(composite)
    D1_Z = max(D1_Z,[],3);
    comp_flag = 'C';
else
    D1_Z = D1_Z(:,:,level);
    comp_flag = 'NC';
end
D1_Z = D1_Z';
D1_Z = flip(D1_Z);
D1_trim = D1_Z(n_idx:s_idx,w_idx:e_idx);


if(composite)
    D2_Z = max(D2_Z,[],3);
    comp_flag = 'C';
else
    D2_Z = D2_Z(:,:,level);
    comp_flag = 'NC';
end
D2_Z = flip(D2_Z);
D2_trim = D2_Z(n_idx:s_idx,w_idx:e_idx);

data_dif = D1_trim - D2_trim;

%%

%data_to_plot = data_dif;
data_to_plot = D2_trim;

if(plot_base_ref)

    exp_name_a = 'AIR';
    exp_name_b = 'AIR-F';

    input_path_br = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020/BASE_REF/';
    
    load(sprintf('%s/%s/%s',input_path_br,exp_name_a,'air_base_d02_001_refl_202002071800.mat'),'wrf_base_refl');
    D1_Z = wrf_base_refl;
    load(sprintf('%s/%s/%s',input_path_br,exp_name_b,'air-f_base_d02_001_refl_202002071800.mat'),'wrf_base_refl');
    D2_Z = wrf_base_refl;
    data_to_plot = D1_Z-D2_Z;
    comp_flag = 'B';

end



%%

% Plot LARGE
f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
h = pcolor(lon_wrf,lat_wrf,data_to_plot); % Plot the data
set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
shading interp; % Smooth out plot from grid-boxes
colorbar('FontSize',axes_font_size); % Make colorbar
colormap(cmap); % Set colors
if(use_clim)
    caxis([clim_lower clim_upper]);
end
% Plot state borders
hold on;
borders('continental us','black','linewidth',1); 
hold off;

if(limit_borders)
    xlim([w_lim e_lim]);
    ylim([s_lim n_lim]);
end

% Apply labels
title(sprintf('A-AF 18z Dif Test: %s(%s)[%s,Z%d]',var_name,units,comp_flag,level),'FontSize',title_font_size); 
xlabel('Longitude','FontSize',label_font_size);
ylabel('Latitude','FontSize',label_font_size);
set(gca,'Fontsize',axes_font_size);
saveas(h,sprintf('%s/%s',output_path,plot_filename));

