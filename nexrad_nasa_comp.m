%% NASA Data comparison with NEXRAD obs
% Jon Seibert
% Last modified: 7 Sep 2021

%% TODO

% Forgoing regrid comparison for now- sticking to image slice
% 1. Take min/max longitude from NASA, min/max latitude for each longitude
% point (small range to approximate point if inconsistent?)
% 2. Take that lat/lon slice from NEXRAD data
% 3. Take composite of NASA data refl, plot both on same lat/lon minmax-
% picture should turn out the same despite dif grid

% For wind- figure out how to get not composite from Seth code, pick a
% level, do the same thing

% Regrid process would be quite intensive given how the flight data isn't on a consistent lat/lon grid- ask if necessary to have a
% dif/RMSE

% Proportional decrease instead of "regridding" the list? slope isn't
% consistent, can't use that

% Temporary solution: plot full NEXRAD domain, mask down to EXRAD by same
% process as GEM mask?

% TODO: create equation based on the top and bottom bounding lines to
% determine what lats are in the box for each NEXRAD lon

% ALTERNATE IDEA: Use linear interpolation (for each NEXRAD lon, grab 2
% nearest neighbor NASA lons on N and S edges of bounded box, determine if
% the point on the line between them is above or below the NEXRAD lon's...
% lat? finalize later



% Brute force workaround: use bilinear interpolation to extract NEXRAD
% value at each specific (non-NaN) point in EXRAD data!
    % Can eventually use this for WRF results too
    % Can also get true grid-point quantitative diff this way
    
% Brute-force and trimming of results successful
% NEXT: Composite the NASA EXRAD data, do quant dif, make pics of composite
% NASA and of dif

%% 0A. Timegroup math
% NASA Timegroups: each pair describes the relevant timestamps contained
% within (1&2, 3&4, etc)
timegroups = [135456,142903,143307,145251,145658,150453,150806,152406,152952,155034,155544,161317,161829,170931];
dategroups = [20200207135456,20200207142903,20200207143307,20200207145251,20200207145658,20200207150453,20200207150806,20200207152406,20200207152952,20200207155034,20200207155544,20200207161317,20200207161829,20200207170931];

date_strings = string(dategroups);
dt = datetime(date_strings,'InputFormat','yyyyMMddHHmmss'); % DateTimes

for n = 1:7
    %n
    window = (dt(n*2) - dt(n*2 - 1))/2;
    midpoint = dt(n*2) - window;
end

midpoints = [141159, 144259, 150055, 151606, 154013, 160430, 164400]; %(All with 20200207 in front)
% Windows = [001703, 000952, 000357, 000800, 001021, 000846, 002531]

% NEXRAD Station List
%stations = ["KBGM","KBOX","KBUF","KCCX","KCLE","KCXX","KDIX","KDTX","KENX","KGYX","KOKX","KTYX"];

% Data types: 
    % NCR = 230km Composite Reflectivity
    % NVW = VAD Wind Profile
    
% WRF variable names: lon = XLON, lat = XLAT, time = XTIME, refl =
% REFL_10CM

%% Colormap

%radar = [255 255 255;  % -Inf < x < -30 (White)
%					230 230 230; 	% -30 <= x < 0 dbz (Light Grey)
%					0 255 255;		% 0 <= x < 5 dbz (Cyan)
%					75 200 255; 	% 5 <= x < 10 dbz (Light Blue)
%					30 144 255; 	% 10 <= x < 15 dbz (Blue)
%					0 0 255;		% 15 <= x < 20 dbz (Indigo)
%			    	0 255 0; 		% 20 <= x < 25 dbz (Lime)
%			    	0 205 0; 		% 25 <= x < 30 dbz (Green)
%			    	0 139 0; 		% 30 <= x < 35 dbz (Dark Green)
%			    	255 255 0; 		% 35 <= x < 40 dbz (Yellow)
%			    	255 215 0;		% 40 <= x < 45 dbz (Gold)
%			    	255 127 0; 		% 45 <= x < 50 dbz (Orange)
%			    	255 0 0;		% 50 <= x < 55 dbz (Scarlet)
%			    	205 0 0; 		% 55 <= x < 60 dbz (Red)
%			    	139 0 0; 		% 60 <= x < 65 dbz (Dark Red)
%			    	255 0 255;		% 65 <= x < 70 dbz (Fuschia)
%			    	145 44 238]; 	% 70 <= x < 75 dbz (Purple)
%			    	255 255 255];	% 75 < x (White)
                    
% Cutting out one of the greys makes the colors line up to the numbers properly on caxis([-30 75])
radar = [255 255 255; % -35 (and lower) < x < -30 (White)
        %230 230 230;  % -30:-25 (Light Grey)  
        230 230 230;  % -25:-20 (Light Grey)
        230 230 230;  % -20:-15 (Light Grey)
        230 230 230;  % -15:-10 (Light Grey)
        230 230 230;  % -10:-5 (Light Grey)
        230 230 230;  % -5 <= x < 0 dbz (Light Grey)
        0 255 255;	  % 0 <= x < 5 dbz (Cyan)
        75 200 255;   % 5 <= x < 10 dbz (Light Blue)
        30 144 255;   % 10 <= x < 15 dbz (Blue)
        0 0 255;	  % 15 <= x < 20 dbz (Indigo)
        0 255 0; 	  % 20 <= x < 25 dbz (Lime)
        0 205 0; 	  % 25 <= x < 30 dbz (Green)
        0 139 0; 	  % 30 <= x < 35 dbz (Dark Green)
        255 255 0; 	  % 35 <= x < 40 dbz (Yellow)
        255 215 0;	  % 40 <= x < 45 dbz (Gold)
        255 127 0; 	  % 45 <= x < 50 dbz (Orange)
        255 0 0;	  % 50 <= x < 55 dbz (Scarlet)
        205 0 0; 	  % 55 <= x < 60 dbz (Red)
        139 0 0; 	  % 60 <= x < 65 dbz (Dark Red)
        255 0 255;	  % 65 <= x < 70 dbz (Fuschia)
        145 44 238];  % 70 <= x < 75 dbz (Purple)
                
radar = radar/255;

%% 0B. Options

% Plot font sizes
title_font_size = 18;
label_font_size = 14;
axes_font_size = 12;

z = 3; % 1.5 km
num_data = 4;
date = '2020-02-07';

base_path = pwd;

nasa_files = string(ls("input/exrad_data/exrad_impacts*"));
nexrad_files = string(ls("nexrad_data/IMPACTS*"));
intermediate_path = "intermediate";

% Nexrad stuff
w_limit = -82; % Western lon plot limit, to match with WRF domain
plot_width = 1600;
plot_height = 900;

% Radar colormap
use_custom_colormap = true;

%% 1. Extract basic dimensions of NASA and NEXRAD data

% NASA EXRAD data
% Read in data
filename = nasa_files(n);
nas_date = extractBetween(filename,15,36);

filepath = sprintf("input/exrad_data/%s",filename);

% Dims: lat x lon x 57 (across x along x height)
lon_nas = ncread(filepath, "longitude"); % 3D [degrees]

% Define dimensions
container_dims = size(lon_nas);

%% 1. Loop: For each timestamp, plot both NASA and NEXRAD data

% Define output format and location
output_path = "output/";
output_name_format = 'refl_%d_%s.png';
output_format = strcat(output_path,output_name_format);
output_format_hist = strcat(output_path,'refl_comp_hist_%s_%d.png');

% Create container for all interpolated NEXRAD data
% refl_nex_fitted_all =
% zeros(container_dims(1),container_dims(2),num_data); %CAN'T GET THIS TO
% WORK RN: HARDCODING TEMPORARILY
refl_nex_fitted_1 = [];
refl_nex_fitted_2 = [];
refl_nex_fitted_3 = [];
refl_nex_fitted_4 = [];

for n = 1:num_data
    
    % NASA EXRAD data
    % Read in data
    filename = nasa_files(n);
    nas_date = extractBetween(filename,15,36);
    
    filepath = sprintf("input/exrad_data/%s",filename);

    % Dims: lat x lon x 57 (across x along x height)
    across_track_grid = ncread(filepath, "across_track_grid"); % 41x1 [km]
    along_track_grid = ncread(filepath, "along_track_grid"); % 517x1 [km]
    vertical_grid = ncread(filepath, "vertical_grid"); % 57x1 [km]
    zonal_wind = ncread(filepath, "zonal_wind"); % 3D [m/s]
    meridional_wind = ncread(filepath, "meridional_wind"); % 3D [m/s]
    vertical_wind = ncread(filepath, "vertical_wind"); % 3D [m/s]
    Ku_band_reflectivity = ncread(filepath, "Ku_band_reflectivity"); % 3D [dBZ]
    lat_nas = ncread(filepath, "latitude"); % 3D [degrees]
    lon_nas = ncread(filepath, "longitude"); % 3D [degrees]


    % Process and plot
    lon_nas = lon_nas - 360; % Convert to negative values
    dims = size(lon_nas);
    lat = dims(1); % This is really just the N-S grid size
    lon = dims(2); % E-W grid size
    height = dims(3);

    % Plot 1: Reflectivity on lat/lon map [then add wind vectors] with colorbar
    % at 1.5 km height (aka z coordinate = 03 on the vertical grid)

    lat_slice = lat_nas(:,:,z);
    lon_slice = lon_nas(:,:,z);
    refl_slice = Ku_band_reflectivity(:,:,z);
    
    % Calculate the bounds of the real data for proper comparison with
    % NEXRAD / count how many NaN rows and columns there are
    [first_data_lat,last_data_lat,first_data_lon,last_data_lon] = find_data_edges_nas(lat,lon,refl_slice);
    
    %zonal_slice = zonal_wind(:,:,z);
    %merid_slice = meridional_wind(:,:,z);
    %wind_x = zeros(size(zonal_slice));
    %wind_x(9:4:(lat-3),5:6:lon) = zonal_slice(9:4:(lat-3),5:6:lon);
    %wind_y = zeros(size(merid_slice));
    %wind_y(9:4:(lat-3),5:6:lon) = merid_slice(9:4:(lat-3),5:6:lon);

    figure('Position',[0 200 2100 620]);
    h = pcolor(lon_slice,lat_slice,refl_slice);
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    if(use_custom_colormap)
        colormap(radar);
        caxis([-30 75]);
    else
        colormap(jet);
        caxis([15 45]);
    end
    %q = quiver(lon_slice,lat_slice,wind_x,wind_y,3,'black','linewidth',1.1,'MaxHeadSize',1); 
    %q.MaxHeadSize = 10;
    ax = gca();
    ax.FontSize = axes_font_size;
    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(strcat("IMPACTS EXRAD Reflectivity (dBZ) | z = 1.5 km | ",nas_date),'FontSize',title_font_size);

    saveas(h,sprintf(output_format,n,'nas'));
    %saveas(h,strcat('output/','refl_',string(n),'_nas.png'));
    
    % --------------------------------------------------------------------
    % Plot vertical slice (lon vs z) of IMPACTS EXRAD data along central
    % spine of latitude/flight path
    mid_lat_idx = floor(lat/2);
    vertical_slice = squeeze(Ku_band_reflectivity(mid_lat_idx,:,:));
    lon_slice_vert = squeeze(lon_nas(mid_lat_idx,:,:));
    
    % Make height grid into 2D grid by copying it along the lon axis
    vertical_grid_wide = zeros(lon,height);
    for lon_idx = 1:lon
        vertical_grid_wide(lon_idx,:) = vertical_grid;
    end
    
    output_name_format_vert = 'refl_zslice_%d_%s.png';
    output_format_vert = strcat(output_path,output_name_format_vert);
    
    figure('Position',[0 200 2100 620]);
    h = pcolor(lon_slice_vert,vertical_grid_wide,vertical_slice);
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    if(use_custom_colormap)
        colormap(radar);
        caxis([-30 75]);
    else
        colormap(jet);
        caxis([15 45]);
    end
    ax = gca();
    ax.FontSize = axes_font_size;
    ylim([1.25 8]);
    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Height (km)','FontSize',label_font_size);
    title(strcat("IMPACTS EXRAD Reflectivity (dBZ) | Vertical Slice along flight path | ",nas_date),'FontSize',title_font_size);

    saveas(h,sprintf(output_format_vert,n,'nas'));
    
    % --------------------------------------------------------------------
    % Now plot NEXRAD data
    % Read in data
    % Use previously defined lat_nas,lon_nas

    % Step 2: Extract NEXRAD data
    input_path_nexrad = "nexrad_data";
    filename = nexrad_files(n);
    filepath = sprintf("%s/%s/%s",base_path,input_path_nexrad,filename);

    lat_nex = ncread(filepath, "latitude"); % 1D [degrees]
    lon_nex = ncread(filepath, "longitude"); % 1D [degrees]
    refl_nex = ncread(filepath, "comp_reflect"); % 2D [dBZ]
    
    
    % -------------------------------------------------------------------
    % Interjection: plot full NEXRAD radar snapshot
    output_name_format_fullnex = 'refl_%d_full_%s.png';
    output_format_fullnex = strcat(output_path,output_name_format_fullnex);
    
    [LON,LAT] = meshgrid(lon_nex,flip(lat_nex));
    
    % Find bounding box of EXRAD data
    s_border = min(min(min((lat_nas))));
    n_border = max(max(max((lat_nas))));
    w_border = min(min(min((lon_nas))));
    e_border = max(max(max((lon_nas))));
    
    bounding_sequence_x = [w_border e_border e_border w_border w_border];
    bounding_sequence_y = [n_border n_border s_border s_border n_border];
 
%%
    figure('Position',[50 50 plot_width plot_height]);
    h = pcolor(LON,LAT,flip(refl_nex'));
    %hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    if(use_custom_colormap)
        colormap(radar);
        caxis([-30 75]);
    else
        colormap(jet);
        caxis([0 45]);
    end
    %q = quiver(lon_slice,lat_slice,wind_x,wind_y,3,'black','linewidth',1.1,'MaxHeadSize',1); 
    %q.MaxHeadSize = 10;
    ax = gca();
    ax.FontSize = axes_font_size;
    hold on;
    borders('continental us','black','linewidth',1);
    hold on;
    plot(bounding_sequence_x,bounding_sequence_y,'r','linewidth',2.2);
    xlim([w_limit max(max(LON))]);
    
    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf("NEXRAD Composite Reflectivity: Full (dBZ) | %s-%s",date,string(midpoints(n))),'FontSize',title_font_size);
%%
    saveas(h,sprintf(output_format_fullnex,n,'nex'));

    % 2. Brute force bilinear interp to make the NEXRAD data match the shape of the NASA data
    % NOTE: Assumes that the NEXRAD data fully encompasses the NASA data, so that
    % edge cases can work as normal!

    % Define empty grid to fit NEX data into
    dims = size(lon_nas);
    lon_slice = lon_nas(:,:,1);
    lat_slice = lat_nas(:,:,1);
    refl_nex_fitted = zeros(dims(1),dims(2)); % Composite anyway, no need for vertical

    % Loop through each NASA point and find corresponding NEXRAD value
    for x = 1:dims(1)
        for y = 1:dims(2)
            [lon1,lon2,lat1,lat2] = find_nn_idx_nex(lon_slice(x,y),lat_slice(x,y),lon_nex,lat_nex); % Get nn index values
            if(lon1 == 0 || lon2 == 0 || lat1 == 0 || lat2 == 0) % If indices could not be found / reference points are NaN
                refl_nex_fitted(x,y) = -30;
            else
                refl_nex_fitted(x,y) = bilinear_interpolate(lon_nex(lon1),lon_nex(lon2),lon_slice(x,y),lat_nex(lat1),lat_nex(lat2),lat_slice(x,y),refl_nex(lon1,lat1),refl_nex(lon2,lat1),refl_nex(lon1,lat2),refl_nex(lon2,lat2));
            end
        end
    end
    
    % Fill in NaNs on the border where they are in EXRAD
    if(first_data_lat > 1)
        refl_nex_fitted(1:first_data_lat,:) = NaN;
    end
    if(last_data_lat < dims(1))
        refl_nex_fitted(last_data_lat:dims(1),:) = NaN;
    end
    if(first_data_lon > 1)
        refl_nex_fitted(:,1:first_data_lon) = NaN;
    end
    if(last_data_lon < dims(2))
        refl_nex_fitted(:,last_data_lon:dims(2)) = NaN;
    end
    
    % We're plotting the fitted NEX data on the NAS grid
    %[lon_slice,lat_slice] = meshgrid(lon_nex,flip(lat_nex));
    lat_slice = lat_nas(:,:,z);
    lon_slice = lon_nas(:,:,z);
    refl_slice = refl_nex_fitted;
    %refl_nex_fitted_all(:,:,n) = refl_nex_fitted; % Store this result for later
    switch n % TEMP method of storing variable size outputs for later
        case 1
            refl_nex_fitted_1 = refl_nex_fitted;
        case 2
            refl_nex_fitted_2 = refl_nex_fitted;
        case 3
            refl_nex_fitted_3 = refl_nex_fitted;
        case 4
            refl_nex_fitted_4 = refl_nex_fitted;
    end
    
    % Save moulded NEXRAD data for ease of later reference
    save(sprintf('Intermediate/refl_nex_fitted_%d',n),'refl_nex_fitted');

    figure('Position',[0 200 2100 620]);
    h = pcolor(lon_slice,lat_slice,refl_slice);
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    if(use_custom_colormap)
        colormap(radar);
        caxis([-30 75]);
    else
        colormap(jet);
        caxis([15 45]);
    end
    %q = quiver(lon_slice,lat_slice,wind_x,wind_y,3,'black','linewidth',1.1,'MaxHeadSize',1); 
    %q.MaxHeadSize = 10;
    
    ax = gca();
    ax.FontSize = axes_font_size;
    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf("NEXRAD Composite Reflectivity on EXRAD Grid (dBZ) | %s-%s",date,string(midpoints(n))),'FontSize',title_font_size);

    saveas(h,sprintf(output_format,n,'nex'));
    
    % --------------------------------------------------------------------
    % Plot histogram/PDF of composite reflectivity values [NEXRAD]
   
    bins = -30:1:75;
    
    figure('Position',[0 200 2100 620]);
    h = histogram(refl_slice,'BinEdges',bins,'Normalization','probability');
    hold on;

    ax = gca();
    ax.FontSize = axes_font_size;
    hold off;

    xlabel('Reflectivity (dBZ)','FontSize',label_font_size);
    ylabel('Probability','FontSize',label_font_size);
   
    title(sprintf("Histogram of NEXRAD Composite Reflectivity on EXRAD Grid (dBZ) | %s-%s",date,string(midpoints(n))),'FontSize',title_font_size);
    saveas(h,sprintf(output_format_hist,'nex',n));
    
    % -------------------------------------------------------
    % 3. Calculate and plot dif between composite NEXRAD and z = 1.5km EXRAD, for lulz
    % Perform pixelwise dif and RMSE
    output_path = "output/";
    output_name_format_dif = 'refl_dif_uneven_%d.png';
    output_format_dif = strcat(output_path,output_name_format_dif); 
    
    refl_nas = Ku_band_reflectivity(:,:,z);
    refl_dif = refl_nas - refl_nex_fitted;
    
    SE = 0;
    for x = 1:dims(1)
        for y = 1:dims(2)
            SE_addition = ((refl_nas(x,y) - refl_nex_fitted(x,y))^2);
            if(~isnan(SE_addition))
                SE = SE + SE_addition;
            end
        end
    end
    MSE = SE/(dims(1)*dims(2));
    RMSE = sqrt(MSE);

    sprintf('RMSE for EXRAD(nc)-NEXRAD(c) dif case %d = %d',n,RMSE)
    
    % We're plotting the dif on the NAS grid
    figure('Position',[0 200 2100 620]);
    h = pcolor(lon_slice,lat_slice,refl_dif);
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    colormap(redblue);
    %q = quiver(lon_slice,lat_slice,wind_x,wind_y,3,'black','linewidth',1.1,'MaxHeadSize',1); 
    %q.MaxHeadSize = 10;
    caxis([-30 30]);
    ax = gca();
    ax.FontSize = axes_font_size;
    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf("EXRAD (z = 1.5km) - NEXRAD Composite Reflectivity Dif (dBZ) | %s-%s",date,string(midpoints(n))),'FontSize',title_font_size);

    saveas(h,sprintf(output_format_dif,n));
    
    close all;
    
end

%% 2. Loop: For each timestamp, composite NASA EXRAD data, take dif with NEXRAD, and plot composite and dif

output_path = "output/";
output_name_format_comp = 'refl_nas_%d_comp.png';
output_format_comp = strcat(output_path,output_name_format_comp);
output_name_format_dif = 'refl_dif_%d.png';
output_format_dif = strcat(output_path,output_name_format_dif);

for n = 1:num_data
    

    
    % NASA EXRAD data
    % Read in data
    filename = nasa_files(n);
    nas_date = extractBetween(filename,15,36);
    
    filepath = sprintf("%s/%s","nasa_data",filename);

    % Dims: lat x lon x 57 (across x along x height)
    across_track_grid = ncread(filepath, "across_track_grid"); % 41x1 [km]
    along_track_grid = ncread(filepath, "along_track_grid"); % 517x1 [km]
    vertical_grid = ncread(filepath, "vertical_grid"); % 57x1 [km]
    Ku_band_reflectivity = ncread(filepath, "Ku_band_reflectivity"); % 3D [dBZ]
    lat_nas = ncread(filepath, "latitude"); % 3D [degrees]
    lon_nas = ncread(filepath, "longitude"); % 3D [degrees]

    % Dimensions
    lon_nas = lon_nas - 360; % Convert to negative values
    dims = size(lon_nas);
    lat = dims(1); % This is really just the N-S grid size
    lon = dims(2); % E-W grid size
    height = dims(3);
    
    % Composite NAS reflectivity
    lat_slice = lat_nas(:,:,1);
    lon_slice = lon_nas(:,:,1);
    refl_comp = max(Ku_band_reflectivity(:,:,2:57),[],3); % Take only heights above 1.25 km!!!!
    
    % Save for exterior use
    refl_ex = Ku_band_reflectivity(:,:,z);
    refl_ex_comp = refl_comp;
    
    save(sprintf('%s/refl_ex_3km_%d',intermediate_path,n),'refl_ex');
    save(sprintf('%s/refl_ex_comp_%d',intermediate_path,n),'refl_ex_comp');
    
    % Create figure
    figure('Position',[0 200 2100 620]);
    h = pcolor(lon_slice,lat_slice,refl_comp);
    hold on;
    %if(n > 0 && n < 5)
    set(h, 'EdgeColor', 'none');
    shading interp;
    %end
    c = colorbar('FontSize',axes_font_size);
    if(use_custom_colormap)
        colormap(radar);
        caxis([-30 75]);
    else
        colormap(jet);
        caxis([15 45]);
    end
    %q = quiver(lon_slice,lat_slice,wind_x,wind_y,3,'black','linewidth',1.1,'MaxHeadSize',1); 
    %q.MaxHeadSize = 10;
    ax = gca();
    ax.FontSize = axes_font_size;
    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(strcat("IMPACTS EXRAD Composite Reflectivity (dBZ) | ",nas_date),'FontSize',title_font_size);

    saveas(h,sprintf(output_format_comp,n));
    
    % --------------------------------------------------------------------
    % Plot histogram/PDF of composite reflectivity values [EXRAD]
    
    bins = -30:1:75;
    
    figure('Position',[0 200 2100 620]);
    h = histogram(refl_comp,'BinEdges',bins,'Normalization','probability');
    hold on;

    ax = gca();
    ax.FontSize = axes_font_size;
    hold off;

    xlabel('Reflectivity (dBZ)','FontSize',label_font_size);
    ylabel('Probability','FontSize',label_font_size);
    title(strcat("Histogram of EXRAD Composite Reflectivity (dBZ) | ",nas_date),'FontSize',title_font_size);
    saveas(h,sprintf(output_format_hist,'ex',n));
    
    % --------------------------------------------------------------------
    % Now dif with fitted NEXRAD data
    % Read in prior data
    % Use previously defined lat_nas,lon_nas
    
    switch n % TEMP method of retrieving stored fit data
        case 1
            refl_nex_fitted = refl_nex_fitted_1;
        case 2
            refl_nex_fitted = refl_nex_fitted_2;
        case 3
            refl_nex_fitted = refl_nex_fitted_3;
        case 4
            refl_nex_fitted = refl_nex_fitted_4;
    end
    
    % Perform pixelwise dif and RMSE
    refl_dif = refl_comp - refl_nex_fitted;
    
    SE = 0;
    for x = 1:dims(1)
        for y = 1:dims(2)
            SE_addition = ((refl_comp(x,y) - refl_nex_fitted(x,y))^2);
            if(~isnan(SE_addition))
                SE = SE + SE_addition;
            end
        end
    end
    MSE = SE/(dims(1)*dims(2));
    RMSE = sqrt(MSE);

    sprintf('RMSE for EXRAD-NEXRAD dif case %d = %d',n,RMSE)
    
    % We're plotting the dif on the NAS grid
    figure('Position',[0 200 2100 620]);
    h = pcolor(lon_slice,lat_slice,refl_dif);
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    colormap(redblue);
    %q = quiver(lon_slice,lat_slice,wind_x,wind_y,3,'black','linewidth',1.1,'MaxHeadSize',1); 
    %q.MaxHeadSize = 10;
    caxis([-30 30]);
    ax = gca();
    ax.FontSize = axes_font_size;
    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf("EXRAD - NEXRAD Composite Reflectivity Dif (dBZ) | %s-%s",date,string(midpoints(n))),'FontSize',title_font_size);

    saveas(h,sprintf(output_format_dif,n));
    %saveas(h,strcat('output/','refl_',string(n),'_nex.png'));
    
    close all;

end
