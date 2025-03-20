% radial_velocity.m
%% Calculate partial-3D radial velocity from 3D wind fields relative to a single radar point
%  For now: KBUF
% 
% Author(s): Jon Seibert
% Last updated: 30 Nov 2021

% CURRENT PROBLEMS:
% - HOW TO PLOT NANs IN A DIFFERENT COLOR (GREY)
% - INTERPOLATING Vr SEEMS TO BE FAILING TO PRODUCE ANYTHING BUT HUGE
% NUMBERS- INVESTIGATE/TRY CREATING THE VALUES FROM SCRATCH PER ANGLE INSTEAD?


%% 0. Options

% Specify paths
base_path = pwd;
input_path = "nasa_data";
intermediate_path = "intermediate";
output_path = "output/radial_velocity";

% File lists
mean_files = string(ls("wrf_data/ncea/wrf*"));
mem_files = string(ls("wrf_data/m1/wrf*"));
ex_files = string(ls("nasa_data/exrad_impacts*"));
radar_data_filename = "NEXRAD_data.mat";

% Plot font sizes
title_font_size = 18;
label_font_size = 14;
axes_font_size = 12;

% Plot dimensions
plot_width = 1600;
plot_height = 900;

radar_idx = 4; % Radar choice index [4 = KBUF]
radar_pad = 0.5; % Amount of +- padding on the radar elevation angle to determine beam coverage
angle_idx = 3; % 1-39 : 0.5-19.5 in 0.5 steps
cmap = 1;
exrad_file_choice = 1; % Which EXRAD file to sample for development

% Specify the list of radar elevation angles
% Try elevations: 0.5*, 10* to start (0.5-19.5 is WSR-88D range, in
% increments of 0.5*)
elevation_angles = single(0.5:0.5:19.5); 

% Specify maximum range of radar
radar_max_dist = 200*1.852*1000; %200 nautical miles * (1.852 km/nm) * (1000 m/km) [FOR WIND- REFL MAX DIST IS 150KM, FROM GEM]
radar_max_dist = radar_max_dist*100000; % TESTING

%% 1. Extract sample of EXRAD 3D wind data to work with

filename = sprintf("%s/%s",input_path,ex_files(exrad_file_choice));

date = extractBetween(filename,15,36);

% Dims: lat x lon x 57 (across x along x height) [41 x 517 x 57]
across_track_grid = ncread(filename, "across_track_grid"); % 41x1 [km]
along_track_grid = ncread(filename, "along_track_grid"); % 517x1 [km]
vertical_grid = ncread(filename, "vertical_grid"); % 57x1 [km]
zonal_wind = ncread(filename, "zonal_wind"); % 3D [m/s]
meridional_wind = ncread(filename, "meridional_wind"); % 3D [m/s]
vertical_wind = ncread(filename, "vertical_wind"); % 3D [m/s]
Ku_band_reflectivity = ncread(filename, "Ku_band_reflectivity"); % 3D [dBZ]
latitude = ncread(filename, "latitude"); % 3D [degrees]
longitude = ncread(filename, "longitude"); % 3D [degrees]
zonal_wind_std = ncread(filename, "zonal_wind_std"); % 3D [m/s]
meridional_wind_std = ncread(filename, "meridional_wind_std"); % 3D [m/s]
vertical_wind_std = ncread(filename, "vertical_wind_std"); % 3D [m/s]

longitude = longitude - 360; % Convert to negative values [ONLY WORKS BC ALL ARE > 180]
dims = size(longitude); % Specify dimensions
lat = dims(1);
lon = dims(2);
height = dims(3);

u_full = zonal_wind;
v_full = meridional_wind;
w_full = vertical_wind;

num_angles = length(elevation_angles);

%% 1. Convert sample to radial velocity relative to KBUF, where in range and in beam

% Radial velocity: Vr
% Elevation angle: theta-e
% Azimuthal angle: phi
% Ground-relative range from the radar: r

% Vr = u*sin(phi)*cos(theta_e) + v*cos(phi)*cos(theta_e) + w*sin(theta_e)

% Radar azimuth angle defined as clockwise from N
% Math angle defined as counter-clockwise from E
% To convert to math form in radians from radar form in degrees:
% phi_math = -1*(phi - 90)*(pi/180)

% HOW TO COMPUTE WHETHER A POINT IS CLOSE ENOUGH TO THE RADAR BEAM TO BE
% INCLUDED AT A GIVEN ELEVATION ANGLE:
%   - Beam width can be considered an extra 0.5* on either side of the
%   central angle: elevation 0.5 includes 0*-1*
%   - therefore, to determine if a point is within the beam at a specific
%   elevation angle, calculate the height of the top and bottom of the beam
%   and compare to height of point

height_list = vertical_grid*1000; % Convert heights from km to m

% Load radar data
load(sprintf("%s/%s",intermediate_path,radar_data_filename),'lons','lats','ids');

id_radar = ids(radar_idx);
lat_radar = lats(radar_idx); %lat_radar = 42.9489; % for KBUF
lon_radar = lons(radar_idx); %lon_radar = 78.7367;

num_angles = length(elevation_angles);

Vr = single(zeros(lat,lon,height,num_angles)); % Empty holding variable for radial velocity

% Loop over all points in the 3D wind data for each angle
% For each point: is it in range of the radar? is it within the beam at the
% current elevation angle? if both are true, compute radial velocity.
% otherwise, NaN
for angle_idx = 1:length(elevation_angles)
    angle = elevation_angles(angle_idx);
    for x = 1:lon
        for y = 1:lat
            % Check if point is within horizontal range (ignores height)
            % Note that lat/lon does not vary with height for EXRAD data
            dist = haversine_distance(longitude(y,x,1),latitude(y,x,1),lon_radar,lat_radar);
            if(dist <= radar_max_dist)

                % Compute azimuth angle
                lat_target = latitude(y,x,1);
                lon_target = longitude(y,x,1);
                az = azimuth('gc',lat_radar,lon_radar,lat_target,lon_target); % Great circle azimuth
                az_math = -1*(az - 90)*(pi/180); % convert azimuth to counterclockwise from East, in radians [math form]

                % Compute height of beam edges
                lower_bound = dist*tand(angle - radar_pad);
                upper_bound = dist*tand(angle + radar_pad);

                % If point falls within beam, compute Vr, else NaN
                for z = 1:height
                    h = height_list(z);
                    % Extract wind values at the chosen point
                    u = u_full(y,x,z);
                    v = v_full(y,x,z);
                    w = w_full(y,x,z);
                    if(h < upper_bound && h > lower_bound)
                        Vr(y,x,z,angle_idx) = u*sin(az_math)*cosd(angle) + v*cos(az_math)*cosd(angle) + w*sind(angle);
                    else
                        Vr(y,x,z,angle_idx) = NaN;
                    end
                end
            else
                Vr(y,x,:,angle_idx) = NaN;
            end
        end
    end
end

save(sprintf("%s/EXRAD_Vr_full_%d",intermediate_path,exrad_file_choice),'Vr');

%% Skip above section- reload saved result

load(sprintf("%s/EXRAD_Vr_full_%d",intermediate_path,exrad_file_choice),'Vr');

%% Reformat Vr to progress along elevation angles instead of height

% Set up blank matrix [x y angle]

Vr_angled = single(zeros(lat,lon,num_angles));

% For each point:
for y = 1:lat
    for x = 1:lon
        % For each angle in the the list of elevation angles:
        for angle_idx = 1:length(elevation_angles)
            angle = elevation_angles(angle_idx);
            % Calculate distance
            dist = haversine_distance(longitude(y,x,1),latitude(y,x,1),lon_radar,lat_radar);
            
            % If point is within range:
            if(dist <= radar_max_dist)
                % Calculate height of beam center
                beam_height = dist*tand(angle);
                
                % Locate nearest data points above and below height of beam
                nn_z = [0 0];
                nn_val = [0 0];
                
                for z = 2:height
                    if(height_list(z) >= beam_height)
                        nn_z = [height_list(z-1) height_list(z)];
                        nn_val = [Vr(y,x,z-1,angle_idx) Vr(y,x,z,angle_idx)];
                        break;
                    end
                end
                
                % If beam falls outside the vertical range or hits its edge s.t.
                % above/below points fall outside, or either vertical neighbor
                % point is NaN: NaN
                if(nn_z(1) <= 0) || any(isnan(nn_val))
                    Vr_angled(y,x,angle_idx) = NaN;
                else
                    % Interpolate to beam
                    % Place value in its slot in the matrix of [x y angle]
                    proportion = (beam_height - nn_z(1))/(nn_z(2) - nn_z(1)); % Proprtional distance between beam and lower data point
                    Vr_angled(y,x,angle_idx) = nn_val(1)*proportion + nn_val(2)*(1 - proportion);
                end
            end
        end
    end
end

Vr_angled_iV = Vr_angled;

%% Plot along elevation angles

for angle_idx = 1:5
%for angle_idx = 1:length(elevation_angles)
    angle = elevation_angles(angle_idx);

    % Plot
    figure('Position',[0 200 2100 620]);
    h = pcolor(longitude(:,:,1),latitude(:,:,1),Vr_angled(:,:,angle_idx));
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    if(cmap == 1)
        cmap_name = "jet";
        colormap(jet);
    else
        cmap_name = "redblue";
        colormap(redblue)
    end
    caxis([-20 20]);

    ax = gca();
    ax.FontSize = axes_font_size;

    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf("Radial Velocity Test: EXRAD | %g* elevation | KBUF",angle),'FontSize',title_font_size);

    saveas(h,sprintf("%s/radial_velocity_iVr_e%d_c%s.png",output_path,angle_idx,cmap_name));
    close all;

end


%% Try creating the per-angle Vr from scratch instead of interpolating
% Instead of computing and interpolating Vr, interpolate 3D wind and use it
% to directly compute angle-relative Vr

Vr_angled = single(zeros(lat,lon,num_angles)); % Empty holding variable for radial velocity

% Loop over all points in the 3D wind data for each angle
% For each point: is it in range of the radar? If true, interpolate the 3D
% winds to the radar beam centerpoint and compute Vr if possible (has above
% and below to interpolate from)
% otherwise, NaN

for x = 1:lon
    for y = 1:lat
        % Check if point is within horizontal range (ignores height)
        % Note that lat/lon does not vary with height for EXRAD data
        dist = haversine_distance(longitude(y,x,1),latitude(y,x,1),lon_radar,lat_radar);
        if(dist <= radar_max_dist)

            % Compute azimuth angle
            lat_target = latitude(y,x,1);
            lon_target = longitude(y,x,1);
            az = azimuth('gc',lat_radar,lon_radar,lat_target,lon_target); % Great circle azimuth
            az_math = -1*(az - 90)*(pi/180); % convert azimuth to counterclockwise from East, in radians [math form]

            for angle_idx = 1:length(elevation_angles)
                angle = elevation_angles(angle_idx);

                % Compute height of beam center
                beam_height = dist*tand(angle);

                % If beam height falls within two interpolate-able heights,
                % compute Vr: otherwise NaN
                Vr_angled(y,x,angle_idx) = NaN;
                for z = 1:height
                    h = height_list(z);
                    
                    if(h > beam_height && z > 1) % If we passed the beam height moving upward
                        lower_bound = height_list(z-1);
                        upper_bound = h;
                        
                        % Compute proportional distance between beam center
                        % and neighbor heights
                        proportion = (beam_height - lower_bound)/(upper_bound - lower_bound);
                        
                        % Extract wind values at the chosen point
                        u = u_full(y,x,z-1)*(1-proportion) + u_full(y,x,z)*(proportion);
                        v = v_full(y,x,z-1)*(1-proportion) + v_full(y,x,z)*(proportion);
                        w = w_full(y,x,z-1)*(1-proportion) + w_full(y,x,z)*(proportion);
                        
                        % Compute Vr
                        Vr_angled(y,x,angle_idx) = u*sin(az_math)*cosd(angle) + v*cos(az_math)*cosd(angle) + w*sind(angle);
                        break; % If found the correct heights and computed, no need to keep looping for this [y x angle] instance
                    end
                end
            end
        else
            Vr_angled(y,x,:) = NaN;
        end
    end
end

save(sprintf('%s/EXRAD_Vr_angled_%d',intermediate_path,exrad_file_choice),'Vr_angled');

%load(sprintf('%s/EXRAD_Vr_angled_%d',intermediate_path,exrad_file_choice),'Vr_angled');

Vr_angled_iW = Vr_angled;


%% Plot along elevation angles

for angle_idx = 1:5
%for angle_idx = 1:length(elevation_angles)
    angle = elevation_angles(angle_idx);

    % Plot
    figure('Position',[0 200 2100 620]);
    h = pcolor(longitude(:,:,1),latitude(:,:,1),Vr_angled(:,:,angle_idx));
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    if(cmap == 1)
        cmap_name = "jet";
        colormap(jet);
    else
        cmap_name = "redblue";
        colormap(redblue)
    end
    caxis([-20 20]);

    ax = gca();
    ax.FontSize = axes_font_size;

    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf("Radial Velocity Test: EXRAD | %g* elevation | KBUF",angle),'FontSize',title_font_size);

    saveas(h,sprintf("%s/radial_velocity_iW_e%d_c%s.png",output_path,angle_idx,cmap_name));
    close all;

end

%% testing: plotting vs. z

angle_idx = 2;
for z = 1:10
    
    plot_height = vertical_grid(z);

    % Plot
    figure('Position',[0 200 2100 620]);
    h = pcolor(longitude(:,:,1),latitude(:,:,1),Vr(:,:,z,1));
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    if(cmap == 1)
        cmap_name = "jet";
        colormap(jet);
    else
        cmap_name = "redblue";
        colormap(redblue)
    end
    caxis([-20 20]);

    %q = quiver(lon_slice,lat_slice,wind_x,wind_y,3,'black','linewidth',1.1,'MaxHeadSize',1); 
    %q.MaxHeadSize = 10;

    ax = gca();
    ax.FontSize = axes_font_size;

    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf("Radial Velocity Test: EXRAD %gkm %g* elevation",plot_height,angle),'FontSize',title_font_size);

    saveas(h,sprintf("%s/radial_velocity_e%d_z%d_c%s.png",output_path,angle_idx,z,cmap_name));
    close all;

end
