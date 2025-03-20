intermediate_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
load(sprintf('%s/%s',intermediate_path,'nexrad_stations.mat')); % wban, station_ids,_station_names,lon_deg,lat_deg,elevations,tower_heights

elevations_meters = elevations*0.3048;
elevations = elevations_meters;

save(sprintf('%s/%s',intermediate_path,'nexrad_stations.mat'),'wban','station_ids','station_names','lon_deg','lat_deg','elevations','tower_heights'); % wban, station_ids,_station_names,lon_deg,lat_deg,elevations,tower_heights
