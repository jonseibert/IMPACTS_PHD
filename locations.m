% Spatial limits of analysis and plot (degrees lat/lon)
w_lim = -77;
e_lim = -69.75;
s_lim = 36;
n_lim = 46;


% Plot locations on map
path_to_borders = './borders';
addpath(path_to_borders);

fig_x = 100;
fig_y = 100;
fig_width = 900; 
fig_height = 900;
fig_width_small = 350;
fig_height_small = 350;


f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
hold on;
h = scatter(loc_lon_list,loc_lat_list);
title(sprintf('Plume Locations | %s',timestamp_title),'FontSize',title_font_size); 
xlabel('Longitude','FontSize',label_font_size);
ylabel('Latitude','FontSize',label_font_size);
set(gca,'Fontsize',axes_font_size);
borders('continental us','black','linewidth',1); 
xlim([w_lim e_lim]);
ylim([s_lim n_lim]);

text(loc_lon_list+0.08,loc_lat_list+0.08,loc_name_list);

saveas(h,'output/Plume_locations.png');
hold off;