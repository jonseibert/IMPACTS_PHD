
% Computes the Haversine distance between two pairs of latitude and
% longitude points
function d = haversine_distance (lon1,lat1,lon2,lat2)
% Inputs must be in degrees!

R = 6371000; % Radius of Earth (m)

a = ((sind((lat2-lat1)/2))^2) + cosd(lat1)*cosd(lat2)*((sind((lon2-lon1)/2))^2);
d = R*2*atan2(sqrt(a),sqrt(1-a)); % Distance in m

end