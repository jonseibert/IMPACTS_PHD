% generate_n0q_list.m
% Purpose: Creates a .txt file of all the IEM Mesonet Archive n0q .png reflectivity
%          composites for a given day, for ease of using wget to download them.
% Author: Jon Seibert
% Date last modified: 27 April 2022
% Notes: 
%   - In current form, creates the .txt file in the directory where this
%     script is run. 
%   - Substitute desired date when using. 
%   - n0q files come in intervals of 5 minutes, so there should be 12 
%     per hour for 288 total per day.
%   - "wget -i n0q_files.txt" will download all the files in the list

%% Settings

year = 2020;
month = 2;
day = 8;

%% Main

% Reference
% path = 'https://mesonet.agron.iastate.edu/archive/data/2020/02/07/GIS/uscomp/';
% base_string = 'n0q_202002070000.png';

fID = fopen("n0q_files.txt",'w');

format = 'https://mesonet.agron.iastate.edu/archive/data/%04.f/%02.f/%02.f/GIS/uscomp/n0q_%04.f%02.f%02.f%02.f%02.f.png\n';

minute_counter = 0;
hour_counter = 0;

for i = 1:288
    
    if(minute_counter == 12)
        minute_counter = 0;
        hour_counter = hour_counter + 1;
    end
    
    min = minute_counter*5;
    hr = hour_counter;
    
    fprintf(fID,format,year,month,day,year,month,day,hr,min);
    
    minute_counter = minute_counter + 1;

end

fclose(fID);