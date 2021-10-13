clear 

% preview the big file:

% ds = datastore('124108_201911_PKSAq_Offshore_data.txt');
% preview(ds) 


% set save name for mat file
savename = 'LS_PKSAq_prelim.mat';

dateformatIN = 'yyyy-mm-dd HH:MM:ss.FFF';
T = readtable('124108_201911_PKSAq_marsh_data.txt');
s124108.time = datenum(T.Time); 
% note: datetime objects don't need a format input (for matlab v2020?)
% see: https://stackoverflow.com/questions/45234957/apply-datenum-to-a-column-within-table
s124108.pdata = T.Pressure;
clear T

T = readtable('124107_201911_PKSAq_offshore_data.txt');
s124107.time = datenum(T.Time);
s124107.pdata = T.Pressure;
clear T

% -- Add sensor elevation data to arrays- - 
% height of transducer above marsh bed:  
s124108.z = 0.085; % 0ffshore
s124107.z = 0.125; % Marsh

% % NAVD elevation of transducer = m_navd  - S
% this was not measured for prelim deployment 
% s124107.sensor_elev_navd = -0.416;  % Offshore
% s124108.sensor_elev_navd = -0.421;   % 0 m 
% s124109.sensor_elev_navd = -0.121;  % 5 m
% s41428.sensor_elev_navd  = 0.142;    % 10 m
% s41429.sensor_elev_navd  = 0.294;    % 20 m 

save(savename,'s124108','s124107');