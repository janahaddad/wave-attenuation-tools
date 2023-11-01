%% discretize elev profile for marsh wave model 
% inputs:
%   elevin = 2 column text file with col 1 = x, col 2 = elevation (navd88,m)
%   dx = x-grid interval
%   v0 = distance of start of veg from the offshore model boundary 
%   Lv = length of veg 
%   Sta = distances of wave stations from offshore boundary
%           (input [] if there are no stations)
%   method = 'pci' or 'lin' for elevation interpolation method 
%   vegmode = 1 for uniform veg, 2 for varying veg 
%   Nvmean = mean stem density for uniform veg
%   bvmean = mean stem diameter for uniform veg
%   lsmean = mean shoot length for uniform veg 
%   vegin = varying veg input (distances from x=0, Nv, bv, ls, and std's)
%         (see text file formats in elev directory)
%% 
function [x,veg,z,xn,idxsta,Nv,bv,ls] = gridgen(elevin,...
                                dx,v0,Lv,Sta,method,vegmode,...
                                Nvmean,bvmean,lsmean,vegin)
x0 = 0; 
z_orig = elevin(:,2);
xmeters = elevin(:,1)-elevin(1,1); % shift in x so x(1) = 0; 
xn = floor(xmeters(end)*10)/10;    % round down to nearest 0.1
x = x0:dx:xn;  
if strcmp(method,'pci')==1
z_navd = interp1(xmeters, z_orig, x,'pchip');
elseif strcmp(method,'lin')==1
z_navd = interp1(xmeters, z_orig, x,'linear','extrap');
end 
z = z_navd + (-min(z_navd));
%StaShift = Sta-elevin(1,1);
% ID stations on grid:
idxsta = zeros(1,length(Sta));
for i = 1:length(Sta)
    %idxsta(i) = find(abs(x-StaShift(i))<0.001);
    idxsta(i) = find(abs(x-Sta(i))<0.001);
end 
% Find grid points with vegetation:
veg = zeros(size(x));
if Lv>0 
   % v0 = v0-elevin(1,1);
   xveg = [v0 v0+Lv];     % location of vegetation on x grid 
   veg(x >= xveg(1) & x < xveg(2)) = 1;
end 
% interpolate vegetation to the grid (vegmode 2) or use constant (vegmode
% 1):
Nv = zeros(1,length(x));
bv = zeros(1,length(x)); 
ls = zeros(1,length(x)); 
if vegmode == 2
    %xveg = 0:dx:Lv;  % 0 here is veg start point
    %Nv = NaN(size(x));

    VegShift = vegin(:,1) - elevin(1,1);
    Nv(x >= xveg(1) & x < xveg(2)) = ...
        interp1(VegShift,vegin(:,2),...
        x(x >= xveg(1) & x < xveg(2)),'linear','extrap');
    bv(x >= xveg(1) & x < xveg(2)) = ...
        interp1(VegShift,vegin(:,4),...
        x(x >= xveg(1) & x < xveg(2)),'linear','extrap');
    ls(x >= xveg(1) & x < xveg(2)) = ...
        interp1(VegShift,vegin(:,6), ...
        x(x >= xveg(1) & x < xveg(2)),'linear','extrap');
elseif vegmode ==1            
    Nv(x >= xveg(1) & x < xveg(2))  = Nvmean;
    bv(x >= xveg(1) & x < xveg(2)) = bvmean; 
    ls(x >= xveg(1) & x < xveg(2)) = lsmean; 
end 


end 