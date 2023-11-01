% Example wave transformation across Deerfield marsh 
clear 

%% figure set up 
dockit = @()set(gcf,'windowstyle','docked'); 
fig = figure; dockit()
%% ---- Water & Wave Parameters 
g     = 9.81; 
rho   = 1025;
BRK   = 4;             % wave breaking : 
                       %     0 = off, 
                       %     1 = TG83 eq 25 w/ bottom slope (after Salmon et al 2015), 
                       %     2 = TG83 eq 26 w/ bottom slope
                       %     3 = TG83 eq 25 w/ constant, 
                       %     4 = TG83 eq 26 w/ constant 
gammac = 0.6;          % wave breaking parameter
deplim = 0.1;          % lower limit for local water depth (m)
T = 1.5;               % peak period of the waves 
Hrms0 = 0.4;           % rms wave height at x = 0 (offshore boundary) 
h0  = 1.4;               % water depth at x = 0
wavemode = 'rand' ;    % should be 'rand', for random waves. 
                       %    ('mono' is for regular waves, not fully implemented)

%% --- Veg Parameters (to turn off veg, vegmode=1 & bv=0 or Lv1=0) 
DRAG = -3;            % DRAG > 0  : constant value
                       %      = -1 : Cd(KC) 
                       %      = -2 : Cd(KC & plant emergence) 
                       %      = -3 : Cd(Re) from AS2014
                       %      = -4 : Cd(KC) from AS2014
                       %      = -5 : Cd(KC) from Jadhav et al 2013...
                       %             need to spec stem height as ls to
                       %             use -5
v0    = 35;            % distance from x = 0 to start of veg canopy (m)  (35 for full, 5 for marsh, 0 for S3S4)
Lv   = 20;             % vegetation length in profile [m]      
vegmode = 2;           % 1 = constant, 2 = spatially varying, linear inter.
E = 8e7;               % Pa, young's mod for flexural rigidity (used for Cd(Ca) options) 
% Enter constant values averaged over area of interest: 
lsmean   = 0.6;           % vegetation height (m) = (alpha*h) 
bvmean = 0.004;        % stem diameter / blade width (m)
Nvmean   = 392;           % density (stems/m2) -- avg over DF18 cross-section 

% load the measured veg params and their location from S0: 
vegin = load('../veg/DF3vegin.txt'); % used if VEGMODE = 2; 
%% --- Bottom Friction Parameters
Cf = 0.003;             % friction coeff
FRIC = 0;              
%% --- Grid Parameters 
elevin = load('../elev/DF18_elev_full.txt');
dx = 0.1;              % choose the model grid size (meters) 
% Locations for validation; need to correspond w elevin (ie, from x=0): 
Sta = [5];   % or [] if there are none 
%% make the grid and related info: 
[x,veg,z,xn,idxsta,Nv,bv,ls] = gridgen(elevin,...
                                dx,v0,Lv,Sta,'pci',...
                                vegmode,Nvmean,bvmean,lsmean,vegin);
%% plot the profile to check gridgen: 
s1 = subplot(2,1,1);
plot ([0 xn], [h0 h0], 'b'); 
hold on 
area(x,z,min(z),'FaceColor',[0.6 0.6 0.6]);

hold on 
for i = 1 : length(x)-1
    if veg(i) == 1
        plot([x(i) x(i)],[z(i) z(i)+ls(i)],'Color',[76 153 0]/255);
    end 
end 
for i = 1:length(idxsta)
    scatter(x(idxsta(i)),z(idxsta(i)),15,'r','filled')
end
%% run the model 
n  = length(x);


[Hrms,h,Dv,Db,Df,gamma,KC,Re,Ca,Cd,sigma,L,k,Ur,S,Er,Uc,a] = ...
    marshwavemodel(dx,n,Hrms0,h0,z,T,veg,ls,g,rho,Nv,bv,E,...
    gammac,Cf,BRK,FRIC,deplim,DRAG,wavemode);



% plot the output rms wave heights across the marsh profile:  
s2 =  subplot(2,1,2);
plot(x,Hrms,'linewidth',1.2);   

%% figure formatting
s1.XLim = [0 xn];
s2.XLim = [0 xn];
s2.YLim = [0 Hrms0];

s1.FontSize = 9;
s1.XTickLabel = [];
s2.FontSize = 9;

ylabel(s1,'elev. (m)');
ylabel(s2,'H_{rms} (m)');
xlabel(s2,'distance from offshore boundary (m)');

s1.Position = [0.1 0.6 0.85 0.35];
s2.Position = [0.1 0.15 0.85 0.35];

s1.YAxis.Label.FontSize = 10;
s2.YAxis.Label.FontSize = 10;

