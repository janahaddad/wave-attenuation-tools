clear
%% set up 
PLOT  = 0;  
SAVE  = 0;   % on/off save variables to mat file.

% name of mat file to save sensor varialbes to:
savename = 'PKSaq_prelim_spectra.mat';   

met=load('./met/clkn7h2019.txt');   % load met file 
load ('./LS_PKSAq_prelim.mat');     % load sensor array pdata

% clip data to an input start and end: 
startdate = datenum(2019,09,28,00,00,0);   
enddate = datenum(2019,10,25,00,00,0);     
% sensor SNs:
SN=[124107]; 
% Wave calc parameters:
rho = 1018; % density of water 
g   = 9.81; % acceleration due to gravity
f_s = 16;   % sampling frequency in Hz
% Format atmospheric pressure data
metdatetime_old=datenum(met(:,1),met(:,2),met(:,3),met(:,4),met(:,5),zeros(length(met),1));
metdatetime=metdatetime_old-4/24; % convert from UTC 
metatmosp=met(:,9)/100;           % hPa to dbar 

%% Correct pressure measurements for atmospheric pressure
for j=1:length(SN)
    eval(['psensor=s',num2str(SN(j)),';'])
% clip pdata and time series to start and end dates :
    istart = find(psensor.time == startdate);
    iend = find(psensor.time == enddate);
    psensor.pdata = psensor.pdata(istart: iend);
    psensor.time = psensor.time(istart: iend);
% interpolate met to sensor's pdata timeseries :
    atmosp_interp = interp1(metdatetime, metatmosp, psensor.time,'linear'); 
% SNxxxxxx.dataadj = pressure adjusted for interpolated atmosphere :
    psensor.dataadj = psensor.pdata-atmosp_interp;     %(pressure in dbars) 
% SNxxxxxx.pcor = pressure in Pascals :
    psensor.pcor = psensor.dataadj*10000; 
% Sensor depths and total water depth:
    psensor.sensor_depth = psensor.pcor/(rho*g);
    psensor.water_depth = psensor.sensor_depth + psensor.z ; 
    %psensor.eta_navd = psensor.sensor_depth + psensor.sensor_elev_navd;
    
% build moving average of data over 5 minute window:
    %psensor.eta_navd_mean_mov = smooth(psensor.eta_navd,(10*60*f_s),'moving');
    %psensor.eta_navd_mean = smooth(psensor.eta_navd,(10*60*f_s),'sgolay');
%   psensor.water_depth_mean = smooth(psensor.water_depth,(10*60*f_s),'sgolay');
% set everything below navd cutoff depth to NaN
    %clipid = find(psensor.eta_navd_mean<0.35); 
    %psensor.eta_navd_mean(clipid) = NaN; 
% write new data to sensor array :
    eval(['s',num2str(SN(j)),'=psensor;'])  
end
clear psensor

%% a way to estimate error via max diff in eta
% % find largest delta in eta across the marsh
% A = nchoosek(SN,2);
% for i = 1:length(A)
%     strsensor1 = ['s' num2str(A(i,1))];
%     psensortemp = eval(strsensor1);
%     eta1 = psensortemp.eta_navd_mean; 
%     strsensor2 = ['s' num2str(A(i,2))];
%     psensortemp = eval(strsensor2);
%     eta2 = psensortemp.eta_navd_mean; 
%     M(:,i) = abs(eta1 - eta2); 
% end 
% clear psensortemp
% maxdelta = max(max(M)); 
%% Plot 5min-mean surface elevations or water depth
switch PLOT
case 1  
cmap = cmocean('haline',7);
figure;
plot(s124107.time, s124107.water_depth_mean, 'Linewidth', 2,'color',cmap(1,:));
hold on
plot(s124108.time, s124108.water_depth_mean, 'LineWidth', 2,'color',cmap(2,:));
legend('s124107 = offshore','s124108 = marsh');
grid on
grid minor
dateaxis('x',6); % datetick('x')
xlabel('day')
ylabel('water depth, m')
end

%% WAVE ANALYSIS

avint=10*60*f_s; % number of samples in a block, for spectrum calculation
nspec=2^10;      % number of samples to use in each fft calculation 
   
for j = 1:length(SN)
tic
    eval(['z=s',num2str(SN(j)),'.z;']);
    eval(['t=s',num2str(SN(j)),'.time;']);
    eval(['p=s',num2str(SN(j)),'.pcor;']);  %sensor pressure in Pa
    N=length(p); % length of entire record
    Nblocks=floor(N/avint); % number of blocks in entire record  
    
% Set the size of new variables before loop 
    t_block=NaN*ones(Nblocks,1); 
    p_mean=NaN*ones(Nblocks,1);
    h_block= NaN*ones(Nblocks,1);
    Paa=NaN*ones(Nblocks,nspec/2+1);
    Ppp= NaN*ones(Nblocks,nspec/2+1);
    p_tilde= NaN*ones(Nblocks,avint);
    badblocks=zeros*ones(Nblocks,1);
    f_c=NaN*ones(Nblocks,1);
    Paa_new=NaN*ones(Nblocks,nspec/2+1);
    Paa_new_smooth=NaN*ones(Nblocks,nspec/2+1);
    sensor_depth=NaN*ones(Nblocks,1);
    Ppp_max = NaN*ones(Nblocks,1);
    Ppp_noise = NaN*ones(Nblocks,1);
    f_peak=NaN*ones(Nblocks,1);
    f_peak_Ppp=NaN*ones(Nblocks,1);
    f_peak_smooth=NaN*ones(Nblocks,1);

% additional variables if the error estimate method is used:
%     h_block_low = NaN*ones(Nblocks,1);
%     h_block_up  = NaN*ones(Nblocks,1);
%     Paa_low = NaN*ones(Nblocks,nspec/2+1);
%     Paa_up = NaN*ones(Nblocks,nspec/2+1);
%     Paa_new_low = NaN*ones(Nblocks,nspec/2+1);
%     Paa_new_up = NaN*ones(Nblocks,nspec/2+1);

% Loop through blocks in time series
    for i=1:Nblocks
        istart=(i-1)*avint+1; % index for start of block i
        iend=i*avint; % index for end of block i 
        t_block(i)=t(istart); % start time for each block
        p_block=p(istart:iend,1); % select ith block of pressure data
        p_mean(i,1)=mean(p_block); % compute mean p for each block (Pa)
        
% if less than 10 cm water (1000 Pa) above the sensor,
% then log the bad data block and set the mean p value to NaN. 
        if p_mean(i,1)<1000 
            badblocks(i)=1;
            p_mean(i,1)=NaN;
% if more than 10 cm of water (1000 Pa) above the sensor... 
        else
            p_tilde0=p_block-p_mean(i,1); % subtract mean from each block to get fluctuations
            p_tilde(i,:)=detrend(p_tilde0); % remove linear trend from pressure measurements in block i
            sensor_depth(i,1)=p_mean(i,1)/(rho*g); % depth of sensor below water surface
            h_block(i,1)=sensor_depth(i,1)+z; % total water depth for each block
% if the error estimate method is used:
%           h_block_up(i,1)=sensor_depth(i,1)+z+maxdelta;
%           h_block_low(i,1)=sensor_depth(i,1)+z-maxdelta;
%
% compute pressure spectrum for each block, segments of length nspec, 50% overlap
           [Ppp(i,:),f]=pwelch(p_tilde(i,:),hanning(nspec),...
               nspec/2,nspec,f_s); 
           f=f';
           [~,k]=wavenumber(2*pi*f,h_block(i)); % compute wavenumber that corresponds to each frequency 
% convert pressure spectrum to wave amplitude spectrum
           Paa(i,:)=Ppp(i,:)/(rho*g)^2.*...
               (cosh(k*h_block(i))./cosh(k*z)).^2;
% if the error estimate method is used:
%            % upper limit: 
%            [~,k]=wavenumber(2*pi*f,h_block_up(i));
%            Paa_up(i,:)=Ppp(i,:)/(rho*g)^2.*...
%              (cosh(k*h_block_up(i))./cosh(k*z)).^2; 
%            % lower limit:
%            [~,k]=wavenumber(2*pi*f,h_block_low(i));
%            Paa_low(i,:)=Ppp(i,:)/(rho*g)^2.*...
%              (cosh(k*h_block_low(i))./cosh(k*z)).^2; 
%
% Cut off at local min in original wave height spectrum, after
% the peak found in the pressure spectrum: 
           Ppp_max(i)=max(Ppp(i,:)); % find peak value of spectrum
           idx_maxPpp=find(Ppp(i,:)==Ppp_max(i)); % id Ppp peak 
           [~,idx_c0] = min(Paa(i,idx_maxPpp:end)); % id local min in Paa
           idx_c=idx_c0+idx_maxPpp-1; % calculate cutoff index 
           f_peak_Ppp(i) = f(idx_maxPpp);
% if cutoff frequency is not sufficiently high relative to wave peak, 
% set this block to NaNs
           if idx_c<1.1*idx_maxPpp % 1.1 factor -> Jones & Monismith (2007)
               badblocks(i)=2;
               Paa_new(i,:)=NaN;
% if the error estimate method is used:
%              Paa_up(i,:)=NaN;
%              Paa_low(i,:) = NaN;
               Paa_new_smooth(i,:) = NaN;
               f_c(i)=NaN;
% if cutoff frequency is too low, set this block to NaNs 
           elseif f(idx_c)<0.5 
               badblocks(i)=3;
               Paa_new(i,:)=NaN;
%                Paa_up(i,:)=NaN;
%                Paa_low(i,:) = NaN;
               Paa_new_smooth(i,:) = NaN;
               f_c(i)=NaN;
% if no value of the spectrum is sufficiently larger
% than the noise floor, set this block to NaNs
           elseif f(idx_maxPpp)<0.1 
               badblocks(i)=4;
               Paa_new(i,:)=NaN;
%                Paa_up(i,:)=NaN;
%                Paa_low(i,:) = NaN;
               Paa_new_smooth(i,:) = NaN;
               f_c(i)=NaN;
% otherwise, cut the spectrum off at the cutoff frequency
% and extrapolate using a f^-5 fit to the spectrum just below the cutoff
           else 
                f_c(i)=f(idx_c); % evaluate and store cutoff frequency
                idx_fit=idx_c-5:idx_c; % indices to use for f^-5 fit 
                hf_fit=fit(f(idx_fit)',Paa(i,idx_fit)',...
                    'a*x^-5','StartPoint',10^-5); 
% merge original spectrum below cutoff with spectral fit above cutoff
                Paa_new(i,:)=[Paa(i,1:idx_c),hf_fit(f(idx_c+1:end))']; 
%                 % cut & fit spectra for uppper:
%                 hf_fit_up=fit(f(idx_fit)',Paa_up(i,idx_fit)',...
%                    'a*x^-5','StartPoint',10^-5);
%                 Paa_new_up(i,:)=[Paa_up(i,1:idx_c),...
%                    hf_fit(f(idx_c+1:end))'];
%                 % cut & fit spectra for lower:
%                 hf_fit_low=fit(f(idx_fit)',Paa_low(i,idx_fit)',...
%                    'a*x^-5','StartPoint',10^-5);
%                 Paa_new_low(i,:)=[Paa_low(i,1:idx_c),...
%                    hf_fit(f(idx_c+1:end))'];
% Smooth out the spectrum to extract a reliable peak frequency: 
% moving average over 5 pt window: 
                span = 5; 
                kernel = ones(1, span) / span;
                % do convolution of spectrum with kernel
	            Paa_new_smooth(i,:) = conv(Paa_new(i,:), kernel, 'same');
                % id peak freq. of Paa and smoothed Paa
                idx_max = find(Paa_new(i,:)==max(Paa_new(i,:)));
                idx_max_smooth=find(Paa_new_smooth(i,:)==...
                    max(Paa_new_smooth(i,:)));
                f_peak(i) = f(idx_max);
                f_peak_smooth(i) = f(idx_max_smooth);          
           end
        end 
    end

% create new time series that starts as day 1 of the year you're in 
    t_block_new=t_block-datenum(2019,1,1,0,0,0);
    
%% compute wave stats from the spectra 
    df=f(2)-f(1); % freq bandwidth 
% freq cutoffs for wave stats
    f_start=find(f>.1,1,'first'); 
    f_end=find(f<4,1,'last');
% compute variance of surface elevation == zero moment of spectrum 
    eta_var=sum(Paa_new(:,f_start:f_end),2)*df; 
% convert to wave height 
    Hrms=sqrt(eta_var)*2*sqrt(2); 
    % Hm0 = 4*sqrt(eta_var); 
    
% if the error estimate method is used:
%     %m0 and Hrms for upper limit:
%     m0_up=sum(Paa_new_up(:,f_start:f_end),2)*df;
%     Hrms_up = sqrt(m0_up)*2*sqrt(2);
%     %m0 and Hrms for lower limit:
%     m0_low=sum(Paa_new_low(:,f_start:f_end),2)*df;
%     Hrms_low = sqrt(m0_low)*2*sqrt(2);
%     Hrms_errorband = max(Hrms_up - Hrms_low); 

% first moment of spectrum:
    temp = zeros(size(Paa_new));
    for i = 1:size(Paa_new,1)
       temp(i,:) = f.*Paa_new(i,:); 
    end
    m1 = sum(temp(:,f_start:f_end),2)*df;
    % Tmean = m0/m1 :
     T_m01 = eta_var./m1;

%% Save wave parameters to array

%     eval(['s',num2str(SN(j)),'.badblockcode=badblocks;'])
    eval(['s',num2str(SN(j)),'.t_block=t_block_new;'])
    eval(['s',num2str(SN(j)),'.Paa=Paa;'])
    eval(['s',num2str(SN(j)),'.Paa_new=Paa_new;'])
    eval(['s',num2str(SN(j)),'.Paa_new_smooth=Paa_new_smooth;'])
%     eval(['s',num2str(SN(j)),'.Ppp=Ppp;'])
%     eval(['s',num2str(SN(j)),'.Ppp_max=Ppp_max;'])
%     eval(['s',num2str(SN(j)),'.Ppp_noise=Ppp_noise;'])
    eval(['s',num2str(SN(j)),'.f=f;'])
    eval(['s',num2str(SN(j)),'.h_block=h_block;'])
    eval(['s',num2str(SN(j)),'.Hrms=Hrms;'])
    eval(['s',num2str(SN(j)),'.f_c=f_c;'])
    eval(['s',num2str(SN(j)),'.f_peak=f_peak;'])
    eval(['s',num2str(SN(j)),'.m0=eta_var;'])
    eval(['s',num2str(SN(j)),'.m1=m1;'])
    eval(['s',num2str(SN(j)),'.T_m01=T_m01;'])
    eval(['s',num2str(SN(j)),'.f_peak_smooth=f_peak_smooth;'])
%     eval(['s',num2str(SN(j)),'.Hrms_up=Hrms_up;'])
%     eval(['s',num2str(SN(j)),'.Hrms_low=Hrms_low;'])
%    eval(['s',num2str(SN(j)),'.Hrms_errorband=Hrms_errorband;'])

toc % ~0.5 to 1 min per sensor (6 day deployment) on my Dell XPS (m2018a)
end

%% save mat
if SAVE == 1
    save(savename, 's124107', 's124108');
end 
    
%% PLOT FIGURES
switch PLOT
case 1  
    
dockit=@()set(gcf,'windowstyle','docked');

for j = 1:length(SN)
     eval(['psensor=s',num2str(SN(j)),';']) % go to next sensor 
     
% Plot wave height spectrograms (spectrum plotted in color with 
% frequency on y axis. x axis is the time of the block)

imaxwvht = find(psensor.Hrms == max(psensor.Hrms));
maxwvhtday = psensor.t_block(imaxwvht);

figure     
dockit()
% Plot pressure spectrogram
subplot 311
pcolor(psensor.t_block',psensor.f',log10(real(psensor.Paa')))
shading flat
caxis([-0.5 5])
ylim([0 3])
hold on
plot(psensor.t_block',psensor.f_c,'k-')
plot([psensor.t_block(imaxwvht) psensor.t_block(imaxwvht)],[0 3],'w--');
ylabel('frequency (Hz)')
title(['Sensor ' num2str(SN(j)) ' Wave height spectrum, log_{10}(P_{aa})'])
colorbar
colormap jet
% xlim([t_block_new(imaxwvht-20) t_block_new(imaxwvht+20)])

% Plot wave height spectrogram
subplot 312
pcolor(psensor.t_block',psensor.f',log10(real(psensor.Paa_new')))
shading flat
caxis([-6 -2])
ylim([0 3])
hold on
plot(psensor.t_block',psensor.f_c,'k-')
hold on
plot([psensor.t_block(imaxwvht) psensor.t_block(imaxwvht)],[0 3],'w--');
ylabel('frequency (Hz)')
title(['Sensor ' num2str(SN(j)) ' Corrected wave height spectrum, log_{10}(P_{aa})'])
colorbar
colormap jet
%     xlim([t_block_new(imaxwvht-20) t_block_new(imaxwvht+20)])

% plot mean wave heights
subplot 313
plot(psensor.t_block,psensor.Hrms);
hold on 
scatter(psensor.t_block(imaxwvht),0.5,60,'r');
ylim([0 0.5]);
title(['Sensor ' num2str(SN(j)) ' Mean wave height (max = ' num2str(max(psensor.Hrms)) 'm)']);
xlabel('day of 2017')
ylabel('[m]')
%     xlim([t_block_new(imaxwvht-20) t_block_new(imaxwvht+20)])
end
end % end the PLOT switch
%% ADDITIONAL TEST PLOTS

%PLOT EXAMPLE SPECTRUM: pressure, original wave, new wave
% figure
% dockit()
% 
% subplot 311
% plot(f(:),log10(real(Paa(imaxwvht,:))));
% hold on
% plot([f_c(imaxwvht) f_c(imaxwvht)],[-8 8]);
% title(['Example original wave height spectrum, day ' num2str(maxwvhtday)])
% ylim([-8 8]);
% xlim([0 8]);
% ylabel('log Paa')
% 
% subplot 312
% plot(f(:),log10(real(Paa_new(imaxwvht,:))));
% hold on
% plot([f_c(imaxwvht) f_c(imaxwvht)],[-8 0]);
% title(['Example wave height spectrum, day ' num2str(maxwvhtday)])
% ylim([-8 0]);
% ylabel('log Paa')
% 
% subplot 313
% plot(f(:),log10(real(Ppp(imaxwvht,:))));
% hold on
% plot([f(1) f(1+ nspec/2)], [log10(Ppp_noise(imaxwvht)) log10(Ppp_noise(imaxwvht))] );
% hold on 
% scatter(f(fidx_max),log10(real(Ppp(imaxwvht,fidx_max))));    % spectral peak 
% scatter(f(fidx_c),log10(real(Ppp(imaxwvht,fidx_c))));    % freq & p spectrum at cutoff 
% title(['Example pressure spectrum, day ' num2str(maxwvhtday)])
% legend('Ppp','Noise Floor', 'spectral peak','cutoff');
% xlabel('freq, Hz')
% ylabel('log Ppp')
% ylim([-1 6]);
% 
% % SCATTER PLOTS: onshore sensors vs. offshore sensor data 
% figure; dockit(); subplot 141
% scatter(s124107.Hrms, s124108.Hrms,10,s124107.h_block);
% refline(1,0);
% colorbar
% colormap jet 
% subplot 142 
% scatter(s124107.Hrms, s124109.Hrms, 10,s124107.h_block);
% refline(1,0);
% colorbar
% colormap jet 
% subplot 143 
% scatter(s124107.Hrms, s41428.Hrms, 10,s124107.h_block);
% refline(1,0);
% colorbar
% colormap jet 
% subplot 144
% scatter(s124107.Hrms, s41429.Hrms, 10,s124107.h_block);
% refline(1,0);
% colorbar
% colormap jet 
% 
% % BREAKER INDEX TIME SERIES (Hrms/h) with offshore water depth series 
% c = parula(7);
% figure; dockit();
% subplot 211
% plot(s124107.t_block,s124107.Hrms./s124107.h_block,'color' ,c(1,:),'linewidth',2);
% hold on 
% plot(s124107.t_block,s124108.Hrms./s124108.h_block,'color' ,c(2,:),'linewidth',2);
% plot(s124107.t_block,s124109.Hrms./s124109.h_block,'color' ,c(3,:),'linewidth',2);
% plot(s124107.t_block,s41428.Hrms./s41428.h_block,'color' ,c(4,:),'linewidth',2);
% plot(s124107.t_block,s41429.Hrms./s41429.h_block,'color' ,c(5,:),'linewidth',2);
% legend('offshore','0m','5m','10m','20m');
% axis tight
% subplot 212
% plot(s124107.t_block,s124107.h_block);
% axis tight
% 
% % Hrms TIME SERIES 
% c = parula(7);
% figure; 
% plot(s124107.t_block-365,s124107.Hrms,'color' ,c(1,:),'linewidth',2);
% hold on 
% plot(s124108.t_block-365,s124108.Hrms,'color' ,c(2,:),'linewidth',2);
% 
% legend('offshore','marsh');
% xlabel('day of 2019');
% ylabel('Hrms [m]');
% axis tight

% 
% % Hm0/h breaker index 
% figure; subplot 211
% plot(s124107.t_block,(4*sqrt(s124109.m0))./s124109.h_block,'color' ,c(1,:),'linewidth',2);
% hold on 
% plot(s124107.t_block,(4*sqrt(s124108.m0))./s124108.h_block,'color' ,c(2,:),'linewidth',2);
% plot(s124107.t_block,(4*sqrt(s41428.m0))./s41428.h_block,'color' ,c(3,:),'linewidth',2);
% plot(s124107.t_block,(4*sqrt(s41429.m0))./s41429.h_block,'color' ,c(4,:),'linewidth',2);
% legend('0m','5m','10m','20m');
% axis tight
% 
% subplot 212
% plot(s124107.t_block,s124107.h_block);
% axis tight

% % Hrms error plots:
% figure;
% plot(s124107.t_block,s124107.Hrms,'k','linewidth',2)
% hold on 
% plot(s124107.t_block,s124107.Hrms_up,'k--')
% plot(s124107.t_block,s124107.Hrms_low,'k--')
% axis tight
% 
% figure;
% plot(s124107.t_block,s124107.Hrms_up-s124107.Hrms_low,'k','linewidth',2)
% axis tight



