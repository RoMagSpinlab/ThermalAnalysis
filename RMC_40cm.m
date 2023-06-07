%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMC Data Analysis for Gamma = 0.5
% 
% Yufan Xu
% RMC is developed from the previous lpMCX. 
% This routine is used to analysis the experimental results from RoMag, 
% interpreting the temperature and magnetic field signals via plots and 
% figures.
% 
% This program uses the latest material properties that appear in Aurnou&
% Bertin(2018). This program also requires Mat_props.m and 
% compile_convection_files.m.
% 
% Matlab Version: R2019b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mise en place
clear all
close all
clc
addpath('core/');

%% Channels Correlation
filename = uigetfile('*.txt');
% filename='xxx.txt';

% Read in data from compiled file
data = dlmread(filename,'\t',0,0);
% fileID = fopen('PathtoRMCheader.txt');

% Select time frame
t0 = 1:length(data(:,1));

%% Rotation and magnetic field 

RPM = 11.9460; %11.9460; %9.5; %12.587; %11.9048; % rotation speed in rpm %RPM=20.4131;9.401
Bfield = 0.0871; %0.0244; %0.0769;% 0.0863; % applied magnetic field

%% Macro
run('Romag_macro.m');
    % Channels Correlation (channel_map.m)
    % Geometry and material properties (constants.m)
    % Heat loss estimate (p_loss_corr.m)
    % Material properties (Mat_props.m)
    % Scales & dimensionless parameters (scale_param.m)
    % Frequency estimates (freq_predict.m)
    
%% titles

% RBC
titlestring0 = [' $\ Ra=$',num2str(Ra,3),'$,\ Nu=$',num2str(Nu,3),...
    '$,\ Ch = 0,\ Ek= \infty$'];

% RMC  
titlestring1 = [' $\Lambda=$',num2str(Elsa,3),'$,\ Ek=$',num2str(Ek,3),...
    '$,\ Ch=$',num2str(Ch,3),'$,\ Ra=$',num2str(Ra,3), ...
    '$,\ Nu=$',num2str(Nu,3),'$,\ Ro_c=$',num2str(Ro_c,3)];

% RC
titlestring2 = [' $\Lambda = 0,\ Ek=$',num2str(Ek,3),...
    '$,\ Ch=0,\ Ra=$',num2str(Ra,3), ...
    '$,\ Nu=$',num2str(Nu,3),'$,\ Ro_c=$',num2str(Ro_c,3)];

% MC
titlestring3 = [' $\Lambda = \infty,\ Ek= \infty,\ Ch=$', num2str(Ch,3),...
    '$,\ Ra=$',num2str(Ra,3),'$,\ Nu=$',num2str(Nu,3),'$,\ N_c=$',num2str(N,3)];

%% Top Mid Bottom Lid Temperature Raw Data
hFig1 = figure(1);
    set(hFig1, 'Position', [100 10 1400 1200])
    set(gca,'OuterPosition',[0 0 1 0.95]);
    timelim = [min(time_tau_kappa),max(time_tau_kappa)*1.1];

subplot(5,1,1)
    plot(time_tau_kappa,Toplid_therm,'LineWidth',1); 
    hold on
    grid on
    xlim (timelim);
    %ylim([35.4 35.6]); % Ra = 3e6 RMC
    %ylim([35.0 36.0]); % Ra = 3e7 RMC
    %ylim ([35 36]);
    ft = legend({'$ 0^{\circ}$','$ 60^{\circ}$',...
        '$ 120^{\circ} $','$ 180^{\circ}$',...
        '$ 240^{\circ}$','$ 300^{\circ}$'},'interpreter',...
        'latex','fontsize',15);    
    xlabel('$t_{\kappa}$','interpreter','latex','fontsize',18) % x-axis label
    ylabel('$T_{H}(\mathrm{^\circ C})$','interpreter','latex','fontsize',18) % y-axis label
    set(ft,'Location','East')
    set(gca,'fontsize',18,'FontName','Times')


    title(titlestring1,'interpreter','latex',...
        'fontsize',20);
     
subplot(5,1,2)
    plot(time_tau_kappa,SW_34,'LineWidth',1); 
    grid on
    hold on
    xlim (timelim);
    %ylim([35.0 35.2]); % Ra = 3e6 RMC
    %ylim([35.3 36.3]); % Ra = 3e7 RMC
    %ylim ([35.3 36.3]);
    %yticks([35.3 35.8 36.3]);
    f34 = legend({'$0^{\circ}$','$60^{\circ}$',...
        '$120^{\circ} $','$180^{\circ}$',...
        '$240^{\circ}$','$300^{\circ}$'},...
        'interpreter','latex','fontsize',15);
    xlabel('$t_{\kappa}$','interpreter','latex','fontsize',18) % x-axis label
    ylabel('$T_{3/4H}(\mathrm{^\circ C})$','interpreter','latex','fontsize',18) % y-axis label
    set(f34,'Location','East')
    set(gca,'fontsize',18, 'FontName','Times')
    
subplot(5,1,3)
    plot(time_tau_kappa,SW_12,'LineWidth',1); 
    grid on
    hold on
    xlim (timelim);
    %ylim([35.0 35.2]); % Ra = 3e6 RMC
    %ylim([35.5 36.5]); % Ra = 3e7 RMC
    %ylim ([35.5 36.5]);
    f12 = legend({'$0^{\circ}$','$30^{\circ}$','$60^{\circ}$',...
        '$90^{\circ}$','$120^{\circ} $','$150^{\circ}$','$180^{\circ}$',...
        '$210^{\circ}$','$240^{\circ}$','$270^{\circ}$','$300^{\circ}$',...
        '$330^{\circ}$'},'interpreter','latex','fontsize',15);
    xlabel('$t_{\kappa}$','interpreter','latex','fontsize',18) % x-axis label
    ylabel('$T_{1/2H}(\mathrm{^\circ C})$','interpreter','latex','fontsize',18) % y-axis label
    set(f12,'Location','East')
    set(gca,'fontsize',18, 'FontName','Times')

subplot(5,1,4)
    plot(time_tau_kappa,SW_14,'LineWidth',1); 
    grid on
    hold on
    xlim(timelim);
    %ylim([35.0 35.2]); % Ra = 3e6 RMC
    %ylim([35.8 36.8]); % Ra = 3e7 RMC
    %ylim ([35.8 36.8]);
    %yticks([35.8 36.3 36.8]);
    f14 = legend({'$0^{\circ}$','$60^{\circ}$',...
        '$120^{\circ} $','$180^{\circ}$',...
        '$240^{\circ}$','$300^{\circ}$'},...
        'interpreter','latex','fontsize',15);
    xlabel('$t_{\kappa}$','interpreter','latex','fontsize',18) % x-axis label
    ylabel('$T_{1/4H}(\mathrm{^\circ C})$','interpreter','latex','fontsize',18) % y-axis label
    set(f14,'Location','East')
    set(gca,'fontsize',18, 'FontName','Times')

subplot(5,1,5)
    plot(time_tau_kappa,Botlid_therm,'LineWidth',1); 
    hold on
    grid on
    xlim (timelim);
    %ylim([35.6 35.8]); % Ra = 3e6 RMC
    %ylim([36.8 37.8]); % Ra = 3e7 RMC
    %ylim ([36.8 37.8]);
    %yticks([36.8 37.3 37.8]);
    fb = legend({'$0^{\circ}$','$60^{\circ}$',...
        '$120^{\circ} $','$180^{\circ}$',...
        '$240^{\circ}$','$300^{\circ}$'},'interpreter',...
        'latex','fontsize',15);    
    xlabel('$t_{\kappa}$','interpreter','latex','fontsize',18) % x-axis label
    ylabel('$T_{0}(\mathrm{^\circ C})$','interpreter','latex','fontsize',18) % y-axis label
    set(fb,'Location','East')
    set(gca,'fontsize',18, 'FontName','Times')
hold off

% saveas(hFig1,[num2str(P,'%.0f') 'W_timeseries'],'epsc');   
% hold on

%% Interiors
% 
hFig2 = figure(2);
    set(hFig2, 'Position', [100 10 1000 600])
    set(gca,'OuterPosition',[0 0 1 0.95]);
plot(time_tau_kappa,Int_therm,'LineWidth',1);
    hold on 
    grid on

plot(time_tau_kappa,Int_therm(:,1),'k-','LineWidth',1);
% Int_therm  % internal thermistors
% intcen_55  % center, 55mm
% int30_16  % 30 degree, 2/3R, 16mm  
% int90_105  % 90 degree, 2/3R, 105mm 
% int150_33  % 150 degree, 2/3R, 33mm
% int210_16 % 210 degree, 2/3R, 16mm

    xlim (timelim);
    %ylim ([35 37]);
    fi = legend({'center, $55 \mathrm{mm}$',...
        '$30^{\circ},\ 2/3 R,\ 16 \mathrm{mm}$',...
        '$90^{\circ},\ 2/3 R,\ 105 \mathrm{mm} $',...
        '$210^{\circ},\ 2/3 R,\ 33 \mathrm{mm}$',...
        '$270^{\circ},\ 2/3 R,\ 16 \mathrm{mm}$'},'interpreter',...
        'latex','fontsize',15);    
    xlabel('$t_{\kappa}$','interpreter','latex','fontsize',18) % x-axis label
    ylabel('$T(\mathrm{^\circ C})$','interpreter','latex','fontsize',18) % y-axis label
    set(fb,'Location','East')
    set(gca,'fontsize',18, 'FontName','Times')

    

%% SideWall Array Plots Mid-Plane
angle1=0:60:360;
angle2=0:30:360;
columnMean = mean(SW_mid13,1);

% normalize the T field
%SW_mid_spec =(SW_mid13-meanT_mean)./mean(DeltaT);
SW_mid_spec =(SW_mid13-T_mean)./(DeltaT);
%SW_34_spec =(SW_347-meanT_mean)./mean(DeltaT);
SW_34_spec =(SW_347-T_mean)./(DeltaT);
%SW_14_spec =(SW_147-meanT_mean)./mean(DeltaT);
SW_14_spec =(SW_147-T_mean)./(DeltaT);
%vert_spec = (vert_array'-meanT_mean)./mean(DeltaT);
vert_spec = (vert_array-T_mean)./(DeltaT);
% normalize the T field

% SW_mid13;
% Toplid_therm7;
% Botlid_therm7;
% SW_347;
% SW_147;
[thermdiffspec1,anglespec1]=meshgrid(angle1,time_tau_kappa);
[thermdiffspec2,anglespec2]=meshgrid(angle2,time_tau_kappa);
[thermdiffspecf,anglespecf]=meshgrid(angle1,time_ff);

%%
hFig3=figure(3);
set(hFig3, 'Position', [100 100 800 1200]);
set(gca,'OuterPosition',[0 0 0.98 1]);
colormap(cmap_YX);

contourf(thermdiffspec2,anglespec2,SW_mid13,20,'LineStyle','none');
%contourf(thermdiffspecf,anglespecf,SW_147,20,'LineStyle','none');
%contourf(thermdiffspecf,anglespecf,SW_347,20,'LineStyle','none');

%t_colormap = text(480,1.7,'$(T-\overline T)/\Delta T$', ...
%    'interpreter','latex','fontsize',20,'Rotation',90);

%text(450,1300,'$T(^\circ \rm C)$','interpreter','latex','FontSize',23, ...
%    'Rotation',270)
% ylim([0 10*tau_kappa/tau_ff]);
% ylim([0 5]);
%xticklabels({'0','60','120','180','240','300','360'});
%ylim([0 0.3])
hold on 

title(titlestring1,'interpreter','latex','fontsize',20);

xticklabels({'0','30','60','90','120','150',...
    '180','210','240','270','300','330','360'});
xl = get(gca,'XTickLabel');
set(gca,'XTickLabel',xl,'FontName','Times','FontSize',18)

ylabel('$t/t_\kappa$','interpreter','latex','FontSize',18);
xlabel('$\phi(^{\circ})$','interpreter','latex','FontSize',18)
set(gca,'fontsize',18,'FontName','Times')

%caxis([55 78])
colorbar
%colorbar('northoutside')
xticks(0:30:360)
hold off

%% vertical hovmoller
hFig4=figure(4);
set(hFig4, 'Position', [100 100 1200 400]);
set(gca,'OuterPosition',[0 0 0.98 1]);

colormap(cmap_SH); 
heights=[1/8 1/4 1/2 3/4 7/8];
%[heightspec,thermdiffspec]=meshgrid(time_ff,heights);
[heightspec,thermdiffspec]=meshgrid(time_tau_kappa,heights);
contourf(heightspec,thermdiffspec,vert_array',20,'LineStyle','none')

hold on
%xlim([0 0.4]);
title(titlestring1,'interpreter','latex','fontsize',20);

ylabel('$H/0.4 \mathrm m$','interpreter','latex','FontSize',18);
xlabel('$t/t_{\kappa}$','interpreter','latex','FontSize',18)
set(gca,'fontsize',18,'FontName','Times')

%caxis([-0.4 0.4])
colorbar
%xlim([0 0.3]);
hold off
% % 

% % 


%% Averaged FFT of the signals with Hann Window 
% FFT of the internal thermistors
% Int_therm: internal thermistors ensemble
% intcen_55: center, 55mm
% int30_16:  30 degree, 2/3R, 16mm  
% int90_105: 90 degree, 2/3R, 105mm 
% int210_33: 210 degree, 2/3R, 33mm
% int270_16: 270 degree, 2/3R, 16mm
%
maxf = zeros(5,1);
num=0;
% angle = [0 60 120 180 240 300];
Fs = 1/(time(3)-time(2));
titlestringFFT = {'Center, 55mm', '30-degree, 2/3R, 16mm', ...
    '90-degree, 2/3R, 105mm', '180-degree, 2/3R, 33mm', '270-degree, 2/3R, 16mm'};
% time = time(end);
% % 
hFig6 = figure(6);
set(hFig6, 'Position', [100 100 1200 1200])
set(gca,'OuterPosition',[0 0 1 0.95]);
% 
% Internal thermistor FFTs
%
for i = 1:5 % count of thermistor
    num=num+1;
    IntFFT = Int_therm(:,i); % averaged SWFFT  
% Sampling frequency
% Window Length of FFT
% nfft = 2^nextpow2(time);     %Transform length
    y = detrend(IntFFT-mean(Int_therm(:,i)));
    %y = detrend((IntFFT-mean(T_mean)));
% y = detrend(SWFFT);
    L = length(y);
% freq = Fs/2*linspace(0,1,nfft/2+1); 
    freq = Fs/2*linspace(0,1,L/2+1); 
% FFT_SW_1=abs(fft(detrend((SWFFT-T_mean)./DeltaT),n)/n);
% Hanning Window
    fk=freq/freq_kappa;
    y_HannWnd = y.*hann(L);
    Ydft_HannWnd = fft(y_HannWnd,L)/L;
    FFT_int_1=abs(Ydft_HannWnd);
    L = L-1;
%  FFT_SW_1=abs(fft(detrend(SWFFT),L)/L);
    FFT_int_2= FFT_int_1(1:L/2+1);
    FFT_int_2(2:end-1) = 2*FFT_int_2(2:end-1);

subplot(3,2,i)

    loglog(fk,2*FFT_int_2,'k','LineWidth',1); % Multiply by 2 because 
%                                                % Hanning wnd Amplitude 
%                                                % Correction Factor = 2  
    grid on
    [pk,MaxFreq] = findpeaks(2*FFT_int_2,'NPeaks',1,'SortStr','descend');%peakfreq
    hold on
    peak=loglog(fk(MaxFreq),pk,'or');
    maxf(num) = fk(MaxFreq);
    text(fk(MaxFreq)*1.5,pk+0.1,...
       ['(' num2str(fk(MaxFreq),3) ' , ' num2str(pk,3) ')' ],'fontsize',20,...
       'FontName','Times');

   
    % axis setting
    ylim([10^-5 10]);
    xlim([0.5 10^4]);
    xticks([1 10 10^2 10^3 10^4])
    yticks([10^-4 10^-2 1])
    xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
    ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
    set(gca,'fontsize',20,'FontName','Times');
    title(titlestringFFT(i),'interpreter','latex','fontsize',22);    
    
     %l1 = loglog([f_MAC/freq_kappa f_MAC/freq_kappa], ylim, 'r--','linewidth',1);
%    l2 = loglog([f_CIA/freq_kappa f_CIA/freq_kappa], ylim, 'b--','linewidth',1);
%    l3 = loglog([f_MC2/freq_kappa f_MC2/freq_kappa], ylim, 'c--','linewidth',1);
%    l4 = loglog([f_MC/freq_kappa f_MC/freq_kappa], ylim, 'g--','linewidth',1);
     l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
     l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
%    l7 = loglog([f_RCwall/freq_kappa f_RCwall/freq_kappa], ylim, 'b:','linewidth',2);
    legend([l5,l6],{'RMC wall mode','RMC oscillatory'},'interpreter','latex',...
        'fontsize',16,'location','SW'); 

end   

% 
%fp=mean(maxf)
% 
% lasttime=time_tau_kappa(end)
% % saveas(hFig9,[num2str(P,'%.0f') 'W_FFT'],'epsc');   
% % 
 
subplot(3,2,6)
    %IntFFT = SW_34(:,2);
    IntFFT = SW_vert0(:,5);
    %IntFFT = ExpTank(:,3);
    %y = detrend(IntFFT);
    y = detrend(IntFFT-mean(SW_vert0(:,5)));
    %y = detrend((IntFFT-mean(T_mean)));
    %y = detrend((IntFFT-meanT_mean)./mean(DeltaT));
    L = length(y);
    freq = Fs/2*linspace(0,1,L/2+1); 
% Hanning Window
    fk=freq/freq_kappa;
    %fk = freq/freq_omega;
    y_HannWnd = y.*hann(L);
    Ydft_HannWnd = fft(y_HannWnd,L)/L;
    FFT_int_1=abs(Ydft_HannWnd);
    L = L-1;
    FFT_int_2= FFT_int_1(1:L/2+1);
    FFT_int_2(2:end-1) = 2*FFT_int_2(2:end-1);
% plot
    loglog(fk,2*FFT_int_2,'k','LineWidth',1); % Multiply by 2 because 
%                                                % Hanning wnd Amplitude 
%                                                % Correction Factor = 2  
    grid on
    [pk,MaxFreq] = findpeaks(2*FFT_int_2,'NPeaks',1,'SortStr','descend');%peakfreq
    hold on
    peak=loglog(fk(MaxFreq),pk,'or');
    maxf(num) = fk(MaxFreq);
%     text(fk(MaxFreq)*1.5, pk+0.1,...
%         ['(' num2str(fk(MaxFreq),'%.3f') ' , ' num2str(pk,3) ')' ],'fontsize',20,...
%         'FontName','Times');
    text(fk(MaxFreq)*1.5, pk+0.1, ['(' num2str(fk(MaxFreq),'%.3f'),')' ],...
        'fontsize',20,'FontName','Times');
    ylim([10^-5 10]);
    xlim([0.5 10^4]);
    xticks([1 10 10^2 10^3 10^4])
    yticks([10^-4 10^-2 1])
    xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
    ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
    set(gca,'fontsize',20,'FontName','Times');
    title('Sidewall, 7/8H, 0-degree','interpreter','latex','fontsize',22);    
     l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
     l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
     legend([l5,l6],{'RMC wall mode','RMC oscillatory'},'interpreter','latex',...
         'fontsize',16,'location','SW'); 
     
     
figure
Nu_t = P.*H/A_fluid./k./DeltaT;
plot(time_tau_kappa, Nu_t, 'k-', 'linewidth',2);
xlim([0 0.5]);
     
%% Print
fprintf('%0s\t %5s\t %5s\t %5s\t %5s\t %5s\t %5s\t %10s\t %8s\t %5s\t %5s\t %5s\n', ...
    'RPM','Bfield','T_mean','DeltaT','Power','Pr','Ra','Ek','Ch','Nu','Ro_c','N');
% 
% % % 
fprintf('%.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3e\t %.3e\t %.3e\t %.3f\t %.3f\t %.3f\n', ...
    RPM,Bfield,meanT_mean,meanDeltaT,P,Pr,Ra,Ek,Ch,Nu,Ro_c,N);

    
    
% hFig7 = figure(7);
% set(hFig7, 'Position', [100 100 800 600])
% set(gca,'OuterPosition',[0 0 1 0.95]);
%     
% histfit(Int_therm(:,1)-mean(Int_therm(:,1)),100,'normal')
% hold on 
% %cent_RMC_Ra9e8 = Int_therm(:,1);
% %cent_RC_Ra8e8 = Int_therm(:,1);
% %cent_MC_Ra8e8 = Int_therm(:,1);
% %R23_10cm_RMC_Ra9e8 = Int_therm(:,3);
% %R23_10cm_RC_Ra8e8 = Int_therm(:,3);
% %R23_10cm_MC_Ra8e8 = Int_therm(:,3);
% %save('cent_RMC_Ra9e8.mat','cent_RMC_Ra9e8')
% %save('cent_RC_Ra8e8.mat','cent_RC_Ra8e8')
% %save('cent_MC_Ra8e8.mat','cent_MC_Ra8e8')
% %save('R23_10cm_RMC_Ra9e8.mat','R23_10cm_RMC_Ra9e8')
% %save('R23_10cm_RC_Ra8e8.mat','R23_10cm_RC_Ra8e8')
% %save('R23_10cm_MC_Ra8e8.mat','R23_10cm_MC_Ra8e8')


%%% %% snapshot of the sidewall
% hFig5=figure(5);
% set(hFig5, 'Position', [100 100 400*pi 800]);
% set(gca,'OuterPosition',[0 0 0.98 1]);
% colormap( [R_1(:), G_1(:), B_1(:)
%            R_2(:), G_2(:), B_2(:)
%            R_3(:), G_3(:), B_3(:)
%            R_4(:), G_4(:), B_4(:)] );
%        
% % make a fancy movie!       
% tstep = 10;
% set(gca,'nextplot','replacechildren'); 
% v = VideoWriter('RMC_steady_Ra3e6.avi');
% %v = VideoWriter('RMC_wallmode_Ra3e7.avi');
% %v = VideoWriter('RMC_wallmode_Ra7e7.avi');
% %v = VideoWriter('RMC_magnetostrophic_9e8.avi');
% %v = VideoWriter('RMC_mix_wm_ms_3e8.avi');
% %v = VideoWriter('RC_geostrophic_Ra8e8.avi');
% %v = VideoWriter('MC_Ra_8e8.avi');
% open(v);
% tic
% for ti = 1:tstep:time(end)
% %ti = 1000; % pick a time,
% Tinfo = [Botlid_therm7(ti,:) SW_vert0(ti,1) SW_vert0(ti,1) SW_147(ti,:)...
%     SW_mid13(ti,:) SW_347(ti,:) SW_vert0(ti,5) SW_vert0(ti,5) Toplid_therm7(ti,:)];
% phi = [0:60:360 0 360 0:60:360 0:30:360 0:60:360 0 360 0:60:360];
% Hinfo = [zeros(1,7) 1/8 1/8 1/4*ones(1,7) 1/2*ones(1,13) 3/4*ones(1,7) 7/8 7/8 ones(1,7)];
% 
% xlin = linspace(0, 360, 360);
% ylin = linspace(min(Hinfo), max(Hinfo), 100);
% 
% [X,Y] = meshgrid(xlin, ylin);
% 
% % Z = griddata(x,y,z,X,Y,'natural');
% % Z = griddata(x,y,z,X,Y,'cubic');
% Z = griddata(phi,Hinfo,Tinfo,X,Y,'natural');
% contourf(X,Y,Z,20,'LineStyle','none')
% caxis([min(SW_therm(:)) max(SW_therm(:))])
% 
% %text(350, 1, ['t = ' num2str(ti*tstep/5) 'sec'],'fontsize',20)
% title([titlestring3 '$,\ t = $' num2str(ti) 'sec'],...
%     'interpreter','latex','fontsize',20);
% ylabel('$H/0.4 \mathrm m$','interpreter','latex','FontSize',18);
% xlabel('$\phi (\mathrm{^{\circ}})$','interpreter','latex','FontSize',18)
% set(gca,'fontsize',18,'FontName','Times')
% xticks(0:30:360)
% yticks([0 1/8 1/4 1/2 3/4 7/8 1])
% 
% colorbar
% 
%    frame = getframe(gcf);
%    writeVideo(v,frame);
% end
% toc
