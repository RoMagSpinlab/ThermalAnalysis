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

%% Channels Correlation
filename = uigetfile('*.txt');
% filename='xxx.txt';

% Read in data from compiled file
data = dlmread(filename,'\t',0,0);
% fileID = fopen('PathtoRMCheader.txt');
% header = textscan(fileID,'%s',length(data(1,:)),'Delimiter','\t');
% fclose(fileID);

% Channels
time = data(:,1)- data(1,1); % read time
begin_time = time(1);
end_time = time(end);
Toplid_therm = data(:,2:7); % top thermistors 0, 60, 120, 180, 240, 300
Botlid_therm = data(:,8:13); % bottom thermistors 0, 60, 120, 180, 240, 300
Toplid_therm7 = [data(:,2:7) data(:,2)];
Botlid_therm7 = [data(:,8:13) data(:,8)];
Int_therm = data(:,14:18); % internal thermistors
intcen_55 = data(:,14); % center, 55mm
int30_16 = data(:,15); % 30 degree, 2/3R, 16mm  
int90_105 = data(:,16); % 90 degree, 2/3R, 105mm 
int150_33 = data(:,17); % 150 degree, 2/3R, 33mm
int210_16 = data(:,18); % 210 degree, 2/3R, 16mm
PMU = data(:,19:20); % Power measurement  
%HallP = data(:,21:22); % 
ExpTank= data(:,23:25); % Expansion tank
SW_therm = data(:,26:51); % sidewall ensemble
SW_34 = data(:,26:31); % 3/4 height
SW_347 = [data(:,26:31) data(:,26)];
SW_12 = data(:,32:43); % 1/2 height
SW_mid13 = [data(:,32:43) data(:,32)];
SW_14 = data(:,44:49); % 1/4 height
SW_147 = [data(:,44:49) data(:,44)];
SW_vert0 = [data(:,50) data(:,44) data(:,32) data(:,26) data(:,51)]; % 1/8, 1/4, 1/2, 3/4, 7/8
vert_array = [data(:,8) SW_vert0 data(:,2)];
SW_loss = [data(:,32) data(:,51)]; % sidewall loss in and out
roomtemp = data(:,53); % room temp

%% Geometry and fixed solid properties
%Circum = 638.18E-3;
R_outer = 0.1016; % outer radius
r_sw = 0.003175; % sidewall thickness
r_insul = 1.0E-2; % insulation thickness
Rfluid = R_outer-r_sw; % fluid radius
A_fluid = pi*(R_outer-r_sw).^2; % horizontal surface area
A_sw = pi*R_outer^2-A_fluid; % SW horizontal surface area
H = 0.4; % tank height
Aspect = 2*(R_outer-r_sw)/H; % aspect ratio
H_insul =  0.19685; % insulation height 
H_cutop = 0.047; % top copper plate thickness
H_cubot = 0.012; % bottom copper plate thicknes
H_cu = 0.002; % Distance between top/bottom thermistors and the surfaces
k_insul = 0.015; % aerogel
% H_altop = 0.019; % 
k_cu = 390; % thermal conductivity of copper
k_sw = 16; % thermal conductivity of stainless steel
k_al = 167; % thermal conductivity of aluminum
g = 9.81; % gravitational acceleration
RPM = 9.401; % rotation speed in rpm
omega = RPM/60*2*pi; % in rad/s
Bfield = 0.0975; % applied magnetic field
freq_omega = RPM/60;
%RPM=20.4131;
%RPS=RPM/60;
%freq_rot=RPS;
%Omega=2*pi*freq_rot;
%Bfield=0.01;
%% Parameters Calculation

% % Sidewall
% SW_Midmean = mean(SW_Mid,2);
% % SW_Lowermean = mean(SW_Lower,2);
% % SW_Uppermean = mean(SW_Upper,2);
% % SW_Midmean = mean(SW_Mid,2);
% Cen_intlmean = mean(Cen_intl,2);
% diff_censw = mean(Cen_intlmean-SW_Midmean);

% Imposed Magnetic Field Strength (Tesla)
%Bfield = 200*mean(HallP(:,1))*10^(-4)

% Radial Induced Field (Gauss)
% HallPr = HallP(:,2);
% for j=1:length(HallP(:,2)) 
%     HallPr(j) = 10.*HallP(j,2)+0.85;
% end
% Vertical Field

% HallPz = HallP(:,1);
% for k=1:length(HallP(:,1)) 
%     HallPz(k) = 200.*HallP(k,1);
% end

% Raw Power
I_pmu = (squeeze(PMU(:,1)))*297; % I from Shunt
V_pmu = (squeeze(PMU(:,2)))*31;  % V from Voltage Divider
Praw = I_pmu.*V_pmu; % raw power output

% Mean Temperature
T_mean = (mean(Botlid_therm,2) + mean(Toplid_therm,2))/2;
meanT_mean = mean(T_mean)

% Radial Power Loss
SW_in = mean(SW_loss(:,1));
SW_out = mean(SW_loss(:,2));
Ploss_insul = (mean(SW_in-SW_out)*2*pi*H_insul*k_insul)/...
    log((R_outer+0.01)/(R_outer));
% Ploss_insul=mean(SW_innertemp-SW_outertemp)/(log((R_outer+0.019)/R_outer)/...
%     (2*pi*H_insul*k_insul)+ ...
%     log((R_outer+0.05)/(R_outer+0.019))/(2*pi*H_insul*0.07));

% Calculating Corrected Power to include losses
Pcorr = Praw-Ploss_insul;

% Delta T
% Correcting for top and bottom thermal blocks
DeltaT_lid_corr = 2*Pcorr.*H_cu./A_fluid./k_cu;

% Calculating DeltaT
DeltaT = mean(Botlid_therm,2) - mean(Toplid_therm,2) - DeltaT_lid_corr;
meanDeltaT = mean(DeltaT)

%% Calculating Power from corrected power and conducted power through the sidewall
P_SWcond = k_sw*DeltaT*A_sw/H;
P = mean(Pcorr-P_SWcond)

 %% Get Material Properties
[rho,alpha,Cp,k,sigma,eta,kappa,nu] = Mat_props(mean(T_mean),2);

%% Dimensionless Time Scales

kappa_av = mean(kappa);
tau_kappa = H^2/mean(kappa);
time_tau_kappa = time./tau_kappa;
freq_kappa = 1/tau_kappa;
tau_nu = H^2/mean(nu);
U_ff = sqrt(alpha*g*mean(DeltaT)*H);
U_tw = alpha*g*mean(DeltaT)/(2*omega);
freq_ff = U_ff/H;
freq_LSC = freq_ff/6;
freq_nu = 1/tau_nu;
tau_ff = H./U_ff;
tau_omega = 1/(2*omega);
time_ff = time./tau_ff;
l = sqrt(U_tw*H/(2*omega)); % local thermal wind scale 
                            % *** assume CIA balance rapid rotating ***
tau_tw = l/U_tw;


%% Dimensionless Numbers
Pr = nu/kappa
Pm = nu/eta
Ra = alpha*g*mean(DeltaT)*H^3/nu/kappa
Nu = mean(P.*H/A_fluid./k./DeltaT)
Ch = mean(sigma.*Bfield^2*H^2./rho./nu)
N = mean(Ch.*sqrt(Pr./Ra))
Rm = mean(U_ff.*H./eta)
Ra_f = alpha*g*P./A_fluid*H^4/nu/kappa/k
Re_ff = U_ff*H/nu
Ro_c = sqrt(alpha*g*mean(DeltaT)/(4*omega^2*H))
Re_tw = Ro_c^2*Re_ff
Ek = nu/(2*omega*H^2)
Elsa = Ch*Ek

%% Nu Error
staterrNu = std(P.*H/A_fluid./k./DeltaT,1);

%totalerr = sqrt(staterrNu^2+Nu^2*((1/P)^2+(0.01/98.6)^2+(2*1e-5/A_fluid)^2+...
%    (0.1/mean(DeltaT))^2))
totalerr = sqrt(Nu^2*((Ploss_insul/P)^2+(0.01/98.6)^2+(2*1e-5/A_fluid)^2+...
    (0.1/mean(DeltaT))^2))

%% Estimates

f_RC_osc = 4.7*(Ek/Pr)^(1/3)*freq_omega %RC overstable oscillatory mode frequency
b_ind = U_tw*Rfluid/eta*Bfield
Conv = Ra*(Nu-1)/Pr^2
%
f_MAC = (Conv*Ek)^(1/2)*nu/H^2 % Starchenko & Jones, 2002
f_CIA = Conv^(2/5)*Ek^(1/5)*nu/H^2 % Aubert et al., 2001
f_MC = (1+0.68*N^0.87)^(-1)*U_ff/H % Zurner et al., 2020
f_MC2 = 0.0324*Ra^(0.9542)*Ch^(-0.6461)/Pr * nu/H^2 % Yan et al., 2019
f_RMC_osc = 0.096*freq_omega % Horn & Aurnou., 2022 
f_RMC_wall = (8.9236*10^(-4))*freq_omega % Horn & Aurnou., 2022 
f_BZF = 0.03*Pr^(-4/3)*Ra*Ek^(5/3)*freq_omega % Zhang et al. 2021
f_RCwall = ((pi^2*sqrt(3)*sqrt(2+sqrt(3)))*Ek/Pr-290.6*Ek^(4/3)/Pr) % Zhang & Liao, 2009

%% title

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
        '$150^{\circ},\ 2/3 R,\ 33 \mathrm{mm}$',...
        '$210^{\circ},\ 2/3 R,\ 16 \mathrm{mm}$'},'interpreter',...
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
SW_mid_spec =(SW_mid13-meanT_mean)./mean(DeltaT);
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

%colorbar;

n = round(100/4);               %// number of colors
R_1 = linspace(0,0,n);  
G_1 = linspace(0.403,0.667,n);
B_1 = linspace(0.603,1,n);
R_2 = linspace(0,1,n);  
G_2 = linspace(0.667,1,n);
B_2 = linspace(1,1,n);
R_3 = linspace(1,1,n);  
G_3 = linspace(1,0,n);
B_3 = linspace(1,0.49,n);
R_4 = linspace(1,0.596,n);  
G_4 = linspace(0,0,n);
B_4 = linspace(0.49,0.298,n);
colormap( [R_1(:), G_1(:), B_1(:)
           R_2(:), G_2(:), B_2(:)
           R_3(:), G_3(:), B_3(:)
           R_4(:), G_4(:), B_4(:)] );  %// create colormap
%set(hFig, 'Position', [0 0 175 3200])
contourf(thermdiffspec2,anglespec2,SW_mid13,20,'LineStyle','none')

%t_colormap = text(480,1.7,'$(T-\overline T)/\Delta T$', ...
%    'interpreter','latex','fontsize',20,'Rotation',90);

title(titlestring1,'interpreter','latex','fontsize',20);
%text(450,1300,'$T(^\circ \rm C)$','interpreter','latex','FontSize',23, ...
%    'Rotation',270)
% ylim([0 10*tau_kappa/tau_ff]);
% ylim([0 5]);
%xticklabels({'0','60','120','180','240','300','360'});
ylim([0 0.1])
xticklabels({'0','30','60','90','120','150',...
    '180','210','240','270','300','330','360'});
xl = get(gca,'XTickLabel');
set(gca,'XTickLabel',xl,'FontName','Times','FontSize',18)

ylabel('$t/t_\kappa$','interpreter','latex','FontSize',18);
xlabel('$\phi(^{\circ})$','interpreter','latex','FontSize',18)
set(gca,'fontsize',18,'FontName','Times')

%caxis([-0.4 0.4])
colorbar
%colorbar('northoutside')
xticks(0:30:360)
hold off
%%
hFig4=figure(4);
set(hFig4, 'Position', [100 100 1200 400]);
set(gca,'OuterPosition',[0 0 0.98 1]);
colormap( [R_1(:), G_1(:), B_1(:)
           R_2(:), G_2(:), B_2(:)
           R_3(:), G_3(:), B_3(:)
           R_4(:), G_4(:), B_4(:)] );  %// create colormap
%set(hFig, 'Position', [0 0 175 3200])
heights=[0 1/8 1/4 1/2 3/4 7/8 1];
[heightspec,thermdiffspec]=meshgrid(time_tau_kappa,heights);

contourf(heightspec,thermdiffspec,vert_array',20,'LineStyle','none')
hold on
xlim([0 0.4]);
title(titlestring1,'interpreter','latex','fontsize',20);
ylabel('$H/0.4 \mathrm m$','interpreter','latex','FontSize',18);
xlabel('$t/t_\kappa$','interpreter','latex','FontSize',18)
set(gca,'fontsize',18,'FontName','Times')

%caxis([-0.4 0.4])
colorbar
hold off

% %% snapshot of the sidewall
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

% % 


%% Averaged FFT of the signals with Hann Window 
% FFT of the internal thermistors
% Int_therm: internal thermistors ensemble
% intcen_55: center, 55mm
% int30_16:  30 degree, 2/3R, 16mm  
% int90_105: 90 degree, 2/3R, 105mm 
% int150_33: 150 degree, 2/3R, 33mm
% int210_16: 210 degree, 2/3R, 16mm
%
maxf = zeros(5,1);
num=0;
% angle = [0 60 120 180 240 300];
Fs = 1/(time(3)-time(2));
titlestringFFT = {'Center, 55mm', '30-degree, 2/3R, 16mm', ...
    '90-degree, 2/3R, 105mm', '150-degree, 2/3R, 33mm', '210-degree, 2/3R, 16mm'};
% time = time(end);
% % 
hFig6 = figure(6);
set(hFig6, 'Position', [100 100 1200 1200])
set(gca,'OuterPosition',[0 0 1 0.95]);
% 
    
for i = 1:5 % count of thermistor
    num=num+1;
    IntFFT = Int_therm(:,i); % averaged SWFFT  
% Sampling frequency
% Window Length of FFT
% nfft = 2^nextpow2(time);     %Transform length
    %y = detrend((IntFFT-mean(Int_therm(:,i))));
    y = detrend((IntFFT-mean(T_mean)));
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
    %text(fk(MaxFreq)*1.5,pk+0.1,...
    %    ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
    %    'FontName','Times');
    
    % axis setting
    ylim([10^-5 10]);
    xlim([0.5 10^4]);
    xticks([1 10 10^2 10^3 10^4])
    yticks([10^-4 10^-2 1])
    
    % plot frequency predictions 
    l1 = loglog([f_MAC/freq_kappa f_MAC/freq_kappa], ylim, 'r--','linewidth',1);
    l2 = loglog([f_CIA/freq_kappa f_CIA/freq_kappa], ylim, 'b--','linewidth',1);
    l3 = loglog([f_MC2/freq_kappa f_MC2/freq_kappa], ylim, 'c--','linewidth',1);
    l4 = loglog([f_MC/freq_kappa f_MC/freq_kappa], ylim, 'g--','linewidth',1);
    l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
    l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
    %l7 = loglog([f_BZF/freq_kappa f_BZF/freq_kappa], ylim, 'b:','linewidth',2);
    l8 = loglog([f_RCwall/freq_kappa f_RCwall/freq_kappa], ylim, 'b:','linewidth',2);    
    xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
    ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
    legend([l1 l2 l3 l4 l5 l6 l8],...
        {'MAC','CIA', 'MC (Yan)', 'MC(Zurner)','RMC wall mode', ...
        'RMC oscillatory', 'RC wallmode' },...
        'fontsize',12,'location', 'northeast');
    set(gca,'fontsize',20,'FontName','Times');
    title(titlestringFFT(i),'interpreter','latex','fontsize',22);    
    

    
    
        
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
    y = detrend((IntFFT-mean(T_mean)));
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
%         ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
%         'FontName','Times');
    ylim([10^-5 10]);
    xlim([0.5 10^4]);
    xticks([1 10 10^2 10^3 10^4])
    yticks([10^-4 10^-2 1])
    
    % plot frequency preditions
    l1 = loglog([f_MAC/freq_kappa f_MAC/freq_kappa], ylim, 'r--','linewidth',1);
    l2 = loglog([f_CIA/freq_kappa f_CIA/freq_kappa], ylim, 'b--','linewidth',1);
    l3 = loglog([f_MC2/freq_kappa f_MC2/freq_kappa], ylim, 'c--','linewidth',1);
    l4 = loglog([f_MC/freq_kappa f_MC/freq_kappa], ylim, 'g--','linewidth',1);
    l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
    l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
    %l7 = loglog([f_BZF/freq_kappa f_BZF/freq_kappa], ylim, 'b:','linewidth',2);
    l8 = loglog([f_RCwall/freq_kappa f_RCwall/freq_kappa], ylim, 'b:','linewidth',2);
        
    xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
    ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
    legend([l1 l2 l3 l4 l5 l6 l8],...
        {'MAC','CIA', 'MC (Yan)', 'MC (Zurner)','RMC wallmode', ...
        'RMC oscillatory', 'RC wallmode' },'fontsize',12, 'location',...
        'northeast');
    set(gca,'fontsize',20,'FontName','Times');
    title('Sidewall, 7/8H, 0-degree','interpreter','latex','fontsize',22);    

    
    
%% PDF of temp variations

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
% % % %% Print
% % % fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', ...
% % %     'Bfield','Tmean','DeltaT','Power','Pr','Pm','Ra','Nu','Ch','N','fp');
% % % % 
% % % 
% % % fprintf('%.4f\t%.3f\t%.3f\t%.3f\t%.4f\t%.2e\t%.3e\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
% % %     Bfield,meanT_mean,meanDeltaT,P,Pr,Pm,Ra,Nu,Q,N,fp);
% % 
% % 
