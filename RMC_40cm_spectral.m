%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMC Data Analysis for Gamma = 0.5 -- Spectral Analysis
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

totalerr = sqrt(Nu^2*((Ploss_insul/P)^2+(0.01/98.6)^2+(2*1e-5/A_fluid)^2+...
    (0.1/mean(DeltaT))^2))

%% Estimates

Ra_w = 1.849e7 %RMC wall mode onset
Ra_osc = 8.936e7 % RMC oscillatory mode onset
Ra_ms = 5.700e7 % RMC magnetostrophic mode

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




%% FFT of the signals with Hann Window 
% FFT of the internal thermistors
% Int_therm: internal thermistors ensemble
% intcen_55: center, 55mm
% int30_16:  30 degree, 2/3R, 16mm  
% int90_105: 90 degree, 2/3R, 105mm 
% int150_33: 150 degree, 2/3R, 33mm
% int210_16: 210 degree, 2/3R, 16mm

maxf = zeros(5,1);
num=0;
% angle = [0 60 120 180 240 300];
Fs = 1/(time(3)-time(2));
titlestringFFT = {'Center, 55mm', '30-degree, 2/3R, 16mm', ...
    '90-degree, 2/3R, 105mm', '150-degree, 2/3R, 33mm', '210-degree, 2/3R, 16mm'};
% time = time(end);
% % 
hFig1 = figure(1);
set(hFig1, 'Position', [100 100 1200 1200])
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

    
    % plot frequency predictions 
    l1 = loglog([f_MAC/freq_kappa f_MAC/freq_kappa], [10^-5 10], 'r-','linewidth',1);
    hold on
    grid on
    ylim([10^-5 10]);
    xlim([0.5 10^4]);
    xticks([1 10 10^2 10^3 10^4])
    yticks([10^-4 10^-2 1])
    l2 = loglog([f_CIA/freq_kappa f_CIA/freq_kappa], ylim, 'b-','linewidth',1);
    l3 = loglog([f_MC2/freq_kappa f_MC2/freq_kappa], ylim, 'c--','linewidth',1);
    l4 = loglog([f_MC/freq_kappa f_MC/freq_kappa], ylim, 'g--','linewidth',1);
    l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
    l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
    l7 = loglog([f_RC_osc/freq_kappa f_RC_osc/freq_kappa], ylim, 'b:','linewidth',2);
    l8 = loglog([f_RCwall/freq_kappa f_RCwall/freq_kappa], ylim, 'b-.','linewidth',1);    
    l9 = loglog([f_BZF/freq_kappa f_BZF/freq_kappa], ylim, 'y:','linewidth',2);
    loglog(fk,2*FFT_int_2,'k','LineWidth',1); % Multiply by 2 because 
%                                                % Hanning wnd Amplitude 
%                                                % Correction Factor = 2  
%    grid on
    [pk,MaxFreq] = findpeaks(2*FFT_int_2,'NPeaks',1,'SortStr','descend');%peakfreq
%    hold on
    peak=loglog(fk(MaxFreq),pk,'or');
    maxf(num) = fk(MaxFreq);
    %text(fk(MaxFreq)*1.5,pk+0.1,...
    %    ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
    %    'FontName','Times');
    
    % axis setting

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
   
    % plot frequency preditions
    l1 = loglog([f_MAC/freq_kappa f_MAC/freq_kappa], [10^-5 10], 'r-','linewidth',1);
    hold on
    grid on 
    ylim([10^-5 10]);
    xlim([0.5 10^4]);
    xticks([1 10 10^2 10^3 10^4])
    yticks([10^-4 10^-2 1]) 
    l2 = loglog([f_CIA/freq_kappa f_CIA/freq_kappa], ylim, 'b-','linewidth',1);
    l3 = loglog([f_MC2/freq_kappa f_MC2/freq_kappa], ylim, 'c--','linewidth',1);
    l4 = loglog([f_MC/freq_kappa f_MC/freq_kappa], ylim, 'g--','linewidth',1);
    l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
    l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
    l7 = loglog([f_RC_osc/freq_kappa f_RC_osc/freq_kappa], ylim, 'b:','linewidth',2);
    %l7 = loglog([f_BZF/freq_kappa f_BZF/freq_kappa], ylim, 'b:','linewidth',2);
    l8 = loglog([f_RCwall/freq_kappa f_RCwall/freq_kappa], ylim, 'b-.','linewidth',1);   
    
    loglog(fk,2*FFT_int_2,'k','LineWidth',1); % Multiply by 2 because 
%                                                % Hanning wnd Amplitude 
%                                                % Correction Factor = 2  
  
    [pk,MaxFreq] = findpeaks(2*FFT_int_2,'NPeaks',1,'SortStr','descend');%peakfreq
  
    peak=loglog(fk(MaxFreq),pk,'or');
    maxf(num) = fk(MaxFreq);
%     text(fk(MaxFreq)*1.5, pk+0.1,...
%         ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
%         'FontName','Times');

        
    xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
    ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
    legend([l1 l2 l3 l4 l5 l6 l8],...
        {'MAC','CIA', 'MC (Yan)', 'MC (Zurner)','RMC wallmode', ...
        'RMC oscillatory', 'RC wallmode' },'fontsize',12, 'location',...
        'northeast');
    set(gca,'fontsize',20,'FontName','Times');
    title('Sidewall, 7/8H, 0-degree','interpreter','latex','fontsize',22);    

    
    %% midplane sidewall 1/2H 
hFig2 = figure(2);
set(hFig2, 'Position', [100 200 800 500])
set(gca,'OuterPosition',[0 0 1 0.95]);
    
    IntFFT = SW_12(:,6); % 150degree
    %IntFFT = SW_vert0(:,5);
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
   
    % plot frequency preditions
    l1 = loglog([f_MAC/freq_kappa f_MAC/freq_kappa], [10^-5 10], 'r-','linewidth',1);
    hold on
    grid on 
    ylim([10^-5 10]);
    xlim([0.5 10^4]);
    xticks([1 10 10^2 10^3 10^4])
    yticks([10^-4 10^-2 1]) 
    l2 = loglog([f_CIA/freq_kappa f_CIA/freq_kappa], ylim, 'b-','linewidth',1);
    %l3 = loglog([f_MC2/freq_kappa f_MC2/freq_kappa], ylim, 'c--','linewidth',1);
    %l4 = loglog([f_MC/freq_kappa f_MC/freq_kappa], ylim, 'g--','linewidth',1);
    l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
    l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
    l7 = loglog([f_RC_osc/freq_kappa f_RC_osc/freq_kappa], ylim, 'b:','linewidth',2);
    %l7 = loglog([f_BZF/freq_kappa f_BZF/freq_kappa], ylim, 'b:','linewidth',2);
    l8 = loglog([f_RCwall/freq_kappa f_RCwall/freq_kappa], ylim, 'b-.','linewidth',1);   
%     for i = 2:4
%     l9 = loglog([i*f_RCwall/freq_kappa i*f_RCwall/freq_kappa], [0.1 10], 'k--','linewidth',1);
%     end
    
    
    loglog(fk,2*FFT_int_2,'k','LineWidth',1); % Multiply by 2 because 
%                                                % Hanning wnd Amplitude 
%                                                % Correction Factor = 2  
  
    [pk,MaxFreq] = findpeaks(2*FFT_int_2,'NPeaks',1,'SortStr','descend');%peakfreq
  
    peak=loglog(fk(MaxFreq),pk,'or');
    maxf(num) = fk(MaxFreq);
%     text(fk(MaxFreq)*1.5, pk+0.1,...
%         ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
%         'FontName','Times');

        
    xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
    ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
    legend([l5 l1 l8 l2 l6 l7],...
        {'RMC wallmode','MAC','RC wallmode','CIA', ...
        'RMC oscillatory', 'RC oscillatory' },...
        'fontsize',15,...
        'location','northeast'); % and 'resonance'
    set(gca,'fontsize',20,'FontName','Times');
    title(['Sidewall, 1/2H, 150-degree, $Ra/Ra_w =\ $', num2str(Ra/Ra_w,3), ...
        ', $Ro_c =\ $', num2str(Ro_c,3), ', $\Lambda = 0 $'],...
        'interpreter','latex','fontsize',22);   
%    RCFFT15W = 2*FFT_int_2;
%   RMCFFT15W = 2*FFT_int_2;
%save('RCFFT15W.mat','RCFFT15W')

% hFig3 = figure(3);
% set(hFig3, 'Position', [100 100 800 500])
% set(gca,'OuterPosition',[0 0 1 0.95]);
% 
%  l1 = loglog([f_MAC/freq_kappa f_MAC/freq_kappa], [10^-5 10], 'r-','linewidth',1);
%     hold on
%     grid on 
%     ylim([10^-5 10]);
%     xlim([0.5 10^4]);
%     xticks([1 10 10^2 10^3 10^4])
%     yticks([10^-4 10^-2 1]) 
%     l2 = loglog([f_CIA/freq_kappa f_CIA/freq_kappa], ylim, 'b-','linewidth',1);
%     %l3 = loglog([f_MC2/freq_kappa f_MC2/freq_kappa], ylim, 'c--','linewidth',1);
%     %l4 = loglog([f_MC/freq_kappa f_MC/freq_kappa], ylim, 'g--','linewidth',1);
%     l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
%     l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
%     l7 = loglog([f_RC_osc/freq_kappa f_RC_osc/freq_kappa], ylim, 'b:','linewidth',2);
%     %l7 = loglog([f_BZF/freq_kappa f_BZF/freq_kappa], ylim, 'b:','linewidth',2);
%     l8 = loglog([f_RCwall/freq_kappa f_RCwall/freq_kappa], ylim, 'b-.','linewidth',1);   
%     for j = 2:4
%     l9 = loglog([j*f_RCwall/freq_kappa j*f_RCwall/freq_kappa], [0.1 10], 'k--','linewidth',1);
%     end
%     
%     
%     loglog(fk,RCFFT15W,'LineWidth',1); % Multiply by 2 because 
% %                                                % Hanning wnd Amplitude 
% %                                                % Correction Factor = 2  
%   
%     [pk,MaxFreq] = findpeaks(2*RCFFT15W,'NPeaks',1,'SortStr','descend');%peakfreq
%   
%     peak=loglog(fk(MaxFreq),pk,'or');
%     maxf(num) = fk(MaxFreq);
% %     text(fk(MaxFreq)*1.5, pk+0.1,...
% %         ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
% %         'FontName','Times');
% %     loglog(fk,RMCFFT50W,'LineWidth',1);
%         
%     xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
%     ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
%     legend([l5 l1 l8 l9 l2 l6 l7],...
%         {'RMC wallmode','MAC','RC wallmode', 'resonance','CIA', ...
%         'RMC oscillatory', 'RC oscillatory' },...
%         'fontsize',15,...
%         'location','northeast');
%     set(gca,'fontsize',20,'FontName','Times');
%     title(['Radial, 150-degree, $Ra/Ra_w =\ $', num2str(Ra/Ra_w,3), ...
%         ', $Ro_c =\ $', num2str(Ro_c,3), ', $\Lambda =\ $', num2str(Elsa,3)],...
%         'interpreter','latex','fontsize',22);    
 
    


% %%
% hFig3 = figure(3);
% set(hFig3, 'Position', [100 100 800 500])
% set(gca,'OuterPosition',[0 0 1 0.95]);
% 
%     IntFFT1 =  SW_12(:,6); 
%     IntFFT2 =  Int_therm(:,4);  
%     IntFFT3 =  Int_therm(:,1); 
%     IntFFT = [IntFFT1 IntFFT2 IntFFT3];
%     %IntFFT = SW_vert0(:,5);
%     %IntFFT = ExpTank(:,3);
%     %y = detrend(IntFFT);
% for i = 1:3
%     y = detrend((IntFFT(:,i)-mean(T_mean)));
% 
%     %y = detrend((IntFFT-meanT_mean)./mean(DeltaT));
%     L = length(y);
%     freq = Fs/2*linspace(0,1,L/2+1); 
% % Hanning Window
%     fk=freq/freq_kappa;
%     %fk = freq/freq_omega;
%     y_HannWnd = y.*hann(L);
%     Ydft_HannWnd = fft(y_HannWnd,L)/L;
%     FFT_int_1=abs(Ydft_HannWnd);
%     L = L-1;
%     FFT_int_2 = FFT_int_1(1:L/2+1);
%     FFT_int_2(2:end-1) = 2*FFT_int_2(2:end-1);
% % plot
%    
%     % plot frequency preditions
%     l1 = loglog([f_MAC/freq_kappa f_MAC/freq_kappa], [10^-5 10], 'r-','linewidth',1);
%     hold on
%     grid on 
%     ylim([10^-5 10]);
%     xlim([0.5 10^4]);
%     xticks([1 10 10^2 10^3 10^4])
%     yticks([10^-4 10^-2 1]) 
%     l2 = loglog([f_CIA/freq_kappa f_CIA/freq_kappa], ylim, 'b-','linewidth',1);
%     %l3 = loglog([f_MC2/freq_kappa f_MC2/freq_kappa], ylim, 'c--','linewidth',1);
%     %l4 = loglog([f_MC/freq_kappa f_MC/freq_kappa], ylim, 'g--','linewidth',1);
%     l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
%     l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
%     l7 = loglog([f_RC_osc/freq_kappa f_RC_osc/freq_kappa], ylim, 'b:','linewidth',2);
%     %l7 = loglog([f_BZF/freq_kappa f_BZF/freq_kappa], ylim, 'b:','linewidth',2);
%     l8 = loglog([f_RCwall/freq_kappa f_RCwall/freq_kappa], ylim, 'b-.','linewidth',1);   
%     for j = 2:4
%     l9 = loglog([j*f_RCwall/freq_kappa j*f_RCwall/freq_kappa], [0.1 10], 'k--','linewidth',1);
%     end
%     
%     
%     loglog(fk,2*FFT_int_2,'LineWidth',1); % Multiply by 2 because 
% %                                                % Hanning wnd Amplitude 
% %                                                % Correction Factor = 2  
%   
%     [pk,MaxFreq] = findpeaks(2*FFT_int_2,'NPeaks',1,'SortStr','descend');%peakfreq
%   
%     peak=loglog(fk(MaxFreq),pk,'or');
%     maxf(num) = fk(MaxFreq);
% %     text(fk(MaxFreq)*1.5, pk+0.1,...
% %         ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
% %         'FontName','Times');
% 
%         
%     xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
%     ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
%     legend([l5 l1 l8 l9 l2 l6 l7],...
%         {'RMC wallmode','MAC','RC wallmode', 'resonance','CIA', ...
%         'RMC oscillatory', 'RC oscillatory' },...
%         'fontsize',15,...
%         'location','northeast');
%     set(gca,'fontsize',20,'FontName','Times');
%     title(['Radial, 150-degree, $Ra/Ra_w =\ $', num2str(Ra/Ra_w,3), ...
%         ', $Ro_c =\ $', num2str(Ro_c,3), ', $\Lambda =\ $', num2str(Elsa,3)],...
%         'interpreter','latex','fontsize',22);    
% end   
%     
% %%    
%     
hFig4 = figure(4);
set(hFig4, 'Position', [100 100 800 500])
set(gca,'OuterPosition',[0 0 1 0.95]);
    
    IntFFT =  Int_therm(:,1); 
    %IntFFT = SW_vert0(:,5);
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
   
    % plot frequency preditions
    l1 = loglog([f_MAC/freq_kappa f_MAC/freq_kappa], [10^-5 10], 'r-','linewidth',1);
    hold on
    grid on 
    ylim([10^-5 10]);
    xlim([0.5 10^4]);
    xticks([1 10 10^2 10^3 10^4])
    yticks([10^-4 10^-2 1]) 
    l2 = loglog([f_CIA/freq_kappa f_CIA/freq_kappa], ylim, 'b-','linewidth',1);
    %l3 = loglog([f_MC2/freq_kappa f_MC2/freq_kappa], ylim, 'c--','linewidth',1);
    %l4 = loglog([f_MC/freq_kappa f_MC/freq_kappa], ylim, 'g--','linewidth',1);
    l5 = loglog([f_RMC_wall/freq_kappa f_RMC_wall/freq_kappa], ylim, 'r-.','linewidth',1);
    l6 = loglog([f_RMC_osc/freq_kappa f_RMC_osc/freq_kappa], ylim, 'r:','linewidth',2);
    l7 = loglog([f_RC_osc/freq_kappa f_RC_osc/freq_kappa], ylim, 'b:','linewidth',2);
    %l7 = loglog([f_BZF/freq_kappa f_BZF/freq_kappa], ylim, 'b:','linewidth',2);
    l8 = loglog([f_RCwall/freq_kappa f_RCwall/freq_kappa], ylim, 'b-.','linewidth',1);   
%     for i = 2:4
%     l9 = loglog([i*f_RCwall/freq_kappa i*f_RCwall/freq_kappa], [0.1 10], 'k--','linewidth',1);
%     end
    
    
    loglog(fk,2*FFT_int_2,'k','LineWidth',1); % Multiply by 2 because 
%                                                % Hanning wnd Amplitude 
%                                                % Correction Factor = 2  
  
    [pk,MaxFreq] = findpeaks(2*FFT_int_2,'NPeaks',1,'SortStr','descend');%peakfreq
  
    peak=loglog(fk(MaxFreq),pk,'or');
    maxf(num) = fk(MaxFreq);
%     text(fk(MaxFreq)*1.5, pk+0.1,...
%         ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
%         'FontName','Times');

        
    xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
    ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
    legend([l5 l1 l8 l2 l6 l7],...
        {'RMC wallmode','MAC','RC wallmode','CIA', ...
        'RMC oscillatory', 'RC oscillatory' },...
        'fontsize',15,...
        'location','northeast');
    set(gca,'fontsize',20,'FontName','Times');
    title(['Internal Center 55mm, $Ra/Ra_w =\ $', num2str(Ra/Ra_w,3), ...
        ', $Ro_c =\ $', num2str(Ro_c,3), ', $\Lambda =\ $', num2str(Elsa,3)],...
        'interpreter','latex','fontsize',22);    
    
    
%% Radial 
% figure
% ra = ([0.000, 0.068, R_outer])';
% pkp =([0.0234, 0.1174, 0.2763])';
% f = fit(ra,pkp,'exp1');
% plot(f,ra,pkp)
% hold on
% plot(ra,pkp,'ko')
% 
% xlabel('Radius (m)','interpreter','latex','fontsize',25); % x-axis label
% ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label  
% set(gca,'fontsize',20,'FontName','Times');
%     

% rocs = [0.037 0.0575 0.0671 0.112 0.148 0.196];
% Racs = []
% f = fit(ra,pkp,'exp1');
% plot(ra,pkp,'ko')
% hold on
% xlabel('Radius (m)','interpreter','latex','fontsize',25); % x-axis label
% ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label  
% set(gca,'fontsize',20,'FontName','Times');