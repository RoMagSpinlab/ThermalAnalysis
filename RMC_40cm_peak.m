%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMC Data Analysis for Gamma = 0.5 --- peak frequencies
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
Ra/Ra_ms

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


% angle = [0 60 120 180 240 300];
Fs = 1/(time(3)-time(2));
titlestringFFT = {'Center, 55mm', '30-degree, 2/3R, 16mm', ...
    '90-degree, 2/3R, 105mm', '150-degree, 2/3R, 33mm', '210-degree, 2/3R, 16mm'};


hFig1 = figure(1);
set(hFig1, 'Position', [500 100 800 500])
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
  
    peak =loglog(fk(MaxFreq),pk,'or');
    maxf = fk(MaxFreq);
    text(fk(MaxFreq)*1.5, pk+0.1,...
        ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
        'FontName','Times');

        
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
    
    %%

% time = time(end);
% % 
maxf = zeros(5,1);
num=0;
hFig2 = figure(2);
set(hFig2, 'Position', [100 100 1200 1200])
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
    text(fk(MaxFreq)*1.5,pk+0.1,...
       ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
       'FontName','Times');
    
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
    text(fk(MaxFreq)*1.5, pk+0.1,...
        ['(' num2str(fk(MaxFreq)) ' , ' num2str(pk) ')' ],'fontsize',20,...
        'FontName','Times');

        
    xlabel('$f/f_{\kappa}$','interpreter','latex','fontsize',25); % x-axis label
    ylabel('$\|FFT\|$','interpreter','latex','fontsize',25); % y-axis label
    legend([l1 l2 l3 l4 l5 l6 l8],...
        {'MAC','CIA', 'MC (Yan)', 'MC (Zurner)','RMC wallmode', ...
        'RMC oscillatory', 'RC wallmode' },'fontsize',12, 'location',...
        'northeast');
    set(gca,'fontsize',20,'FontName','Times');
    title('Sidewall, 7/8H, 0-degree','interpreter','latex','fontsize',22);    

    
%% f/f_kappa vs. Ra/Ra_ms

% Ra/Ra_ms for all the RMC cases
% 15W, 50W, 70W, 320W, 800W, 2000W | 
% 100W, 200W, 400W, 1200W, 1600W | from induction tests
Ras = [2.93E+07 7.13E+07 9.77E+07 2.81E+08 5.05E+08 9.19E+08 ...
       1.29E+08 2.13E+08 3.43E+08 7.11E+08 8.62E+08];

Nus = [2.41 4.87 5.06 8.45 12 17.1 ...
       5.4901 6.8106 8.691 12.944 14.4438];

Rocs = [0.037 0.058 0.067 0.112 0.148 0.196 ... 
        0.077 0.0983 0.1229 0.1748 0.1907];
    
RaRams = [0.5139 1.2514 1.7146 4.9224 8.8548 16.1187 ...
          2.2577 3.7317 6.0144 12.4782 15.1159];

fpk_c_kappa = [15.1729 19.1042 30.9522 43.0896 116.6629 28.2613 ...
               27.7262 32.5552 43.68 48.5477 42.8247]; 
           
fpk_SW78H_kappa = [15.1729 19.1042 15.4761 19.8875 6.8357 4.915 ...
                   16.4304 14.08 1.477 17.4772 14.2749];



%% Nu vs. Ra

hFig3 = figure(3);
set(hFig3, 'Position', [400 100 800 500])
set(gca,'OuterPosition',[0 0 1 0.95]);

% plot central thermistor FFTs 
f3a = loglog(Ras,Nus, 'ko', 'markersize',10);
hold on 
grid on

x0 = logspace(7,10,100);

%RBC trend
    
p1 = polyfit(log(Ras),log(Nus),1);
p2 = polyfit(log(Ras),log(Nus-1),1);
p1n = num2str(p1(1),3);
f3b = loglog(x0,exp(p1(2)).* x0.^p1(1), 'k:','LineWidth',2);

text(6e7,3,[' Nu = ' num2str(exp(p1(2)),'%.3e') 'Ra^{' p1n '}'],'fontsize',20);
xlim([10^7 2*10^9]);
ylim([1 20]);
xlabel('$Ra$','interpreter','latex','fontsize',25); % x-axis label
ylabel('$Nu$','interpreter','latex','fontsize',25); % y-axis label

% legend([f3a f3b f3c f3d f3e], ...
%     {'Central $7/8H$', 'Sidewall $7/8H$', 'MAC', 'RC wall mode', 'CIA'},...
%     'interpreter','latex','fontsize',20, 'location','northwest');
set(gca,'fontsize',20,'FontName','Times');


Nufit = exp(p2(2)).* x0.^p2(1);

f_MAC_kappa = (Ras.*(Nus-1)/0.027^2.*(1e-6)).^(0.5)*(nu/H^2)./freq_kappa;

f_CIA_kappa = ((Ras.*(Nus-1)/0.027^2).^(2/5)*(1e-6).^(1/5)*(nu/H^2))./freq_kappa;

f_MAC_kappa2 = (x0.*Nufit/0.027^2.*(1e-6)).^(0.5)*(nu/H^2)./freq_kappa;

f_CIA_kappa2 = ((x0.*Nufit/0.027^2).^(2/5)*(1e-6).^(1/5)*(nu/H^2))./freq_kappa;
%%

hFig4 = figure(4);
set(hFig4, 'Position', [500 100 700 500])
set(gca,'OuterPosition',[0 0 1 0.95]);

ax1 = gca();

% plot central thermistor FFTs 
f4a = loglog(RaRams,fpk_c_kappa, 'b^', 'markersize',15,  'linewidth', 2);
hold on 
grid on

% plot MAC predictions
f4b = loglog(RaRams,fpk_SW78H_kappa, 'ko', 'markersize',10, 'linewidth', 1);

%f3c = loglog(sort(RaRams),sort(f_MAC_kappa), 'k--', 'linewidth', 1);

f4c = loglog(x0/Ra_ms,f_MAC_kappa2, 'b--', 'linewidth', 2);

f4d = loglog([0.3 20], [f_RCwall/freq_kappa f_RCwall/freq_kappa], 'k-.', 'linewidth', 1);
%loglog(RaRams, 10.*RaRams.^(0.5), 'k-', 'markersize',10);

f4e = loglog(x0/Ra_ms,f_CIA_kappa2, 'k--', 'linewidth', 1);

loglog([1 1], [5 300], 'k-', 'linewidth', 1);

xlim([0.3 20]);
ylim([5 300]);
xlabel('$Ra/Ra_{ms}$','interpreter','latex','fontsize',25); % x-axis label
ylabel('$f_{pk}/f_{\kappa}$','interpreter','latex','fontsize',25); % y-axis label

legend([f4a f4b f4e f4c f4d], ...
    {'Central $7/8H$', 'Sidewall $7/8H$','CIA (RC Bulk)', 'MAC (RMC Bulk)', 'RC wall mode'},...
    'interpreter','latex','fontsize',18, 'location','northwest');
set(gca,'fontsize',25,'FontName','Times');

% handle second X-axis
ax2 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'YAxisLocation','right', ...
    'Color','none');

ax2.XLim = sqrt([0.3 20].*Ra_ms.*(1e-6)^2/0.027);
xlabel('$Ro_c$', 'interpreter','latex', 'fontsize',25);
ax2.XTick = [0.03 0.05 0.1 0.2];
%ax2.XTickLabel = {'$Ro_c$', 'interpreter','latex' 'fontsize',25};
ax2.XScale = 'log';
ax2.YAxis.Visible = 'off';
set(gca,'fontsize',25,'FontName','Times');
% title('Peak Frequency vs. Magnetostrophic Supercriticality',...
%     'interpreter','latex','fontsize',22);    
%ax1.Box = 'off';
ax2.Box = 'off';