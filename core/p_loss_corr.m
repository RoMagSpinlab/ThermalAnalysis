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
    log((R_outer+0.02)/(R_outer));
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

% Calculating Power from corrected power and conducted power through the sidewall
P_SWcond = k_sw*DeltaT*A_sw/H;
P = mean(Pcorr-P_SWcond)