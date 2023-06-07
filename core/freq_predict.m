%% Frequency predictions

f_RC_osc = 4.7*(Ek/Pr)^(1/3)*freq_omega; %RC overstable oscillatory mode frequency
b_ind = U_tw*Rfluid/eta*Bfield;
Conv = Ra*(Nu-1)/Pr^2;
f_MAC = (Conv*Ek)^(1/2)*nu/H^2; % Starchenko & Jones, 2002
f_CIA = Conv^(2/5)*Ek^(1/5)*nu/H^2; % Aubert et al., 2001
f_MC = (1+0.68*N^0.87)^(-1)*U_ff/H; % Zurner et al., 2020
f_MC2 = 0.0324*Ra^(0.9542)*Ch^(-0.6461)/Pr * nu/H^2; % Yan et al., 2019
f_RMC_osc = 0.096*freq_omega; % Horn & Aurnou., 2022 
f_RMC_wall = (8.9236*10^(-4))*freq_omega; % Horn & Aurnou., 2022 
f_BZF = 0.03*Pr^(-4/3)*Ra*Ek^(5/3)*freq_omega; % Zhang et al. 2021
f_RCwall =(132.1*Ek/Pr-1465.5*Ek^(4/3)/Pr)*freq_omega; % Zhang & Liao, 2009