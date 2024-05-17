%%%
%
%   An 1-D thermal conduction profile for Cu-Ga-Cu
%
%%%
%% mise en place
clear all
close all
clc

%%
hFig1 = figure(1);
    set(hFig1, 'Position', [100 300 600 600])
    set(gca,'OuterPosition',[0 0 1 1]);

k_Cu = 390; % W/(m K)
k_Ga = 31.4; % W/(m K)
P = 20; % W
d1 = 0.047; % m
d2 = 0.012; % m
H = 0.4; % m
T0 = 35; % C
R = 0.0984; % m
A = (R^2)*pi; % m^2

T1 = P/A * d1 / k_Cu + T0;

T2 = P/A * H / k_Ga + T1;

T3 = P/A * d2 / k_Cu + T2;
     
plot([T3 T2],[0 d2],'r-','LineWidth',2);   
    hold on 
    grid on
plot([T2 T1],[d2 H+d2],'r-','LineWidth',2);   
plot([T1 T0],[H+d2 H+d1+d2],'r-','LineWidth',2);   

text(36,0.43,'Cu Top Lid','fontsize',22, 'FontName','Times');
text(40,0.4,'$T_1 = 35.08 ^{\circ} \mathrm C$','interpreter','latex','fontsize',22);
text(36,0.2,'Ga (subcritical)','fontsize',22, 'FontName','Times');
text(40,0.025,'$T_2 = 43.45 ^{\circ} \mathrm C$','interpreter','latex','fontsize',22)
text(36, 0.025, 'Cu Bottom Lid','fontsize',22, 'FontName','Times');

%     plot(time_tau_kappa,int90_105,...
%         'color', [0 0.4470 0.7410],'LineWidth',1);
xlim ([floor(T0)-1, floor(T3)+1]);
ylim ([0, H+d1+d2]);

plot(xlim, [d2 d2], 'k--', 'LineWidth',1);
plot(xlim, [H+d2 H+d2], 'k--', 'LineWidth',1);
xlabel('$T(\mathrm{^\circ C})$','interpreter','latex','fontsize',22) % x-axis label
ylabel('$H(\mathrm{m})$','interpreter','latex','fontsize',22) % y-axis label
set(gca,'fontsize',22, 'FontName','Times');
    
    
    
    
    
    
    