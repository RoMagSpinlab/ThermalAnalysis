%% mise en place
clear all
close all
clc

% Ra/Ra_ms for all the RMC cases
% 15W, 50W, 70W, 320W, 800W, 2000W | 
% 100W, 200W, 400W, 1200W, 1600W | from induction tests
Ras = [2.93E+07 7.13E+07 9.77E+07 2.81E+08 5.05E+08 9.19E+08 ...
       1.29E+08 2.13E+08 3.43E+08 7.11E+08 8.62E+08];

Nus = [2.41 4.87 5.06 8.45 12 17.1 ...
       5.4901 6.8106 8.691 12.944 14.4438];
%2.93E+07 
RaRC = [1.11E+08 4.44E+08 7.95E+08];
%3.4208 
NuRC = [3.1224 10.1423 18.5051];

RaMC = [2.00E+08 7.88E+08];

NuMC = [9.1228 20.0457];

Rocs = [0.037 0.058 0.067 0.112 0.148 0.196 ... 
        0.077 0.0983 0.1229 0.1748 0.1907];
    
RaRams = [0.5139 1.2514 1.7146 4.9224 8.8548 16.1187 ...
          2.2577 3.7317 6.0144 12.4782 15.1159];

               
% King & Aurnou 2013 Ek = 1e-6
kiRCRa = [1.74e7 3.34e7 4.22e7 5.54e7 7.3e7 1.05e8 1.38e8 1.68e8];

kiRCNu = [0.972 1.05 1.24 1.29 1.51 1.82 2.13 2.39];

%kiRMCRa = [1.68e8 1.65e8 1.52e8 1.43e8]; 
kiRMCRa = 1.43e8; 

%kiRMCNu = [2.41 2.45 2.63 2.79];
kiRMCNu = 2.79;

hFig1 = figure(1);
set(hFig1, 'Position', [400 100 800 500])
set(gca,'OuterPosition',[0 0 1 0.95]);

f1a = loglog(Ras,Nus, 'ko', 'markersize',10);
hold on 
grid on


%f1a_RC = loglog(RaRC, NuRC, 'k^-', 'markersize',10);

%f1a_MC = loglog(RaMC, NuMC,'kd', 'markersize',10');


f1b = loglog(kiRCRa, kiRCNu, 'b^', 'markersize',10);

f1c = loglog(kiRMCRa, kiRMCNu, 'go', 'markersize',10);


x0 = logspace(7,10,100);

%RBC trend
    
p1 = polyfit(log(Ras),log(Nus),1);
p2 = polyfit(log(Ras),log(Nus-1),1);
p1n = num2str(p1(1),3);
f1d = loglog(x0,exp(p1(2)).* x0.^p1(1), 'k:','LineWidth',2);

text(3e8,6,['$ Nu = ' num2str(exp(p1(2)),'%.3e') 'Ra^{' p1n '}$'],...
    'interpreter','latex','fontsize',20);
xlim([10^7 2*10^9]);
ylim([1 30]);
xlabel('$Ra$','interpreter','latex','fontsize',25); % x-axis label
ylabel('$Nu$','interpreter','latex','fontsize',25); % y-axis label

% legend([f1a f1a_RC f1a_MC f1b f1c], ...
%      {'RMC, $\Lambda = 3,\ Ek = 10^{-6}$', ...
%      'RC, $Ek = 10^{-6}$', ...
%      'MC, $Ch = 3\times 10^{-6}$', ...
%      'KA15 RC, $Ek = 10^{-6}$', ...
%      'KA15 RMC, $\Lambda = 0.97,\ Ek \approx 8.75\times 10^{-7}$'},...
%      'interpreter','latex','fontsize',16, 'location','northwest');

legend([f1a  f1b f1c], ...
     {'RMC, $\Lambda = 3,\ Ek = 10^{-6}$', ...
     'KA15 RC, $Ek = 10^{-6}$', ...
     'KA15 RMC, $\Lambda = 0.97,\ Ek \approx 8.75\times 10^{-7}$'},...
     'interpreter','latex','fontsize',16, 'location','northwest');
set(gca,'fontsize',20,'FontName','Times');
