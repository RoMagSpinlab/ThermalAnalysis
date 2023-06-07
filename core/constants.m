%% constants

% Circum = 638.18E-3;
R_outer = 0.1016; % outer radius
r_sw = 0.003175; % sidewall thickness
r_insul = 1.0E-2; % insulation thickness
Rfluid = R_outer-r_sw; % fluid radius
A_fluid = pi*(R_outer-r_sw).^2; % horizontal surface area
A_sw = pi*R_outer^2-A_fluid; % SW horizontal surface area
H = 0.4; % tank height
Aspect = 2*(R_outer-r_sw)/H; % aspect ratio
H_insul =  0.4; % insulation height 
H_cutop = 0.047; % top copper plate thickness
H_cubot = 0.012; % bottom copper plate thicknes
H_cu = 0.002; % Distance between top/bottom thermistors and the surfaces
k_insul = 0.015; % aerogel
% H_altop = 0.019; % 
k_cu = 390; % thermal conductivity of copper
k_sw = 16; % thermal conductivity of stainless steel
k_al = 167; % thermal conductivity of aluminum
g = 9.81; % gravitational acceleration