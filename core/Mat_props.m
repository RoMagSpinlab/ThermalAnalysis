function [rho,alpha,Cp,k,sigma,eta,kappa,nu] = Mat_props(Tmean,fluid)
% Determine fluid properties given the the mean temp
%   Detailed explanation goes here

if fluid==1 %% Water
    % King Thesis
    rho=999.8 + 0.1041.*Tmean - 9.718E-3.*Tmean.^2 + 5.184E-5.*Tmean.^3;
    alpha=-6.82E-5 + 1.70E-5*Tmean^2 - 1.82E-7*Tmean^2 + 1.05E-9*Tmean^3;
    kappa=1.312E-7 + 6.972E-10*Tmean - 5.631E-12*Tmean^2 + 2.633E-14^Tmean^3;
    k=0.5529 + 2.662E-3*Tmean - 2.374E-5*Tmean^2 + 1.108E-7*Tmean^3;
    Cp=k/rho/kappa;
    nu=6.581/rho;
    sigma=20E-1;
    eta=1/(4*pi*1E-7*sigma);
elseif fluid==2 % Gallium
    % Combination of Sources
    %Thermal Expansivity
    alpha=1.25E-4;
    % Density
    rho=6.09E3.*(1-alpha.*(Tmean - 29.8));
    % Specific Heat
    Cp=397.6;
    % Thermal Conductivity
    k=31.4;
    % Viscous Diffusivity
    nu=0.46e-3.*exp(4000./8.3144./(Tmean+273.15))./rho;
    % Thermal Diffusivity
    kappa=k./rho./Cp;
    % Electrical conductivity
    sigma=3.88E6;
    % Magnetic Diffusivity
    eta=1/(4*pi*1E-7*sigma);  
end


end

