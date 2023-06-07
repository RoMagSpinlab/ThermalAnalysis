%% scales
omega = RPM/60*2*pi;        % in rad/s
freq_omega = RPM/60; 
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

totalerr = sqrt(Nu^2*((Ploss_insul/P)^2+(1/400)^2+(2*1e-5/A_fluid)^2+...
    (0.1/mean(DeltaT))^2))