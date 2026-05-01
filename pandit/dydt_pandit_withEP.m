function deriv = dydt_pandit_withEP(t,statevar,i_Stim)

%This file contains the ordinary differential equations that are called by
%the main file using the ode15s solver.

global a b g_t tau_Ca_inact E_Ca_L g_Ca_L K_NaCa d_NaCa gamma_NaCa K_SR
global K_fb K_rb N_fb N_rb Vmaxf Vmaxr k_a_minus k_a_plus k_b_minus k_b_plus
global k_c_minus k_c_plus v1 g_B_Ca g_B_K g_B_Na f_Na g_f HTRPN_tot% E_Ca
global LTRPN_tot k_htrpn_minus k_htrpn_plus k_ltrpn_minus k_ltrpn_plus
global tau_tr tau_xfer CMDN_tot CSQN_tot K_mCMDN K_mCSQN
global V_JSR V_NSR V_SS V_myo g_K1 Cm F R T i_Ca_P_max g_Na K_m_K K_m_Na
global i_NaK_max Ca_o K_o Na_o g_ss
global A_p K_di K_dsr
global cell_layer
global fluo3_conc fluo3_kon fluo3_koff
global model_version caff_flag
global N_p

if ~isreal(statevar)
    deriv = zeros(28,1);
else
  
%The formulations of Ito gating differ between endocardial and epicardial
%models.  This flag determines which is used by the solver.
if strcmp(cell_layer, 'endo')
    endoflag = 1;
else
    endoflag = 0;
end

%The original and adjusted formulations have different SERCA formulations.
%This flag determines which is used by the solver.
if strcmp(model_version,'adjusted')
    adjflag = 1;
elseif strcmp(model_version,'original')
    adjflag = 0;
else
    error('Error: model_version improperly specified.')
end


statevarcell = num2cell(statevar) ;
[Y(1), Y(2), Y(3), Y(4), Y(5), Y(6), Y(7), Y(8), Y(9), Y(10), Y(11), Y(12),...
    Y(13), Y(14), Y(15), Y(16), Y(17), Y(18), Y(19), Y(20), Y(21), Y(22),...
    Y(23), Y(24), Y(25), Y(26), Y(27), Y(28)] = deal(statevarcell{:}) ;

% time (second)
%Equilibrium potentials:
E_K = R*T/F*log(K_o/Y(19)); %Potassium membrane equilibrium potential
E_Na = R*T/F*log(Na_o/Y(20)); %Sodium membrane equilibrium potential
E_Ca = 0.5*R*T/F*log(Ca_o/Y(17)); %Calcium membrane equilibrium potential, used for background calcium current only.

%Calcium-independent transient outward K+ current, I_to
i_t = g_t*Y(1)*(a*Y(2)+b*Y(3))*(Y(21)-E_K); %Transient-outward K+ current, pA/pF
r_max = 1.0/(1.0+exp((Y(21)+10.6)/-11.42)); %I_to activation gate maximum
tau_r = 1.0/(45.16*exp(0.03577*(Y(21)+50.0))+98.9*exp(-0.1*(Y(21)+38.0))); %I_to activation gate rate constant
dY(1, 1) = (r_max-Y(1))/tau_r; %dr: change in I_to activation gate rate constant
s_max = 1.0/(1.0+exp((Y(21)+45.3)/6.8841)); %I_to fast inactivation gate maximum
if endoflag == 0 %Epicardial cell
    tau_s = 0.35*exp(-((Y(21)+70.0)/15.0)^2.0)+0.035; %I_to fast inactivation gate rate constant, epicardial cell
elseif endoflag == 1 %Endocardial cell
    tau_s = 0.55*exp(-((Y(21)+70.0)/25.0)^2.0)+0.049; %I_to fast inactivation gate rate constant, endocardial cell
else
    error('Error: value for "choose_endo" not specified')
end
dY(2, 1) = (s_max-Y(2))/tau_s; %ds: change in I_to fast inactivation gate
s_slow_max = 1.0/(1.0+exp((Y(21)+45.3)/6.8841)); %I_to slow inactivation gate maximum
if endoflag == 0 %Epicardial cell
    tau_s_slow = 3.7*exp(-((Y(21)+70.0)/30.0)^2.0)+0.035; %I_to slow inactivation gate rate constant for epicardial cell
elseif endoflag == 1 %Endocardial cell
    tau_s_slow = 3.3*exp(-((Y(21)+70.0)/30.0)^2.0)+0.049; %I_to slow inactivation gate rate constant for endocardial cell
else
    error('Error: value for "choose_endo" not specified')
end
dY(3, 1) = (s_slow_max-Y(3))/tau_s_slow; %ds_slow: change in I_to slow inactivation gate rate constant

%L-type Ca2+ current, I_CaL
i_Ca_L = g_Ca_L*Y(5)*((0.9+Y(4)/10.0)*Y(6)+(0.1-Y(4)/10.0)*Y(7))*(Y(21)-E_Ca_L); %L-type Ca2+ current
Ca_inact_max = 1.0/(1.0+Y(18)/0.01); %I_CaL calcium-dependent inactivation gate maximum
dY(4, 1) = (Ca_inact_max-Y(4))/tau_Ca_inact; %dCa_inact: Change in I_CaL calcium-dependent inactivation gate
d_max = 1.0/(1.0+exp((Y(21)+15.3)/-5.0)); %I_CaL V-dependent activation gate maximum
tau_d = 0.00305*exp(-0.0045*(Y(21)+7.0)^2.0)+0.00105*exp(-0.002*(Y(21)-18.0)^2.0)+0.00025; %I_CaL V-dependent activation gate rate constant
dY(5, 1) = (d_max-Y(5))/tau_d; %dd: Change in I_CaL V-dependent activation gate
f_11_max = 1.0/(1.0+exp((Y(21)+26.7)/5.4)); %I_CaL fast V-dependent inactivation gate maximum
tau_f_11 = 0.105*exp(-((Y(21)+45.0)/12.0)^2.0)+0.04/(1.0+exp((-Y(21)+25.0)/25.0))+0.015/(1.0+exp((Y(21)+75.0)/25.0))+0.0017; %I_CaL fast V-dependent inactivation gate rate constant
dY(6, 1) = (f_11_max-Y(6))/tau_f_11; %df_11: Change in I_CaL fast V-dependent inactivation gate
f_12_max = 1.0/(1.0+exp((Y(21)+26.7)/5.4)); %I_CaL slow V-dependent inactivation gate maximum
tau_f_12 = 0.041*exp(-((Y(21)+47.0)/12.0)^2.0)+0.08/(1.0+exp((Y(21)+55.0)/-5.0))+0.015/(1.0+exp((Y(21)+75.0)/25.0))+0.0017; %I_CaL slow V-dependent inactivation gate rate constant
dY(7, 1) = (f_12_max-Y(7))/tau_f_12; %df_12: Change in I_CaL slow V-dependent inactivation gate

%Sodium-calcium exchanger, I_Na_Ca
i_NaCa = K_NaCa*(Y(20)^3.0*Ca_o*exp(0.03743*Y(21)*gamma_NaCa)-Na_o^3.0*Y(17)*exp(0.03743*Y(21)*(gamma_NaCa-1.0)))/(1.0+d_NaCa*(Y(17)*Na_o^3.0+Ca_o*Y(20)^3.0)); %Sodium-calcium exchanger current

%SERCA2a Ca2+ pump
if adjflag == 0;
    fb = (Y(17)/K_fb)^N_fb;
    rb = (Y(16)/K_rb)^N_rb;
    J_up = K_SR*(Vmaxf*fb-Vmaxr*rb)/(1.0+fb+rb);
elseif adjflag == 1;
    Ca_NSR = Y(16)*1e-3; %converted to M from mM. Although K_d,sr has units of mM, multiplication by 1e-3 suggests that Ca_sr is in M
    Ca_i = Y(17)*1e-3; %converted to M from mM. Although K_d,sr has units of mM, multiplication by 1e-3 suggests that Ca_sr is in M
    Ki = (Ca_i / (1e-3*K_di) ) ^ 2; %unitless
    Ksr = (Ca_NSR / (1e-3*K_dsr) ) ^ 2; %unitless
    D_cycle = 0.104217 + 17.923*Ksr +Ki*(1.75583e6 + 7.61673e6*Ksr)+ ...
        Ki^2*(6.08463e11 + 4.50544e11*Ksr) ; %Denominator for v_cycle. Unitless.
    N_cycle = 3.24873e12*Ki^2 + Ki*(9.17846e6 - 11478.2*Ksr) - 0.329904*Ksr; %Numerator for v_cycle. Unitless
    v_cycle = N_cycle / D_cycle ; %Cycling rate per molecule (s-1)
    J_up = 2 * v_cycle * A_p ; %uM/s, since A_p is in uM;
    J_up = K_SR * J_up * 1e-3; %Convert to mM/s, since Ca_i is in mM in this model
end

%Calcium release channel in sarcoplasmic reticulum (RyR)
dY(8, 1) = -k_a_plus*Y(18)^4.0*Y(8)+k_a_minus*Y(10); %dP_C1: Change in proportion of RyRs in closed state 1
dY(10, 1) = k_a_plus*Y(18)^4.0*Y(8)-(k_a_minus*Y(10)+k_b_plus*Y(18)^3.0*Y(10)+k_c_plus*Y(10))+k_b_minus*Y(11)+k_c_minus*Y(9); %dP_O1: Change in proportion of RyRs in open state 1
dY(11, 1) = k_b_plus*Y(18)^3.0*Y(10)-k_b_minus*Y(11);%dP_O2: Change in proportion of RyRs in open state 2
dY(9, 1) = k_c_plus*Y(10)-k_c_minus*Y(9);%dP_C2: Change in proportion of RyRs in closed state 2
if caff_flag
    J_rel = v1*1*(Y(15)-Y(18)); %If simulating caffeine for this section, J_rel is calculated assuming fully open RyRs
else
    J_rel = v1*(Y(10)+Y(11))*(Y(15)-Y(18)); %Otherwise, calculated using normal RyR gating
end

%Background currents:
i_B_Na = g_B_Na*(Y(21)-E_Na); %Background Na+ current
i_B_Ca = g_B_Ca*(Y(21)-E_Ca); %Background Ca2+ current
i_B_K = g_B_K*(Y(21)-E_K); %Background K+ current
i_B = i_B_Na+i_B_Ca+i_B_K; %Total background current

%Hyperpolarization-activated current (funny current), I_f:
f_K = 1.0-f_Na;
i_f_Na = g_f*Y(12)*f_Na*(Y(21)-E_Na);
i_f_K = g_f*Y(12)*f_K*(Y(21)-E_K);
i_f = i_f_Na+i_f_K; %Total funny current
y_infinity = 1.0/(1.0+exp((Y(21)+138.6)/10.48)); %Hyperpolarization-activated current y gate rate constant maximum
tau_y = 1.0/(0.11885*exp((Y(21)+80.0)/28.37)+0.5623*exp((Y(21)+80.0)/-14.19)); %Hyperpolarization-activated current y gate rate constant
dY(12, 1) = (y_infinity-Y(12))/tau_y; %dy: change in hyperolarization-activated current y gate

%Intracellular and SR Ca2+ fluxes
J_tr = (Y(16)-Y(15))/tau_tr;
J_xfer = (Y(18)-Y(17))/tau_xfer;
J_HTRPNCa = k_htrpn_plus*Y(17)*(HTRPN_tot-Y(13))-k_htrpn_minus*Y(13);
dY(13, 1) = J_HTRPNCa; %dHTRPNCa: Change in high affinity Ca2+-binding site on troponin
J_LTRPNCa = k_ltrpn_plus*Y(17)*(LTRPN_tot-Y(14))-k_ltrpn_minus*Y(14);
dY(14, 1) = J_LTRPNCa; %dLTRPNCa: Change in low affinity Ca2+-binding site on troponin
J_trpn = J_HTRPNCa+J_LTRPNCa;

%Calcium buffering factors
beta_i = 1.0/(1.0+CMDN_tot*K_mCMDN/(K_mCMDN+Y(17))^2.0);
beta_SS = 1.0/(1.0+CMDN_tot*K_mCMDN/(K_mCMDN+Y(18))^2.0);
beta_JSR = 1.0/(1.0+CSQN_tot*K_mCSQN/(K_mCSQN+Y(15))^2.0); %Y(15): JSR free calcium

%Sarcolemmal Ca2+ pump current, ICaP
i_Ca_P = i_Ca_P_max*Y(17)/(Y(17)+0.0004); %%Sarcolemmal Ca2+ pump current

%Voltage-dependent Na+ current, I_Na (cont below)
i_Na = g_Na*Y(24)^3.0*Y(22)*Y(23)*(Y(21)-E_Na); %Voltage-dependent sodium current

%Sodium potassium ATPase, I_NaK
sigma = (exp(Na_o/67.3)-1.0)/7.0;
i_NaK = i_NaK_max/(1.0+0.1245*exp(-0.1*Y(21)*F/(R*T))+0.0365*sigma*exp(-Y(21)*F/(R*T)))*K_o/(K_o+K_m_K)/(1.0+(K_m_Na/Y(20))^1.5); %NKA current

%Steady-state outward K+ current, I_ss (cont below)
i_ss = g_ss*Y(25)*Y(26)*(Y(21)-E_K);

%Inward rectifier K+ current, I_K1
i_K1 = (48.0/(exp((Y(21)+37.0)/25.0)+exp((Y(21)+37.0)/-25.0))+10.0)*0.001/(1.0+exp((Y(21)-(E_K+76.77))/-17.0))+g_K1*(Y(21)-(E_K+1.73))/((1.0+exp(1.613*F*(Y(21)-(E_K+1.73))/(R*T)))*(1.0+exp((K_o-0.9988)/-0.124)));

%Dynamic fluo-3 binding
J_F3_i = fluo3_kon * Y(17)*(fluo3_conc - Y(27)) - fluo3_koff*Y(27); %mM/s. Binding of free fluo-3 to Cai
J_F3_ss = fluo3_kon * Y(18)*(fluo3_conc - Y(28)) - fluo3_koff*Y(28); %mM/s. Binding of free fluo-3 to Cass

%% Electroporation currents -----------------------------------------------

% Define constants
z_Na = 1; % Sodium ion valence
z_K = 1; % Potassium ion valence
z_Cl = -1; % Chloride ion valence
z_Ca = 2; % Calcium ion valence
D_Na = 1.334e-9*(1 + 0.021*(T-298)); % Sodium ion diffusion coefficient [m^2/s]
D_K =  1.957e-9*(1 + 0.021*(T-298)); % Potassium ion diffusion coefficient [m^2/s]
D_Cl = 2.032e-9*(1 + 0.021*(T-298)) ; % Chloride ion diffusion coefficient [m^2/s]
D_Ca = 0.792e-9*(1 + 0.021*(T-298)); % Calcium ion diffusion coefficient [m^2/s]
dm = 5e-9; % [m]
rp = 0.76e-9; % [m]

% Extracellular and intracellular conductivity
Cl_o = Na_o + K_o * 2*Ca_o;
Na_i = Y(20); K_i = Y(19); Ca_i = Y(17);
Cl_i = Na_i + K_i + 2*Ca_i;
sigma_e = F^2 / (R*1e-3 * T) * (z_Na * D_Na * Na_o + z_K * D_K * K_o + z_Cl * D_Cl * Cl_o + z_Ca * D_Ca * Ca_o); % in S/m
sigma_i = F^2 / (R*1e-3 * T) * (z_Na * D_Na * Na_i + z_K * D_K * K_i + z_Cl * D_Cl * Cl_i + z_Ca * D_Ca * Ca_i); % in S/m
chi = sigma_e / sigma_i;

% Molar flux across the pores [mol/s] - outward direction is positive
um = Y(21) * F / (R * T); % nondimensional transmembrane voltage, i.e. no units
J_ep_Na = (N_p * pi * rp^2) * D_Na * (z_Na*um + log(chi))/dm * (chi - 1)/log(chi) * (Na_o - Na_i*exp(z_Na*um))/(1 - chi * exp(z_Na*um));
J_ep_K  = (N_p * pi * rp^2) * D_K  * (z_K *um + log(chi))/dm * (chi - 1)/log(chi) * ( K_o -  K_i*exp(z_K *um))/(1 - chi * exp(z_K *um));
J_ep_Cl = (N_p * pi * rp^2) * D_Cl * (z_Cl*um + log(chi))/dm * (chi - 1)/log(chi) * (Cl_o - Cl_i*exp(z_Cl*um))/(1 - chi * exp(z_Cl*um));
J_ep_Ca = (N_p * pi * rp^2) * D_Ca * (z_Ca*um + log(chi))/dm * (chi - 1)/log(chi) * (Ca_o - Ca_i*exp(z_Ca*um))/(1 - chi * exp(z_Ca*um));

% Electric current across the pores [A * 1e9 = nA] - outward direction is positive
i_ep_Na = J_ep_Na*z_Na*F * 1e9; 
i_ep_K  = J_ep_K *z_K *F * 1e9; 
i_ep_Cl = J_ep_Cl*z_Cl*F * 1e9; 
i_ep_Ca = J_ep_Ca*z_Ca*F * 1e9; 
i_ep_tot = i_ep_Na + i_ep_K + i_ep_Cl + i_ep_Ca;

% If you want to neglect the contribution of the of electroporation current
% on the intracellular concentrations. 
% i_ep_Na = 0; i_ep_K = 0; i_ep_Cl = 0; i_ep_Ca = 0;

%% ------------------------------------------------------------------------
%Ion fluxes
dY(17, 1) = beta_i*(J_xfer-J_F3_i - (J_up+J_trpn+(i_B_Ca+i_ep_Ca-2.0*i_NaCa+i_Ca_P)/(2.0*V_myo*F))); %Change in intracellular Ca2+
dY(20, 1) = -(i_Na+i_B_Na+i_ep_Na+i_NaCa*3.0+i_NaK*3.0+i_f_Na)/(V_myo*F); %dNa_i: change in intracellular sodium
dY(19, 1) = -(i_Stim+i_ss+i_B_K+i_ep_K+i_t+i_K1+i_f_K+i_NaK*-2.0)/(V_myo*F); %dK_i: change in intracellular K+
dY(18, 1) = beta_SS*(-1*J_F3_ss + J_rel*V_JSR/V_SS-J_xfer*V_myo/V_SS-i_Ca_L/(2.0*V_SS*F)); %dCa_ss: change in subspace calcium
dY(15, 1) = beta_JSR*(J_tr-J_rel); %dCa_JSR: Change in junctional SR calcium
dY(16, 1) = J_up*V_myo/V_NSR-J_tr*V_JSR/V_NSR; %dCa_NSR: Change in network SR calcium
dY(27, 1) = J_F3_i;
dY(28, 1) = J_F3_ss;

%Membrane potential:
dY(21, 1) = -(i_Na+i_Ca_L+i_t+i_ss+i_f+i_K1+i_B+i_ep_tot+i_NaK+i_NaCa+i_Ca_P+i_Stim)/Cm; %dV: change in membrane potential

%Voltage-dependent Na+ current, I_Na (began above)
h_max = 1.0/(1.0+exp((Y(21)+76.1)/6.07)); %I_Na fast inactivation gate maximum
if (Y(21) >= -40.0)
    tau_h = 0.0004537*(1.0+exp(-(Y(21)+10.66)/11.1));
else
    tau_h = 0.00349/(0.135*exp(-(Y(21)+80.0)/6.8)+3.56*exp(0.079*Y(21))+310000.0*exp(0.35*Y(21)));
end; %I_Na fast inactivation gate rate constant
dY(22, 1) = (h_max-Y(22))/tau_h; %dh: change in I_Na fast inactivation gate rate constant
j_max = 1.0/(1.0+exp((Y(21)+76.1)/6.07));  %I_Na slow inactivation gate maximum
if (Y(21) >= -40.0)
    tau_j = 0.01163*(1.0+exp(-0.1*(Y(21)+32.0)))/exp(-0.0000002535*Y(21));
else
    tau_j = 0.00349/((Y(21)+37.78)/(1.0+exp(0.311*(Y(21)+79.23)))*(-127140.0*exp(0.2444*Y(21))-0.00003474*exp(-0.04391*Y(21)))+0.1212*exp(-0.01052*Y(21))/(1.0+exp(-0.1378*(Y(21)+40.14))));
end; %I_Na slow inactivation gate rate constant
dY(23, 1) = (j_max-Y(23))/tau_j; %dj: change in slow inactivation gate
m_max = 1.0/(1.0+exp((Y(21)+45.0)/-6.5)); %I_Na fast activation gate maximum
tau_m = 0.00136/(0.32*(Y(21)+47.13)/(1.0-exp(-0.1*(Y(21)+47.13)))+0.08*exp(-Y(21)/11.0)); %I_Na fast activation gate rate constant
dY(24, 1) = (m_max-Y(24))/tau_m; %dm1: Change in I_Na fast activation gate

%Steady-state outward K+ current, I_ss (begun above)
r_ss_max = 1.0/(1.0+exp((Y(21)+11.5)/-11.82)); %I_ss inactivation gate, maximum
tau_r_ss = 10.0/(45.16*exp(0.03577*(Y(21)+50.0))+98.9*exp(-0.1*(Y(21)+38.0)));%I_ss inactivation gate, rate constant
dY(25, 1) = (r_ss_max-Y(25))/tau_r_ss; %dr_ss: change in I_ss inactivation gate
s_ss_max = 1.0/(1.0+exp((Y(21)+87.5)/10.3));%I_ss activation gate, maximum
tau_s_ss = 2.1; %I_ss activation gate, rate constant
dY(26, 1) = (s_ss_max-Y(26))/tau_s_ss; %ds_ss: change in I_ss activation gate

deriv = dY ;
end
