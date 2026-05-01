clear all

global a b g_t tau_Ca_inact E_Ca_L g_Ca_L K_NaCa d_NaCa gamma_NaCa K_SR
global K_fb K_rb N_fb N_rb Vmaxf Vmaxr k_a_minus k_a_plus k_b_minus k_b_plus
global k_c_minus k_c_plus v1 g_B_Ca g_B_K g_B_Na f_Na g_f HTRPN_tot
global LTRPN_tot k_htrpn_minus k_htrpn_plus k_ltrpn_minus k_ltrpn_plus
global tau_tr tau_xfer CMDN_tot CSQN_tot K_mCMDN K_mCSQN
global V_JSR V_NSR V_SS V_myo g_K1 Cm F R T i_Ca_P_max g_Na K_m_K K_m_Na
global i_NaK_max Ca_o K_o Na_o g_ss
global caff_flag
global A_p K_di K_dsr
global cell_layer model_version
global fluo3_conc fluo3_kon fluo3_koff
global N_p

cell_layer = 'endo'; % 'epi' or 'endo'
model_version = 'adjusted'; % 'original' or 'adjusted'
fluo3_conc = 5e-3;

g_Na = 0.8;

if strcmp(cell_layer, 'endo')
    g_Na = 1.33 * g_Na;
end

g_K1 = 0.024;
g_B_Na = 0.00008015;
g_B_K = 0.000138;
g_f = 0.00145;
i_NaK_max = 0.08;
K_SR = 1.0;

if strcmp(model_version, 'original')

    g_Ca_L = 0.031;
    g_t = 0.035;

    if strcmp(cell_layer, 'endo')
        g_t = 0.4647 * g_t;
    end

    g_ss = 0.007;
    g_B_Ca = 0.0000324;
    i_Ca_P_max = 0.004;
    v1 = 1800;
    K_NaCa = 0.000009984;

elseif strcmp(model_version, 'adjusted')

    K_NaCa = 2.29632e-6;
    i_Ca_P_max = 0.002;
    g_Ca_L = 0.0124;
    g_t = 0.0175;

    if strcmp(cell_layer, 'endo')
        g_t = 0.4647 * g_t;
    end

    g_ss = 0.0049;
    v1 = 1530;
    g_B_Ca = 8.1e-6;

end

Cm = 0.0001;
F = 96487.0;
T = 295.0;
R = 8314.5;

Ca_o = 1.0;
K_o = 5.4;
Na_o = 140.0;

V_JSR = 5.6e-8;
V_NSR = 5.04e-7;
V_SS = 1.2e-9;
V_myo = 9.36e-6;

if strcmp(model_version, 'original')
    HTRPN_tot = 0.14;
    LTRPN_tot = 0.07;
elseif strcmp(model_version, 'adjusted')
    HTRPN_tot = 0.21;
    LTRPN_tot = 0.105;
end

k_htrpn_minus = 0.066;
k_htrpn_plus = 200000.0;
k_ltrpn_minus = 40.0;
k_ltrpn_plus = 40000.0;

if strcmp(model_version, 'original')
    CMDN_tot = 0.05;
elseif strcmp(model_version, 'adjusted')
    CMDN_tot = 0.075;
end

CSQN_tot = 15.0;
K_mCMDN = 0.00238;
K_mCSQN = 0.8;

tau_tr = 0.0005747;
tau_xfer = 0.0267;

a = 0.886;
b = 0.114;

if strcmp(cell_layer, 'endo')
    a = 0.583;
    b = 0.417;
end

tau_Ca_inact = 0.009;
d_NaCa = 0.0001;
gamma_NaCa = 0.5;

k_a_minus = 576.0;
k_a_plus = 12.15e12;
k_b_plus = 4.05e9;

if strcmp(model_version, 'original')
    k_b_minus = 1930;
    k_c_minus = 0.8;
    k_c_plus = 100;
elseif strcmp(model_version, 'adjusted')
    k_b_minus = 5790;
    k_c_minus = 0.8 / 3;
    k_c_plus = 300;
end

f_Na = 0.2;
E_Ca_L = 65.0;

if strcmp(model_version, 'original')

    N_fb = 1.2;
    N_rb = 1.0;
    K_fb = 0.000168;
    K_rb = 3.29;
    Vmaxf = 0.04;
    Vmaxr = 0.9;

elseif strcmp(model_version, 'adjusted')

    A_p = 150;
    K_di = 0.91;
    K_dsr = 2.24;

end

K_m_K = 1.5;
K_m_Na = 10.0;

fluo3_kon = 80e3;
fluo3_koff = 90;

stimdelay = 1e-12;
stimdur = 0.5e-3;
stim_amp = 10;
PCL = 1;

n_stimuli_ctrl = 10;
n_stimuli = 10;

[ctrl_start_times, ctrl_stim_starts, ctrl_stim_ends, ctrl_end_times, ctrl_intervals, ctrl_Istim] = ...
    establish_stim_protocol(PCL, stimdelay, stimdur, stim_amp, n_stimuli_ctrl, 0);

[start_times, stim_starts, stim_ends, end_times, intervals, Istim] = ...
    establish_stim_protocol(PCL, stimdelay, stimdur, stim_amp, n_stimuli, n_stimuli_ctrl * PCL);

sv_begin = load_statevars(cell_layer, model_version, fluo3_conc > 0);

t_ctrl = 0;
sv_ctrl = sv_begin;
N_p = 0;

for iiii = 1:length(ctrl_intervals)

    [t_temp, sv_temp] = ode15s(@dydt_pandit_withEP, ctrl_intervals{iiii}, sv_begin, [], ctrl_Istim(iiii));

    t_ctrl = [t_ctrl; t_temp(2:end)];
    sv_ctrl = [sv_ctrl; sv_temp(2:end,:)];

    sv_begin = sv_temp(end,:);

end

Npores = [0 5 8 10 12 15];

for jj = 1:length(Npores)

    N_p = Npores(jj);
    sv_begin = sv_ctrl(end,:);

    t_ep{jj} = 0;
    sv_ep{jj} = sv_begin;

    for iiii = 1:length(intervals)

        [t_temp, sv_temp] = ode15s(@dydt_pandit_withEP, intervals{iiii}, sv_begin, [], Istim(iiii));

        t_ep{jj} = [t_ep{jj}; t_temp(2:end)];
        sv_ep{jj} = [sv_ep{jj}; sv_temp(2:end,:)];

        sv_begin = sv_temp(end,:);

    end

end

close all

for jj = 1:length(Npores)

    figure;

    width = 8;
    height = 6;

    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [0, 0, width, height]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0, 0, width, height]);

    t_ep_temp = t_ep{jj};
    sv_ep_temp = sv_ep{jj};

    t = [t_ctrl; t_ep_temp];
    Um = [sv_ctrl(:,21); sv_ep_temp(:,21)];
    Ca = [sv_ctrl(:,22); sv_ep_temp(:,17)];

    subplot(2,1,1)
    hold on
    box on
    xlabel('Time (s)');
    ylabel('U (mV)');
    set(gca, 'FontSize', 15)
    xlim([10 20])
    ylim([-100 80])
    yticks(-100:20:100)
    plot(t, Um)

    subplot(2,1,2)
    hold on
    box on
    xlabel('Time (s)');
    ylabel('[Ca]i (\muM)');
    set(gca, 'FontSize', 15)
    xlim([0 10])
    ylim([0 10])
    plot(t, Ca, 'r')

    saveas(gcf, ['Results_for_Npores_' num2str(Npores(jj)) '.png'])

end
