v_init = -84.5286; % resting voltage
m_init = 0.0017; % sodium current activation gate
h_init = 0.9832; % sodium current fast inactivation
J_init = 0.995484; % slow inactivation
d_init = 0.000003; % Calcium activation gate
f_init = 1.0000; % Calcium inactivation gate
X_init = 0.0057; % activation gate
Ca_init = 0.0002; % Calcium

data.F = 96.5; % Faraday constant, coulombs/mmol
data.R = 8.314; % gas constant, J/K
data.T = 273 + 37; % absolute temperature, K
data.RTF = (data.R * data.T / data.F);

data.K_o = 5.4; % mM
data.K_i = 145; % mM
data.Na_o = 140; % mM
data.Na_i = 18; % mM
data.Ca_o = 1.8;
data.sqrt = sqrt(data.K_o / 5.4);

x0 = [v_init h_init m_init J_init d_init f_init X_init Ca_init];
x0_original = x0;

data.PK = 1.66e-6;
data.PNa_K = 0.01833;

data.GNa_max = 23; % mS/cm^2
data.Gsi_max = 0.09; % mS/cm^2
data.GK_max = 0.282 * data.sqrt; % mS/cm^2
data.GK1_max = 0.6047 * data.sqrt; % mS/cm^2
data.GKp_max = 0.0183; % mS/cm^2
data.Gb_max = 0.03921;
data.Gep = 0; % mS/cm^2

% Equilibrium state simulation
data.Ns = 2;
data.Is = 0;
data.ts = 0.5;
data.tp = 1000;

State_equil = [];
time_equil = [];

for ii = 1:data.Ns
    [t, X] = ode15s('funct', [0 data.tp], x0, [], data);
    State_equil = [State_equil; X];
    time_equil = [time_equil; t + (ii - 1) * data.tp];
    x0 = X(end, :);
end

x0_equil = State_equil(end, :);

figure;
subplot(2,1,1)
plot(time_equil, State_equil(:,1), 'LineWidth', 1.5)
xlabel('Time (ms)');
ylabel('U (mV)');
set(gca, 'FontSize', 15)

subplot(2,1,2)
plot(time_equil, State_equil(:,8) * 1e3, 'LineWidth', 1.5)
xlabel('Time (ms)');
ylabel('[Ca]i (\muM)');
set(gca, 'FontSize', 15)

set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 8, 5]);
saveas(gcf, 'Equilibrium_state_results.png');

fprintf('Equilibrium State Values:\n');
fprintf('%.4f ', x0_equil);
fprintf('\n');

% Response to stimulus simulation
data.Ns = 10;
data.Is = 80;
data.ts = 0.5;
data.tp = 1000;

x0 = x0_equil;
State_stimulus = [];
time_stimulus = [];

for ii = 1:data.Ns
    [t, X] = ode15s('funct', [0 data.tp], x0, [], data);
    State_stimulus = [State_stimulus; X];
    time_stimulus = [time_stimulus; t + (ii - 1) * data.tp];
    x0 = X(end, :);
end

% Electroporation parameters
sigma_e = 1.5;
sigma_i = 0.5;
dm = 5e-9;
rp = 0.76e-9;

sigma_p = (sigma_e - sigma_i) / log(sigma_e / sigma_i);
Gp = 2 * sigma_p * pi * rp^2 / (pi * rp + 2 * dm); % S

num_pores = [0 5 12 15 18];

A_cardiomyocyte = 150e-12 / (0.01 * 1e-4);
Gep_values = num_pores * Gp * 1e3 / A_cardiomyocyte;

fprintf('Calculated Gep values:\n');
for i = 1:length(Gep_values)
    fprintf('%.4f mS/cm^2 \n', Gep_values(i));
end

% Electroporation simulations
for idx = 1:length(Gep_values)

    Gep = Gep_values(idx);
    data.Gep = 0;

    x0 = x0_equil;
    State_electroporation = [];
    time_electroporation = [];

    for ii = 1:data.Ns

        if ii == 4
            data.Gep = Gep;
        end

        [t, X] = ode15s('funct', [0 data.tp], x0, [], data);

        State_electroporation = [State_electroporation; X];
        time_electroporation = [time_electroporation; t + (ii - 1) * data.tp];

        x0 = X(end, :);
    end

    figure;

    width = 8;
    height = 6;
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [0, 0, width, height]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0, 0, width, height]);

    subplot(2,1,1)
    plot(time_electroporation * 1e-3, State_electroporation(:,1), 'LineWidth', 1.5)
    xlabel('Time (s)');
    ylabel('U (mV)');
    set(gca, 'FontSize', 15)
    ylim([-100 100])
    yticks(-100:30:100)
    title(['Gep = ', num2str(Gep), ' mS/cm^2'])

    subplot(2,1,2)
    plot(time_electroporation * 1e-3, State_electroporation(:,8) * 1e3, 'LineWidth', 1.5)
    xlabel('Time (s)');
    ylabel('[Ca]i (\muM)');
    set(gca, 'FontSize', 15)
    ylim([0 10])

    saveas(gcf, ['Gep_', num2str(Gep), '_results.png']);

end
