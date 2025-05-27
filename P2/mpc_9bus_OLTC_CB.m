% === MPC for 9-Bus System with PVs, EVs, CBs, OLTC, and DistFlow Constraints ===
clear; clc;
yalmip('clear');

%% === Parameters ===
W = 6; N = 9;
V0 = 1.0; Vmin = 0.98; Vmax = 1.02;
Sbase = 100; dt = 1;
Cbat = 1.0; P_EV_max = 1.0;
S_PV = 1.2;
Q_CB_fixed = 0.0625;         % Reactive power per CB unit
CB_OP_LIMIT = 6;           % Max switching operations per bank
CB_BUNDLE = 3;             % Number of CB units per bus
Tset = -5:1:5; DV = 0.00625;
tap_values = (1 + Tset * DV).^2; nTaps = length(tap_values);

%% === Network Topology and Line Parameters ===
parent = [0 1 2 2 3 4 5 6 8];  % 0 = slack bus
r = [0; 0.01; 0.015; 0.02; 0.01; 0.01; 0.01; 0.015; 0.02];
x = [0; 0.03; 0.04; 0.05; 0.03; 0.03; 0.03; 0.04; 0.05];

%% === Load and DER Data ===
Pload = repmat([0;0.3;0.3;0.2;0.2;0.1;0.1;0.25;0.2],1,W);
Qload = 0.4 * Pload;

PV_buses = [3, 8]; nPV = length(PV_buses);
PPV = repmat(0.1, nPV, W);

EV_buses = [6, 9]; nEV = length(EV_buses);
Psch_EV = repmat(0.1, nEV, W);
SOC_init = [0.3; 0.4]; SOC_req = [0.7; 0.7];

CB_buses = [4, 7]; nCB = length(CB_buses);

%% === Decision Variables ===
P = sdpvar(N, W, 'full');
Q = sdpvar(N, W, 'full');
v = sdpvar(N, W, 'full');
I2 = sdpvar(N, W, 'full');

Q_PV = sdpvar(nPV, W);
dP_EV = sdpvar(nEV, W);
SOC = sdpvar(nEV, W+1);
absdP = sdpvar(nEV, W);
absQ = sdpvar(nPV, W);
slack_soc = sdpvar(nEV, 1);

CB_on = binvar(nCB, CB_BUNDLE, W, 'full');
CB_switch = binvar(nCB, CB_BUNDLE, W-1, 'full');
Q_CB = sdpvar(nCB, W);

d = binvar(nTaps, W);             % OLTC tap selector
Ttap_val = sdpvar(1, W);          % Actual tap value (squared ratio)

%% === Constraints ===
constraints = [];
constraints = [constraints, SOC(:,1) == SOC_init];

for t = 1:W
    % OLTC tap selection and substation voltage setting
    constraints = [constraints, sum(d(:,t)) == 1];
    constraints = [constraints, Ttap_val(t) == tap_values * d(:,t)];
    constraints = [constraints, v(1,t) == V0^2 * Ttap_val(t)];

    for i = 2:N
        k = parent(i);
        if k == 0, continue; end

        % DistFlow voltage drop
        constraints = [constraints, ...
            v(i,t) == v(k,t) - 2*(r(i)*P(i,t) + x(i)*Q(i,t)) + ...
                       (r(i)^2 + x(i)^2)*I2(i,t)];

        % SOCP relaxation
        constraints = [constraints, ...
            P(i,t)^2 + Q(i,t)^2 <= I2(i,t) * v(k,t)];

        % Power balance with EVs, PVs, CBs
        Pinj = 0; Qinj = 0;
        idx_ev = find(EV_buses == i);
        if ~isempty(idx_ev)
            Pinj = Psch_EV(idx_ev,t) + dP_EV(idx_ev,t);
        end
        idx_pv = find(PV_buses == i);
        if ~isempty(idx_pv)
            Qinj = Qinj + Q_PV(idx_pv,t);
        end
        idx_cb = find(CB_buses == i);
        if ~isempty(idx_cb)
            Qinj = Qinj + Q_CB(idx_cb,t);
        end

        constraints = [constraints, ...
            P(i,t) == P(k,t) - Pload(i,t) + Pinj - r(i)*I2(i,t), ...
            Q(i,t) == Q(k,t) - Qload(i,t) + Qinj - x(i)*I2(i,t)];
    end

    % Voltage magnitude limits
    constraints = [constraints, Vmin^2 <= v(:,t) <= Vmax^2];

    % EV constraints
    for e = 1:nEV
        constraints = [constraints, -P_EV_max <= dP_EV(e,t) <= P_EV_max];
        P_actual = Psch_EV(e,t) + dP_EV(e,t);
        constraints = [constraints, ...
            SOC(e,t+1) == SOC(e,t) + dt * P_actual / Cbat, ...
            0.2 <= SOC(e,t+1) <= 1.0, ...
            absdP(e,t) >= dP_EV(e,t), absdP(e,t) >= -dP_EV(e,t)];
    end

    % PV constraints
    for p = 1:nPV
        P_now = PPV(p,t);
        Qmax = sqrt(S_PV^2 - P_now^2);
        constraints = [constraints, ...
            -Qmax <= Q_PV(p,t) <= Qmax, ...
            absQ(p,t) >= Q_PV(p,t), absQ(p,t) >= -Q_PV(p,t)];
    end

    % CB constraints
    for c = 1:nCB
        constraints = [constraints, Q_CB(c,t) == sum(CB_on(c,:,t)) * Q_CB_fixed];
        if t < W
            for j = 1:CB_BUNDLE
                constraints = [constraints, ...
                    CB_switch(c,j,t) >= CB_on(c,j,t+1) - CB_on(c,j,t), ...
                    CB_switch(c,j,t) >= CB_on(c,j,t) - CB_on(c,j,t+1)];
            end
        end
    end
end

% CB switching limit
for c = 1:nCB
    for j = 1:CB_BUNDLE
        constraints = [constraints, sum(CB_switch(c,j,:)) <= CB_OP_LIMIT];
    end
end

% Soft SOC satisfaction
constraints = [constraints, SOC(:,W+1) + slack_soc >= SOC_req, slack_soc >= 0];

%% === Objective ===
objective = 0.5 * sum(absdP(:)) + 0.5 * sum(absQ(:)) + 1000 * sum(slack_soc);

%% === Solve ===
opts = sdpsettings('solver','gurobi','verbose',1);
diagnostics = optimize(constraints, objective, opts);

if diagnostics.problem == 0
    disp('✅ MPC solved successfully');
    disp('Final Voltages:'); disp(value(v(:,W)).^0.5);
    disp('Final SOC:'); disp(value(SOC(:,W+1)));
else
    disp('❌ MPC failed');
    disp(diagnostics.info);
    check(constraints)
end


if diagnostics.problem == 0
    disp('✅ MPC solved successfully');
    disp('Final Voltages:'); disp(value(v(:,W)).^0.5);
    disp('Final SOC:'); disp(value(SOC(:,W+1)));

    %% === Plotting ===
    time = 1:W;

    % 1. Bus Voltage Magnitudes (sqrt of v)
    figure;
    plot(time, sqrt(value(v))');
    title('Bus Voltage Magnitudes');
    xlabel('Time Step'); ylabel('Voltage (p.u.)');
    legend(arrayfun(@(i) sprintf('Bus %d', i), 1:N, 'UniformOutput', false));
    grid on;

    % 2. EV Charging Deviation
    figure;
    plot(time, value(dP_EV)');
    title('EV Charging Deviation');
    xlabel('Time Step'); ylabel('Deviation (p.u.)');
    legend(arrayfun(@(e) sprintf('EV at Bus %d', EV_buses(e)), 1:nEV, 'UniformOutput', false));
    grid on;

    % 3. PV Reactive Power
    figure;
    plot(time, value(Q_PV)');
    title('PV Reactive Power Injection');
    xlabel('Time Step'); ylabel('Q (p.u.)');
    legend(arrayfun(@(p) sprintf('PV at Bus %d', PV_buses(p)), 1:nPV, 'UniformOutput', false));
    grid on;

    % 4. CB Reactive Power
    figure;
    plot(time, value(Q_CB)');
    title('Capacitor Bank Reactive Power');
    xlabel('Time Step'); ylabel('Q (p.u.)');
    legend(arrayfun(@(c) sprintf('CB at Bus %d', CB_buses(c)), 1:nCB, 'UniformOutput', false));
    grid on;

    % 5. OLTC Tap Position
    taps = value(d)' * Tset';
    figure;
    plot(time, taps, '-o');
    title('OLTC Tap Position');
    xlabel('Time Step'); ylabel('Tap Position (Steps)');
    grid on;
end
