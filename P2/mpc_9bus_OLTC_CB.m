% === MPC with Clustering and Sliding Window for 33-Bus System (Feasibility-Adjusted) ===
clear; clc;
yalmip('clear');

%% === Parameters ===
W = 24; N = 33;
V0 = 1.0; Vmin = 0.96; Vmax = 1.04;
Sbase = 100; dt = 1;
Cbat = 1.5; P_EV_max = 0.75;
S_PV = 1.5;
Q_CB_fixed = 0.1;
CB_OP_LIMIT = 8;
CB_BUNDLE = 3;
Tset = -2:1:2; DV = 0.00625;
tap_values = (1 + Tset * DV).^2; nTaps = length(tap_values);
OLTC_OP_LIMIT = 5;

%% Network topology (parent of each node; 0 = slack)
parent = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 2 19 20 21 3 23 24 6 26 27 28 29 30 31 32 33];  % Node i's parent
r = [0; 0.01; 0.012; 0.011; 0.01; 0.015; 0.015; 0.017; 0.018; 0.013; 0.015; 0.015; 0.017; 0.018; 0.013; 0.015; 0.01; 0.012; 0.011; 0.009; 0.01; 0.012; 0.011; 0.01; 0.015; 0.015; 0.017; 0.018; 0.013; 0.015; 0.015; 0.017; 0.018;]/10;  % Ohm
x = [0; 0.03; 0.025; 0.02; 0.03; 0.035; 0.03; 0.04; 0.045; 0.033; 0.035; 0.03; 0.04; 0.045; 0.033; 0.03; 0.025; 0.02; 0.03; 0.022; 0.03; 0.025; 0.02; 0.03; 0.035; 0.03; 0.04; 0.045; 0.033; 0.035; 0.03; 0.04; 0.045;]/10;     % Ohm

%% Load Data
base_loads = [0; 50; 120; 110; 0; 45; 30; 0; 25; 20; 0; 15; 20; 30; 0; 25; 0; 0; 20; 15; 20; 0; 50; 0; 0; 120; 0; 110; 0; 0; 45; 0; 15];  % kW
pf = 0.95; theta = acos(pf); tan_phi = tan(theta);
load_profile = [ ...
    0.3 0.3 0.35 0.4 0.5 0.7 0.9 1.1 ...
    1.2 1.25 1.2 1.15 1.1 1.0 0.95 0.9 ...
    0.85 0.75 0.6 0.5 0.45 0.4 0.35 0.3];

% Real and reactive power loads (pu)
Pload = zeros(N, W);
Qload = zeros(N, W);
for t = 1:W
    Pload(:, t) = base_loads * load_profile(t) / Sbase;
    Qload(:, t) = Pload(:, t) * tan_phi;
end

PV_buses = [6, 15, 33]; 
nPV = length(PV_buses);
PV_profile = 0.9 * sin(pi*(1:W)/W); 
PV_profile(PV_profile<0) = 0;
PPV = repmat(0.6*PV_profile, nPV, 1);

EV_buses = [10, 20, 25, 30]; nEV = length(EV_buses);
EV_pattern = [0.05 0.1 0.1 0.15 0.1 0.05 0.05 0.1 0.1 0.15 0.15 0.1 0.05 0.05 0.1 0.1 0.15 0.15 0.1 0.05 0.05 0.1 0.1 0.1];
Psch_EV = repmat(EV_pattern, nEV, 1);
SOC_init = 0.2 + 0.3 * rand(nEV,1);
SOC_req = 0.7 * ones(nEV,1);

CB_buses = [5, 17, 28]; nCB = length(CB_buses);

%% === Cluster Definitions ===
clusters = {1:11, 12:22, 23:33};

%% === Local Control Weights ===
w_ev = 0.3 * ones(nEV, 1);
w_pv = 0.002 * ones(nPV, 1);

%% === MPC Sliding Window Simulation ===
window_size = 6;
num_steps = W - window_size + 1;

SOC_total = zeros(nEV, W + 1);
SOC_total(:,1) = SOC_init;
all_v = zeros(N, W);
all_Ttap = zeros(1, W);
Q_PV_full = zeros(nPV, W);
Q_CB_full = zeros(nCB, W);

OLTC_switch = binvar(W - window_size, 1); % OLTC tap switching tracker
prev_tap_idx = sdpvar(1,1);

for t0 = 1:num_steps
    Wm = window_size;
    t_idx = t0:(t0 + Wm - 1);

    dP_EV = sdpvar(nEV, Wm);
    SOC = sdpvar(nEV, Wm+1);
    SOC(:,1) = SOC_total(:,t0);
    absdP = sdpvar(nEV, Wm);
    absQ = sdpvar(nPV, Wm);
    slack_soc = sdpvar(nEV, 1);
    Q_PV = sdpvar(nPV, Wm);
    v = sdpvar(N, Wm);
    P = sdpvar(N, Wm);
    Q = sdpvar(N, Wm);
    I2 = sdpvar(N, Wm);
    d = binvar(nTaps, Wm);
    bt = binvar(W-1, 1);                    % OLTC Switching Limit (not used yet)
    Ttap_val = sdpvar(1, Wm);
    Q_CB = sdpvar(nCB, Wm);
    CB_on = binvar(nCB, CB_BUNDLE, Wm, 'full');
    CB_switch = binvar(nCB, CB_BUNDLE, Wm-1, 'full');

    constraints = [];
    constraints = [constraints, SOC(:,1) == SOC_total(:,t0)];

    for t = 1:Wm-1
        constraints = [constraints,
            Ttap_val(1,t+1) - Ttap_val(1,t) <= 2*2 * bt(t), ...
            Ttap_val(1,t) - Ttap_val(1,t+1) <= 2*(-2) * bt(t)
            ];
    end
    constraints = [constraints, sum(bt) <= OLTC_OP_LIMIT];

    for t = 1:Wm
        constraints = [constraints, sum(d(:,t)) == 1];
        constraints = [constraints, Ttap_val(t) == tap_values * d(:,t)];
        constraints = [constraints, v(1,t) == V0^2 * Ttap_val(t)];

        for i = 2:N
            k = parent(i);
            if k == 0, continue; end
            constraints = [constraints,
                v(i,t) == v(k,t) - 2*(r(i)*P(i,t) + x(i)*Q(i,t)) + (r(i)^2 + x(i)^2)*I2(i,t),
                P(i,t)^2 + Q(i,t)^2 <= I2(i,t) * v(k,t)];

            Pinj = 0; Qinj = 0;
            idx_ev = find(EV_buses == i);
            if ~isempty(idx_ev)
                Pinj = Psch_EV(idx_ev,t_idx(t)) + dP_EV(idx_ev,t);
            end
            idx_pv = find(PV_buses == i);
            if ~isempty(idx_pv)
                Qinj = Qinj + Q_PV(idx_pv,t);
            end
            idx_cb = find(CB_buses == i);
            if ~isempty(idx_cb)
                Qinj = Qinj + Q_CB(idx_cb,t);
            end

            constraints = [constraints,
                P(i,t) == P(k,t) - Pload(i,t_idx(t)) + Pinj - r(i)*I2(i,t),
                Q(i,t) == Q(k,t) - Qload(i,t_idx(t)) + Qinj - x(i)*I2(i,t)];
        end

        constraints = [constraints, Vmin^2 <= v(:,t) <= Vmax^2];

        for e = 1:nEV
            constraints = [constraints,
                -P_EV_max <= dP_EV(e,t) <= P_EV_max];
            P_actual = Psch_EV(e,t_idx(t)) + dP_EV(e,t);
            constraints = [constraints,
                SOC(e,t+1) == SOC(e,t) + dt * P_actual / Cbat,
                0.2 <= SOC(e,t+1) <= 1.0,
                absdP(e,t) >= dP_EV(e,t), absdP(e,t) >= -dP_EV(e,t)];
        end

        for p = 1:nPV
            P_now = PPV(p,t_idx(t));
            Qmax = sqrt(S_PV^2 - P_now^2);
            constraints = [constraints,
                -Qmax <= Q_PV(p,t) <= Qmax,
                absQ(p,t) >= Q_PV(p,t), absQ(p,t) >= -Q_PV(p,t)];
        end
    end

    constraints = [constraints, SOC(:,Wm+1) + slack_soc >= SOC_req, slack_soc >= 0];

    for c = 1:nCB
        
        for t = 1:Wm
            constraints = [constraints, Q_CB(c,t) == sum(CB_on(c,:,t)) * Q_CB_fixed];
        end
        for j = 1:CB_BUNDLE
            for t = 1:(Wm-1)
                constraints = [constraints, ...
                    CB_switch(c,j,t) >= CB_on(c,j,t+1) - CB_on(c,j,t), ...
                    CB_switch(c,j,t) >= CB_on(c,j,t) - CB_on(c,j,t+1)];
            end
            constraints = [constraints, sum(CB_switch(c,j,:)) <= CB_OP_LIMIT];
        end
    end

    objective = 0;
    for cl = 1:length(clusters)
        idx = clusters{cl};
        for t = 1:Wm
            for i = 1:length(idx)
                b = idx(i);
                if ismember(b, EV_buses)
                    ev_idx = find(EV_buses == b);
                    objective = objective + w_ev(ev_idx) * absdP(ev_idx, t);
                end
                if ismember(b, PV_buses)
                    pv_idx = find(PV_buses == b);
                    objective = objective + w_pv(pv_idx) * absQ(pv_idx, t);
                end
            end
        end
    end
    objective = objective + 1000 * sum(slack_soc);

    opts = sdpsettings('solver','gurobi','verbose',0);
    diagnostics = optimize(constraints, objective, opts);

    if diagnostics.problem == 0
        [~, tap_idx] = max(value(d));
        prev_tap_idx = tap_idx(1);

        SOC_total(:,t0+1) = value(SOC(:,2));
        all_v(:,t0) = value(v(:,1));
        all_Ttap(t0) = value(Ttap_val(1));
        Q_PV_full(:,t0) = value(Q_PV(:,1));
        Q_CB_full(:,t0) = value(Q_CB(:,1));
    else
        disp(['❌ MPC failed at step ', num2str(t0)]);
        disp(diagnostics.info);
        break;
    end
end

%% === Visualization ===
time_axis = 1:(W - window_size + 1);

figure;
plot(time_axis, sqrt(all_v(:,1:length(time_axis))'));
title('Bus Voltage Magnitudes Over 24 Hours');
xlabel('Hour'); ylabel('Voltage (p.u.)');
grid on;

figure;
plot(time_axis, all_Ttap(1:length(time_axis)), '-o');
title('OLTC Tap Ratio Over Time');
xlabel('Hour'); ylabel('Tap Ratio');
grid on;

figure;
plot(0:W, SOC_total');
title('EV State of Charge Over 24 Hours');
xlabel('Hour'); ylabel('SOC');
grid on;

figure;
plot(time_axis, Q_PV_full(:,1:length(time_axis))');
title('PV Reactive Power Over Time');
xlabel('Hour'); ylabel('Q_{PV} (p.u.)');
grid on;

figure;
plot(time_axis, Q_CB_full(:,1:length(time_axis))', 'LineWidth', 1.5);
title('Capacitor Bank Reactive Power Over Time');
xlabel('Hour'); ylabel('Q_{CB} (p.u.)');
grid on;
