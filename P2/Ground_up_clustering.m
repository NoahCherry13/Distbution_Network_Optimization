% Format: [from_bus, to_bus, R (p.u.), X (p.u.)]
W = 6;
PV_buses = [5, 12, 24];  % buses with PV inverters
EVA_buses = [10, 22];    % buses with EV aggregators
% Initialization
nClusters = length(clusters);
boundary_values.P = zeros(nBuses, W);
boundary_values.Q = zeros(nBuses, W);
boundary_values.V2 = ones(nBuses, W); % assume 1.0^2 to start

cluster_states = cell(nClusters, 1);


lines = [
    1  2  0.0922  0.0470;
    2  3  0.4930  0.2511;
    3  4  0.3660  0.1864;
    4  5  0.3811  0.1941;
    5  6  0.8190  0.7070;
    6  7  0.1872  0.6188;
    7  8  1.7114  1.2351;
    8  9  1.0300  0.7400;
    9 10  1.0440  0.7400;
    10 11  0.1966  0.0650;
    11 12  0.3744  0.1238;
    12 13  1.4680  1.1550;
    13 14  0.5416  0.7129;
    14 15  0.5910  0.5260;
    15 16  0.7463  0.5450;
    16 17  1.2890  1.7210;
    17 18  0.7320  0.5740;
    2 19  0.1640  0.1565;
    19 20  1.5042  1.3554;
    20 21  0.4095  0.4784;
    21 22  0.7089  0.9373;
    3 23  0.4512  0.3083;
    23 24  0.8980  0.7091;
    24 25  0.8960  0.7011;
    6 26  0.2030  0.1034;
    26 27  0.2842  0.1447;
    27 28  1.0590  0.9337;
    28 29  0.8042  0.7006;
    29 30 0.5075  0.2585;
    30 31 0.9744  0.9630;
    31 32 0.3105  0.3619;
    32 33 0.3410  0.5302;
    ];

nBuses = length(lines);
parents = buildParentArray(lines, nBuses);

c = louvainCluster(lines, PV_buses, EVA_buses);

clusters = defineClusters(c, PV_buses, EVA_buses);

[P_load, Q_load, P_PV, S_PV, P_EV_base, V_base] = generateTestData(nBuses, PV_buses, EVA_buses, 24);

[Q_ctrl1, P_ctrl1, V_result1] = mpcController(clusters(1), lines, P_PV, S_PV, P_EV_base, P_load, Q_load, V_base, parents, W);
[Q_ctrl2, P_ctrl2, V_result2] = mpcController(clusters(2), lines, P_PV, S_PV, P_EV_base, P_load, Q_load, V_base, parents, W);
[Q_ctrl3, P_ctrl3, V_result3] = mpcController(clusters(3), lines, P_PV, S_PV, P_EV_base, P_load, Q_load, V_base, parents, W);

%{
for i = 1:length(clusters)
    disp("here")
    [Q_ctrl(i), P_ctrl(i), V_result(i)] = mpcController(clusters(i), lines, P_PV, S_PV, P_EV_base, P_load, Q_load, V_base, parents, W);
end
%}
function parent = buildParentArray(lines, nBuses)
% Build parent array from line topology
% lines: [from, to, R, X]
% nBuses: total number of buses

parent = zeros(nBuses, 1);  % initialize with 0 (slack bus)

for i = 1:size(lines,1)
    from = lines(i,1);
    to = lines(i,2);
    parent(to) = from;
end
end

%% Implements modified Louvain clustering for 33-bus VVC system
function [cluster] = louvainCluster(lines, PV_buses, EVA_buses)
nBuses = length(lines) + 1;

%% Voltage Sensitivity
V0 = 1.0;  % Base voltage
SVP = zeros(nBuses); SVQ = zeros(nBuses);
alpha1 = 0.5; alpha2 = 0.5; alpha3 = 0.5; alpha4 = 0.5;

P_PV = rand(length(PV_buses), 1);     % example PV active power at time t
Srated_PV = ones(length(PV_buses), 1);% rated power
Rup = 0.5 * ones(length(EVA_buses), 1); % example upward flexibility
DV = rand(nBuses, 1) * 0.1;            % voltage deviation placeholder

%% --- Step 1: Build Adjacency Matrix ---
A = zeros(nBuses);
R = zeros(nBuses); X = zeros(nBuses);
for i = 1:size(lines,1)
    f = lines(i,1); t = lines(i,2);
    A(f,t) = 1; A(t,f) = 1;
    R(f,t) = lines(i,3); R(t,f) = lines(i,3);
    X(f,t) = lines(i,4); X(t,f) = lines(i,4);
end

%% --- Compute Voltage Sensitivity ---
for i = 1:nBuses
    for n = 1:nBuses
        path = shortestpath(graph(A), 1, n);
        branch_indices = sub2ind(size(R), path(1:end-1), path(2:end));
        SVP(i,n) = sum(R(branch_indices)) / V0;
        SVQ(i,n) = sum(X(branch_indices)) / V0;
    end
end

%% --- Compute Regulation Capabilities ---
Qmax = sqrt(Srated_PV.^2 - P_PV.^2);
SVQ_avg = mean(SVQ(PV_buses,:), 1);
SVP_avg = mean(SVP(EVA_buses,:), 1);

pv_support = sum(Qmax) * SVQ_avg;  % sum(Qmax) is scalar, SVQ_avg is 1 x nBuses
psi_Q = min(1, pv_support ./ max(abs(DV'), 1e-4));
ev_support = sum(Rup) * SVP_avg;
psi_P = min(1, ev_support ./ max(abs(DV'), 1e-4));
psi_C = min(psi_Q + psi_P, 1);

%% --- Step 4: Compute Modified Modularity ---
degree = sum(A, 2);
m = sum(degree) / 2;
cluster = 1:nBuses;


rho_mod = -Inf; improvement = true;
while improvement
    improvement = false;
    for i = 1:nBuses
        neighbors = find(A(i,:) == 1);
        for j = neighbors
            if cluster(i) ~= cluster(j)
                trial_cluster = cluster;
                trial_cluster(trial_cluster == cluster(j)) = cluster(i);
                % Compute new modularity
                rho_topo = 0;
                for p = 1:nBuses
                    for q = 1:nBuses
                        if trial_cluster(p) == trial_cluster(q)
                            rho_topo = rho_topo + (A(p,q) - degree(p)*degree(q)/(2*m));
                        end
                    end
                end
                rho_topo = rho_topo / (2*m);
                rho_new = alpha1 * rho_topo + alpha2 * mean(alpha3*SVP_avg + alpha4*SVQ_avg) * mean(psi_C);
                if rho_new > rho_mod
                    cluster = trial_cluster;
                    rho_mod = rho_new;
                    improvement = true;
                end
            end
        end
    end
end
end

function [clusters] = defineClusters(cluster_assignments, PV_buses, EVA_buses)
clusters = struct();
unique_clusters = unique(cluster_assignments);
num_clusters = length(unique_clusters);
for i = 1:num_clusters
    members = cluster_assignments == unique_clusters(i);
    clusters(i).id = unique_clusters(i);
    clusters(i).members = find(members);
    clusters(i).PVs = intersect(clusters(i).members, PV_buses);
    clusters(i).PVs = find(ismember(PV_buses, clusters(i).PVs));
    clusters(i).EVAs = intersect(clusters(i).members, EVA_buses);
    clusters(i).EVAs = find(ismember(EVA_buses, clusters(i).EVAs));
    clusters(i).size = length(clusters(i).members);
end
end

function [P_load, Q_load, P_PV, S_PV, P_EV_base, V_base] = generateTestData(nBuses, PV_buses, EVA_buses, W)
% Generates synthetic load, PV, and EV data for testing MPC on 33-bus system

% Base values
P_base = 0.5;  % average per-unit active load
Q_base = 0.2;  % average per-unit reactive load
S_PV_base = 1.0; % per-unit PV apparent power
P_PV_peak = 0.8; % max per-unit PV output
P_EV_mean = 0.6; % base charging load per EV

% Load profiles (random variation around base)
P_load = P_base + 0.1*randn(nBuses, W);
Q_load = Q_base + 0.05*randn(nBuses, W);

% Ensure no negative loads
P_load = max(P_load, 0);
Q_load = max(Q_load, 0);

% PV forecast
nPV = length(PV_buses);
P_PV = zeros(nPV, W);
S_PV = S_PV_base * ones(nPV, 1);

for i = 1:nPV
    profile = linspace(0.3, P_PV_peak, W) + 0.05*randn(1,W); % morning ramp-up
    P_PV(i,:) = max(0, min(profile, S_PV_base));
end

% EV base charging profiles
nEV = length(EVA_buses);
P_EV_base = zeros(nEV, W);
for i = 1:nEV
    profile = P_EV_mean + 0.05*randn(1, W);
    P_EV_base(i,:) = max(0, profile);
end

% Voltage base profile (flat initially)
V_base = ones(nBuses, W);
end

function [Q_ctrl, P_ctrl, V_result, cluster_state] = mpcController(cluster, lines, P_PV_full, S_PV_full, P_EV_base_full, P_load_full, Q_load_full, V_base_full, parent, W, boundary_values)
% MPC controller with inter-cluster coordination using passed boundary values
% Inputs:
% - cluster: struct with .nodes, .PVs, .EVAs
% - lines: [from, to, R, X] line matrix
% - parent: parent vector indexed by global bus numbers
% - boundary_values: struct with fields 'P', 'Q', 'V2' indexed by global bus numbers

nodes = cluster.members;
nodes_map = containers.Map(nodes, 1:length(nodes));
nBuses = length(nodes);

% Extract per-cluster data
P_load = P_load_full(nodes, 1:W);
Q_load = Q_load_full(nodes, 1:W);

% Line resistance and reactance for each bus relative to parent
R = zeros(nBuses, 1);
X = zeros(nBuses, 1);
for i = 2:nBuses
    k = parent(nodes(i));
    idx = find((lines(:,1) == nodes(i) & lines(:,2) == k) | (lines(:,1) == k & lines(:,2) == nodes(i)), 1);
    if ~isempty(idx)
        R(i) = lines(idx,3)/100;
        X(i) = lines(idx,4)/100;
    end
end

% Create YALMIP variables
V2 = sdpvar(nBuses, W);
P = sdpvar(nBuses, W);
Q = sdpvar(nBuses, W);

Vmin = 0.95;
Vmax = 1.05;

constraints = [];

% Fix voltage at root node if it exists
for t = 1:W
    if nodes(1) == 1  % slack bus present in cluster
        constraints = [constraints, V2(1,t) == 1.0^2];
    end
    constraints = [constraints, Vmin^2 <= V2(:,t) <= Vmax^2];
end

for t = 1:W
    for i = 1:nBuses
        global_idx = nodes(i);
        k_global = parent(global_idx);

        if k_global == 0
            continue;
        end

        if isKey(nodes_map, k_global)
            k = nodes_map(k_global);
            constraints = [constraints,
                P(i,t) == P(k,t) - P_load(i,t) - R(i)*(P(k,t)^2 + Q(k,t)^2),
                Q(i,t) == Q(k,t) - Q_load(i,t) - X(i)*(P(k,t)^2 + Q(k,t)^2),
                V2(i,t) == V2(k,t) - 2*(R(i)*P(k,t) + X(i)*Q(k,t))];
        else
            % Use passed boundary values
            P_k = boundary_values.P(k_global,t);
            Q_k = boundary_values.Q(k_global,t);
            V2_k = boundary_values.V2(k_global,t);

            constraints = [constraints,
                P(i,t) == P_k - P_load(i,t) - R(i)*(P_k^2 + Q_k^2),
                Q(i,t) == Q_k - Q_load(i,t) - X(i)*(P_k^2 + Q_k^2),
                V2(i,t) == V2_k - 2*(R(i)*P_k + X(i)*Q_k)];
        end
    end
end

objective = sum(sum(P)) + sum(sum(Q));

ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
diagnostics = optimize(constraints, objective, ops);

Q_ctrl = [];
P_ctrl = [];
if diagnostics.problem == 0
    V_result = value(V2);
    cluster_state.P = value(P);
    cluster_state.Q = value(Q);
    cluster_state.V2 = value(V2);
    cluster_state.nodes = nodes;
else
    warning('DistFlow feasibility problem.');
    V_result = NaN(nBuses, W);
    cluster_state = struct('P', [], 'Q', [], 'V2', [], 'nodes', nodes);
end
end
