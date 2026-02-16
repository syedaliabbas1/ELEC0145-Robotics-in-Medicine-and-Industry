%% =========================================================
%  ELEC0145 - Assignment 1 - Task 3: Surgical Planning
%  Adaptive Surgical Planner: MVEE Ellipsoid Slicing + OBB
%  =========================================================
%  PIPELINE:
%    VerticesUnique (from provided starter code)
%       -> PCA  (orientation + elongation metric)
%       -> MVEE (minimum-volume enclosing ellipsoid + fit ratio)
%       -> Adaptive selection: Ellipsoid Slicing OR OBB
%       -> DOF 1-6 waypoints + DOF 7 depth schedule + tip selection
%
%  TWO CASES DEMONSTRATED:
%    Case A: Loaded X/Y/Z data (roughly ellipsoidal tumour)
%            -> Algorithm selects ELLIPSOID SLICING
%    Case B: Synthetic elongated tumour (a=30, b=5, c=5)
%            -> Algorithm selects OBB
%  =========================================================

clear; close all; clc;

%% =========================================================
%  SECTION 0: PARAMETERS
%  All tunable values defined once here.
%  =========================================================
params.margin               = 5;      % mm  - surgical safety margin (consistent with Task 2)
params.w_bulk               = 7.0;    % mm  - serrated bulk-removal tip width (Task 2)
params.w_fine               = 5.5;    % mm  - smooth fine-boundary tip width  (Task 2)
params.standoff             = 30;     % mm  - robot TCP standoff before bone contact
params.DOF7_max_stroke      = 50;     % mm  - physical depth slide stroke limit (Task 2)
params.khachiyan_tol        = 1e-6;   % convergence tolerance for MVEE algorithm
params.khachiyan_maxiter    = 2000;   % max iterations for MVEE
params.fitRatio_threshold   = 2.0;    % F >= this -> use OBB (ellipsoid too conservative)
params.elongation_threshold = 5.0;   % E >= this -> use OBB (tumour too needle-like)

fprintf('============================================================\n');
fprintf('  ELEC0145 Task 3: Adaptive Surgical Planner\n');
fprintf('============================================================\n\n');

%% =========================================================
%  SECTION 1A: CASE A - Load provided data (starter code)
%  =========================================================
fprintf('--- CASE A: Loading provided tumour data ---\n');

Choice = 2;

if Choice == 1
    X = []; Y = []; Z = [];
    eul = [pi/3 pi/6 pi/4];
    rotmZYX = eul2rotm(eul);
    xc = 25; yc = 25; zc = 25;
    a_gen = 20; b_gen = 15; c_gen = 10;
    for i = 1:2:50
        for j = 1:2:50
            for k = 1:2:50
                if ((i-xc)^2/a_gen^2 + (j-yc)^2/b_gen^2 + (k-zc)^2/c_gen^2 + randn/5) < 1
                    RotTrans = rotmZYX * [i-xc; j-yc; k-zc] + [xc; yc; zc];
                    X = [X; RotTrans(1)];
                    Y = [Y; RotTrans(2)];
                    Z = [Z; RotTrans(3)];
                end
            end
        end
    end
else
    load X.mat; load Y.mat; load Z.mat;
end

% Reproduce starter code exactly
[k2_A, V_tumour_A] = convhull(X, Y, Z, 'Simplify', true);

Vertices_A = [X(k2_A(:,1)), Y(k2_A(:,1)), Z(k2_A(:,1));
              X(k2_A(:,2)), Y(k2_A(:,2)), Z(k2_A(:,2));
              X(k2_A(:,3)), Y(k2_A(:,3)), Z(k2_A(:,3))];
VerticesUnique_A = unique(Vertices_A, 'rows');

% Tumour surface area
surfArea_A = compute_surface_area(k2_A, X, Y, Z);

fprintf('  Vertices: %d points\n', size(VerticesUnique_A, 1));
fprintf('  Volume:   %.2f mm^3\n', V_tumour_A);
fprintf('  Surface:  %.2f mm^2\n\n', surfArea_A);

%% =========================================================
%  SECTION 1B: CASE B - Synthetic elongated tumour
%  =========================================================
fprintf('--- CASE B: Generating elongated tumour (a=30, b=5, c=5) ---\n');

X_B = []; Y_B = []; Z_B = [];
xc_B = 25; yc_B = 25; zc_B = 25;
a_B = 30; b_B = 5; c_B = 5;

for i = 1:1:50
    for j = 1:1:50
        for k = 1:1:50
            if ((i-xc_B)^2/a_B^2 + (j-yc_B)^2/b_B^2 + (k-zc_B)^2/c_B^2) < 1
                X_B = [X_B; i]; %#ok<AGROW>
                Y_B = [Y_B; j]; %#ok<AGROW>
                Z_B = [Z_B; k]; %#ok<AGROW>
            end
        end
    end
end

[k2_B, V_tumour_B] = convhull(X_B, Y_B, Z_B, 'Simplify', true);

Vertices_B = [X_B(k2_B(:,1)), Y_B(k2_B(:,1)), Z_B(k2_B(:,1));
              X_B(k2_B(:,2)), Y_B(k2_B(:,2)), Z_B(k2_B(:,2));
              X_B(k2_B(:,3)), Y_B(k2_B(:,3)), Z_B(k2_B(:,3))];
VerticesUnique_B = unique(Vertices_B, 'rows');

surfArea_B = compute_surface_area(k2_B, X_B, Y_B, Z_B);

fprintf('  Vertices: %d points\n', size(VerticesUnique_B, 1));
fprintf('  Volume:   %.2f mm^3\n', V_tumour_B);
fprintf('  Surface:  %.2f mm^2\n\n', surfArea_B);

%% =========================================================
%  RUN FULL PLANNING PIPELINE FOR BOTH CASES
%  =========================================================
fprintf('============================================================\n');
fprintf('  Running planning pipeline...\n');
fprintf('============================================================\n\n');

results_A = run_surgical_planner(VerticesUnique_A, V_tumour_A, surfArea_A, params, 'Case A');
results_B = run_surgical_planner(VerticesUnique_B, V_tumour_B, surfArea_B, params, 'Case B');

%% =========================================================
%  FIGURES
%  =========================================================

%% Figure 1a: Case A - Tumour + resection volume + PCA axes
figure('Name', 'Fig 1a - Case A Tumour and Resection Volume', ...
       'Color', 'w', 'Position', [50 50 700 600]);
plot_resection_volume(VerticesUnique_A, results_A, 'Case A (Ellipsoid Slicing)');

%% Figure 1b: Case B - Tumour + resection volume + PCA axes
figure('Name', 'Fig 1b - Case B Tumour and Resection Volume', ...
       'Color', 'w', 'Position', [800 50 700 600]);
plot_resection_volume(VerticesUnique_B, results_B, 'Case B (OBB)');

%% Figure 2: Slice planes through ellipsoid (Case A)
if strcmp(results_A.method, 'ellipsoid')
    figure('Name', 'Fig 2 - Ellipsoid Slice Planes (Case A)', ...
           'Color', 'w', 'Position', [50 50 700 600]);

    plot_slice_planes(VerticesUnique_A, k2_A, X, Y, Z, results_A);
    title('Fig 2: Ellipsoid Slice Planes Coloured by Depth (Case A)', ...
          'FontSize', 11, 'FontWeight', 'bold');
end

%% Figure 3: Representative layer cross-sections (Case A)
if strcmp(results_A.method, 'ellipsoid')
    figure('Name', 'Fig 3 - Layer Cross-Sections (Case A)', ...
           'Color', 'w', 'Position', [50 50 1300 430]);

    plot_layer_crosssections(results_A);
    sgtitle('Fig 3: Representative Layer Cross-Sections -- Raster (green) + Boundary (orange)', ...
            'FontSize', 11, 'FontWeight', 'bold');
end

%% Figure 4: Method comparison
figure('Name', 'Fig 4 - Method Comparison', ...
       'Color', 'w', 'Position', [50 50 900 650]);

plot_comparison(results_A, results_B);
sgtitle('Fig 4: Adaptive Method Comparison — Case A vs Case B', ...
        'FontSize', 12, 'FontWeight', 'bold');

%% =========================================================
%  FINAL SUMMARY TABLE
%  =========================================================
print_summary(results_A, results_B, params);


%% =========================================================
%  =================== FUNCTIONS ==========================
%  =========================================================

%% ---------------------------------------------------------
function results = run_surgical_planner(V, V_tumour, surfArea, params, caseName)
%  Core planning pipeline. Takes VerticesUnique as input.
%  Returns struct with all planning outputs.
% ---------------------------------------------------------

fprintf('------------------------------------------------------------\n');
fprintf('  %s\n', caseName);
fprintf('------------------------------------------------------------\n');

results.caseName  = caseName;
results.V_tumour  = V_tumour;
results.surfArea  = surfArea;
results.V         = V;

%% STEP 1: PCA — tumour orientation frame
centroid = mean(V, 1);
V0 = V - centroid;
[coeff, score, latent] = pca(V0);

e1 = coeff(:,1)';  % longest axis
e2 = coeff(:,2)';  % second axis
e3 = coeff(:,3)';  % shortest axis — depth axis (DOF 7)
R_pca = coeff;

lambda1 = latent(1);
lambda3 = latent(3);
elongation = lambda1 / lambda3;

results.centroid   = centroid;
results.e1 = e1; results.e2 = e2; results.e3 = e3;
results.R_pca      = R_pca;
results.latent     = latent;
results.elongation = elongation;

fprintf('  PCA: elongation E = lambda1/lambda3 = %.3f\n', elongation);

%% STEP 2: MVEE via Khachiyan algorithm
[A_mvee, c_mvee, Q_mvee, r_axes, V_MVEE] = fit_mvee(V, ...
    params.khachiyan_tol, params.khachiyan_maxiter);

r1 = r_axes(1); r2 = r_axes(2); r3 = r_axes(3);
fitRatio = V_MVEE / V_tumour;

results.A_mvee   = A_mvee;
results.c_mvee   = c_mvee;
results.Q_mvee   = Q_mvee;
results.r1 = r1; results.r2 = r2; results.r3 = r3;
results.V_MVEE   = V_MVEE;
results.fitRatio = fitRatio;

fprintf('  MVEE: semi-axes = [%.2f, %.2f, %.2f] mm\n', r1, r2, r3);
fprintf('  MVEE: V_MVEE = %.2f mm^3,  fit ratio F = %.3f\n', V_MVEE, fitRatio);

%% STEP 3: Adaptive method selection
fprintf('\n  --- ADAPTIVE METHOD SELECTION ---\n');
fprintf('  Fit ratio  F = %.3f  (threshold %.1f)\n', fitRatio, params.fitRatio_threshold);
fprintf('  Elongation E = %.3f  (threshold %.1f)\n', elongation, params.elongation_threshold);

if fitRatio < params.fitRatio_threshold && elongation < params.elongation_threshold
    method = 'ellipsoid';
    fprintf('  SELECTED: Ellipsoid Slicing\n');
    fprintf('  Reason: MVEE fits tightly and tumour is compact\n\n');
else
    method = 'obb';
    fprintf('  SELECTED: OBB (6-plane resection)\n');
    if fitRatio >= params.fitRatio_threshold
        fprintf('  Reason: F=%.2f exceeds threshold (ellipsoid too large)\n\n', fitRatio);
    else
        fprintf('  Reason: E=%.2f exceeds threshold (tumour too elongated)\n\n', elongation);
    end
end
results.method = method;

%% STEP 4A: Ellipsoid slicing
if strcmp(method, 'ellipsoid')
    results = plan_ellipsoid_slicing(results, params);
end

%% STEP 4B: OBB
if strcmp(method, 'obb')
    results = plan_obb(results, score, params);
end

end

%% ---------------------------------------------------------
function [A_mvee, c_mvee, Q_mvee, r_axes, V_MVEE] = fit_mvee(V, tol, maxiter)
%  Khachiyan (Todd-Yildirim) iterative MVEE algorithm.
%
%  Solves: min_{A,c} -log det(A)
%          s.t. (v_i - c)^T A (v_i - c) <= 1  for all i
%
%  INPUT:  V       - Nx3 surface vertices
%  OUTPUT: A_mvee  - 3x3 ellipsoid matrix
%          c_mvee  - 1x3 ellipsoid centre
%          Q_mvee  - 3x3 orientation (cols = ellipsoid axes)
%          r_axes  - [r1 r2 r3] semi-axes, r1>=r2>=r3
%          V_MVEE  - ellipsoid volume (mm^3)
% ---------------------------------------------------------

P = V';          % 3xN — each column is a point
[d, n] = size(P);

% Lift to (d+1)-dimensional space for MVEE with free centre
Q_lift = [P; ones(1, n)];                       % (d+1)xN

% Initialise uniform weights
u = ones(n, 1) / n;

for iter = 1:maxiter
    % Weighted scatter matrix in lifted space
    X_mat = Q_lift * diag(u) * Q_lift';         % (d+1)x(d+1)

    % Mahalanobis distances for all points
    M = diag(Q_lift' * (X_mat \ Q_lift));       % Nx1

    % Maximum distance and its index
    [max_val, j] = max(M);

    % Convergence check
    if (max_val - d - 1) < tol * (d + 1)
        break;
    end

    % Step size (optimal for this j)
    step = (max_val - d - 1) / ((d + 1) * (max_val - 1));

    % Update weights
    u = (1 - step) * u;
    u(j) = u(j) + step;
end

% Ellipsoid centre and shape matrix (using centred weighted scatter)
c_col  = P * u;                                            % 3x1 centre
c_mvee = c_col';                                           % 1x3 centre
S = P * diag(u) * P' - c_col * c_col';                    % 3x3 centred scatter
A_mvee = inv(S) / d;                                      % 3x3 ellipsoid matrix

% Extract semi-axes via eigendecomposition of A_mvee
% A = Q * Lambda * Q^T  =>  r_k = lambda_k^{-1/2}
[Q_raw, Lambda_raw] = eig(A_mvee);
lambda_vals = diag(Lambda_raw);

% Sort ascending so that r_k = 1/sqrt(lambda_k) gives r1>=r2>=r3
[lambda_sorted, idx] = sort(lambda_vals, 'ascend');
Q_mvee = Q_raw(:, idx);
r_axes = 1 ./ sqrt(lambda_sorted);             % [r1 r2 r3], r1>=r2>=r3

V_MVEE = (4/3) * pi * r_axes(1) * r_axes(2) * r_axes(3);

end

%% ---------------------------------------------------------
function results = plan_ellipsoid_slicing(results, params)
%  Compute full resection plan for ellipsoid slicing method.
%  Slices along shortest MVEE axis (r3 direction = depth axis).
% ---------------------------------------------------------

m       = params.margin;
w_bulk  = params.w_bulk;
w_fine  = params.w_fine;

r1 = results.r1; r2 = results.r2; r3 = results.r3;
c_mvee = results.c_mvee;
Q_mvee = results.Q_mvee;

% Margin-expanded semi-axes
% r'_k = r_k + m  (semi-axis expansion; see report for limitation note)
r1p = r1 + m;
r2p = r2 + m;
r3p = r3 + m;

% Resection volume
V_resection = (4/3) * pi * r1p * r2p * r3p;
V_excess    = V_resection - results.V_tumour;
pct_excess  = (V_excess / results.V_tumour) * 100;

% AABB baseline for comparison (what OBB would require for same tumour)
V = results.V;
V_aabb = (max(V(:,1))-min(V(:,1)) + 2*m) * ...
         (max(V(:,2))-min(V(:,2)) + 2*m) * ...
         (max(V(:,3))-min(V(:,3)) + 2*m);
pct_excess_aabb = ((V_aabb - results.V_tumour) / results.V_tumour) * 100;
volume_saved_pct = ((V_aabb - V_resection) / V_aabb) * 100;

% Layer parameterisation
% Slices along e3 (depth axis), spaced by w_bulk
% z_k = layer centre depth from ellipsoid centre along depth axis
% Range: [-r3p, +r3p]
n_layers = ceil(2 * r3p / w_bulk);

% Layer centre depths, symmetric about centre
% Clamped to avoid sqrt(negative) at very tips
z_centres = linspace(-r3p + w_bulk/2, r3p - w_bulk/2, n_layers);
z_centres = max(min(z_centres, r3p - 1e-4), -r3p + 1e-4);

% DOF 7 depth: distance from tool entry face (+r3p) to layer centre
DOF7_vals = r3p + z_centres;   % ranges from ~0 (top layer) to ~2*r3p (bottom)
DOF7_max  = 2 * r3p;

fprintf('  Ellipsoid slicing:\n');
fprintf('    Margin-expanded semi-axes: [%.2f, %.2f, %.2f] mm\n', r1p, r2p, r3p);
fprintf('    Number of layers: %d\n', n_layers);
fprintf('    DOF7 max required: %.2f mm  (stroke limit: %.0f mm)  -- %s\n', ...
    DOF7_max, params.DOF7_max_stroke, ...
    ternary(DOF7_max <= params.DOF7_max_stroke, 'FEASIBLE', 'EXCEEDS STROKE'));
fprintf('    Resection volume: %.2f mm^3  (+%.1f%% over tumour)\n', ...
    V_resection, pct_excess);
fprintf('    AABB baseline:    %.2f mm^3  (+%.1f%% over tumour)\n', ...
    V_aabb, pct_excess_aabb);
fprintf('    Volume saved vs AABB: %.1f%%\n\n', volume_saved_pct);

% Per-layer path generation — pre-allocate ALL elements with defaults
% so degenerate tip layers that hit 'continue' still have valid fields
empty_layer = struct('z', 0, 'a_k', 0, 'b_k', 0, 'DOF7', 0, ...
    'bulk_path', [], 'fine_path', [], ...
    'n_bulk', 0, 'n_fine', 0, 'n_chords', 0);
layers = repmat(empty_layer, 1, n_layers);

total_bulk = 0; total_fine = 0;

for k = 1:n_layers
    z_k = z_centres(k);

    % Cross-section semi-axes at this depth
    scale_k = sqrt(max(0, 1 - (z_k / r3p)^2));
    a_k = r1p * scale_k;
    b_k = r2p * scale_k;

    layers(k).z     = z_k;
    layers(k).a_k   = a_k;
    layers(k).b_k   = b_k;
    layers(k).DOF7  = r3p + z_k;

    % Skip degenerate layers (at exact tips)
    if a_k < 1e-3 || b_k < 1e-3
        layers(k).bulk_path = [];
        layers(k).fine_path = [];
        continue;
    end

    %% Interior raster (bulk tip, spacing w_bulk)
    % Raster lines parallel to e1 (r1 direction), stepping along e2 (r2 direction)
    % Exclude outermost w_fine strip — reserved for boundary pass
    y_start = -b_k + w_fine;
    y_end   =  b_k - w_fine;

    if y_start > y_end
        % Very small layer — only boundary pass
        y_raster = [];
    else
        y_raster = y_start : w_bulk : y_end;
    end

    bulk_path = zeros(2 * length(y_raster), 3);
    for idx = 1:length(y_raster)
        y_j = y_raster(idx);
        x_half = a_k * sqrt(max(0, 1 - (y_j / b_k)^2));

        % Alternating pass direction (boustrophedon)
        if mod(idx, 2) == 1
            pt_A_local = [-x_half;  y_j;  z_k];
            pt_B_local = [ x_half;  y_j;  z_k];
        else
            pt_A_local = [ x_half;  y_j;  z_k];
            pt_B_local = [-x_half;  y_j;  z_k];
        end
        bulk_path(2*idx-1, :) = (c_mvee' + Q_mvee * pt_A_local)';
        bulk_path(2*idx,   :) = (c_mvee' + Q_mvee * pt_B_local)';
    end

    %% Boundary polygon (fine tip, chord approximation)
    % Ramanujan ellipse perimeter approximation
    h_ram = ((a_k - b_k) / (a_k + b_k))^2;
    P_k   = pi * (a_k + b_k) * (1 + 3*h_ram / (10 + sqrt(4 - 3*h_ram)));
    n_chords = ceil(P_k / w_fine);

    theta = linspace(0, 2*pi, n_chords + 1);
    fine_path = zeros(2 * n_chords, 3);

    for jj = 1:n_chords
        pt_A_local = [a_k * cos(theta(jj));   b_k * sin(theta(jj));   z_k];
        pt_B_local = [a_k * cos(theta(jj+1)); b_k * sin(theta(jj+1)); z_k];
        fine_path(2*jj-1, :) = (c_mvee' + Q_mvee * pt_A_local)';
        fine_path(2*jj,   :) = (c_mvee' + Q_mvee * pt_B_local)';
    end

    layers(k).bulk_path = bulk_path;
    layers(k).fine_path = fine_path;
    layers(k).n_bulk    = length(y_raster);
    layers(k).n_chords  = n_chords;
    layers(k).n_fine    = n_chords;

    total_bulk = total_bulk + length(y_raster);
    total_fine = total_fine + n_chords;
end

% Store results
results.layers        = layers;
results.n_layers      = n_layers;
results.z_centres     = z_centres;
results.DOF7_vals     = DOF7_vals;
results.DOF7_max      = DOF7_max;
results.r1p = r1p; results.r2p = r2p; results.r3p = r3p;
results.V_resection   = V_resection;
results.V_excess      = V_excess;
results.pct_excess    = pct_excess;
results.V_aabb        = V_aabb;
results.pct_excess_aabb = pct_excess_aabb;
results.volume_saved_pct = volume_saved_pct;
results.total_bulk    = total_bulk;
results.total_fine    = total_fine;
results.n_cuts        = n_layers;
results.n_reorient    = 1;   % one approach direction for all layers

end

%% ---------------------------------------------------------
function results = plan_obb(results, score, params)
%  Compute OBB resection plan using PCA-aligned bounding box.
%  6 cutting planes with raster paths and tip selection.
% ---------------------------------------------------------

m      = params.margin;
w_bulk = params.w_bulk;
w_fine = params.w_fine;

centroid = results.centroid;
R_pca    = results.R_pca;
e1 = results.e1; e2 = results.e2; e3 = results.e3;

% Bounding box in PCA frame with safety margin
minP = min(score, [], 1) - m;
maxP = max(score, [], 1) + m;
boxSize = maxP - minP;

L1 = boxSize(1); L2 = boxSize(2); L3 = boxSize(3);
h1 = L1/2; h2 = L2/2; h3 = L3/2;

% Box centre in world frame
centrePCA   = (minP + maxP) / 2;
centreWorld = centroid + (R_pca * centrePCA')';

V_OBB      = L1 * L2 * L3;
V_excess   = V_OBB - results.V_tumour;
pct_excess = (V_excess / results.V_tumour) * 100;

fprintf('  OBB:\n');
fprintf('    Dimensions: L1=%.2f  L2=%.2f  L3=%.2f mm\n', L1, L2, L3);
fprintf('    DOF7 required: %.2f mm  (stroke limit: %.0f mm)  -- %s\n', ...
    h3, params.DOF7_max_stroke, ...
    ternary(h3 <= params.DOF7_max_stroke, 'FEASIBLE', 'EXCEEDS STROKE'));
fprintf('    Resection volume: %.2f mm^3  (+%.1f%% over tumour)\n\n', ...
    V_OBB, pct_excess);

% Define 6 cutting planes
faceData = {
    '+e1',  e1,  centreWorld + h1*e1,  e2, e3, L2, L3;
    '-e1', -e1,  centreWorld - h1*e1,  e2, e3, L2, L3;
    '+e2',  e2,  centreWorld + h2*e2,  e1, e3, L1, L3;
    '-e2', -e2,  centreWorld - h2*e2,  e1, e3, L1, L3;
    '+e3',  e3,  centreWorld + h3*e3,  e1, e2, L1, L2;
    '-e3', -e3,  centreWorld - h3*e3,  e1, e2, L1, L2;
};

planes(6) = struct('name','','normal',[],'point',[],'uAxis',[],'vAxis',[],'uLen',0,'vLen',0);
motion(6) = struct('approach',[],'bulk_path',[],'fine_path',[],'DOF7',0,'n_bulk',0,'n_fine',0);

total_bulk = 0; total_fine = 0;

for i = 1:6
    planes(i).name   = faceData{i,1};
    planes(i).normal = faceData{i,2};
    planes(i).point  = faceData{i,3};
    planes(i).uAxis  = faceData{i,4};
    planes(i).vAxis  = faceData{i,5};
    planes(i).uLen   = faceData{i,6};
    planes(i).vLen   = faceData{i,7};

    n_vec = planes(i).normal;
    p0    = planes(i).point;
    u_ax  = planes(i).uAxis;
    v_ax  = planes(i).vAxis;
    hU    = planes(i).uLen / 2;
    hV    = planes(i).vLen / 2;

    motion(i).approach = p0 + params.standoff * n_vec;

    % DOF 7: faces perpendicular to e3 need no extra depth extension
    if i <= 4
        motion(i).DOF7 = h3;
    else
        motion(i).DOF7 = 0;
    end

    % Interior raster (bulk tip)
    y_vals = (-hV + w_fine) : w_bulk : (hV - w_fine);
    bulk_path = zeros(2 * length(y_vals), 3);
    for idx = 1:length(y_vals)
        vv = y_vals(idx);
        if mod(idx, 2) == 1
            ptA = p0 - hU*u_ax + vv*v_ax;
            ptB = p0 + hU*u_ax + vv*v_ax;
        else
            ptA = p0 + hU*u_ax + vv*v_ax;
            ptB = p0 - hU*u_ax + vv*v_ax;
        end
        bulk_path(2*idx-1, :) = ptA;
        bulk_path(2*idx,   :) = ptB;
    end

    % Boundary passes (fine tip) — outermost strip each side
    fine_pts = {p0 - hU*u_ax - hV*v_ax, p0 + hU*u_ax - hV*v_ax;   % bottom edge
                p0 - hU*u_ax + hV*v_ax, p0 + hU*u_ax + hV*v_ax;   % top edge
                p0 - hU*u_ax - hV*v_ax, p0 - hU*u_ax + hV*v_ax;   % left edge
                p0 + hU*u_ax - hV*v_ax, p0 + hU*u_ax + hV*v_ax};  % right edge
    fine_path = zeros(8, 3);
    for fi = 1:4
        fine_path(2*fi-1, :) = fine_pts{fi, 1};
        fine_path(2*fi,   :) = fine_pts{fi, 2};
    end

    motion(i).bulk_path = bulk_path;
    motion(i).fine_path = fine_path;
    motion(i).n_bulk    = length(y_vals);
    motion(i).n_fine    = 4;  % 4 boundary edges per face

    total_bulk = total_bulk + length(y_vals);
    total_fine = total_fine + 4;
end

% OBB corners for visualisation
cornersP = [minP(1) minP(2) minP(3);
            maxP(1) minP(2) minP(3);
            maxP(1) maxP(2) minP(3);
            minP(1) maxP(2) minP(3);
            minP(1) minP(2) maxP(3);
            maxP(1) minP(2) maxP(3);
            maxP(1) maxP(2) maxP(3);
            minP(1) maxP(2) maxP(3)];
cornersW = centroid + (R_pca * cornersP')';

results.planes       = planes;
results.motion       = motion;
results.cornersW     = cornersW;
results.centreWorld  = centreWorld;
results.L1 = L1; results.L2 = L2; results.L3 = L3;
results.h3           = h3;
results.V_resection  = V_OBB;
results.V_excess     = V_excess;
results.pct_excess   = pct_excess;
results.total_bulk   = total_bulk;
results.total_fine   = total_fine;
results.n_cuts       = 6;
results.DOF7_max     = h3;
results.n_reorient   = 6;

end

%% ---------------------------------------------------------
function plot_resection_volume(V, results, titleStr)
%  Plot tumour vertices + resection volume + PCA axes.
% ---------------------------------------------------------

scatter3(V(:,1), V(:,2), V(:,3), 20, [0.8 0.2 0.2], 'filled', ...
         'DisplayName', 'Tumour vertices'); hold on;

centroid = results.centroid;
e1 = results.e1; e2 = results.e2; e3 = results.e3;
scale = 15;

if strcmp(results.method, 'ellipsoid')
    r1p = results.r1p; r2p = results.r2p; r3p = results.r3p;
    c_mvee = results.c_mvee;
    Q_mvee = results.Q_mvee;

    % Original MVEE surface
    [xs, ys, zs] = ellipsoid(0, 0, 0, results.r1, results.r2, results.r3, 25);
    sz = size(xs);
    pts_local = [xs(:)'; ys(:)'; zs(:)'];
    pts_world = c_mvee' + Q_mvee * pts_local;
    xs_w = reshape(pts_world(1,:), sz);
    ys_w = reshape(pts_world(2,:), sz);
    zs_w = reshape(pts_world(3,:), sz);
    surf(xs_w, ys_w, zs_w, 'FaceColor', [0.2 0.8 0.9], 'FaceAlpha', 0.15, ...
         'EdgeColor', 'none', 'DisplayName', 'MVEE (no margin)');

    % Margin-expanded ellipsoid
    [xs2, ys2, zs2] = ellipsoid(0, 0, 0, r1p, r2p, r3p, 25);
    pts_local2 = [xs2(:)'; ys2(:)'; zs2(:)'];
    pts_world2 = c_mvee' + Q_mvee * pts_local2;
    xs_w2 = reshape(pts_world2(1,:), sz);
    ys_w2 = reshape(pts_world2(2,:), sz);
    zs_w2 = reshape(pts_world2(3,:), sz);
    surf(xs_w2, ys_w2, zs_w2, 'FaceColor', [0.2 0.4 0.9], 'FaceAlpha', 0.20, ...
         'EdgeColor', [0.1 0.2 0.7], 'EdgeAlpha', 0.3, ...
         'DisplayName', 'Resection volume (+5mm margin)');

else
    % OBB edges
    cornersW = results.cornersW;
    edgeList = [1 2;2 3;3 4;4 1; 5 6;6 7;7 8;8 5; 1 5;2 6;3 7;4 8];
    firstOBB = true;
    for ei = 1:size(edgeList, 1)
        pA = cornersW(edgeList(ei,1),:);
        pB = cornersW(edgeList(ei,2),:);
        if firstOBB
            plot3([pA(1) pB(1)],[pA(2) pB(2)],[pA(3) pB(3)], ...
                  'b-', 'LineWidth', 2, 'DisplayName', 'OBB resection volume');
            firstOBB = false;
        else
            plot3([pA(1) pB(1)],[pA(2) pB(2)],[pA(3) pB(3)], ...
                  'b-', 'LineWidth', 2, 'HandleVisibility', 'off');
        end
    end

    % Show raster path on largest face
    bp = results.motion(5).bulk_path;
    fp = results.motion(5).fine_path;
    if ~isempty(bp)
        for pi_idx = 1:size(bp,1)/2
            plot3([bp(2*pi_idx-1,1) bp(2*pi_idx,1)], ...
                  [bp(2*pi_idx-1,2) bp(2*pi_idx,2)], ...
                  [bp(2*pi_idx-1,3) bp(2*pi_idx,3)], ...
                  '-', 'Color', [0.1 0.7 0.1], 'LineWidth', 1.2, ...
                  'HandleVisibility', ternary(pi_idx==1,'on','off'), ...
                  'DisplayName', 'Bulk raster (+e3 face)');
        end
    end
end

% PCA axes
quiver3(centroid(1),centroid(2),centroid(3), scale*e1(1),scale*e1(2),scale*e1(3), ...
        0,'r','LineWidth',2.5,'MaxHeadSize',0.6,'DisplayName','e_1 (PC1)');
quiver3(centroid(1),centroid(2),centroid(3), scale*e2(1),scale*e2(2),scale*e2(3), ...
        0,'g','LineWidth',2.5,'MaxHeadSize',0.6,'DisplayName','e_2 (PC2)');
quiver3(centroid(1),centroid(2),centroid(3), scale*e3(1),scale*e3(2),scale*e3(3), ...
        0,'b','LineWidth',2.5,'MaxHeadSize',0.6,'DisplayName','e_3 (depth)');
scatter3(centroid(1),centroid(2),centroid(3),80,'k','filled','DisplayName','Centroid');

axis equal; grid on;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title(titleStr, 'FontSize', 10, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 7);
set(gca, 'FontSize', 9);

end

%% ---------------------------------------------------------
function plot_slice_planes(V, k2, X, Y, Z, results)
%  Plot ellipsoid slice planes coloured by depth (Case A).
% ---------------------------------------------------------

% Tumour hull (transparent)
trisurf(k2, X, Y, Z, 'FaceColor', 'cyan', 'FaceAlpha', 0.15, ...
        'EdgeColor', 'none', 'DisplayName', 'Tumour hull'); hold on;
scatter3(V(:,1), V(:,2), V(:,3), 10, [0.8 0.2 0.2], 'filled', ...
         'DisplayName', 'Tumour vertices');

layers  = results.layers;
n_layers = results.n_layers;
cmap    = jet(n_layers);

for k = 1:n_layers
    a_k = layers(k).a_k;
    b_k = layers(k).b_k;
    z_k = layers(k).z;
    if a_k < 1e-3; continue; end

    % Ellipse boundary for this layer
    theta = linspace(0, 2*pi, 60);
    pts_local = [a_k * cos(theta); b_k * sin(theta); z_k * ones(1,60)];
    pts_world = results.c_mvee' + results.Q_mvee * pts_local;

    plot3(pts_world(1,:), pts_world(2,:), pts_world(3,:), ...
          '-', 'Color', cmap(k,:), 'LineWidth', 2, 'HandleVisibility', 'off');
end

% Colorbar for depth
colormap(jet); caxis([results.z_centres(1), results.z_centres(end)]);
cb = colorbar; cb.Label.String = 'Layer depth z_k (mm)'; cb.FontSize = 9;

axis equal; grid on;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
legend('Tumour hull', 'Tumour vertices', 'Location', 'best', 'FontSize', 8);
set(gca, 'FontSize', 9);

end

%% ---------------------------------------------------------
function plot_layer_crosssections(results)
%  Plot 3 representative layer cross-sections (top, mid, bottom).
% ---------------------------------------------------------

layers   = results.layers;
n_layers = results.n_layers;
indices  = [1, round(n_layers/2), n_layers];
labels   = {'Top layer', 'Middle layer', 'Bottom layer'};

for sp = 1:3
    k = indices(sp);
    subplot(1, 3, sp);

    a_k = layers(k).a_k;
    b_k = layers(k).b_k;
    z_k = layers(k).z;

    if a_k < 1e-3
        title(sprintf('%s\n(degenerate)', labels{sp}), 'FontSize', 9);
        continue;
    end

    hold on;

    % Ellipse boundary (black dashed)
    theta = linspace(0, 2*pi, 200);
    plot(a_k*cos(theta), b_k*sin(theta), 'k--', 'LineWidth', 1.5, ...
         'DisplayName', 'Ellipse boundary');

    % Bulk raster paths (green) — shown in local PCA frame for clarity
    bp = layers(k).bulk_path;
    if ~isempty(bp)
        % Project back to local frame for 2D plot
        for pi_idx = 1:size(bp,1)/2
            pA_w = bp(2*pi_idx-1,:)';
            pB_w = bp(2*pi_idx,:)';
            pA_l = results.Q_mvee' * (pA_w - results.c_mvee');
            pB_l = results.Q_mvee' * (pB_w - results.c_mvee');
            if pi_idx == 1
                plot([pA_l(1) pB_l(1)],[pA_l(2) pB_l(2)], ...
                     '-', 'Color',[0.1 0.7 0.1],'LineWidth',1.5,'DisplayName','Bulk tip raster');
            else
                plot([pA_l(1) pB_l(1)],[pA_l(2) pB_l(2)], ...
                     '-', 'Color',[0.1 0.7 0.1],'LineWidth',1.5,'HandleVisibility','off');
            end
        end
    end

    % Fine boundary polygon (orange)
    fp = layers(k).fine_path;
    if ~isempty(fp)
        for fi = 1:size(fp,1)/2
            pA_w = fp(2*fi-1,:)';
            pB_w = fp(2*fi,:)';
            pA_l = results.Q_mvee' * (pA_w - results.c_mvee');
            pB_l = results.Q_mvee' * (pB_w - results.c_mvee');
            if fi == 1
                plot([pA_l(1) pB_l(1)],[pA_l(2) pB_l(2)], ...
                     '-', 'Color',[0.9 0.5 0.0],'LineWidth',2,'DisplayName','Fine-boundary tip');
            else
                plot([pA_l(1) pB_l(1)],[pA_l(2) pB_l(2)], ...
                     '-', 'Color',[0.9 0.5 0.0],'LineWidth',2,'HandleVisibility','off');
            end
        end
    end

    axis equal; grid on;
    xlabel('e_1 direction (mm)', 'FontSize', 9);
    ylabel('e_2 direction (mm)', 'FontSize', 9);
    title(sprintf('%s\nz_k=%.1fmm  a_k=%.1f  b_k=%.1f\nDOF7=%.1fmm  %d bulk + %d chords', ...
        labels{sp}, z_k, a_k, b_k, layers(k).DOF7, layers(k).n_bulk, layers(k).n_chords), ...
        'FontSize', 8);
    legend('Location', 'best', 'FontSize', 7);
    set(gca, 'FontSize', 9);
end

end

%% ---------------------------------------------------------
function plot_comparison(rA, rB)
%  Side-by-side comparison bar charts for both cases.
% ---------------------------------------------------------

categories = {rA.caseName, rB.caseName};
colors_A = [0.2 0.5 0.8];
colors_B = [0.8 0.4 0.2];

metrics = {
    'Cutting operations',  rA.n_cuts,       rB.n_cuts;
    'Robot reorientations', rA.n_reorient,   rB.n_reorient;
    'Excess volume (%)',    rA.pct_excess,   rB.pct_excess;
    'DOF7 max (mm)',        rA.DOF7_max,     rB.DOF7_max;
};

titles_m = {'Cutting Operations', 'Robot Reorientations', ...
            'Excess Bone Removed (%)', 'DOF 7 Max Extension (mm)'};
ylabels_m = {'No. of layers / planes', 'No. of reorientations', ...
             'Excess volume (%)', 'Extension (mm)'};

for sp = 1:4
    subplot(2, 2, sp);
    vals = [metrics{sp,2}, metrics{sp,3}];
    b = bar(vals, 0.5);
    b.FaceColor = 'flat';
    b.CData(1,:) = colors_A;
    b.CData(2,:) = colors_B;
    set(gca, 'XTickLabel', categories, 'FontSize', 9);
    ylabel(ylabels_m{sp}, 'FontSize', 9);
    title(titles_m{sp}, 'FontSize', 10, 'FontWeight', 'bold');
    grid on;
    % Value labels on bars
    for bi = 1:2
        text(bi, vals(bi) * 1.04, sprintf('%.1f', vals(bi)), ...
             'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    end
    % Method label
    text(1, vals(1)*0.5, rA.method, 'HorizontalAlignment','center',...
         'FontSize',8,'Color','w','FontWeight','bold');
    text(2, vals(2)*0.5, rB.method, 'HorizontalAlignment','center',...
         'FontSize',8,'Color','w','FontWeight','bold');
end

end

%% ---------------------------------------------------------
function print_summary(rA, rB, params)
%  Print final summary table to console.
% ---------------------------------------------------------

fprintf('\n');
fprintf('============================================================\n');
fprintf('         SURGICAL PLANNING — FINAL SUMMARY\n');
fprintf('============================================================\n');
fprintf('%-32s  %-14s  %-14s\n', '', rA.caseName, rB.caseName);
fprintf('%s\n', repmat('-', 1, 64));
fprintf('%-32s  %-14s  %-14s\n', 'Method selected', rA.method, rB.method);
fprintf('%-32s  %-14.2f  %-14.2f\n', 'Tumour volume (mm^3)', rA.V_tumour, rB.V_tumour);
fprintf('%-32s  %-14.2f  %-14.2f\n', 'Tumour surface area (mm^2)', rA.surfArea, rB.surfArea);
fprintf('%-32s  %-14.3f  %-14.3f\n', 'MVEE fit ratio F', rA.fitRatio, rB.fitRatio);
fprintf('%-32s  %-14.3f  %-14.3f\n', 'PCA elongation E', rA.elongation, rB.elongation);
fprintf('%-32s  %-14.2f  %-14.2f\n', 'Semi-axis r1 (mm)', rA.r1, rB.r1);
fprintf('%-32s  %-14.2f  %-14.2f\n', 'Semi-axis r2 (mm)', rA.r2, rB.r2);
fprintf('%-32s  %-14.2f  %-14.2f\n', 'Semi-axis r3 depth (mm)', rA.r3, rB.r3);
fprintf('%-32s  %-14.2f  %-14.2f\n', 'Safety margin (mm)', params.margin, params.margin);
fprintf('%-32s  %-14.2f  %-14.2f\n', 'Resection volume (mm^3)', rA.V_resection, rB.V_resection);
fprintf('%-32s  %-14.1f  %-14.1f\n', 'Excess volume (%)', rA.pct_excess, rB.pct_excess);
fprintf('%-32s  %-14d  %-14d\n', 'Cutting operations', rA.n_cuts, rB.n_cuts);
fprintf('%-32s  %-14d  %-14d\n', 'Robot reorientations', rA.n_reorient, rB.n_reorient);

dof7A = rA.DOF7_max;
dof7B = rB.DOF7_max;
fprintf('%-32s  %-14.2f  %-14.2f\n', 'DOF7 max required (mm)', dof7A, dof7B);
fprintf('%-32s  %-14s  %-14s\n', 'DOF7 feasible? (<=50mm)', ...
    ternary(dof7A<=50,'YES','NO'), ternary(dof7B<=50,'YES','NO'));
fprintf('%-32s  %-14d  %-14d\n', 'Total bulk tip passes', rA.total_bulk, rB.total_bulk);
fprintf('%-32s  %-14d  %-14d\n', 'Total fine tip passes', rA.total_fine, rB.total_fine);
fprintf('%s\n', repmat('-', 1, 64));

fprintf('\nDOF Assignment Summary:\n');
fprintf('  DOF 1-6  (robot joints)  : Position + orientation of tool frame\n');
fprintf('  DOF 7    (depth slide)   : Per-layer/per-plane depth extension\n');
fprintf('  DOF 8    (tip selector)  : Bulk tip (interior) / Fine tip (boundary)\n');
fprintf('============================================================\n');

end

%% ---------------------------------------------------------
function area = compute_surface_area(k2, X, Y, Z)
%  Compute total surface area from convex hull triangulation.
% ---------------------------------------------------------
area = 0;
for i = 1:size(k2, 1)
    v1 = [X(k2(i,1)), Y(k2(i,1)), Z(k2(i,1))];
    v2 = [X(k2(i,2)), Y(k2(i,2)), Z(k2(i,2))];
    v3 = [X(k2(i,3)), Y(k2(i,3)), Z(k2(i,3))];
    area = area + 0.5 * norm(cross(v2-v1, v3-v1));
end
end

%% ---------------------------------------------------------
function out = ternary(cond, a, b)
%  Inline conditional helper.
% ---------------------------------------------------------
if cond; out = a; else; out = b; end
end