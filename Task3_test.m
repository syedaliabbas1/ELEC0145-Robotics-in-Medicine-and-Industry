%% =========================================================
%  Task3_test.m — Validation suite for the adaptive surgical planner
%  Tests the algorithm with 6 different tumour geometries
%  to verify correct method selection and numerical outputs.
%  =========================================================

clear; close all; clc;

%% Parameters (same as Task3.m)
params.margin               = 5;
params.w_bulk               = 7.0;
params.w_fine               = 5.5;
params.standoff             = 30;
params.DOF7_max_stroke      = 50;
params.khachiyan_tol        = 1e-6;
params.khachiyan_maxiter    = 2000;
params.fitRatio_threshold   = 2.0;
params.elongation_threshold = 5.0;

%% Define test cases
% Each row: {name, centre, semi-axes [a b c], expected_method}
test_cases = {
    'T1: Sphere (r=15)',             [25 25 25], [15 15 15],   'ellipsoid';
    'T2: Compact ellipsoid',         [25 25 25], [20 15 10],   'ellipsoid';
    'T3: Moderate elongation',       [25 25 25], [25 10 10],   'obb';
    'T4: Needle (a=30, b=5, c=5)',   [25 25 25], [30  5  5],   'obb';
    'T5: Flat disc (a=20, b=20, c=3)', [25 25 25], [20 20  3], 'obb';
    'T6: Asymmetric (a=25, b=15, c=5)', [25 25 25], [25 15  5], 'obb';
};

n_tests = size(test_cases, 1);
all_results = cell(n_tests, 1);
all_pass = true;

fprintf('============================================================\n');
fprintf('  TASK 3 ALGORITHM VALIDATION SUITE\n');
fprintf('  Testing %d synthetic tumour geometries\n', n_tests);
fprintf('  Thresholds: F < %.1f AND E < %.1f -> ellipsoid, else OBB\n', ...
    params.fitRatio_threshold, params.elongation_threshold);
fprintf('============================================================\n\n');

for t = 1:n_tests
    name     = test_cases{t, 1};
    centre   = test_cases{t, 2};
    abc      = test_cases{t, 3};
    expected = test_cases{t, 4};

    a = abc(1); b = abc(2); c = abc(3);
    xc = centre(1); yc = centre(2); zc = centre(3);

    fprintf('############################################################\n');
    fprintf('  %s\n', name);
    fprintf('  True semi-axes: a=%.0f, b=%.0f, c=%.0f mm\n', a, b, c);
    fprintf('  Expected method: %s\n', expected);
    fprintf('############################################################\n\n');

    % Generate integer grid points inside ellipsoid
    X_t = []; Y_t = []; Z_t = [];
    for i = floor(xc - a):ceil(xc + a)
        for j = floor(yc - b):ceil(yc + b)
            for k = floor(zc - c):ceil(zc + c)
                if ((i-xc)^2/a^2 + (j-yc)^2/b^2 + (k-zc)^2/c^2) < 1
                    X_t = [X_t; i]; %#ok<AGROW>
                    Y_t = [Y_t; j]; %#ok<AGROW>
                    Z_t = [Z_t; k]; %#ok<AGROW>
                end
            end
        end
    end

    if length(X_t) < 4
        fprintf('  [SKIP] Too few interior points (%d) for convex hull\n\n', length(X_t));
        continue;
    end

    [k2, V_tumour] = convhull(X_t, Y_t, Z_t);
    Vertices = [X_t(k2(:,1)), Y_t(k2(:,1)), Z_t(k2(:,1));
                X_t(k2(:,2)), Y_t(k2(:,2)), Z_t(k2(:,2));
                X_t(k2(:,3)), Y_t(k2(:,3)), Z_t(k2(:,3))];
    VerticesUnique = unique(Vertices, 'rows');
    surfArea = compute_surface_area(k2, X_t, Y_t, Z_t);

    V_true_ellipsoid = (4/3) * pi * a * b * c;

    fprintf('  --- Geometry ---\n');
    fprintf('  Interior grid points: %d\n', length(X_t));
    fprintf('  Convex hull vertices: %d\n', size(VerticesUnique, 1));
    fprintf('  True ellipsoid volume:  %.2f mm^3\n', V_true_ellipsoid);
    fprintf('  Convex hull volume:     %.2f mm^3  (%.1f%% of true)\n', ...
        V_tumour, V_tumour/V_true_ellipsoid*100);
    fprintf('  Surface area: %.2f mm^2\n\n', surfArea);

    % Run planner
    results = run_surgical_planner(VerticesUnique, V_tumour, surfArea, params, name);
    all_results{t} = results;

    % Verification checks
    fprintf('\n  === VERIFICATION CHECKS ===\n');
    n_fail = 0;

    % Check 1: Method selection matches expected
    if strcmp(results.method, expected)
        fprintf('  [PASS] Method selection: %s (expected %s)\n', results.method, expected);
    else
        fprintf('  [FAIL] Method selection: %s (expected %s)\n', results.method, expected);
        n_fail = n_fail + 1;
    end

    % Check 2: MVEE semi-axes should be close to true semi-axes (sorted descending)
    true_sorted = sort(abc, 'descend');
    mvee_sorted = [results.r1, results.r2, results.r3];
    axis_error = abs(mvee_sorted - true_sorted);
    max_axis_err = max(axis_error);
    if max_axis_err < 2.0
        fprintf('  [PASS] MVEE semi-axes: [%.2f, %.2f, %.2f] vs true [%.0f, %.0f, %.0f]  (max err=%.2f mm)\n', ...
            mvee_sorted(1), mvee_sorted(2), mvee_sorted(3), ...
            true_sorted(1), true_sorted(2), true_sorted(3), max_axis_err);
    else
        fprintf('  [WARN] MVEE semi-axes: [%.2f, %.2f, %.2f] vs true [%.0f, %.0f, %.0f]  (max err=%.2f mm)\n', ...
            mvee_sorted(1), mvee_sorted(2), mvee_sorted(3), ...
            true_sorted(1), true_sorted(2), true_sorted(3), max_axis_err);
    end

    % Check 3: MVEE should enclose all vertices
    A_mvee = results.A_mvee;
    c_mvee = results.c_mvee;
    max_maha = 0;
    for vi = 1:size(VerticesUnique, 1)
        d_vec = (VerticesUnique(vi,:) - c_mvee)';
        maha = d_vec' * A_mvee * d_vec;
        if maha > max_maha
            max_maha = maha;
        end
    end
    if max_maha <= 1.01
        fprintf('  [PASS] MVEE enclosure: max Mahalanobis = %.6f (<= 1.0)\n', max_maha);
    else
        fprintf('  [FAIL] MVEE enclosure: max Mahalanobis = %.6f (> 1.0, vertex outside ellipsoid!)\n', max_maha);
        n_fail = n_fail + 1;
    end

    % Check 4: Fit ratio should be >= 1.0 (MVEE >= tumour volume)
    if results.fitRatio >= 1.0
        fprintf('  [PASS] Fit ratio F = %.3f (>= 1.0)\n', results.fitRatio);
    else
        fprintf('  [FAIL] Fit ratio F = %.3f (< 1.0, MVEE smaller than tumour!)\n', results.fitRatio);
        n_fail = n_fail + 1;
    end

    % Check 5: Semi-axes ordering r1 >= r2 >= r3
    if results.r1 >= results.r2 - 1e-6 && results.r2 >= results.r3 - 1e-6
        fprintf('  [PASS] Semi-axis ordering: r1=%.2f >= r2=%.2f >= r3=%.2f\n', ...
            results.r1, results.r2, results.r3);
    else
        fprintf('  [FAIL] Semi-axis ordering violated: r1=%.2f, r2=%.2f, r3=%.2f\n', ...
            results.r1, results.r2, results.r3);
        n_fail = n_fail + 1;
    end

    % Check 6: DOF7 feasibility
    if results.DOF7_max <= params.DOF7_max_stroke
        fprintf('  [PASS] DOF7 feasible: %.2f mm <= %.0f mm stroke\n', ...
            results.DOF7_max, params.DOF7_max_stroke);
    else
        fprintf('  [WARN] DOF7 exceeds stroke: %.2f mm > %.0f mm\n', ...
            results.DOF7_max, params.DOF7_max_stroke);
    end

    % Check 7: Resection volume > tumour volume
    if results.V_resection > results.V_tumour
        fprintf('  [PASS] Resection volume %.2f > tumour volume %.2f\n', ...
            results.V_resection, results.V_tumour);
    else
        fprintf('  [FAIL] Resection volume %.2f <= tumour volume %.2f\n', ...
            results.V_resection, results.V_tumour);
        n_fail = n_fail + 1;
    end

    % Check 8: Method-specific checks
    if strcmp(results.method, 'ellipsoid')
        % Verify layer count
        r3p = results.r3 + params.margin;
        expected_layers = ceil(2 * r3p / params.w_bulk);
        if results.n_layers == expected_layers
            fprintf('  [PASS] Layer count: %d (expected %d)\n', results.n_layers, expected_layers);
        else
            fprintf('  [FAIL] Layer count: %d (expected %d)\n', results.n_layers, expected_layers);
            n_fail = n_fail + 1;
        end

        % Verify resection volume = (4/3)*pi*r1p*r2p*r3p
        r1p = results.r1 + params.margin;
        r2p = results.r2 + params.margin;
        expected_Vres = (4/3) * pi * r1p * r2p * r3p;
        vol_err = abs(results.V_resection - expected_Vres);
        if vol_err < 1.0
            fprintf('  [PASS] Resection volume: %.2f (expected %.2f, err=%.4f)\n', ...
                results.V_resection, expected_Vres, vol_err);
        else
            fprintf('  [FAIL] Resection volume: %.2f (expected %.2f, err=%.4f)\n', ...
                results.V_resection, expected_Vres, vol_err);
            n_fail = n_fail + 1;
        end

        % Verify all layers have valid cross-sections
        degenerate = 0;
        for lk = 1:results.n_layers
            if results.layers(lk).a_k < 1e-3
                degenerate = degenerate + 1;
            end
        end
        fprintf('  [INFO] Layers: %d valid, %d degenerate\n', ...
            results.n_layers - degenerate, degenerate);

    else
        % OBB checks
        if results.n_cuts == 6
            fprintf('  [PASS] OBB cutting planes: 6\n');
        else
            fprintf('  [FAIL] OBB cutting planes: %d (expected 6)\n', results.n_cuts);
            n_fail = n_fail + 1;
        end

        if results.n_reorient == 6
            fprintf('  [PASS] OBB reorientations: 6\n');
        else
            fprintf('  [FAIL] OBB reorientations: %d (expected 6)\n', results.n_reorient);
            n_fail = n_fail + 1;
        end

        % Verify OBB volume = L1*L2*L3
        expected_Vobb = results.L1 * results.L2 * results.L3;
        vol_err = abs(results.V_resection - expected_Vobb);
        if vol_err < 1.0
            fprintf('  [PASS] OBB volume: %.2f (expected L1*L2*L3 = %.2f)\n', ...
                results.V_resection, expected_Vobb);
        else
            fprintf('  [FAIL] OBB volume: %.2f (expected L1*L2*L3 = %.2f)\n', ...
                results.V_resection, expected_Vobb);
            n_fail = n_fail + 1;
        end
    end

    % Check 9: PCA elongation sanity (compare to true axis ratio)
    true_elongation_ratio = max(abc)^2 / min(abc)^2;
    fprintf('  [INFO] PCA elongation E = %.3f  (true axis ratio squared = %.3f)\n', ...
        results.elongation, true_elongation_ratio);

    % Summary for this test
    if n_fail == 0
        fprintf('\n  >> TEST %d PASSED (all checks)\n\n', t);
    else
        fprintf('\n  >> TEST %d: %d FAILURES\n\n', t, n_fail);
        all_pass = false;
    end
end

%% Final summary table
fprintf('\n');
fprintf('============================================================\n');
fprintf('  VALIDATION SUMMARY\n');
fprintf('============================================================\n');
fprintf('%-35s  %-10s  %-8s  %-8s  %-7s  %-7s  %-10s  %-8s\n', ...
    'Test', 'Method', 'F', 'E', 'r1', 'r3', 'Excess%', 'DOF7');
fprintf('%s\n', repmat('-', 1, 105));

for t = 1:n_tests
    r = all_results{t};
    if isempty(r); continue; end
    expected = test_cases{t, 4};
    match = ternary(strcmp(r.method, expected), '', ' !!MISMATCH');
    fprintf('%-35s  %-10s  %-8.3f  %-8.3f  %-7.2f  %-7.2f  %-10.1f  %-8.2f%s\n', ...
        test_cases{t,1}, r.method, r.fitRatio, r.elongation, ...
        r.r1, r.r3, r.pct_excess, r.DOF7_max, match);
end

fprintf('%s\n', repmat('-', 1, 105));
if all_pass
    fprintf('\n  ALL TESTS PASSED\n');
else
    fprintf('\n  SOME TESTS FAILED -- review above\n');
end
fprintf('============================================================\n');


%% =========================================================
%  FUNCTIONS (copied from Task3.m for standalone execution)
%  =========================================================

function results = run_surgical_planner(V, V_tumour, surfArea, params, caseName)
fprintf('------------------------------------------------------------\n');
fprintf('  %s\n', caseName);
fprintf('------------------------------------------------------------\n');

results.caseName  = caseName;
results.V_tumour  = V_tumour;
results.surfArea  = surfArea;
results.V         = V;

% PCA
centroid = mean(V, 1);
V0 = V - centroid;
% Manual PCA via SVD (avoids Statistics Toolbox dependency)
[~, S_svd, coeff] = svd(V0, 'econ');
score  = V0 * coeff;
latent = (diag(S_svd).^2) / (size(V0, 1) - 1);

e1 = coeff(:,1)'; e2 = coeff(:,2)'; e3 = coeff(:,3)';
R_pca = coeff;
elongation = latent(1) / latent(3);

results.centroid   = centroid;
results.e1 = e1; results.e2 = e2; results.e3 = e3;
results.R_pca      = R_pca;
results.latent     = latent;
results.elongation = elongation;

fprintf('  PCA: eigenvalues = [%.3f, %.3f, %.3f]\n', latent(1), latent(2), latent(3));
fprintf('  PCA: elongation E = lambda1/lambda3 = %.3f\n', elongation);

% MVEE
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

fprintf('  MVEE: centre = [%.2f, %.2f, %.2f]\n', c_mvee(1), c_mvee(2), c_mvee(3));
fprintf('  MVEE: semi-axes = [%.2f, %.2f, %.2f] mm\n', r1, r2, r3);
fprintf('  MVEE: V_MVEE = %.2f mm^3,  fit ratio F = %.3f\n', V_MVEE, fitRatio);

% Adaptive selection
fprintf('\n  --- ADAPTIVE METHOD SELECTION ---\n');
fprintf('  Fit ratio  F = %.3f  (threshold %.1f)  %s\n', fitRatio, params.fitRatio_threshold, ...
    ternary(fitRatio < params.fitRatio_threshold, 'OK', 'EXCEEDS'));
fprintf('  Elongation E = %.3f  (threshold %.1f)  %s\n', elongation, params.elongation_threshold, ...
    ternary(elongation < params.elongation_threshold, 'OK', 'EXCEEDS'));

if fitRatio < params.fitRatio_threshold && elongation < params.elongation_threshold
    method = 'ellipsoid';
    fprintf('  SELECTED: Ellipsoid Slicing\n\n');
else
    method = 'obb';
    fprintf('  SELECTED: OBB (6-plane resection)\n');
    if fitRatio >= params.fitRatio_threshold
        fprintf('  Reason: F=%.2f exceeds threshold\n\n', fitRatio);
    else
        fprintf('  Reason: E=%.2f exceeds threshold\n\n', elongation);
    end
end
results.method = method;

if strcmp(method, 'ellipsoid')
    results = plan_ellipsoid_slicing(results, params);
else
    results = plan_obb(results, score, params);
end

end

function [A_mvee, c_mvee, Q_mvee, r_axes, V_MVEE] = fit_mvee(V, tol, maxiter)
P = V';
[d, n] = size(P);
Q_lift = [P; ones(1, n)];
u = ones(n, 1) / n;

for iter = 1:maxiter
    X_mat = Q_lift * diag(u) * Q_lift';
    M = diag(Q_lift' * (X_mat \ Q_lift));
    [max_val, j] = max(M);
    if (max_val - d - 1) < tol * (d + 1)
        break;
    end
    step = (max_val - d - 1) / ((d + 1) * (max_val - 1));
    u = (1 - step) * u;
    u(j) = u(j) + step;
end

c_col  = P * u;
c_mvee = c_col';
S = P * diag(u) * P' - c_col * c_col';
A_mvee = inv(S) / d;

[Q_raw, Lambda_raw] = eig(A_mvee);
lambda_vals = diag(Lambda_raw);
[lambda_sorted, idx] = sort(lambda_vals, 'ascend');
Q_mvee = Q_raw(:, idx);
r_axes = 1 ./ sqrt(lambda_sorted);

V_MVEE = (4/3) * pi * r_axes(1) * r_axes(2) * r_axes(3);
end

function results = plan_ellipsoid_slicing(results, params)
m      = params.margin;
w_bulk = params.w_bulk;
w_fine = params.w_fine;

r1 = results.r1; r2 = results.r2; r3 = results.r3;
c_mvee = results.c_mvee;
Q_mvee = results.Q_mvee;

r1p = r1 + m; r2p = r2 + m; r3p = r3 + m;

V_resection = (4/3) * pi * r1p * r2p * r3p;
V_excess    = V_resection - results.V_tumour;
pct_excess  = (V_excess / results.V_tumour) * 100;

V = results.V;
V_aabb = (max(V(:,1))-min(V(:,1)) + 2*m) * ...
         (max(V(:,2))-min(V(:,2)) + 2*m) * ...
         (max(V(:,3))-min(V(:,3)) + 2*m);
pct_excess_aabb = ((V_aabb - results.V_tumour) / results.V_tumour) * 100;
volume_saved_pct = ((V_aabb - V_resection) / V_aabb) * 100;

n_layers = ceil(2 * r3p / w_bulk);
z_centres = linspace(-r3p + w_bulk/2, r3p - w_bulk/2, n_layers);
z_centres = max(min(z_centres, r3p - 1e-4), -r3p + 1e-4);

DOF7_vals = r3p + z_centres;
DOF7_max  = 2 * r3p;

fprintf('  Ellipsoid slicing:\n');
fprintf('    Margin-expanded semi-axes: [%.2f, %.2f, %.2f] mm\n', r1p, r2p, r3p);
fprintf('    Number of layers: %d\n', n_layers);
fprintf('    DOF7 max required: %.2f mm  (stroke limit: %.0f mm)  -- %s\n', ...
    DOF7_max, params.DOF7_max_stroke, ...
    ternary(DOF7_max <= params.DOF7_max_stroke, 'FEASIBLE', 'EXCEEDS STROKE'));
fprintf('    Resection volume: %.2f mm^3  (+%.1f%% over tumour)\n', V_resection, pct_excess);
fprintf('    AABB baseline:    %.2f mm^3  (+%.1f%% over tumour)\n', V_aabb, pct_excess_aabb);
fprintf('    Volume saved vs AABB: %.1f%%\n', volume_saved_pct);

empty_layer = struct('z', 0, 'a_k', 0, 'b_k', 0, 'DOF7', 0, ...
    'bulk_path', [], 'fine_path', [], 'n_bulk', 0, 'n_fine', 0, 'n_chords', 0);
layers = repmat(empty_layer, 1, n_layers);
total_bulk = 0; total_fine = 0;

for k = 1:n_layers
    z_k = z_centres(k);
    scale_k = sqrt(max(0, 1 - (z_k / r3p)^2));
    a_k = r1p * scale_k;
    b_k = r2p * scale_k;

    layers(k).z    = z_k;
    layers(k).a_k  = a_k;
    layers(k).b_k  = b_k;
    layers(k).DOF7 = r3p + z_k;

    if a_k < 1e-3 || b_k < 1e-3
        continue;
    end

    % Interior raster
    y_start = -b_k + w_fine;
    y_end   =  b_k - w_fine;
    if y_start > y_end
        y_raster = [];
    else
        y_raster = y_start : w_bulk : y_end;
    end

    bulk_path = zeros(2 * length(y_raster), 3);
    for idx = 1:length(y_raster)
        y_j = y_raster(idx);
        x_half = a_k * sqrt(max(0, 1 - (y_j / b_k)^2));
        if mod(idx, 2) == 1
            pt_A_local = [-x_half; y_j; z_k];
            pt_B_local = [ x_half; y_j; z_k];
        else
            pt_A_local = [ x_half; y_j; z_k];
            pt_B_local = [-x_half; y_j; z_k];
        end
        bulk_path(2*idx-1, :) = (c_mvee' + Q_mvee * pt_A_local)';
        bulk_path(2*idx,   :) = (c_mvee' + Q_mvee * pt_B_local)';
    end

    % Boundary polygon
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

% Per-layer detail log
fprintf('\n    --- Layer details ---\n');
for k = 1:n_layers
    fprintf('    Layer %d: z=%.2f  a_k=%.2f  b_k=%.2f  DOF7=%.2f  bulk=%d  fine=%d\n', ...
        k, layers(k).z, layers(k).a_k, layers(k).b_k, layers(k).DOF7, ...
        layers(k).n_bulk, layers(k).n_chords);
end
fprintf('    TOTALS: %d bulk + %d fine = %d passes\n\n', total_bulk, total_fine, total_bulk+total_fine);

results.layers       = layers;
results.n_layers     = n_layers;
results.z_centres    = z_centres;
results.DOF7_vals    = DOF7_vals;
results.DOF7_max     = DOF7_max;
results.r1p = r1p; results.r2p = r2p; results.r3p = r3p;
results.V_resection  = V_resection;
results.V_excess     = V_excess;
results.pct_excess   = pct_excess;
results.V_aabb       = V_aabb;
results.pct_excess_aabb = pct_excess_aabb;
results.volume_saved_pct = volume_saved_pct;
results.total_bulk   = total_bulk;
results.total_fine   = total_fine;
results.n_cuts       = n_layers;
results.n_reorient   = 1;
end

function results = plan_obb(results, score, params)
m      = params.margin;
w_bulk = params.w_bulk;
w_fine = params.w_fine;

centroid = results.centroid;
R_pca    = results.R_pca;
e1 = results.e1; e2 = results.e2; e3 = results.e3;

minP = min(score, [], 1) - m;
maxP = max(score, [], 1) + m;
boxSize = maxP - minP;

L1 = boxSize(1); L2 = boxSize(2); L3 = boxSize(3);
h1 = L1/2; h2 = L2/2; h3 = L3/2;

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
fprintf('    Resection volume: %.2f mm^3  (+%.1f%% over tumour)\n', V_OBB, pct_excess);

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

fprintf('\n    --- Face details ---\n');
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

    if i <= 4
        motion(i).DOF7 = h3;
    else
        motion(i).DOF7 = 0;
    end

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

    fine_pts = {p0 - hU*u_ax - hV*v_ax, p0 + hU*u_ax - hV*v_ax;
                p0 - hU*u_ax + hV*v_ax, p0 + hU*u_ax + hV*v_ax;
                p0 - hU*u_ax - hV*v_ax, p0 - hU*u_ax + hV*v_ax;
                p0 + hU*u_ax - hV*v_ax, p0 + hU*u_ax + hV*v_ax};
    fine_path = zeros(8, 3);
    for fi = 1:4
        fine_path(2*fi-1, :) = fine_pts{fi, 1};
        fine_path(2*fi,   :) = fine_pts{fi, 2};
    end

    motion(i).bulk_path = bulk_path;
    motion(i).fine_path = fine_path;
    motion(i).n_bulk    = length(y_vals);
    motion(i).n_fine    = 4;

    fprintf('    Face %s: dims=%.1fx%.1f  DOF7=%.2f  bulk=%d  fine=%d\n', ...
        planes(i).name, planes(i).uLen, planes(i).vLen, motion(i).DOF7, ...
        motion(i).n_bulk, motion(i).n_fine);

    total_bulk = total_bulk + length(y_vals);
    total_fine = total_fine + 4;
end
fprintf('    TOTALS: %d bulk + %d fine = %d passes\n\n', total_bulk, total_fine, total_bulk+total_fine);

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

function area = compute_surface_area(k2, X, Y, Z)
area = 0;
for i = 1:size(k2, 1)
    v1 = [X(k2(i,1)), Y(k2(i,1)), Z(k2(i,1))];
    v2 = [X(k2(i,2)), Y(k2(i,2)), Z(k2(i,2))];
    v3 = [X(k2(i,3)), Y(k2(i,3)), Z(k2(i,3))];
    area = area + 0.5 * norm(cross(v2-v1, v3-v1));
end
end

function out = ternary(cond, a, b)
if cond; out = a; else; out = b; end
end
