=============================================================
ELEC0145 — Assignment 1 — Task 3: Surgical Planning
COMPLETE IMPLEMENTATION PLAN
=============================================================

=============================================================
OVERVIEW
=============================================================

One unified MATLAB script: Task3_SurgicalPlanner.m

Two tumour cases demonstrated:
  Case A — Loaded X/Y/Z data (roughly ellipsoidal)
            -> Algorithm selects ELLIPSOID SLICING
  Case B — Synthetic elongated tumour (a=30, b=5, c=5)
            -> Algorithm selects OBB (6 planes)

Both cases run through the same code path via a reusable
function, proving the adaptive selection works on contrasting
geometries. This directly satisfies the assignment requirement
to "show how the tool is used for each individual case."

Output:
  - 8 figures (4 per case)
  - Full printed summary table
  - All robot waypoints in world frame
  - DOF 7 depth schedule per case
  - Tip selection sequence per case

=============================================================
SECTION 0 — PARAMETERS (all tunable values in one block)
=============================================================

PURPOSE:
  All numerical constants defined once at the top.
  Examiner can see all design choices immediately.
  No magic numbers buried in code.

VARIABLES TO DEFINE:
  margin              = 5;      % mm — surgical safety margin
  w_bulk              = 7.0;    % mm — serrated bulk-removal tip width
  w_fine              = 5.5;    % mm — smooth fine-boundary tip width
  standoff            = 30;     % mm — robot TCP standoff before bone
  DOF7_max_stroke     = 50;     % mm — physical depth slide limit
  khachiyan_tol       = 1e-6;   % convergence tolerance for MVEE
  khachiyan_maxiter   = 1000;   % max iterations for MVEE
  fitRatio_threshold  = 2.0;    % above this -> use OBB
  elongation_threshold = 4.0;   % above this -> use OBB
  n_chord_boundary    = [];     % computed per layer from Ramanujan

JUSTIFICATION TO STATE IN REPORT:
  margin: consistent with Task 2 swept-volume analysis
  w_bulk / w_fine: directly from Task 2 dual-tip dimensions
  fitRatio_threshold: if MVEE volume > 2x tumour, ellipsoid
    approximation is too conservative; OBB is tighter
  elongation_threshold: ratio lambda1/lambda3 > 4 indicates
    needle-like tumour where ellipsoid performs poorly

=============================================================
SECTION 1 — DATA LOADING AND CONVEX HULL
=============================================================

PURPOSE:
  Reproduce starter code exactly for Case A.
  For Case B, generate synthetic elongated tumour inline.

CASE A — IMPLEMENTATION:
  load X.mat; load Y.mat; load Z.mat;
  X = X(:); Y = Y(:); Z = Z(:);
  [k2, V_tumour] = convhull(X, Y, Z, 'Simplify', true);
  
  % Surface area from hull triangles
  for each triangle in k2:
    area_i = 0.5 * norm(cross(v2-v1, v3-v1));
  tumourSurfaceArea = sum(area_i);
  
  % Extract unique surface vertices
  Vertices = [X(k2(:,1)), Y(k2(:,1)), Z(k2(:,1)); ...x3]
  VerticesUnique = unique(Vertices, 'rows');
  V = VerticesUnique;   % Nx3 working matrix

CASE B — IMPLEMENTATION:
  % Generate elongated ellipsoid, no rotation for clarity
  a=30; b=5; c_r=5;
  xc=25; yc=25; zc=25;
  for i=1:2:50, j=1:2:50, k=1:2:50:
    if ((i-xc)^2/a^2 + (j-yc)^2/b^2 + (k-zc)^2/c_r^2) < 1:
      store point
  [k2, V_tumour] = convhull(X, Y, Z, 'Simplify', true);
  % same surface area + VerticesUnique extraction as Case A

OUTPUTS FROM SECTION 1:
  V              — Nx3 unique surface vertices
  k2             — convex hull triangulation
  V_tumour       — tumour volume (mm^3)
  tumourSurfaceArea — tumour surface area (mm^2)

=============================================================
SECTION 2 — PCA
=============================================================

PURPOSE:
  Compute tumour orientation frame.
  Shared by both methods — OBB uses axes for box alignment,
  ellipsoid uses eigenvalues for elongation check.

IMPLEMENTATION:
  centroid = mean(V, 1);             % 1x3
  V0 = V - centroid;                 % mean-centred Nx3
  [coeff, score, latent] = pca(V0);  % base MATLAB, no toolbox

  e1 = coeff(:,1)';   % longest axis
  e2 = coeff(:,2)';   % second axis
  e3 = coeff(:,3)';   % shortest axis — DEPTH AXIS for DOF 7
  R_pca = coeff;      % 3x3 rotation matrix (cols = axes)

  lambda1 = latent(1);
  lambda3 = latent(3);
  elongation = lambda1 / lambda3;

EQUATIONS FOR REPORT:
  centroid: x_bar = (1/N) * sum(v_i)
  SVD: X_0 = U * Sigma * V^T
  elongation: E = lambda_1 / lambda_3

OUTPUTS FROM SECTION 2:
  centroid       — 1x3 tumour centroid in world frame
  e1, e2, e3     — principal axes (unit vectors, row vectors)
  R_pca          — 3x3 rotation matrix
  latent         — [lambda1, lambda2, lambda3] variances
  elongation     — lambda1 / lambda3

=============================================================
SECTION 3 — MVEE VIA KHACHIYAN ALGORITHM
=============================================================

PURPOSE:
  Fit the minimum-volume enclosing ellipsoid to VerticesUnique.
  Extract semi-axes r1, r2, r3 and orientation matrix Q.
  Compute V_MVEE and fit ratio F.

ALGORITHM — TODD-YILDIRIM / KHACHIYAN ITERATIVE:

  INPUT: P = V' (3 x N matrix, each column is a point)
  
  Initialise:
    n = number of points
    d = 3 (dimension)
    u = ones(n,1) / n     % uniform weights
    
  Iterate until convergence:
    % Weighted scatter matrix
    X = P * diag(u) * P'              % 3x3
    M = diag(P' * inv(X) * P)         % nx1, Mahalanobis distances
    
    % Find point with maximum Mahalanobis distance
    [max_val, j] = max(M)
    
    % Step size
    step = (max_val - d - 1) / ((d + 1) * (max_val - 1))
    
    % Update weights
    u = (1 - step) * u
    u(j) = u(j) + step
    
    % Check convergence
    if max_val - d - 1 < khachiyan_tol * (d + 1): break

  RESULT:
    A_mvee = inv(P * diag(u) * P') / d    % ellipsoid matrix
    c_mvee = P * u                          % ellipsoid centre

  SEMI-AXIS EXTRACTION:
    [Q_mvee, Lambda_mvee] = eig(A_mvee)
    % Sort eigenvalues descending
    r_k = 1 / sqrt(lambda_k)   for k = 1,2,3
    r1 >= r2 >= r3

  VOLUME AND FIT RATIO:
    V_MVEE = (4/3) * pi * r1 * r2 * r3
    F = V_MVEE / V_tumour

IMPLEMENTATION NOTES:
  - Run on VerticesUnique (hull vertices only, not all X/Y/Z)
    This is critical: Khachiyan on thousands of interior points
    is slow. Hull vertices are typically 30-80 points — fast.
  - Q_mvee columns are the ellipsoid orientation axes
  - These may differ from PCA axes — MVEE is tighter
  - Use Q_mvee (not R_pca) for the world-frame transformation
    of ellipsoid slice waypoints

EQUATIONS FOR REPORT:
  Optimisation problem:
    min_{A,c}  -log det(A)
    s.t.  (v_i - c)^T A (v_i - c) <= 1  for all i

  Semi-axis extraction:
    A = Q * Lambda * Q^T
    r_k = lambda_k^{-1/2}

  Fit ratio:
    F = (4/3 * pi * r1 * r2 * r3) / V_tumour

OUTPUTS FROM SECTION 3:
  A_mvee         — 3x3 ellipsoid matrix
  c_mvee         — 1x3 ellipsoid centre in world frame
  Q_mvee         — 3x3 orientation matrix (cols = ellipsoid axes)
  r1, r2, r3     — semi-axes (mm), r1 >= r2 >= r3
  V_MVEE         — MVEE volume (mm^3)
  F              — fit ratio

=============================================================
SECTION 4 — ADAPTIVE METHOD SELECTION
=============================================================

PURPOSE:
  Apply decision thresholds.
  Print clear rationale to console.
  Set method flag for branching in Section 5.

IMPLEMENTATION:
  fprintf('--- ADAPTIVE METHOD SELECTION ---\n');
  fprintf('Fit ratio F = %.3f  (threshold: %.1f)\n', F, fitRatio_threshold);
  fprintf('Elongation E = %.3f  (threshold: %.1f)\n', elongation, elongation_threshold);

  if F < fitRatio_threshold && elongation < elongation_threshold
      method = 'ellipsoid';
      fprintf('SELECTED: Ellipsoid Slicing\n');
      fprintf('Reason: MVEE fits tumour tightly, shape is compact\n');
  else
      method = 'obb';
      fprintf('SELECTED: OBB (6-plane resection)\n');
      if F >= fitRatio_threshold
          fprintf('Reason: MVEE too large (F=%.2f >= %.1f)\n', F, fitRatio_threshold);
      else
          fprintf('Reason: Tumour too elongated (E=%.2f >= %.1f)\n', elongation, elongation_threshold);
      end
  end

OUTPUTS FROM SECTION 4:
  method         — string: 'ellipsoid' or 'obb'

=============================================================
SECTION 5A — ELLIPSOID SLICING [runs if method == 'ellipsoid']
=============================================================

PURPOSE:
  Compute full resection plan using margin-expanded ellipsoid.
  Generate all robot waypoints in world frame.
  Assign DOF 7 depth per layer.
  Assign tip selection per pass.

STEP 1 — MARGIN EXPANSION:
  r1p = r1 + margin;    % r' = r + 5mm for each semi-axis
  r2p = r2 + margin;
  r3p = r3 + margin;    % r3' is the DEPTH semi-axis (maps to DOF 7)
  
  V_resection = (4/3) * pi * r1p * r2p * r3p;
  V_excess = V_resection - V_tumour;
  pct_excess = (V_excess / V_tumour) * 100;

  NOTE FOR REPORT:
    Semi-axis expansion is NOT a perfectly uniform offset shell.
    The physical margin varies from ~5mm at axis tips to slightly
    more elsewhere. Acknowledged as limitation; true offset surface
    would require numerical computation beyond assignment scope.

STEP 2 — LAYER PARAMETERISATION:
  Slice along the shortest ellipsoid axis (r3p direction = depth axis).
  Layer spacing = w_bulk = 7mm (bulk tip width).
  
  n_layers = ceil(2 * r3p / w_bulk);
  
  % Layer centre depths from ellipsoid centre along depth axis
  % Symmetric about centre, covering full range [-r3p, +r3p]
  z_k = -r3p + (k - 0.5) * w_bulk   for k = 1, ..., n_layers
  
  % Clamp z_k to [-r3p, r3p] to avoid sqrt of negative
  z_k = max(min(z_k, r3p - 1e-6), -r3p + 1e-6);
  
  % Cross-section semi-axes at each layer
  scale_k = sqrt(1 - (z_k / r3p)^2);
  a_k = r1p * scale_k;
  b_k = r2p * scale_k;
  
  % DOF 7 depth for each layer
  % Tool enters from top of ellipsoid (+r3p face)
  % DOF 7 extension = distance from top surface to layer centre
  DOF7_k = r3p - z_k   (measured from entry face inward)
  
  % Feasibility check
  DOF7_max = 2 * r3p    (full depth through ellipsoid)
  if DOF7_max > DOF7_max_stroke: flag warning

STEP 3 — PER-LAYER INTERIOR RASTER (BULK TIP):
  For each layer k:
    % Raster lines parallel to e1 (r1p axis), spaced w_bulk apart
    % Covering the elliptical cross-section interior
    
    y_vals = -b_k + w_fine : w_bulk : b_k - w_fine;
    % (exclude outermost w_fine strip — reserved for boundary pass)
    
    For each y_j in y_vals:
      % x extent of ellipse at this y
      x_half = a_k * sqrt(1 - (y_j / b_k)^2);
      
      % Raster line endpoints in MVEE local frame
      pt_A_local = [-x_half;  y_j;  z_k];
      pt_B_local = [+x_half;  y_j;  z_k];
      
      % Transform to world frame
      pt_A_world = c_mvee' + Q_mvee * pt_A_local;
      pt_B_world = c_mvee' + Q_mvee * pt_B_local;
      
      Store as bulk-tip waypoint pair
      Tip = 'bulk'
      DOF7 = DOF7_k

STEP 4 — PER-LAYER BOUNDARY POLYGON (FINE TIP):
  For each layer k:
    % Approximate ellipse boundary by inscribed polygon
    % Chord length <= w_fine ensures cutting coverage
    
    % Ramanujan ellipse perimeter approximation
    h_ram = ((a_k - b_k) / (a_k + b_k))^2;
    P_k = pi * (a_k + b_k) * (1 + 3*h_ram / (10 + sqrt(4 - 3*h_ram)));
    
    % Number of chords
    n_chords_k = ceil(P_k / w_fine);
    
    % Parametric ellipse points at chord vertices
    theta = linspace(0, 2*pi, n_chords_k + 1);
    
    For each j = 1 to n_chords_k:
      % Chord endpoints in MVEE local frame
      pt_A_local = [a_k*cos(theta(j));   b_k*sin(theta(j));   z_k];
      pt_B_local = [a_k*cos(theta(j+1)); b_k*sin(theta(j+1)); z_k];
      
      % Transform to world frame
      pt_A_world = c_mvee' + Q_mvee * pt_A_local;
      pt_B_world = c_mvee' + Q_mvee * pt_B_local;
      
      Store as fine-tip waypoint pair
      Tip = 'fine'
      DOF7 = DOF7_k

STEP 5 — LAYER CUTTING SEQUENCE:
  For each layer k = 1 to n_layers:
    1. Robot moves to approach point (above ellipsoid entry face
       along depth axis direction, standoff = 30mm)
    2. DOF 7 advances to DOF7_k
    3. Execute all BULK raster waypoints for this layer
    4. Execute all FINE boundary chord waypoints for this layer
    5. DOF 7 retracts
    6. Robot advances to next layer approach

  Total robot reorientations: 1 (one approach direction for all layers)
  This is the key advantage over OBB (6) and convex hull (20-40)

OUTPUTS FROM SECTION 5A:
  n_layers           — total number of cutting layers
  layers(k).z        — depth of layer k from ellipsoid centre
  layers(k).a        — cross-section semi-axis along r1 direction
  layers(k).b        — cross-section semi-axis along r2 direction
  layers(k).DOF7     — depth slide extension for layer k (mm)
  layers(k).bulk_waypoints  — Mx3 world-frame raster waypoints
  layers(k).fine_waypoints  — Lx3 world-frame boundary waypoints
  layers(k).n_bulk_passes   — number of raster lines
  layers(k).n_chords        — number of boundary chord segments
  V_resection        — total resection volume (mm^3)
  V_excess           — excess over tumour (mm^3)
  pct_excess         — excess as percentage of tumour volume

=============================================================
SECTION 5B — OBB [runs if method == 'obb']
=============================================================

PURPOSE:
  Compute resection plan using PCA-aligned bounding box.
  Already developed in Task3_SurgicalPlanning.m (Code 1).
  Reproduced here for completeness in the unified script.

IMPLEMENTATION SUMMARY (full detail already in Code 1):

STEP 1 — BOUNDING BOX:
  Project VerticesUnique onto PCA axes (score from pca())
  minP = min(score) - margin;
  maxP = max(score) + margin;
  boxSize = maxP - minP;   % [L1, L2, L3]
  
  L1 = boxSize(1);   % along e1
  L2 = boxSize(2);   % along e2
  L3 = boxSize(3);   % along e3 — DEPTH AXIS
  
  centreWorld = centroid + R_pca * ((minP + maxP)/2)';

STEP 2 — 6 CUTTING PLANES:
  For each of 6 faces (+/-e1, +/-e2, +/-e3):
    normal = +/- e_k
    point  = centreWorld +/- h_k * e_k
    uAxis, vAxis = the two in-plane axes
    uLen, vLen   = face dimensions

STEP 3 — RASTER PATHS WITH TIP SELECTION:
  For each plane:
    Interior passes (bulk tip, 7mm spacing):
      y_vals from w_fine to face_edge - w_fine, step w_bulk
    Boundary passes (fine tip, 5.5mm):
      outermost strip each side
    Alternating pass direction (boustrophedon raster)

STEP 4 — DOF 7 DEPTH:
  For faces perpendicular to e3:  DOF7 = 0
  For faces parallel to e3:       DOF7 = h3 (half box depth)

OUTPUTS FROM SECTION 5B:
  planes(1..6)       — normal, point, uAxis, vAxis, uLen, vLen
  motion(1..6)       — approach, sweepPath, tipSequence, DOF7_depth
  V_OBB              — OBB resection volume (mm^3)
  pct_excess_OBB     — excess volume percentage

=============================================================
SECTION 6 — FIGURES
=============================================================

FIGURE 1 — Tumour + Resection Volume + PCA Axes
  subplot or full figure
  
  scatter3(V) in red — tumour surface vertices
  
  IF ellipsoid method:
    Draw margin-expanded ellipsoid surface mesh:
      [xs,ys,zs] = ellipsoid(0,0,0, r1p, r2p, r3p, 30)
      transform to world frame via Q_mvee and c_mvee
      surf(xs_world, ys_world, zs_world) transparent blue
    Draw original MVEE surface (no margin) transparent cyan
    
  IF OBB method:
    Draw OBB box edges in blue (12 edges from 8 corners)
    
  quiver3 for e1, e2, e3 from centroid (red, green, blue)
  scatter3 centroid in black
  
  Labels: 'X (mm)', 'Y (mm)', 'Z (mm)'
  Title:  'Tumour + [Method] Resection Volume + PCA Axes'
  Legend: tumour vertices, resection volume, MVEE (if shown), axes

FIGURE 2 — Slice Planes / DOF 7 Schedule
  
  IF ellipsoid method:
    For each layer k, draw the elliptical cross-section
    in the plane at depth z_k, coloured by depth (blue=shallow,
    red=deep) using a colormap
    Overlay tumour hull trisurf in transparent cyan
    colorbar labelled 'Layer depth (mm)'
    
  IF OBB method:
    Draw all 6 cutting planes as filled rectangles
    Colour each face differently
    Show normal quivers from each plane centre
    
  axis equal, grid on, view(3)
  Title: 'Cutting Layers / DOF 7 Schedule'

FIGURE 3 — Representative Layer Cross-Sections (3 layers)
  
  Only for ellipsoid method (most informative figure)
  
  3 subplots side by side:
    Layer 1 (top, small ellipse)
    Layer ceil(n_layers/2) (middle, largest ellipse)
    Layer n_layers (bottom, small ellipse)
  
  For each subplot:
    plot bulk raster lines in green
    plot fine boundary polygon in orange
    draw ellipse outline in black dashed
    annotate: a_k, b_k, DOF7_k, n_bulk, n_chords
  
  xlabel('e1 direction (mm)'), ylabel('e2 direction (mm)')
  Title per subplot: 'Layer k: z = %.1f mm, DOF7 = %.1f mm'
  
  FOR OBB METHOD: substitute with 3 representative face
  raster path subplots (one face each from e1, e2, e3 directions)

FIGURE 4 — Method Comparison (Case A vs Case B)
  
  2x2 subplot layout:
  
  Top-left: bar chart — number of cuts
    Case A ellipsoid: n_layers (typically 6-8)
    Case B OBB: 6
    ylabel 'Number of cutting operations'
    
  Top-right: bar chart — total robot reorientations
    Case A: 1
    Case B: 6
    ylabel 'Robot reorientations'
    
  Bottom-left: bar chart — resection volume (mm^3)
    Case A: V_resection_A
    Case B: V_OBB_B
    ylabel 'Resection volume (mm^3)'
    
  Bottom-right: bar chart — excess volume %
    Case A: pct_excess_A
    Case B: pct_excess_B
    ylabel 'Excess bone removed (%)'
  
  sgtitle('Method Comparison: Case A (Ellipsoid) vs Case B (OBB)')
  All bars labelled with values

ADDITIONAL FIGURES FOR CASE B (OBB):
  Replicate Fig 1 and Fig 2 for Case B geometry
  Label clearly as 'Case B — Elongated Tumour'

=============================================================
SECTION 7 — CASE B EXECUTION
=============================================================

PURPOSE:
  Run the complete algorithm on Case B (elongated tumour).
  Demonstrate OBB selection and output.
  Produce Case B versions of Figures 1 and 2.

IMPLEMENTATION:
  Wrap all of Sections 1-6 into a helper function:
  
  function results = run_surgical_planner(X, Y, Z, params)
    % Sections 1-6 inside
    % Returns struct with all outputs
  end
  
  Then main script calls:
    results_A = run_surgical_planner(X_A, Y_A, Z_A, params);
    results_B = run_surgical_planner(X_B, Y_B, Z_B, params);
  
  Figure 4 uses results_A and results_B together.

ALTERNATIVE (simpler for submission):
  If wrapping into a function is complex, simply run the code
  twice with a clear comment block separating Case A and Case B,
  with variable names suffixed _A and _B.

=============================================================
SECTION 8 — SUMMARY PRINTOUT
=============================================================

FORMAT:
  fprintf table with two columns: Case A | Case B

  =====================================================
  SURGICAL PLANNING SUMMARY
  =====================================================
  
                              Case A          Case B
                           (Ellipsoid)        (OBB)
  ---------------------------------------------------
  Tumour volume (mm^3)       XXX.X           XXX.X
  Tumour surface area (mm^2) XXX.X           XXX.X
  MVEE semi-axes (mm)     r1 r2 r3          r1 r2 r3
  Fit ratio F                X.XX            X.XX
  Elongation E               X.XX            X.XX
  Method selected          Ellipsoid          OBB
  Safety margin (mm)           5               5
  Resection volume (mm^3)    XXX.X           XXX.X
  Excess volume (%)          XX.X            XX.X
  Num cutting operations     XX               6
  Robot reorientations        1               6
  DOF7 max extension (mm)   XX.X            XX.X
  DOF7 feasible?            YES/NO          YES/NO
  Bulk tip passes            XXX             XXX
  Fine tip passes            XXX             XXX
  =====================================================

=============================================================
COMPLETE EQUATIONS LIST FOR REPORT (MATHEMATICS SECTION)
=============================================================

All equations below must appear explicitly in the report
with variable definitions. The examiner awards marks for
the mathematical derivation, not just the code.

1. CENTROID:
   x_bar = (1/N) * sum_{i=1}^{N} v_i

2. PCA VIA SVD:
   X_0 = U * Sigma * V^T    (mean-centred data matrix)
   R_pca = V                 (rotation to PCA frame)

3. ELONGATION:
   E = lambda_1 / lambda_3

4. MVEE OPTIMISATION PROBLEM:
   min_{A,c}  -log det(A)
   s.t.  (v_i - c)^T A (v_i - c) <= 1   for all i = 1..N

5. MVEE KHACHIYAN UPDATE RULE:
   M_j = p_j^T A^{-1} p_j                (Mahalanobis distance)
   step = (M_max - d - 1) / ((d+1)(M_max - 1))
   u <- (1 - step) * u;  u_j <- u_j + step
   
   where d = 3 (dimension), j = argmax(M)

6. SEMI-AXIS EXTRACTION:
   A = Q * Lambda * Q^T     (eigendecomposition)
   r_k = lambda_k^{-1/2}   for k = 1, 2, 3

7. MVEE VOLUME AND FIT RATIO:
   V_MVEE = (4/3) * pi * r1 * r2 * r3
   F = V_MVEE / V_tumour

8. FIT RATIO DECISION:
   if F < F_thresh  AND  E < E_thresh:
       use ellipsoid slicing
   else:
       use OBB

9. MARGIN EXPANSION:
   r'_k = r_k + m,   m = 5 mm

   Note: this is semi-axis expansion, not a true parallel
   offset surface. The physical margin varies between
   r'_min - r_min = 5mm (at tips) and slightly more
   elsewhere due to ellipsoid curvature.

10. LAYER DEPTH VALUES:
    n_layers = ceil(2 * r'_3 / w_bulk)
    z_k = -r'_3 + (k - 0.5) * w_bulk,   k = 1..n_layers

11. CROSS-SECTION SEMI-AXES:
    scale_k = sqrt(1 - (z_k / r'_3)^2)
    a_k = r'_1 * scale_k
    b_k = r'_2 * scale_k

12. DOF 7 DEPTH PER LAYER:
    DOF7_k = r'_3 + z_k
    (distance from tool entry face to layer centre)
    Feasibility: DOF7_max = 2 * r'_3 <= 50 mm (stroke limit)

13. RAMANUJAN ELLIPSE PERIMETER:
    h = ((a_k - b_k) / (a_k + b_k))^2
    P_k = pi * (a_k + b_k) * (1 + 3h / (10 + sqrt(4 - 3h)))

14. CHORD COUNT PER LAYER:
    n_chords_k = ceil(P_k / w_fine)

15. WORLD-FRAME WAYPOINT TRANSFORMATION:
    p_world = c_mvee + Q_mvee * [a_k*cos(theta_j);
                                  b_k*sin(theta_j);
                                  z_k             ]

16. OBB DIMENSIONS (if selected):
    L_k = (max_i(v_i . e_k) - min_i(v_i . e_k)) + 2m
    V_OBB = L1 * L2 * L3

=============================================================
REPORT STRUCTURE FOR TASK 3 (5 pages)
=============================================================

PAGE 1 — Algorithm Overview
  - One paragraph: motivation for adaptive approach
  - Full flowchart (hand-drawn or MATLAB-generated):
      VerticesUnique -> PCA -> MVEE -> F and E ->
      Decision -> [Ellipsoid slicing | OBB] ->
      DOF 1-6 trajectory + DOF 7 schedule + tip selection
  - Table: decision thresholds with justification

PAGE 2 — Mathematical Derivation
  - PCA equations (eq 1-3)
  - MVEE formulation and Khachiyan update rule (eq 4-7)
  - Decision metrics (eq 8)
  - Ellipsoid slicing equations (eq 9-15)
  - OBB equations (eq 16)
  - This page is dense with equations — no figures

PAGE 3 — Case A Results (Ellipsoid)
  - Figure 1: tumour + MVEE + margin ellipsoid + PCA axes
  - Figure 3: 3 representative layer cross-sections
  - Short paragraph: F value, E value, why ellipsoid selected
  - Table of key metrics for Case A

PAGE 4 — Case B Results (OBB) + Comparison
  - Figure 1 (Case B): tumour + OBB
  - Figure 4: side-by-side comparison bar charts
  - Short paragraph: why OBB selected for elongated tumour
  - Table of key metrics for Case B

PAGE 5 — Tool Integration and Conclusion
  - Table: surgical parameter -> DOF assignment
      centroid + MVEE orientation -> DOF 1-6 (robot pose)
      DOF7_k per layer            -> DOF 7 (depth slide)
      bulk/fine flag per pass     -> DOF 8 (tip selector)
  - Short paragraph: how this connects to Task 2 design
  - 2-3 sentence conclusion: adaptive planner selects
    optimal geometry based on tumour shape metrics,
    minimising healthy bone removal while satisfying
    all surgical constraints from Task 1 and Task 2

=============================================================
FIGURE PRODUCTION CHECKLIST
=============================================================

Before submitting, verify each figure has:
  [ ] Title
  [ ] xlabel, ylabel, zlabel (where applicable)
  [ ] Legend with descriptive labels
  [ ] axis equal (for 3D geometry plots)
  [ ] FontSize >= 10 on all text
  [ ] Curves distinguishable in black and white
      (use different line styles: -, --, -., : as well as colors)
  [ ] All annotated values correct and matching printed summary

=============================================================
CODE QUALITY CHECKLIST
=============================================================

  [ ] No variable name conflicts (c, a, b used for both
      ellipsoid radii and other things — use c_mvee, r1, r2, r3)
  [ ] No hardcoded paths — only load X.mat, Y.mat, Z.mat
  [ ] All %#ok<AGROW> suppressions in place for growing arrays
  [ ] Khachiyan loop has max iteration safeguard
  [ ] z_k clamped before sqrt to prevent complex numbers
  [ ] DOF7 feasibility check with clear warning message
  [ ] Section headers clearly labelled
  [ ] fprintf outputs match report figures exactly

