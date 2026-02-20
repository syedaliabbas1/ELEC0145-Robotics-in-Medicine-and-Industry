    Task 3: Surgical Planning

    Algorithm Overview

    Given a set of tumour surface vertices (VerticesUnique) obtained from pre-operative imaging data—assuming segmentation has already been performed and only tumour data is provided (TumourData.zip)—the surgical planner must compute all cut parameters for the 8-DOF piezosurgical tool designed in Task 2. These parameters comprise the position and orientation of the tool (DOF 1–6 via the 6-DOF robot), the depth of each cut layer (DOF 7 via the linear slide), and the tip selection per pass -- the serrated bulk-removal tip (7.0 mm, Table 2-1) for interior material excavation or the smooth fine-boundary tip (5.5 mm, Table 2-1) for conformal margin refinement (DOF 8 via the dual-tip selector). The approach is adaptive: it analyses tumour geometry using principal component analysis (PCA) [67] and a minimum-volume enclosing ellipsoid (MVEE) [68] to automatically select between ellipsoid slicing for compact tumours and oriented bounding box (OBB) resection for elongated tumours, thereby minimising healthy bone removal consistent with the conformal resection objective established in Task 1. The complete pipeline is shown in Figure 3-1, and the design parameters are summarised in Table 3-1.

    Figure 3-1 – Adaptive surgical planning pipeline. Tumour vertices are analysed via convex hull, PCA, and MVEE to compute shape metrics F (fit ratio) and E (elongation). Decision logic selects ellipsoid slicing (F < 2.0 and E < 5.0) or OBB resection (otherwise). Both paths produce complete DOF 1–8 parameterisation: tool pose waypoints, depth schedule, and tip selection sequence.

    [INSERT FIGURE 3-1: Flowchart showing: VerticesUnique (Nx3) -> Convex Hull (volume, surface area) -> PCA (centroid, e1/e2/e3, elongation E) -> MVEE via Khachiyan algorithm (semi-axes r1>=r2>=r3, fit ratio F) -> Decision Logic (F < 2.0 AND E < 5.0?) -> YES branch: Ellipsoid Slicing -> Layer params (z_k, a_k, b_k) -> Per-layer boustrophedon raster + boundary chords / NO branch: OBB (6 planes) -> Face params (normal, point, dims) -> Per-face raster + boundary edges -> Both converge to: DOF 1-6 tool pose, DOF 7 depth per layer/plane, DOF 8 bulk (interior) / fine (boundary)]

    Table 3-1 lists the design parameters and decision thresholds used throughout the planning algorithm. All tool-specific values are inherited directly from the Task 2 mechanical design.

    | Parameter | Symbol | Value | Justification |
    |---|---|---|---|
    | Safety margin | m | 5 mm | Consistent with Task 2 swept-volume analysis (Section: Conformality and Access) |
    | Bulk tip width | w_bulk | 7.0 mm | Serrated tip dimension from Task 2 (Table 2-1) |
    | Fine tip width | w_fine | 5.5 mm | Smooth tip dimension from Task 2 (Table 2-1) |
    | Fit ratio threshold | F_thresh | 2.0 | If MVEE volume exceeds 2x tumour volume, ellipsoid approximation is too conservative |
    | Elongation threshold | E_thresh | 5.0 | PCA eigenvalue ratio lambda_1/lambda_3 > 5 indicates needle-like or flat tumour geometry |
    | DOF 7 stroke limit | d_max | 50 mm | Physical constraint from Task 2 linear slide mechanism |

    Table 3-1 – Design parameters and decision thresholds for the adaptive surgical planner.

    Mathematical Derivation

    Centroid and Principal Component Analysis

    The tumour centroid is computed as the arithmetic mean of the N surface vertices:

    x_bar = (1/N) * sum_{i=1}^{N} v_i     (2)

    The mean-centred data matrix is formed by subtracting the centroid from each vertex:

    X_0 = V - 1 * x_bar^T     (N x 3)     (3)

    Singular value decomposition of the centred data yields the principal axes [67]:

    X_0 = U * Sigma * R_pca^T     (4)

    where R_pca = [e_1 | e_2 | e_3] contains the principal axes as columns, ordered by decreasing variance. The eigenvalues (variances) lambda_1 >= lambda_2 >= lambda_3 define the elongation metric:

    E = lambda_1 / lambda_3     (5)

    A high E indicates a strongly anisotropic tumour shape, where either one axis is much longer (needle-like) or one axis is much shorter (flat disc) than the others.

    Minimum-Volume Enclosing Ellipsoid

    The MVEE is the smallest ellipsoid that contains all surface vertices, obtained by solving the convex optimisation problem [68]:

    min_{A,c}  -log det(A)     s.t.  (v_i - c)^T A (v_i - c) <= 1   for all i = 1..N     (6)

    To handle the free centre variable, the data is lifted to (d+1)-dimensional space by appending a row of ones to the data matrix P (3 x N):

    Q = [P; 1^T]     ((d+1) x N)     (7)

    The Khachiyan iterative algorithm [68] solves this problem through the following update procedure. Uniform weights u_i = 1/N are initialised for all points. At each iteration, the lifted scatter matrix is computed as:

    X = Q * diag(u) * Q^T     (8)

    Mahalanobis distances are evaluated for all points:

    M_i = Q_i^T * X^{-1} * Q_i     (9)

    The maximally violating point j = argmax(M_i) is identified, and the step size is computed as:

    alpha = (M_j - d - 1) / ((d+1)(M_j - 1))     (10)

    Weights are updated: u <- (1 - alpha) * u, u_j <- u_j + alpha. The algorithm converges when M_j - d - 1 < tol * (d+1), where tol = 10^{-6}.

    The ellipsoid centre is the weighted mean of the data points:

    c = P * u     (11)

    The centred weighted scatter matrix and ellipsoid shape matrix are computed as:

    S = P * diag(u) * P^T - c * c^T     (12)

    A = S^{-1} / d     (13)

    Semi-axes are extracted via eigendecomposition of A:

    A = Q_mvee * Lambda * Q_mvee^T     (14)

    r_k = lambda_k^{-1/2},   k = 1, 2, 3     with r_1 >= r_2 >= r_3     (15)

    where Q_mvee contains the ellipsoid orientation axes as columns. The MVEE volume and fit ratio are:

    V_MVEE = (4/3) * pi * r_1 * r_2 * r_3     (16)

    F = V_MVEE / V_tumour     (17)

    Adaptive Method Selection

    The resection method is selected based on the quantitative shape metrics F and E:

    if F < F_thresh AND E < E_thresh:    Ellipsoid Slicing     (18)
    else:    OBB (6-plane resection)

    If the MVEE fits the tumour tightly (low F) and the tumour shape is compact (low E), ellipsoid slicing produces the most conformal resection. Otherwise, the extreme axis ratio makes ellipsoid cross-sections impractically thin, and OBB resection provides more reliable planar cuts.

    Ellipsoid Slicing Parameterisation

    Each MVEE semi-axis is expanded by the safety margin m to form the resection ellipsoid:

    r'_k = r_k + m     (19)

    This semi-axis expansion is not a true parallel offset surface; the physical margin equals exactly m at the axis tips but varies slightly elsewhere due to ellipsoid curvature. This is acknowledged as a limitation, though the variation is small for the geometries considered here.

    The number of cutting layers along the shortest axis (the r_3 direction, which maps to the DOF 7 depth axis) and their centre depths are:

    n_layers = ceil(2 * r'_3 / w_bulk)     (20)

    z_k = -r'_3 + (k - 0.5) * w_bulk,   k = 1 .. n_layers     (21)

    At each depth z_k, the elliptical cross-section semi-axes are:

    s_k = sqrt(1 - (z_k / r'_3)^2)     (22)

    a_k = r'_1 * s_k,   b_k = r'_2 * s_k     (23)

    DOF 7 (linear slide extension, as designed in Task 2) for each layer is the distance from the tool entry face to the layer centre:

    d_k = r'_3 + z_k     (24)

    The maximum DOF 7 extension across all layers is 2 * r'_3, which must not exceed the 50 mm stroke limit of the Task 2 slide mechanism.

    Interior material is removed using a boustrophedon (serpentine) raster pattern with the bulk serrated tip, with passes parallel to e_1 spaced by w_bulk along e_2. The outermost strip of width w_fine is excluded from bulk raster and reserved for boundary finishing with the smooth fine tip. The number of boundary chords per layer is determined using the Ramanujan ellipse perimeter approximation [69]:

    h = ((a_k - b_k) / (a_k + b_k))^2     (25)

    P_k = pi * (a_k + b_k) * (1 + 3h / (10 + sqrt(4 - 3h)))     (26)

    n_chords_k = ceil(P_k / w_fine)     (27)

    All waypoints are transformed from the MVEE local frame to the world frame using the MVEE orientation matrix Q_mvee and centre c:

    p_world = c + Q_mvee * [a_k * cos(theta_j);  b_k * sin(theta_j);  z_k]     (28)

    OBB Parameterisation

    When OBB is selected, the PCA-aligned bounding box dimensions are computed from the projections of VerticesUnique onto each principal axis, with the safety margin applied:

    L_k = (max_i(v_i . e_k) - min_i(v_i . e_k)) + 2m     (29)

    V_OBB = L_1 * L_2 * L_3     (30)

    Six cutting planes are defined, one per box face, each specified by its outward normal n_f, a point on the face p_f, two in-plane raster axes u_f and v_f, and the face dimensions. Each face is cut with the same boustrophedon raster (bulk tip interior) and boundary edge (fine tip) strategy used in ellipsoid slicing. DOF 7 extends to h_3 = L_3/2 for the four lateral faces (perpendicular to the depth axis) and remains at zero for the top and bottom faces.

    Case A Results: Ellipsoid Slicing

    The algorithm was applied to the provided tumour data (TumourData.zip). The convex hull of the loaded X, Y, Z coordinates yielded VerticesUnique, which was processed through the full pipeline. Figure 3-2 shows the tumour surface vertices together with the fitted MVEE and the margin-expanded resection ellipsoid.

    [INSERT FIGURE 3-2: Export from Task3.m Fig 1a -- 3D scatter plot showing red tumour surface vertices, cyan transparent MVEE surface (no margin), blue transparent margin-expanded resection ellipsoid (+5 mm), PCA principal axes as red/green/blue quivers from centroid, and black centroid marker. Axes labelled X (mm), Y (mm), Z (mm).]

    Figure 3-2 – Case A tumour surface vertices with MVEE (cyan) and margin-expanded resection ellipsoid (blue, +5 mm). PCA principal axes shown from centroid. The MVEE conforms tightly to the tumour geometry, yielding a low fit ratio that triggers ellipsoid slicing.

    The MVEE semi-axes were [r_1, r_2, r_3] = [23.22, 17.95, 11.74] mm, giving a fit ratio F = 1.319 and a tumour volume of 15545.33 mm^3 with surface area 3278.47 mm^2. The PCA elongation was E = 4.498. Since both F < 2.0 and E < 5.0, the algorithm selected ellipsoid slicing as the resection method. The margin-expanded semi-axes [r'_1, r'_2, r'_3] = [28.22, 22.95, 16.74] mm produced 5 cutting layers with a maximum DOF 7 extension of 33.48 mm, well within the 50 mm stroke limit of the Task 2 slide mechanism. The total resection volume was 45428.45 mm^3, representing +192.2% excess over the tumour volume. For comparison, an axis-aligned bounding box (AABB) baseline would require 87122.56 mm^3 (+460.4%). Applying the healthy-bone savings metric from Task 2, equation (1), the ellipsoid approach reduces excess healthy bone removal by 58.3% relative to the AABB baseline, confirming the conformal advantage of the ellipsoid resection geometry.

    Figure 3-3 shows three representative layer cross-sections at the top, middle, and bottom of the ellipsoid. The interior raster pattern (green, bulk serrated tip) adapts to the varying elliptical cross-section at each depth, while the boundary polygon (orange, smooth fine tip) traces the ellipse perimeter with chord lengths not exceeding w_fine = 5.5 mm. The decreasing cross-section size toward the ellipsoid tips results in fewer raster passes and shorter boundary perimeters, reducing cutting time at those layers.

    [INSERT FIGURE 3-3: Export from Task3.m Fig 3 -- Three side-by-side subplots showing layer cross-sections. Each subplot shows: black dashed ellipse boundary, green interior raster lines (bulk tip, boustrophedon pattern), orange boundary chord segments (fine tip). Each subplot annotated with z_k, a_k, b_k, DOF7 value, number of bulk passes and boundary chords. Left: top layer. Centre: middle layer (largest cross-section). Right: bottom layer.]

    Figure 3-3 – Representative cross-sections at three depths (top, middle, bottom layers). Green lines: serrated bulk tip interior raster. Orange segments: smooth fine tip boundary chords. Black dashed: elliptical cross-section boundary. The boustrophedon pattern adapts to the ellipse dimensions at each layer depth.

    Case B Results: OBB Resection

    To demonstrate the adaptive selection on a contrasting geometry, a synthetic elongated tumour was generated as a discrete ellipsoid with semi-axes a = 30, b = 5, c = 5 mm centred at (25, 25, 25) mm. This geometry has a highly anisotropic shape that is not well-suited to ellipsoid slicing due to the extreme axis ratio. Figure 3-4 shows the tumour vertices with the PCA-aligned OBB.

    [INSERT FIGURE 3-4: Export from Task3.m Fig 1b -- 3D scatter plot showing red elongated tumour vertices, blue OBB wireframe (12 edges connecting 8 corners), green raster lines on one representative face (+e3), PCA principal axes as red/green/blue quivers from centroid, and black centroid marker. Axes labelled X (mm), Y (mm), Z (mm).]

    Figure 3-4 – Case B elongated tumour (a = 30, b = 5, c = 5 mm) with PCA-aligned oriented bounding box. Green lines show the bulk raster path on one representative face. The extreme axis ratio triggers OBB selection.

    The MVEE fit ratio was F = 1.212 and the PCA elongation was E = 16.706, which exceeds the threshold E_thresh = 5.0. The algorithm therefore selected the OBB method. The OBB dimensions were L_1 = 68.00, L_2 = 18.72, L_3 = 18.72 mm, producing 6 cutting planes with 6 robot reorientations. DOF 7 was set to h_3 = 9.36 mm for the four lateral faces and 0 for the top and bottom faces. The resection volume was 23838.43 mm^3, representing +838.0% excess over the tumour volume.

    Method Comparison

    Figure 3-5 compares the two methods across four surgical planning metrics.

    [INSERT FIGURE 3-5: Export from Task3.m Fig 4 -- 2x2 bar chart subplot. Top-left: Cutting operations (Case A: n_layers, Case B: 6). Top-right: Robot reorientations (Case A: 1, Case B: 6). Bottom-left: Excess volume %. Bottom-right: DOF 7 max extension (mm). Each bar labelled with its numerical value and the method name (ellipsoid/obb).]

    Figure 3-5 – Comparison of ellipsoid slicing (Case A) versus OBB (Case B) across four surgical planning metrics: cutting operations, robot reorientations, excess bone volume, and DOF 7 maximum extension.

    Ellipsoid slicing requires only 1 robot reorientation compared to 6 for OBB, because all layers are accessed from a single approach direction along the depth axis. This reduces procedure time and limits the accumulation of registration error from repeated tool repositioning—a practical consideration for the robot integration described in Task 2. For compact tumour geometries such as Case A, ellipsoid slicing also produces a lower excess volume percentage, since the ellipsoidal resection envelope conforms more closely to the actual tumour shape than a rectangular bounding box, directly supporting the conformal resection objective from Task 1.

    OBB resection is necessary for elongated tumours where the extreme eigenvalue ratio makes ellipsoid cross-sections impractically thin in one direction. In Case B, the elongation E exceeds the threshold even though the fit ratio F may remain low—the MVEE fits tightly around the elongated shape, but the issue is the geometric impracticality of slicing along a very elongated axis. The adaptive selection based on quantitative shape metrics (F and E) removes the need for manual geometric judgement from the surgeon, consistent with the automation rationale for robotic surgery.

    DOF Integration and Cutting Sequence

    Table 3-2 summarises how each surgical parameter maps to the 8-DOF tool system from Task 2.

    | Surgical Parameter | DOF | Source | Value |
    |---|---|---|---|
    | Tool position (x, y, z) | DOF 1–3 | Centroid + Q_mvee orientation | World-frame waypoints from Eq. (28) |
    | Tool orientation (roll, pitch, yaw) | DOF 4–6 | Q_mvee columns (ellipsoid) or plane normals (OBB) | Aligned to resection geometry |
    | Cutting depth | DOF 7 | d_k = r'_3 + z_k (ellipsoid) or h_3 (OBB) | Per-layer or per-face extension via linear slide |
    | Tip selection | DOF 8 | Interior pass: bulk serrated; boundary pass: smooth fine | Binary selection per pass via lateral carriage |

    Table 3-2 – Mapping of surgical planning parameters to the 8-DOF piezosurgical tool system (Task 2).

    For ellipsoid slicing (Case A), the cutting sequence proceeds as follows. The robot positions the tool at an approach point along the depth axis (e_3 direction), with a standoff of 30 mm from the entry face. For each layer k = 1 to n_layers: (a) DOF 7 advances to depth d_k via the closed-loop linear slide; (b) the bulk serrated tip executes all interior raster passes in a boustrophedon pattern; (c) DOF 8 switches to the smooth fine tip via the lateral selector carriage, with the position sensor confirming engagement before the ultrasonic generator is enabled (Table 2-2, failure mode: tip selector mid-transition); (d) all boundary chord passes are executed; (e) DOF 7 retracts. A single approach direction is used for all layers, requiring only 1 robot reorientation in total. The mechanical hard stop on the depth slide (Task 2) prevents over-penetration beyond the planned resection boundary regardless of control faults.

    For OBB resection (Case B), 6 separate approach directions are required (one per box face), each with its own raster and boundary passes. DOF 7 extends to h_3 for the four lateral faces and remains at zero for the top and bottom faces.

    Summary of Results

    Table 3-3 presents the key metrics for both cases.

    | Metric | Case A (Ellipsoid) | Case B (OBB) |
    |---|---|---|
    | Tumour volume (mm^3) | 15545.33 | 2541.33 |
    | Tumour surface area (mm^2) | 3278.47 | 1366.66 |
    | MVEE semi-axes r_1, r_2, r_3 (mm) | 23.22, 17.95, 11.74 | 29.80, 4.97, 4.97 |
    | MVEE fit ratio F | 1.319 | 1.212 |
    | PCA elongation E | 4.498 | 16.706 |
    | Method selected | Ellipsoid slicing | OBB |
    | Resection volume (mm^3) | 45428.45 | 23838.43 |
    | Excess volume (%) | 192.2 | 838.0 |
    | Cutting operations | 5 | 6 |
    | Robot reorientations | 1 | 6 |
    | DOF 7 max (mm) | 33.48 | 9.36 |
    | DOF 7 feasible (<=50 mm) | YES | YES |
    | Total bulk tip passes | 21 | 12 |
    | Total fine tip passes | 120 | 24 |

    Table 3-3 – Surgical planning summary for Case A and Case B.

    Limitations

    The semi-axis expansion r'_k = r_k + m is not a true parallel offset surface; the physical margin varies slightly due to ellipsoid curvature, being exactly m at the axis tips but marginally more elsewhere. A true offset surface would require numerical computation beyond the scope of this work. The MVEE is computed on convex hull vertices only, assuming the tumour boundary is convex; non-convex tumours would require alpha-shape preprocessing prior to ellipsoid fitting. The tool path assumes rigid bone fixation with no compliance or deflection during cutting, consistent with the rigid fixation assumption in the Task 2 RoboDK simulation. The elongation metric E = lambda_1/lambda_3 cannot distinguish needle-like shapes (one large eigenvalue) from flat disc shapes (one small eigenvalue), as both produce high E values; a more refined metric using lambda_1/lambda_2 and lambda_2/lambda_3 separately could differentiate these cases for further optimisation of the method selection.

    Summary

    The adaptive surgical planner automatically selects the resection geometry that best matches the tumour morphology, parameterises all 8 DOF of the piezosurgical tool designed in Task 2, and generates complete robot waypoints in the world frame. The ellipsoid slicing method minimises healthy bone removal for compact tumours while the OBB method handles elongated or irregular geometries where a tight ellipsoid fit is not geometrically practical for cutting. Both methods produce a fully determined cutting plan comprising DOF 1–6 tool poses, a DOF 7 depth schedule, and DOF 8 tip selection per pass, ready for execution by the robotic system. The biological safety constraints established in Task 1—thermal control via integrated irrigation, en-bloc integrity via progressive ultrasonic material removal, and soft-tissue selectivity inherent to piezosurgery—are preserved throughout, as the planning algorithm generates tool paths compatible with the co-axial irrigation and suction system and the force-limited operating regime of the Task 2 tool.

    [NOTE: The following references should be added to the global reference list at the end of the report, continuing from [66]:]

    [67] I. T. Jolliffe, Principal Component Analysis, 2nd ed. New York, NY, USA: Springer, 2002.

    [68] M. J. Todd and E. A. Yildirim, "On Khachiyan's algorithm for the computation of minimum-volume enclosing ellipsoids," Discrete Applied Mathematics, vol. 155, no. 13, pp. 1731-1744, 2007.

    [69] S. Ramanujan, "Modular equations and approximations to pi," Quarterly Journal of Mathematics, vol. 45, pp. 350-372, 1914.
