# Task 3: Surgical Planning -- Report Plan

## Constraints
- 5 pages including all images and figures
- Font: 11pt Calibri or Arial, justified
- Figures: captions, axes labels, legends, distinguishable in B&W
- IEEE citation style
- Marking criteria: "Clear explanation of how parameterization of cuts are done, supported by mathematics" (12%)

---

## PAGE 1 -- Algorithm Overview and Problem Statement

### Opening paragraph (~6 lines)
State the problem: given a set of tumour surface vertices (VerticesUnique) from
pre-operative imaging, the surgical planner must compute all cut parameters
for the 8-DOF tool designed in Task 2. This includes:
- Position and orientation of the tool (DOF 1-6 via robot)
- Depth of each cut layer (DOF 7 via linear slide)
- Tip selection per pass -- bulk serrated or fine smooth (DOF 8 via selector)

Mention that the approach is adaptive: it analyses tumour geometry to select
between ellipsoid slicing (for compact tumours) and oriented bounding box
resection (for elongated tumours), minimising healthy bone removal.

### Flowchart (Figure 5 in report, ~40% of page)
Full-width flowchart showing the pipeline:

```
VerticesUnique (Nx3)
       |
       v
  Convex Hull
  (volume, surface area)
       |
       v
     PCA
  (centroid, e1/e2/e3, elongation E)
       |
       v
    MVEE
  (Khachiyan algorithm)
  (semi-axes r1>=r2>=r3, fit ratio F)
       |
       v
  Decision Logic
  F < 2.0 AND E < 5.0?
      / \
    YES   NO
     |     |
     v     v
  Ellipsoid   OBB
  Slicing   (6 planes)
     |     |
     v     v
  Layer params   Face params
  (z_k, a_k, b_k)  (normal, point, dims)
     |     |
     v     v
  Per-layer raster  Per-face raster
  + boundary chords + boundary edges
     \     /
      \   /
       v v
  DOF 1-6: tool pose (world frame waypoints)
  DOF 7:   depth per layer/plane
  DOF 8:   bulk (interior) / fine (boundary)
```

### Decision thresholds table (~15% of page)

| Parameter | Symbol | Value | Justification |
|-----------|--------|-------|---------------|
| Safety margin | m | 5 mm | Consistent with Task 2 swept-volume analysis |
| Bulk tip width | w_bulk | 7.0 mm | From Task 2 serrated tip dimension |
| Fine tip width | w_fine | 5.5 mm | From Task 2 smooth tip dimension |
| Fit ratio threshold | F_thresh | 2.0 | If MVEE volume > 2x tumour, ellipsoid too conservative |
| Elongation threshold | E_thresh | 5.0 | lambda1/lambda3 > 5 indicates elongated/flat shape |
| DOF 7 stroke limit | d_max | 50 mm | Physical constraint from Task 2 slide mechanism |

---

## PAGE 2 -- Mathematical Derivation

This page is equation-dense with minimal prose. Each equation block gets
a 1-2 sentence introduction. Number all equations sequentially.

### 2.1 Centroid and PCA (~20% of page)

Centroid:
  x_bar = (1/N) * sum_{i=1}^{N} v_i                              (1)

Mean-centred data matrix:
  X_0 = V - 1*x_bar^T    (N x 3)                                 (2)

SVD / PCA:
  X_0 = U * Sigma * R_pca^T                                       (3)
  R_pca = [e1 | e2 | e3]   (columns = principal axes)

Elongation metric:
  E = lambda_1 / lambda_3                                          (4)
  where lambda_k are the eigenvalues (variances) from PCA

### 2.2 Minimum-Volume Enclosing Ellipsoid (~35% of page)

Optimisation problem:
  min_{A,c}  -log det(A)                                           (5)
  s.t.  (v_i - c)^T A (v_i - c) <= 1   for all i = 1..N

Lift data to (d+1)-dimensional space to handle free centre:
  Q = [P; 1^T]    ((d+1) x N, appending row of ones)              (6)

Khachiyan iterative update (present as Algorithm 1 box):
  1. Initialise: u_i = 1/N for all i
  2. Compute lifted scatter: X = Q * diag(u) * Q^T                 (7)
  3. Mahalanobis distances: M_i = Q_i^T * X^{-1} * Q_i             (8)
  4. Find j = argmax(M_i)
  5. Step size: alpha = (M_j - d - 1) / ((d+1)(M_j - 1))          (9)
  6. Update: u <- (1-alpha)*u;  u_j <- u_j + alpha
  7. Converge when M_j - d - 1 < tol*(d+1)

Ellipsoid centre (weighted mean of points):
  c = P * u                                                        (10)

Centred weighted scatter and ellipsoid matrix:
  S = P * diag(u) * P^T - c * c^T                                 (11)
  A = S^{-1} / d                                                   (12)

Semi-axis extraction via eigendecomposition:
  A = Q * Lambda * Q^T                                             (13)
  r_k = lambda_k^{-1/2},  k = 1,2,3                               (14)
  with r_1 >= r_2 >= r_3

Volume and fit ratio:
  V_MVEE = (4/3) * pi * r_1 * r_2 * r_3                           (15)
  F = V_MVEE / V_tumour                                            (16)

### 2.3 Adaptive Selection Rule (~5% of page)

  if F < F_thresh AND E < E_thresh:  Ellipsoid Slicing             (17)
  else:  OBB (6-plane resection)

### 2.4 Ellipsoid Slicing Parameterisation (~25% of page)

Margin expansion:
  r'_k = r_k + m                                                   (18)

Layer count and depths:
  n_layers = ceil(2*r'_3 / w_bulk)                                 (19)
  z_k = -r'_3 + (k - 0.5) * w_bulk,  k = 1..n_layers             (20)

Cross-section semi-axes at depth z_k:
  s_k = sqrt(1 - (z_k / r'_3)^2)                                  (21)
  a_k = r'_1 * s_k,   b_k = r'_2 * s_k                           (22)

DOF 7 depth per layer:
  d_k = r'_3 + z_k                                                 (23)

Interior raster: boustrophedon (serpentine) passes parallel to e1,
stepping along e2 with spacing w_bulk, excluding outermost w_fine strip.

Boundary chord count (Ramanujan perimeter approximation):
  h = ((a_k - b_k)/(a_k + b_k))^2                                 (24)
  P_k = pi*(a_k + b_k)*(1 + 3h/(10 + sqrt(4 - 3h)))              (25)
  n_chords_k = ceil(P_k / w_fine)                                  (26)

World-frame waypoint transformation:
  p_world = c + Q * [a_k*cos(theta_j);  b_k*sin(theta_j);  z_k]  (27)

### 2.5 OBB Parameterisation (~15% of page)

Box dimensions from PCA projections:
  L_k = (max_i(v_i . e_k) - min_i(v_i . e_k)) + 2m               (28)
  V_OBB = L_1 * L_2 * L_3                                         (29)

6 cutting planes defined by:
  normal n_f, point p_f, in-plane axes u_f / v_f, dimensions

---

## PAGE 3 -- Case A Results (Ellipsoid Slicing)

### Figure 6: Tumour + MVEE + Resection Ellipsoid + PCA Axes (~45% of page)
- Source: Fig 1a from Task3.m (Case A separate figure window)
- Shows: red tumour vertices, cyan transparent MVEE, blue transparent
  margin-expanded ellipsoid, PCA axes (red/green/blue quivers), black centroid
- Caption: "Fig. 6. Case A tumour surface vertices with MVEE (cyan) and
  margin-expanded resection ellipsoid (blue, +5mm). PCA principal axes shown
  from centroid."

### Results paragraph (~15% of page)
Report the numerical outputs from Case A:
- Tumour volume, surface area
- MVEE semi-axes [r1, r2, r3]
- Fit ratio F (expected < 2.0)
- Elongation E (expected < 5.0)
- Method selected: Ellipsoid Slicing
- Number of layers, DOF7 max, feasibility check
- Resection volume and excess percentage

### Figure 7: Three Representative Layer Cross-Sections (~40% of page)
- Source: Fig 3 from Task3.m (3 subplots side by side)
- Shows: top/middle/bottom layers with ellipse boundary (black dashed),
  bulk raster lines (green), fine boundary chords (orange)
- Each subplot annotated with z_k, a_k, b_k, DOF7, pass counts
- Caption: "Fig. 7. Representative cross-sections at three depths. Green
  lines: bulk tip interior raster. Orange segments: fine tip boundary
  polygon. Black dashed: elliptical cross-section boundary."

---

## PAGE 4 -- Case B Results (OBB) and Comparison

### Figure 8: Tumour + OBB (~30% of page)
- Source: Fig 1b from Task3.m (Case B separate figure window)
- Shows: red tumour vertices (elongated), blue OBB wireframe edges,
  green raster lines on +e3 face, PCA axes, centroid
- Caption: "Fig. 8. Case B elongated tumour (a=30, b=5, c=5 mm) with
  PCA-aligned oriented bounding box. Green lines show bulk raster path
  on one representative face."

### Results paragraph for Case B (~10% of page)
Same metrics as Case A but for OBB:
- F value and E value (E >= 5.0 triggers OBB for elongated tumours)
- Method selected: OBB
- 6 cutting planes, 6 robot reorientations
- DOF7 = h3 for lateral faces, 0 for top/bottom
- Resection volume and excess percentage

### Figure 9: Method Comparison Bar Charts (~35% of page)
- Source: Fig 4 from Task3.m (2x2 subplot)
- Shows: cutting operations, robot reorientations, excess volume %,
  DOF7 max extension -- side by side for Case A vs Case B
- Caption: "Fig. 9. Comparison of ellipsoid slicing (Case A) versus OBB
  (Case B) across four surgical planning metrics."

### Comparison discussion (~25% of page)
Key points to make:
- Ellipsoid slicing requires only 1 robot reorientation vs 6 for OBB
  (faster procedure, less registration error accumulation)
- Ellipsoid slicing produces lower excess volume for compact tumours
  (more conformal to actual tumour shape)
- OBB is necessary for elongated tumours where the extreme axis ratio
  makes ellipsoid cross-sections impractically elongated for cutting
  (E exceeds threshold, even though F may remain low)
- The adaptive selection automatically picks the better method based
  on quantitative shape metrics (F and E), not manual judgement

---

## PAGE 5 -- DOF Integration and Conclusion

### DOF assignment table (~25% of page)

| Surgical Parameter | DOF | Source | Value |
|--------------------|-----|--------|-------|
| Tool position (x,y,z) | DOF 1-3 | Centroid + Q_mvee orientation | World-frame waypoints |
| Tool orientation (roll,pitch,yaw) | DOF 4-6 | Q_mvee columns / plane normals | Aligned to resection geometry |
| Cutting depth | DOF 7 | d_k = r'_3 + z_k (ellipsoid) or h_3 (OBB) | Per-layer/per-face extension |
| Tip selection | DOF 8 | Interior pass -> bulk; boundary pass -> fine | Binary selection per pass |

### Cutting sequence description (~20% of page)
For ellipsoid slicing (Case A), describe the execution order:
1. Robot positions tool at approach point along depth axis (e3 direction),
   standoff = 30mm from entry face
2. For each layer k = 1 to n_layers:
   a. DOF 7 advances to d_k
   b. Execute all bulk tip interior raster passes
   c. Switch to fine tip (DOF 8)
   d. Execute all boundary chord passes
   e. DOF 7 retracts
3. Single approach direction for all layers (1 reorientation total)

For OBB (Case B), note 6 separate approach directions (one per face),
each with its own raster + boundary passes.

### Summary table of both cases (~25% of page)

Print the key rows from the MATLAB summary output:

|  | Case A | Case B |
|--|--------|--------|
| Method | Ellipsoid | OBB |
| Tumour volume | xxx mm^3 | xxx mm^3 |
| Resection volume | xxx mm^3 | xxx mm^3 |
| Excess volume | xx.x% | xx.x% |
| Cutting operations | n_layers | 6 |
| Robot reorientations | 1 | 6 |
| DOF7 max | xx.x mm | xx.x mm |
| DOF7 feasible | YES/NO | YES/NO |

(Fill actual values after running Task3.m)

### Limitations paragraph (~15% of page)
- Semi-axis expansion (r'_k = r_k + m) is not a true parallel offset
  surface; the physical margin varies slightly due to ellipsoid curvature.
  A true offset would require numerical computation beyond scope.
- The MVEE is computed on convex hull vertices only, assuming the tumour
  boundary is convex. Non-convex tumours would require concave hull or
  alpha-shape preprocessing.
- The tool path assumes rigid bone fixation (no compliance or deflection
  during cutting).
- The elongation metric E = lambda1/lambda3 cannot distinguish needle-like
  shapes (one large eigenvalue) from flat disc shapes (one small eigenvalue);
  both produce high E values. A more refined metric using lambda1/lambda2
  and lambda2/lambda3 separately could differentiate these cases.

### Closing statement (~15% of page)
The adaptive surgical planner automatically selects the resection geometry
that best matches the tumour morphology, parameterises all 8 DOF of the
tool system from Task 2, and generates complete robot waypoints in the
world frame. The ellipsoid slicing method minimises healthy bone removal
for compact tumours while the OBB method handles elongated or irregular
geometries where a tight ellipsoid fit is not achievable.

---

## Figures Summary

| Report Fig | Source | Content |
|------------|--------|---------|
| Fig 5 | Hand-drawn or generated | Algorithm flowchart |
| Fig 6 | Task3.m Fig 1a | Case A: tumour + MVEE + resection ellipsoid |
| Fig 7 | Task3.m Fig 3 | Case A: 3 layer cross-sections |
| Fig 8 | Task3.m Fig 1b | Case B: tumour + OBB wireframe |
| Fig 9 | Task3.m Fig 4 | 2x2 comparison bar charts |

Total: 5 figures across 5 pages (1 per page average).
Figure numbering continues from Task 2 (assumes Task 2 used Figs 1-4).

---

## Equations Summary

29 numbered equations across Page 2, covering:
- PCA: equations 1-4
- MVEE: equations 5-16 (including lifted coordinates and centred scatter)
- Decision: equation 17
- Ellipsoid slicing: equations 18-27 (including boustrophedon raster)
- OBB: equations 28-29

---

## Writing Checklist

- [ ] All figure numbers referenced in text before they appear
- [ ] All equation numbers referenced in text
- [ ] IEEE citations for: PCA, Khachiyan/MVEE algorithm, Ramanujan approximation
- [ ] Cross-references to Task 1 (technology choice) and Task 2 (tool design, DOF)
- [ ] Numerical values in tables match Task3.m console output exactly
- [ ] No identifying information (anonymous marking)
- [ ] All figures have captions below them
- [ ] All axes labelled with units (mm, mm^2, mm^3)
- [ ] Curves distinguishable in black and white (line styles vary)
