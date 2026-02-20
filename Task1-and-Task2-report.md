Task 1: Literature Survey

Surgical bone tumour resection requires high-precision cutting tools capable of removing diseased tissue while minimising damage to healthy bone and adjacent soft tissue [1]. This review evaluates four cutting technologies—mechanical burrs, ultrasonic cutting (piezosurgery), laser osteotomy, and waterjet cutting—for automated conformal bone resection across five criteria: conformal cutting capability, biological safety, robotic integration readiness, technical limitations, and practical deployment feasibility [2].

Figure 1-1 – Original schematic comparison of four osteotomy mechanisms. (a) Mechanical burr (b) Piezosurgery (c) Laser (d) Waterjet. Created using BioRender.

Mechanical Cutting Tools (High-Speed Burrs)
Mechanism

Mechanical cutting represents the oldest and most widely adopted approach to bone resection, serving as the clinical baseline [3]. Material removal occurs when contact stresses exceed bone fracture strength, producing chips evacuated by irrigation. While oscillating saws dominate fast planar osteotomies and milling cutters represent high-precision extremes, high-speed burrs occupy the most clinically relevant regime through continuous grinding and shearing, enabling curved, concave, and internal cuts not achievable with planar saws [4], [5].

Surgical Performance & Automation

Burrs provide moderate cutting speed with high geometric flexibility, though cutting forces remain sensitive to feed rate and require steady motion [6]. Thermal damage stays localised with irrigation and continuous motion, but excessive force rapidly elevates temperature. High debris generation introduces dispersion and tumour cell dissemination risks without effective suction [7], [8], [9].

Burrs mount easily to robotic end-effectors due to compact size and simple actuation, and robotic bone milling is well documented experimentally [10], [11], [12]. However, no large-scale clinical trials of fully autonomous burr-based tumour resection exist, highlighting translational and regulatory barriers.

Tool wear degrades accuracy, and performance varies between cortical and cancellous bone, necessitating adaptive control [13], [14]. Despite recoverable failure modes (thermal buildup, chatter, soft tissue contact), burrs remain widely deployed, inexpensive, and sterilization compatible [15]–[19].

Ultrasonic Bone Cutting (Piezosurgery)
Mechanism

High-frequency vibrations (25–35 kHz) selectively cut mineralised tissue through micro-fracture while soft tissues deform elastically without damage [20], [21].

Surgical Performance & Automation

This selective mechanism is the principal clinical advantage, particularly near neurovascular structures where experimental studies demonstrate preserved soft tissue integrity even under direct contact, contrasting sharply with rotary instruments [22].

Thermal performance is comparable to conventional tools with adequate irrigation (≥30 ml/min) and low applied force (1.5–2.0 N) [20], [23], though excessive loading (>900 g) can raise temperatures above 48°C, risking osteonecrosis [24].

Cutting speed is 28–65% slower than rotary instruments [25], [26], creating a throughput limitation for large-volume resections.

Instruments are compact and readily mounted as robotic end-effectors. The stable force profile (~1.5 N) offers potential for force-based boundary detection, though fully autonomous implementations remain limited [22], [27].

Cutting performance is strongest in dense cortical bone and degrades in cancellous regions [28], necessitating feedback control. The dominant failure mode—tip breakage—is recoverable, and selective cutting provides a biological safety margin [29]. Commercial systems exhibit high technology readiness (TRL 8–9) with moderate costs [30], [31].

Clinically, piezosurgery is most established in dental and maxillofacial surgery [32], [33]. Limited large orthopaedic adoption reflects throughput constraints rather than biological unsuitability [34], [35]. When scaled through adapted geometries and robotic feed control, these properties can suit large osteotomies for margin control and critical-structure protection.

Laser Osteotomy
Mechanism

Laser osteotomy employs optical energy absorption and rapid heating for micro-explosive material removal without physical contact [36].

Surgical Performance & Automation

Robotic laser systems achieve narrow kerf widths and accurate trajectory control, outperforming conventional tools in geometric deviation metrics [37].

However, ablation rate and resection time remain limiting for thick cortical bone due to thermal constraints, and debris production introduces dispersion challenges in oncologic settings [38].

The non-contact nature aligns well with robotic path following and image-derived trajectories, though line-of-sight beam delivery restricts access in confined anatomy.

Bone optical and thermal properties vary across cortical and cancellous regions, altering ablation rates unless compensated [39].

While lasers avoid blade wear, optical contamination affects repeatability. Major safety risks include thermal injury and overcut; OCT-gated approaches mitigate this by stopping energy at detected boundaries.

Laser osteotomy for complex orthopaedic tasks is emerging rather than routine—broader adoption is limited by hardware cost, OR integration complexity, and workflow demands [40], [41]. Current evidence supports laser osteotomy as a precision end-effector rather than a universal replacement [42].

Waterjet Cutting
Mechanism

Waterjet cutting employs high-pressure fluid to fracture tissue. For bone tumour resection, pulsed waterjets are most relevant, offering controlled energy delivery with reduced collateral damage [43], [44].

Surgical Performance & Automation

Pulsed jets have been experimentally validated for bone but lack extensive clinical deployment in oncologic resection [45].

The technology demonstrates robotic compatibility [46], [47], yet significant constraints exist:

Large sterile saline volumes required

Mechanical complexity of pulsed jet systems

Residual fluid management challenges [48]–[50]

While waterjets avoid high-temperature thermal damage, pressure refinement is critical to prevent excessive bone boundary damage [51].

Clinical Constraints for Conformal Resection

Automated bone tumour resection must satisfy three critical constraints:

Peak bone temperature below 47°C

Kerf width under 5 mm

En-bloc resection prohibiting tumour fragmentation

[53]–[56]

These constraints prioritise technologies with predictable thermal behaviour and stable robotic execution.

Comparative Analysis & Technology Selection
Weighted Decision Matrix
Criteria	Weight	Piezosurgery	Laser	Burrs	Waterjet
Conformal Cutting	5	20	25	20	15
Biological Safety	7	35	21	28	21
Robotic Integration	5	20	20	20	15
Practical Feasibility	3	9	6	9	6
TOTAL		84	67	77	62

Table 1 – Weighted decision matrix.

Piezosurgery achieves the highest overall score due to exceptional biological safety. Hybrid systems remain underexplored [56], [57].

Task 2: Conformal Tool Design
Design Rationale

Piezosurgery was selected based on its soft-tissue selectivity and robotic readiness. However, ultrasonic tips operate 28–65% slower than rotary instruments [58], [59].

The proposed solution introduces:

Dual-tip architecture

Motorised sliding depth mechanism

Co-axial irrigation and suction

The system operates as an 8-DOF active end-effector (6 robot DOF + depth slide + tip selector).

System Architecture

Figure 2-1 – 8-DOF conformal resection system.

DOF Allocation
DOF	Description	Actuator	Control Mode
1–6	Robot joints	Servo motors	Pre-planned trajectory
7	Depth slide	Linear actuator	Closed-loop position
8	Tip selector	Linear actuator	Position-confirmed switching

Table 2-1 – DOF allocation.

Mechanical Design
Dual Tip Module

Serrated bulk-removal tip (7 mm)

Smooth fine-boundary tip (5.5 mm)

Material: Ti-6Al-4V

Maximum stroke: 50 mm

Sliding Depth Mechanism

Telescoping linear slide

Lead screw with encoder feedback

Mechanical hard stop (independent of software)

Irrigation & Suction Integration

Minimum flow rate: 30 ml/min

Annular irrigation channel

Co-axial suction lumen

Conformality Analysis

The conformal volume reduces unnecessary healthy bone removal compared to a bounding box approach.

𝑋
%
=
(
𝑉
𝐵
−
𝑉
𝑇
)
−
(
𝑉
𝐶
−
𝑉
𝑇
)
𝑉
𝐵
−
𝑉
𝑇
×
100
%
X%=
V
B
	​

−V
T
	​

(V
B
	​

−V
T
	​

)−(V
C
	​

−V
T
	​

)
	​

×100%

Where:

𝑉
𝐵
V
B
	​

 = baseline bounding box volume

𝑉
𝐶
V
C
	​

 = conformal resection volume

𝑉
𝑇
V
T
	​

 = tumour volume

Reduction achieved: 26.7%

Structural & Thermal Feasibility
Structural FEA

Load: 2 N lateral

Peak stress: 0.647 MPa

Safety factor: 15

Max deflection: 
1.66
×
10
−
4
1.66×10
−4
 mm

Thermal Analysis

Heat source: 0.5 W

Peak temperature: 46.1°C

Below 47°C threshold

Failure Modes & Mitigation
Failure Mode	Mitigation
Tip fracture	Force monitoring + shutdown
Over-penetration	Mechanical hard stop
Thermal damage	Flow sensor interlock
Tip mis-selection	Position-confirmed activation
Collision	Pre-op simulation + force halt
Implementation Feasibility

Biocompatible materials

Steam sterilisation compatible

CNC manufacturable

Compatible with approved ultrasonic generators

Summary

The proposed 8-DOF piezosurgical end-effector resolves throughput limitations while preserving biological safety. Structural and thermal validation confirm feasibility. RoboDK simulation verifies conformal access across all bone_5 cases, supporting clinical translation potential.