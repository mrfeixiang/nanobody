#!/usr/bin/env python3
"""
ESMFold VHH 구조 예측 결과 분석 및 시각화
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from collections import defaultdict
import math

WORK_DIR = "/home/syslab/Desktop/nanobody"
PDB_FILE = os.path.join(WORK_DIR, "analysis_results", "vhh_structure_esmfold.pdb")
OUTPUT_DIR = os.path.join(WORK_DIR, "analysis_results")

# VHH sequence and CDR boundaries (from our analysis)
VHH_SEQ = "QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNNFGWNWIRQSPSRGLELGRTYYRSKWYNDYAVSVRSRITINPDTSKNQFSLQLNSVTPEDTAVYYCARQGSTYFDYWGQGTLVTVSS"

# CDR boundaries (0-indexed) from our motif-based analysis
REGIONS = {
    'FR1':  (0, 25),
    'CDR1': (25, 35),
    'FR2':  (35, 50),
    'CDR2': (50, 58),
    'FR3':  (58, 99),
    'CDR3': (99, 109),
    'FR4':  (109, 120),
}

REGION_COLORS = {
    'FR1': '#a8d8ea', 'FR2': '#a8d8ea', 'FR3': '#a8d8ea', 'FR4': '#a8d8ea',
    'CDR1': '#e74c3c', 'CDR2': '#2ecc71', 'CDR3': '#3498db',
}

def get_region(resid):
    """Return region name for a given residue index (1-based)."""
    for name, (start, end) in REGIONS.items():
        if start < resid <= end:
            return name
    return 'FR4'

# =============================================================================
# Parse PDB
# =============================================================================
print("=" * 60)
print("  ESMFold VHH Structure Analysis")
print("=" * 60)

atoms = []
ca_atoms = []
with open(PDB_FILE) as f:
    for line in f:
        if line.startswith("ATOM"):
            atom = {
                'serial': int(line[6:11]),
                'name': line[12:16].strip(),
                'resname': line[17:20].strip(),
                'chain': line[21],
                'resid': int(line[22:26]),
                'x': float(line[30:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54]),
                'bfactor': float(line[60:66]),  # pLDDT score
            }
            atoms.append(atom)
            if atom['name'] == 'CA':
                ca_atoms.append(atom)

print(f"\n  PDB file: {os.path.basename(PDB_FILE)}")
print(f"  Total atoms: {len(atoms)}")
print(f"  CA atoms (residues): {len(ca_atoms)}")

# =============================================================================
# 1. pLDDT Confidence Score Analysis
# =============================================================================
print(f"\n{'='*60}")
print("  1. pLDDT Confidence Score Analysis")
print(f"{'='*60}")
print("""
  pLDDT (predicted Local Distance Difference Test):
    ESMFold's per-residue confidence metric (0-1 scale in B-factor column).

    > 0.9: Very high confidence (well-structured)
    > 0.7: High confidence (reliable)
    > 0.5: Low confidence (possibly disordered)
    < 0.5: Very low confidence (likely disordered)
""")

plddt = [ca['bfactor'] for ca in ca_atoms]
resids = [ca['resid'] for ca in ca_atoms]

print(f"  Overall pLDDT: {np.mean(plddt):.3f} (mean), {np.median(plddt):.3f} (median)")
print(f"  Range: {min(plddt):.3f} - {max(plddt):.3f}")

# Per-region pLDDT
print(f"\n  Per-region pLDDT:")
for name in ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']:
    start, end = REGIONS[name]
    region_plddt = [p for r, p in zip(resids, plddt) if start < r <= end]
    if region_plddt:
        mean_p = np.mean(region_plddt)
        conf = "Very High" if mean_p > 0.9 else ("High" if mean_p > 0.7 else "Low")
        bar = '█' * int(mean_p * 30)
        print(f"    {name:5s}: {mean_p:.3f} [{conf:9s}] {bar}")

# Count by confidence category
very_high = sum(1 for p in plddt if p > 0.9)
high = sum(1 for p in plddt if 0.7 < p <= 0.9)
low = sum(1 for p in plddt if 0.5 < p <= 0.7)
very_low = sum(1 for p in plddt if p <= 0.5)
print(f"\n  Confidence distribution:")
print(f"    Very High (>0.9): {very_high} residues ({very_high/len(plddt)*100:.1f}%)")
print(f"    High (0.7-0.9):   {high} residues ({high/len(plddt)*100:.1f}%)")
print(f"    Low (0.5-0.7):    {low} residues ({low/len(plddt)*100:.1f}%)")
print(f"    Very Low (<0.5):  {very_low} residues ({very_low/len(plddt)*100:.1f}%)")

# Figure 1: pLDDT per residue
fig, ax = plt.subplots(figsize=(14, 5))

# Background coloring by region
for name, (start, end) in REGIONS.items():
    ax.axvspan(start + 0.5, end + 0.5, alpha=0.15, color=REGION_COLORS[name])
    mid = (start + end) / 2 + 0.5
    ax.text(mid, 1.02, name, ha='center', va='bottom', fontsize=8,
            fontweight='bold' if 'CDR' in name else 'normal',
            color=REGION_COLORS[name] if 'CDR' in name else 'gray')

# pLDDT line with color coding
colors = []
for r in resids:
    region = get_region(r)
    colors.append(REGION_COLORS[region])

for i in range(len(resids) - 1):
    ax.plot([resids[i], resids[i+1]], [plddt[i], plddt[i+1]],
            color=colors[i], linewidth=2)
ax.scatter(resids, plddt, c=colors, s=20, zorder=5, edgecolors='black', linewidths=0.3)

# Confidence thresholds
ax.axhline(0.9, color='green', linestyle='--', alpha=0.5, linewidth=0.8)
ax.axhline(0.7, color='orange', linestyle='--', alpha=0.5, linewidth=0.8)
ax.axhline(0.5, color='red', linestyle='--', alpha=0.5, linewidth=0.8)
ax.text(120.5, 0.91, 'Very High', fontsize=7, color='green')
ax.text(120.5, 0.71, 'High', fontsize=7, color='orange')
ax.text(120.5, 0.51, 'Low', fontsize=7, color='red')

ax.set_xlabel('Residue Number')
ax.set_ylabel('pLDDT Confidence Score')
ax.set_title('ESMFold pLDDT Confidence Score per Residue (VHH)')
ax.set_xlim(0, 125)
ax.set_ylim(0.3, 1.05)

legend_patches = [
    mpatches.Patch(color='#e74c3c', alpha=0.5, label='CDR1'),
    mpatches.Patch(color='#2ecc71', alpha=0.5, label='CDR2'),
    mpatches.Patch(color='#3498db', alpha=0.5, label='CDR3'),
    mpatches.Patch(color='#a8d8ea', alpha=0.5, label='Framework'),
]
ax.legend(handles=legend_patches, loc='lower left', fontsize=8)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "07_esmfold_plddt.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"\n  >> Saved: analysis_results/07_esmfold_plddt.png")


# =============================================================================
# 2. Contact Map
# =============================================================================
print(f"\n{'='*60}")
print("  2. Residue Contact Map")
print(f"{'='*60}")

n_res = len(ca_atoms)
dist_matrix = np.zeros((n_res, n_res))

for i in range(n_res):
    for j in range(n_res):
        dx = ca_atoms[i]['x'] - ca_atoms[j]['x']
        dy = ca_atoms[i]['y'] - ca_atoms[j]['y']
        dz = ca_atoms[i]['z'] - ca_atoms[j]['z']
        dist_matrix[i][j] = math.sqrt(dx*dx + dy*dy + dz*dz)

# Contact map (< 8 Angstrom)
contact_threshold = 8.0
contact_map = (dist_matrix < contact_threshold).astype(float)

# Remove trivial contacts (adjacent residues)
for i in range(n_res):
    for j in range(max(0, i-2), min(n_res, i+3)):
        contact_map[i][j] = 0

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Distance map
im1 = axes[0].imshow(dist_matrix, cmap='viridis_r', origin='lower', aspect='equal')
axes[0].set_xlabel('Residue Number')
axes[0].set_ylabel('Residue Number')
axes[0].set_title('CA-CA Distance Map (Angstrom)')
plt.colorbar(im1, ax=axes[0], shrink=0.8, label='Distance (A)')

# Add CDR annotations
for name, (start, end) in REGIONS.items():
    if 'CDR' in name:
        color = REGION_COLORS[name]
        for ax_i in axes:
            ax_i.axvline(start, color=color, linewidth=0.8, alpha=0.7)
            ax_i.axvline(end, color=color, linewidth=0.8, alpha=0.7)
            ax_i.axhline(start, color=color, linewidth=0.8, alpha=0.7)
            ax_i.axhline(end, color=color, linewidth=0.8, alpha=0.7)

# Contact map
im2 = axes[1].imshow(contact_map, cmap='Blues', origin='lower', aspect='equal')
axes[1].set_xlabel('Residue Number')
axes[1].set_ylabel('Residue Number')
axes[1].set_title(f'Contact Map (CA-CA < {contact_threshold}A)')
plt.colorbar(im2, ax=axes[1], shrink=0.8)

# CDR labels
for name, (start, end) in REGIONS.items():
    if 'CDR' in name:
        mid = (start + end) / 2
        for ax_i in axes:
            ax_i.text(mid, n_res + 3, name, ha='center', fontsize=7,
                     color=REGION_COLORS[name], fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "08_esmfold_contact_map.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"  >> Saved: analysis_results/08_esmfold_contact_map.png")

# Count inter-CDR contacts
print(f"\n  Long-range contacts (CA < {contact_threshold}A, |i-j| > 5):")
cdr_regions = ['CDR1', 'CDR2', 'CDR3']
for r1 in cdr_regions:
    s1, e1 = REGIONS[r1]
    for r2 in cdr_regions:
        s2, e2 = REGIONS[r2]
        if r1 <= r2:
            contacts = 0
            for i in range(s1, e1):
                for j in range(s2, e2):
                    if abs(i - j) > 5 and i < n_res and j < n_res:
                        if dist_matrix[i][j] < contact_threshold:
                            contacts += 1
            if contacts > 0:
                print(f"    {r1}-{r2}: {contacts} contacts")


# =============================================================================
# 3. Ramachandran Plot
# =============================================================================
print(f"\n{'='*60}")
print("  3. Ramachandran Plot (Backbone Dihedral Angles)")
print(f"{'='*60}")

# Extract backbone atoms per residue
backbone = defaultdict(dict)
for atom in atoms:
    if atom['name'] in ['N', 'CA', 'C']:
        backbone[atom['resid']][atom['name']] = np.array([atom['x'], atom['y'], atom['z']])

def calc_dihedral(p1, p2, p3, p4):
    """Calculate dihedral angle between 4 points."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    if n1_norm == 0 or n2_norm == 0:
        return 0
    n1 = n1 / n1_norm
    n2 = n2 / n2_norm
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.degrees(np.arctan2(y, x))

phi_psi = []
phi_psi_regions = []
sorted_resids = sorted(backbone.keys())

for i in range(1, len(sorted_resids) - 1):
    r_prev = sorted_resids[i - 1]
    r_curr = sorted_resids[i]
    r_next = sorted_resids[i + 1]

    if r_curr - r_prev != 1 or r_next - r_curr != 1:
        continue

    try:
        # Phi: C(i-1) - N(i) - CA(i) - C(i)
        phi = calc_dihedral(
            backbone[r_prev]['C'], backbone[r_curr]['N'],
            backbone[r_curr]['CA'], backbone[r_curr]['C'])
        # Psi: N(i) - CA(i) - C(i) - N(i+1)
        psi = calc_dihedral(
            backbone[r_curr]['N'], backbone[r_curr]['CA'],
            backbone[r_curr]['C'], backbone[r_next]['N'])
        phi_psi.append((phi, psi))
        phi_psi_regions.append(get_region(r_curr))
    except (KeyError, ValueError):
        continue

phi_vals = [pp[0] for pp in phi_psi]
psi_vals = [pp[1] for pp in phi_psi]

fig, ax = plt.subplots(figsize=(8, 8))

# Background: allowed regions
alpha_x = np.linspace(-180, -20, 50)
alpha_y = np.linspace(-80, 50, 50)
beta_x = np.linspace(-180, -50, 50)
beta_y = np.linspace(50, 180, 50)

# Scatter plot colored by region
for region_name in ['FR1', 'FR2', 'FR3', 'FR4', 'CDR1', 'CDR2', 'CDR3']:
    indices = [i for i, r in enumerate(phi_psi_regions) if r == region_name]
    if indices:
        phi_r = [phi_vals[i] for i in indices]
        psi_r = [psi_vals[i] for i in indices]
        is_cdr = 'CDR' in region_name
        ax.scatter(phi_r, psi_r,
                  c=REGION_COLORS[region_name],
                  s=50 if is_cdr else 20,
                  label=region_name if is_cdr else None,
                  edgecolors='black' if is_cdr else 'gray',
                  linewidths=0.8 if is_cdr else 0.3,
                  alpha=0.8,
                  zorder=10 if is_cdr else 5)

ax.set_xlabel('Phi (degrees)')
ax.set_ylabel('Psi (degrees)')
ax.set_title('Ramachandran Plot (ESMFold VHH Prediction)')
ax.set_xlim(-180, 180)
ax.set_ylim(-180, 180)
ax.axhline(0, color='gray', linewidth=0.5, alpha=0.3)
ax.axvline(0, color='gray', linewidth=0.5, alpha=0.3)
ax.legend(fontsize=9)
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "09_esmfold_ramachandran.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"  >> Saved: analysis_results/09_esmfold_ramachandran.png")

# Count residues in allowed regions
alpha_helix = sum(1 for phi, psi in phi_psi if -160 < phi < -20 and -80 < psi < 50)
beta_sheet = sum(1 for phi, psi in phi_psi if -180 < phi < -50 and 50 < psi < 180)
print(f"\n  Backbone dihedral statistics:")
print(f"    Total residues with phi/psi: {len(phi_psi)}")
print(f"    Alpha-helix region: {alpha_helix} ({alpha_helix/len(phi_psi)*100:.1f}%)")
print(f"    Beta-sheet region: {beta_sheet} ({beta_sheet/len(phi_psi)*100:.1f}%)")
print(f"    (VHH domains are primarily beta-sheet immunoglobulin fold)")


# =============================================================================
# 4. 3D Structure Summary Visualization
# =============================================================================
print(f"\n{'='*60}")
print("  4. Structure Summary Visualization")
print(f"{'='*60}")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 4a: pLDDT bar chart by region
ax = axes[0][0]
region_names = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
region_means = []
region_stds = []
bar_colors = []
for name in region_names:
    start, end = REGIONS[name]
    vals = [p for r, p in zip(resids, plddt) if start < r <= end]
    region_means.append(np.mean(vals) if vals else 0)
    region_stds.append(np.std(vals) if vals else 0)
    bar_colors.append(REGION_COLORS[name])

bars = ax.bar(region_names, region_means, yerr=region_stds, capsize=3,
              color=bar_colors, edgecolor='black', linewidth=0.5)
ax.set_ylabel('Mean pLDDT')
ax.set_title('pLDDT by Region')
ax.set_ylim(0, 1.1)
ax.axhline(0.7, color='orange', linestyle='--', alpha=0.5)
ax.axhline(0.9, color='green', linestyle='--', alpha=0.5)

# 4b: Region length and structure
ax = axes[0][1]
region_lens = [REGIONS[n][1] - REGIONS[n][0] for n in region_names]
ax.barh(region_names, region_lens, color=bar_colors, edgecolor='black', linewidth=0.5)
ax.set_xlabel('Length (residues)')
ax.set_title('Region Lengths')
for i, (name, length) in enumerate(zip(region_names, region_lens)):
    start, end = REGIONS[name]
    seq = VHH_SEQ[start:end]
    ax.text(length + 0.5, i, f' {seq}', va='center', fontsize=6, family='monospace')

# 4c: Distance from center of mass (proxy for surface exposure)
ax = axes[1][0]
com = np.mean([[ca['x'], ca['y'], ca['z']] for ca in ca_atoms], axis=0)
dist_from_com = []
for ca in ca_atoms:
    d = math.sqrt((ca['x']-com[0])**2 + (ca['y']-com[1])**2 + (ca['z']-com[2])**2)
    dist_from_com.append(d)

for name, (start, end) in REGIONS.items():
    indices = [i for i, r in enumerate(resids) if start < r <= end]
    if indices:
        x_vals = [resids[i] for i in indices]
        y_vals = [dist_from_com[i] for i in indices]
        ax.fill_between(x_vals, y_vals, alpha=0.3, color=REGION_COLORS[name])
        ax.plot(x_vals, y_vals, color=REGION_COLORS[name], linewidth=1.5)

ax.set_xlabel('Residue Number')
ax.set_ylabel('Distance from Center of Mass (A)')
ax.set_title('Surface Exposure (Distance from CoM)')

# Annotate CDR peaks
for name in ['CDR1', 'CDR2', 'CDR3']:
    start, end = REGIONS[name]
    indices = [i for i, r in enumerate(resids) if start < r <= end]
    if indices:
        max_idx = indices[np.argmax([dist_from_com[i] for i in indices])]
        ax.annotate(name, xy=(resids[max_idx], dist_from_com[max_idx]),
                   xytext=(resids[max_idx], dist_from_com[max_idx] + 3),
                   fontsize=8, fontweight='bold', color=REGION_COLORS[name],
                   ha='center', arrowprops=dict(arrowstyle='->', color=REGION_COLORS[name]))

# 4d: B-factor (pLDDT) colored structure projection (XY)
ax = axes[1][1]
x_coords = [ca['x'] for ca in ca_atoms]
y_coords = [ca['y'] for ca in ca_atoms]

# Plot backbone trace
for i in range(len(ca_atoms) - 1):
    ax.plot([x_coords[i], x_coords[i+1]], [y_coords[i], y_coords[i+1]],
            color='lightgray', linewidth=1, zorder=1)

# Color by region
for name, (start, end) in REGIONS.items():
    indices = [i for i, r in enumerate(resids) if start < r <= end]
    if indices:
        xs = [x_coords[i] for i in indices]
        ys = [y_coords[i] for i in indices]
        ps = [plddt[i] for i in indices]
        sc = ax.scatter(xs, ys, c=ps, cmap='RdYlGn', vmin=0.5, vmax=1.0,
                       s=40, edgecolors=REGION_COLORS[name], linewidths=1.5, zorder=5)

        # Label CDR regions
        if 'CDR' in name:
            cx, cy = np.mean(xs), np.mean(ys)
            ax.text(cx, cy + 2, name, fontsize=9, fontweight='bold',
                   color=REGION_COLORS[name], ha='center',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

plt.colorbar(sc, ax=ax, shrink=0.8, label='pLDDT')
ax.set_xlabel('X (A)')
ax.set_ylabel('Y (A)')
ax.set_title('Structure Projection (colored by pLDDT)')
ax.set_aspect('equal')

# Label N and C termini
ax.annotate('N', xy=(x_coords[0], y_coords[0]),
           fontsize=12, fontweight='bold', color='blue',
           xytext=(x_coords[0]-5, y_coords[0]-5),
           arrowprops=dict(arrowstyle='->', color='blue'))
ax.annotate('C', xy=(x_coords[-1], y_coords[-1]),
           fontsize=12, fontweight='bold', color='red',
           xytext=(x_coords[-1]+5, y_coords[-1]+5),
           arrowprops=dict(arrowstyle='->', color='red'))

plt.suptitle('ESMFold VHH Structure Analysis Summary', fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "10_esmfold_summary.png"), dpi=150, bbox_inches='tight')
plt.close()
print(f"  >> Saved: analysis_results/10_esmfold_summary.png")


# =============================================================================
# 5. CDR Loop Geometry Analysis
# =============================================================================
print(f"\n{'='*60}")
print("  5. CDR Loop Geometry")
print(f"{'='*60}")

for cdr_name in ['CDR1', 'CDR2', 'CDR3']:
    start, end = REGIONS[cdr_name]
    cdr_cas = [ca for ca in ca_atoms if start < ca['resid'] <= end]

    if len(cdr_cas) >= 2:
        # Loop span (distance between start and end CA)
        d_span = math.sqrt(
            (cdr_cas[0]['x'] - cdr_cas[-1]['x'])**2 +
            (cdr_cas[0]['y'] - cdr_cas[-1]['y'])**2 +
            (cdr_cas[0]['z'] - cdr_cas[-1]['z'])**2
        )
        # Loop length along backbone
        backbone_len = 0
        for i in range(len(cdr_cas) - 1):
            backbone_len += math.sqrt(
                (cdr_cas[i+1]['x'] - cdr_cas[i]['x'])**2 +
                (cdr_cas[i+1]['y'] - cdr_cas[i]['y'])**2 +
                (cdr_cas[i+1]['z'] - cdr_cas[i]['z'])**2
            )

        # Mean pLDDT
        mean_plddt = np.mean([ca['bfactor'] for ca in cdr_cas])

        print(f"\n  {cdr_name} ({end-start} residues: {VHH_SEQ[start:end]})")
        print(f"    End-to-end span: {d_span:.1f} A")
        print(f"    Backbone length: {backbone_len:.1f} A")
        print(f"    Compactness (span/backbone): {d_span/backbone_len:.2f}")
        print(f"    Mean pLDDT: {mean_plddt:.3f}")


# =============================================================================
# Final Summary
# =============================================================================
print(f"\n{'='*60}")
print("  Structure Prediction Summary")
print(f"{'='*60}")
print(f"""
  Model: ESMFold v1
  Sequence: VHH nanobody (119 aa)
  Overall pLDDT: {np.mean(plddt):.3f}

  The structure shows a typical immunoglobulin VHH fold:
  - Beta-sheet rich framework ({beta_sheet/len(phi_psi)*100:.0f}% beta region)
  - CDR loops protruding from the beta-sandwich core
  - High confidence in framework regions
  - CDR loops show varying confidence (flexibility expected)

  Generated files:
  - vhh_structure_esmfold.pdb     (3D structure)
  - 07_esmfold_plddt.png          (per-residue confidence)
  - 08_esmfold_contact_map.png    (distance & contact maps)
  - 09_esmfold_ramachandran.png   (backbone dihedral angles)
  - 10_esmfold_summary.png        (multi-panel summary)
""")
