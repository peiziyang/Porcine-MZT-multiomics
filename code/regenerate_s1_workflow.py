#!/usr/bin/env python
"""Generate Fig S1: Multi-omics workflow overview (pure schematic).

Clean, evenly spaced layout:
  - Row 1 (4 inputs) and Row 2 (6 parallel analyses) share identical left/right edges
  - equal inter-row spacing and uniform within-row gaps
  - larger fonts so text fills most of each box
Data-source arrows follow the manuscript (verified separately):
  GSE235729 -> MOFA+ ;  GSE235729 (44 CGmap) -> CpG reanalysis
  GSE168106 -> Atlas projection + SCENIC ;  (in vivo backbone) -> Waddington OT
  GSE44183  -> Cross-species comparison
  SRP301735/GSE164812 -> PA/IVF module reanalysis ; -> Waddington OT (IVF/PA)
"""
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.family': 'Arial', 'svg.fonttype': 'none', 'savefig.dpi': 300})
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT = r'E:/Workbuddy/2026-07-27-11-58-27/manuscript/submission/figures'
PAD = 0.12

# ============ geometry (top-down, uniform gaps) ============
X0, X1 = 0.6, 19.6              # identical horizontal span for all rows
WIDTH  = X1 - X0
GAP    = 1.5                    # uniform spacing between row bands

OUT_H  = 1.8
AN_H   = 2.5
IN_H   = 1.7
OUT_Y  = 0.5
AN_Y   = OUT_Y + OUT_H + GAP            # = 3.8
IN_Y   = AN_Y + AN_H + GAP              # = 7.8  (top 9.5)
YTOP   = IN_Y + IN_H                     # 9.5

# row 1: 4 inputs, 3 inner gaps of 1.0
IN_GAP = 1.0
IN_W   = (WIDTH - 3*IN_GAP) / 4.0
# row 2: 6 analyses, 5 inner gaps of 0.44
AN_GAP = 0.44
AN_W   = (WIDTH - 5*AN_GAP) / 6.0

fig, ax = plt.subplots(figsize=(14, 8.5))
ax.set_xlim(0, 20); ax.set_ylim(0, 10.4); ax.axis('off')

def box(x, y, w, h, text, fc='#EAF2FA', ec='black', fs=8, bold=False, lw=1.3):
    ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.12',
                                         facecolor=fc, edgecolor=ec, lw=lw))
    ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=fs,
            fontweight='bold' if bold else 'normal', wrap=True, linespacing=1.45)

def arrow(x1, y1, x2, y2, color='gray'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='-|>', lw=1.4, color=color,
                                mutation_scale=14))

fig.suptitle('Multi-omics workflow overview', fontsize=14, fontweight='bold', y=0.99)

# ============ Row 1: Input data ============
in_xs = [X0 + i*(IN_W + IN_GAP) for i in range(4)]
box(in_xs[0], IN_Y, IN_W, IN_H,
    'GSE235729\n32 GV oocytes, matched\nRNA-seq + WGBS (44 CGmap)', '#EAF2FA', fs=8.5)
box(in_xs[1], IN_Y, IN_W, IN_H,
    'GSE168106\npreimplantation atlas\nE0\u2013E10 (in vivo backbone)', '#EAF2FA', fs=8.5)
box(in_xs[2], IN_Y, IN_W, IN_H,
    'GSE44183\nhuman / mouse\npreimplantation RNA-seq', '#EAF2FA', fs=8.5)
box(in_xs[3], IN_Y, IN_W, IN_H,
    'SRP301735 / GSE164812\n42-sample PA / IVF\nbulk RNA-seq', '#EAF2FA', fs=8.5)

# ============ Row 2: Analyses (6 parallel) ============
an_xs = [X0 + i*(AN_W + AN_GAP) for i in range(6)]
box(an_xs[0], AN_Y, AN_W, AN_H,
    'MOFA+ integration\n(RNA + METH, 32 GV,\nF1\u2013F7 factors)\nGSE235729', '#FFF2CC', fs=8)
box(an_xs[1], AN_Y, AN_W, AN_H,
    'Atlas projection\n+ SCENIC regulons\n(E0\u2013E10 dynamics;\nperturbation)\nGSE168106', '#FFF2CC', fs=7.6)
box(an_xs[2], AN_Y, AN_W, AN_H,
    'Cross-species\ncomparison\n(F1 gene set;\nhuman / mouse)\nGSE44183', '#FCE4D6', fs=8)
box(an_xs[3], AN_Y, AN_W, AN_H,
    'CpG reanalysis\n(F1 promoter CpG;\n44 CGmap files)\nGSE235729', '#FCE4D6', fs=8)
box(an_xs[4], AN_Y, AN_W, AN_H,
    'Waddington OT\n(atlas backbone +\nIVF/PA displacement)\nGSE168106\n+ IVF/PA', '#FCE4D6', fs=7.6)
box(an_xs[5], AN_Y, AN_W, AN_H,
    'PA/IVF module\nreanalysis\n(F1/F3/F4/F6 vs PA,\nstage-adjusted)\nSRP301735', '#FFF2CC', fs=7.6)
an_cx  = [x + AN_W/2 for x in an_xs]
an_top = AN_Y + AN_H + PAD
an_bot = AN_Y

# ============ Row 3: Outcome (same full span) ============
box(X0, OUT_Y, WIDTH, OUT_H,
    'Candidate maternal-associated modules (F1\u2013F7)\nwith module-level PA vs IVF differences\n(F1 heterogeneous; F3 stage-declining set reduced in PA)',
    '#DDEBD2', fs=9.5)

# ============ Arrows: inputs -> analyses ============
CX0 = [x + IN_W/2 for x in in_xs]       # input centres
# GSE235729 -> MOFA+ ;  GSE235729 CGmap -> CpG reanalysis
arrow(CX0[0], IN_Y, an_cx[0], an_top)
arrow(CX0[0]+1.1, IN_Y-0.25, an_cx[3], an_top)
# GSE168106 -> Atlas/SCENIC ;  in vivo backbone -> Waddington OT
arrow(CX0[1], IN_Y, an_cx[1], an_top)
arrow(CX0[1]+1.1, IN_Y-0.25, an_cx[4], an_top)
# GSE44183 -> Cross-species
arrow(CX0[2], IN_Y, an_cx[2], an_top)
# SRP301735/GSE164812 -> PA/IVF ;  IVF/PA -> Waddington OT
arrow(CX0[3], IN_Y, an_cx[5], an_top)
arrow(CX0[3]-1.1, IN_Y-0.25, an_cx[4], an_top)

# ============ Arrows: analyses -> outcome ============
for cx in an_cx:
    arrow(cx, an_bot, cx, OUT_Y + OUT_H + PAD)

plt.tight_layout()
fig.savefig(f'{OUT}/Fig_S1_overview.png', dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT}/Fig_S1_overview.svg', format='svg', bbox_inches='tight')
print(f"S1 saved")
print(f"  rows share span X {X0}-{X1}")
print(f"  IN_Y={IN_Y}..{IN_Y+IN_H} | AN_Y={AN_Y}..{AN_Y+AN_H} | OUT_Y={OUT_Y}..{OUT_Y+OUT_H}")
print(f"  gap row1->row2 = {IN_Y-(AN_Y+AN_H):.2f} | gap row2->out = {AN_Y-(OUT_Y+OUT_H):.2f}")
