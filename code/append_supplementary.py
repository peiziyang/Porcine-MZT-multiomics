# -*- coding: utf-8 -*-
"""Append a Supplementary Material section (Supplementary Fig. S1 = the two remaining
Fig. 1 UMAPs: by source / by stage-group) to the v2 manuscript, producing v3."""

SRC = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/manuscript_full_with_figures_v6.md"
OUT = "E:/Workbuddy/2026-07-27-11-58-27/manuscript/manuscript_full_with_figures_v6.md"

BY_SOURCE = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/m2_umap_by_source.png"
BY_STAGE_GROUP = "E:/Workbuddy/2026-07-27-11-58-27/data/processed/m2_umap_by_stage_group.png"

SUPP = r"""# SUPPLEMENTARY MATERIAL

## Supplementary Figure S1. Additional UMAP embeddings of the porcine preimplantation reference atlas.
(A) UMAP coloured by study-of-origin (GSE168106, GSE112380, GSE139512, GSE164812, GSE234116, GSE160334). Despite independent generation, batches separate less cleanly than developmental stage, confirming that stage — not platform — is the primary axis of variation (cf. Fig. 1B–C). (B) UMAP coloured by coarse stage-group (oocyte/GV → cleavage → blastocyst → post-implantation/pgEpiSC), illustrating the progression of the preimplantation programme as a near-continuous trajectory. Both panels derive from the same integrated 1,955-cell, 8,764-gene coordinate system used in Fig. 1.
![Supplementary Fig. S1A. UMAP coloured by study-of-origin.](E:/Workbuddy/2026-07-27-11-58-27/data/processed/m2_umap_by_source.png)
![Supplementary Fig. S1B. UMAP coloured by coarse stage-group.](E:/Workbuddy/2026-07-27-11-58-27/data/processed/m2_umap_by_stage_group.png)
"""

with open(SRC, encoding="utf-8") as f:
    body = f.read().rstrip("\n")

out_text = body + "\n\n" + SUPP

with open(OUT, "w", encoding="utf-8") as f:
    f.write(out_text)

print("Wrote", OUT, "| total chars:", len(out_text))
