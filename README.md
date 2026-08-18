# Tendon Repair Trajectory Analysis

This repository contains the analysis scripts used to investigate mechanotransduction-associated transcriptional programs during tendon repair. The scripts are organized into separate **mouse** and **human** analysis workflows, followed by cross-species comparison.

### Mouse Analysis
The mouse workflow includes preprocessing, cell-type annotation, UCell gene-signature scoring, condition and cell-type analysis, and trajectory inference to characterize mechanotransduction-associated transcriptional dynamics across tendon healing.

### Human Analysis
The human workflow focuses on fibroblast populations from healthy and ruptured tendon samples and evaluates mechanotransduction-associated gene signatures for comparison with the mouse repair-associated populations.

### Data Sources

**Mouse scRNA-seq:**  NCBI Gene Expression Omnibus (GEO), [GSE288443](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288443)

**Human snRNA-seq:**  [CZ CELLxGENE collection](https://cellxgene.cziscience.com/collections/579203e2-182f-47bc-8230-7aa47247e2a4)

The mouse dataset contains samples from WT, 1D, 7D, and 30D following tendon injury, while the human dataset contains healthy and tendon rupture samples. :contentReference[oaicite:0]{index=0}

### Software

The analyses were performed in **R**; package and version numbers can be found in [sessionInfo](/scripts/human/sessionInfo.txt).
