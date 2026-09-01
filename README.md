# Framework_v1

## A Multi-Omics Evidence Framework for Prioritizing Variant–Drug Relationships in Coronary Artery Disease

This repository contains the analysis workflow and supporting code for a multi-layer computational framework developed to prioritize **variant–drug relationships relevant to coronary artery disease (CAD)**.

The framework integrates disease genetics, metabolite-associated genetic variation, Mendelian randomization, human pharmacometabolomic evidence, functional variant annotation, and external pharmacogenomic resources.

The primary objective is not to infer clinical actionability or therapeutic efficacy directly, but to reduce a large multi-omics evidence space into a focused set of **testable pharmacogenomic hypotheses supported by convergent genetic, metabolic, functional, and pharmacological evidence**.

---

## Overview

The analysis follows the general workflow:

**CAD-associated metabolite colocalization**  
→ **metabolite-specific mGWAS instruments**  
→ **LD clumping and harmonization with CAD GWAS**  
→ **two-sample Mendelian randomization**  
→ **drug–metabolite pharmacometabolomic integration**  
→ **directional evidence prioritization**  
→ **variant functional annotation and gene assignment**  
→ **gene–drug-category module construction**  
→ **external benchmarking and tiered candidate prioritization**  
→ **post-prioritization locus and literature evaluation**

The framework was developed as a proof-of-concept in CAD but is conceptually adaptable to other complex diseases for which compatible genetic and molecular datasets are available.

---

## Analysis Workflow

### 1. CAD-associated metabolic traits

Precomputed metabolite–disease colocalization results were obtained from the data release accompanying the work of **Tambets et al.**

Zenodo record:

https://zenodo.org/records/17945143

The following archive and dataset were used:

```text
Suzuki_Aragam_coloc.tar.gz
└── UKBB_EUR_Suzuki_Aragam_coloc_results.tsv
