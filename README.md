# Framework_v1

## Multi-Omics Evidence Framework for Prioritizing Variant–Drug-Class Relationships in Coronary Artery Disease

Framework_v1 is a modular computational workflow that integrates coronary artery disease (CAD) GWAS data, metabolite GWAS (mGWAS), colocalization, Mendelian randomization (MR), population-scale pharmacometabolomic associations, variant-level genetic evidence, and functional variant annotation.

The current implementation uses coronary artery disease as the disease application and the Nature Medicine pharmacometabolomics resource as the drug–metabolite evidence layer.

> **Important:** The current Framework_v1 implementation does **not** use the earlier BioWalk/HMDB drug-network layer. That component was removed and replaced with the Nature Medicine drug–metabolite association atlas.

---

## Current workflow

```text
CAD GWAS
   |
   v
CAD–metabolite colocalization
   |
   v
Nature Medicine metabolite mapping
   |
   v
mGWAS download
   |
   v
Genome-wide significant instruments
   |
   v
LD clumping
   |
   v
Harmonization
   |
   v
Mendelian randomization
   |
   v
Nature Medicine drug–metabolite associations
   |
   v
Evidence matrix and pharmacometabolomic candidate tiers
   |
   v
Variant-level evidence propagation
   |
   v
VEP functional annotation
   |
   v
Variant–gene–drug-class evidence integration
```

The workflow is designed so that each evidence layer remains traceable. Functional annotation is used for interpretation and does not overwrite the pharmacometabolomic evidence tier.

---

## Repository structure

A typical Framework_v1 layout is:

```text
Framework_v1/
├── config/
│   └── framework_config.R
│
├── R/
│   └── reusable R modules
│
├── scripts/
│   ├── 00_prepare_natmed_inputs.R
│   ├── 00_download_exposure_mGWAS.R
│   ├── 00_reset_generated_outputs.R
│   ├── 01_MR_pipeline.R
│   ├── 02_PharmacoMetabolomics.R
│   ├── 03_EvidenceMatrix.R
│   ├── 04_VariantEvidence.R
│   ├── 05_prepare_VEP_coordinates.R
│   ├── 06_parse_VEP_annotation.R
│   ├── 07_variant_annotation_pipeline.R
│   ├── 08_integrate_variant_gene_drug.R
│   └── 10_smoke_test_natmed.R
│
├── data/
│   ├── raw/
│   │   ├── NatMed/
│   │   │   └── 41591_2019_722_MOESM3_ESM.xlsx
│   │   ├── outcome_GWAS/
│   │   │   └── CAD_Aragam2022_GCST90132314_GRCh37.tsv.gz
│   │   ├── exposure_mGWAS/
│   │   └── LD_reference/
│   │       ├── EUR.bed
│   │       ├── EUR.bim
│   │       └── EUR.fam
│   │
│   └── processed/
│
├── data/intermediate/
├── results/
│   └── evidence/
├── logs/
└── README.md
```

The central analysis paths and thresholds are defined in:

```text
config/framework_config.R
```

---

## Main input datasets

### 1. CAD GWAS

The current outcome dataset is:

```text
CAD_Aragam2022_GCST90132314_GRCh37.tsv.gz
```

GWAS Catalog accession:

```text
GCST90132314
```

This dataset is used as the CAD outcome for downstream MR analyses.

---

### 2. CAD–metabolite colocalization

The current colocalization input is:

```text
CAD_PP_H4_099_coloc_with_UKBB_mGWAS_paths.tsv
```

The current implementation starts from the high-confidence CAD–metabolite colocalization set used to identify candidate metabolite traits and their corresponding UK Biobank mGWAS resources.

---

### 3. Nature Medicine pharmacometabolomics workbook

The pharmacometabolomic evidence source is:

```text
data/raw/NatMed/41591_2019_722_MOESM3_ESM.xlsx
```

The workflow uses the following sheets:

```text
Suppl_T1   drug-category information
Suppl_T2   Nightingale metabolite definitions
Suppl_T4   drug–metabolite association results
```

The current pipeline uses the **M2 model as the primary pharmacometabolomic model**.

M2 corresponds to the M1 model additionally adjusted for BMI and smoking.

Whenever available, the row-specific multiple-testing threshold reported in `Suppl_T4` is used to determine M2 significance.

---

## Pipeline stages

### 0. Prepare Nature Medicine inputs

```bash
Rscript scripts/00_prepare_natmed_inputs.R
```

This script:

- reads `Suppl_T1`, `Suppl_T2`, and `Suppl_T4`;
- harmonizes CAD/Tambets metabolite trait names with Nature Medicine Nightingale abbreviations;
- handles dot/underscore naming differences such as `XXL_VLDL_P` vs `XXL.VLDL.P`;
- uses exact reported-trait/explanation matches;
- applies a small explicit alias table only where required;
- creates the Nature Medicine-compatible candidate metabolite manifest.

Principal outputs include:

```text
NatMedMetaboliteMapping.tsv
CandidateMetabolites_NatMed.tsv
NatMedDrugCategories.tsv
NatMedM2Associations.tsv
NatMedMetaboliteMappingSummary.tsv
```

---

### 0b. Smoke test

Before downloading large mGWAS files, validate the source workbook and trait mapping:

```bash
Rscript scripts/10_smoke_test_natmed.R
```

For the source files used in the current implementation, the expected counts are:

```text
CAD PP.H4 >= 0.99 traits                         177
Nature Medicine Nightingale traits              150
matched traits                                  102
Nature Medicine drug categories                  87
Suppl_T4 M1/M2 rows                           13,050
rows for the 102 overlapping traits            8,874
M2 multiple-testing-significant rows           1,260
drug categories represented among these rows      42
```

The earlier provisional overlap count of 103 was incorrect. The source workbook contains **acetate (Ace)** but not **acetone**, so acetone must not be forced to match.

---

### 0c. Download exposure mGWAS files

```bash
Rscript scripts/00_download_exposure_mGWAS.R
```

This script downloads only the mGWAS traits that:

1. are present in the CAD colocalization set, and
2. overlap the Nature Medicine metabolite panel.

Downloads use accession-based filenames and completed files are skipped, allowing safe resumption.

The generated manifest is:

```text
mGWAS_download_manifest.tsv
```

---

### 1. Mendelian randomization

```bash
Rscript scripts/01_MR_pipeline.R
```

The MR methodology remains unchanged from the earlier Framework_v1 implementation.

The main steps are:

```text
genome-wide significant SNP selection
        ↓
LD clumping
        ↓
exposure/outcome harmonization
        ↓
MR analysis
        ↓
sensitivity and QC analyses
```

Current core parameters are:

```text
Genome-wide significance      P < 5e-8
LD clumping r²                0.001
LD clumping window            10,000 kb
Strong-instrument reference   F > 10
Minimum SNP count             3
MR-PRESSO                     enabled
MR-PRESSO simulations         1,000
```

Principal output:

```text
data/processed/MR_summary.tsv
```

The MR pipeline also produces harmonized instruments and intermediate sensitivity/QC outputs.

MR results should not be interpreted in isolation; they are integrated with colocalization and pharmacometabolomic evidence downstream.

---

### 2. Pharmacometabolomic evidence

```bash
Rscript scripts/02_PharmacoMetabolomics.R
```

This script joins the MR results to the Nature Medicine M1/M2 drug–metabolite association data.

Principal output:

```text
data/processed/PharmacoMetabolomicEvidence.tsv
```

Summary output:

```text
data/processed/PharmacoMetabolomicEvidence_summary.tsv
```

Important interpretation:

> The reported drug–metabolite beta is an **observational association coefficient** from the Nature Medicine analysis. It must not be interpreted as a causal drug effect.

M2 is the primary model.

Where available, row-specific multiple-testing thresholds from `Suppl_T4` are used. A fallback threshold is defined in the configuration for cases where the row-specific threshold is unavailable.

---

### 3. Evidence matrix and candidate selection

```bash
Rscript scripts/03_EvidenceMatrix.R
```

The current implementation deliberately does **not** use the older arbitrary weighted `FinalPriorityScore`.

Primary candidate selection requires:

```text
1. Significant M2 drug–metabolite association
2. Opposite signs for MR beta and drug–metabolite beta
3. MR FDR < 0.05
```

The sign criterion represents a mechanistic prioritization rule:

```text
MR beta > 0 and drug beta < 0
or
MR beta < 0 and drug beta > 0
```

That is, a drug class is prioritized when its observed metabolite association is opposite to the direction in which the metabolite is associated with CAD in the MR analysis.

This is a prioritization heuristic and **not proof that the drug modifies CAD risk through that metabolite**.

### Tier definition

M1/M2 robustness is used to distinguish candidate tiers rather than being applied as a universal hard filter.

In the current implementation:

```text
Tier 1
    M2 significant
    MR/drug beta directions opposite
    MR FDR < 0.05
    significant in both M1 and M2
    M1 and M2 beta directions consistent

Tier 2
    passes the primary candidate criteria
    but does not meet the additional Tier 1 robustness requirement
```

Current configuration:

```text
MR FDR threshold                  0.05
MR nominal P threshold            0.05
Top candidates per metabolite       20
Top candidates overall             100
```

Principal outputs include:

```text
EvidenceMatrix.tsv
DrugClassCandidateNetwork.tsv
```

---

## Why there is no mandatory final score

Framework_v1 is intended to preserve the evidence layers rather than collapse them into a single opaque value.

The current implementation therefore favors:

```text
MR evidence
+
colocalization evidence
+
Nature Medicine M1/M2 evidence
+
directional relationship
+
candidate tier
+
variant-level support
+
functional annotation
```

over a mandatory weighted score.

The evidence matrix is the primary scientific output.

---

### 4. Variant-level evidence

```bash
Rscript scripts/04_VariantEvidence.R
```

This stage propagates candidate drug-class evidence to the harmonized instrumental SNPs.

Principal outputs:

```text
InstrumentEvidence.tsv
MetaboliteInstrumentSummary.tsv
VariantDrugClassEvidence.tsv
EvidenceMatrix_withVariants.tsv
DrugClassCandidateNetwork_withVariants.tsv
```

This step connects the metabolite-level causal analysis to the individual genetic variants used as instruments.

---

### 5. Prepare VEP coordinates

```bash
Rscript scripts/05_prepare_VEP_coordinates.R
```

This stage creates the variant coordinate file required for Ensembl Variant Effect Predictor.

The current workflow uses GRCh37 coordinates.

Run VEP externally after this step.

The exact VEP command depends on the local VEP installation, cache, genome FASTA, and selected annotation options and should be recorded for reproducibility.

---

### 6. Parse VEP annotations

After VEP has completed:

```bash
Rscript scripts/06_parse_VEP_annotation.R
```

This stage parses the VEP results and constructs variant–gene annotations.

---

### 7. Variant annotation pipeline

```bash
Rscript scripts/07_variant_annotation_pipeline.R
```

This stage consolidates and QC-checks the functional annotation evidence.

VEP annotations are an **interpretation layer**. A VEP-assigned gene must not automatically be treated as a causal disease gene or a drug target.

---

### 8. Variant–gene–drug-class integration

```bash
Rscript scripts/08_integrate_variant_gene_drug.R
```

Final output:

```text
results/evidence/VariantGeneDrugClassEvidence.tsv
```

QC output:

```text
results/evidence/VariantGeneDrugClassIntegration_QC.tsv
```

Variant and gene annotations do not alter the original pharmacometabolomic candidate tier.

They provide additional biological context for interpreting the evidence.

---

## Clean restart

To deliberately remove previous generated outputs and previously downloaded exposure mGWAS files:

```bash
Rscript scripts/00_reset_generated_outputs.R --yes --delete-mgwas
```

Then run:

```bash
Rscript scripts/00_prepare_natmed_inputs.R
Rscript scripts/10_smoke_test_natmed.R
Rscript scripts/00_download_exposure_mGWAS.R
Rscript scripts/01_MR_pipeline.R
Rscript scripts/02_PharmacoMetabolomics.R
Rscript scripts/03_EvidenceMatrix.R
Rscript scripts/04_VariantEvidence.R
Rscript scripts/05_prepare_VEP_coordinates.R
```

Run VEP externally.

Then continue with:

```bash
Rscript scripts/06_parse_VEP_annotation.R
Rscript scripts/07_variant_annotation_pipeline.R
Rscript scripts/08_integrate_variant_gene_drug.R
```

---

## Configuration summary

The main configuration file is:

```text
config/framework_config.R
```

Current important settings include:

```text
CAD outcome
CAD_Aragam_2022_GCST90132314

Colocalization input
CAD_PP_H4_099_coloc_with_UKBB_mGWAS_paths.tsv

Nature Medicine workbook
41591_2019_722_MOESM3_ESM.xlsx

Primary pharmacometabolomic model
M2

MR P threshold
5e-8

LD r²
0.001

LD window
10,000 kb

F-statistic reference
10

Minimum SNPs
3

MR-PRESSO
TRUE

MR-PRESSO simulations
1,000

MR FDR candidate threshold
0.05

Directional concordance required
TRUE

MR FDR required
TRUE

Top candidates per metabolite
20

Top candidates overall
100
```

---

## Scientific interpretation safeguards

Framework_v1 intentionally keeps several evidence types separate.

### MR

MR supports genetic evidence consistent with a causal metabolite–CAD relationship under the assumptions of instrumental-variable analysis.

It does not prove that pharmacologically modifying the metabolite will reproduce the genetically predicted effect.

### Drug–metabolite association

Nature Medicine M1/M2 betas are observational drug–metabolite associations.

They are not randomized drug effects and should not be interpreted as causal pharmacological estimates.

### Directional prioritization

Opposite MR and drug–metabolite beta signs are used as a prioritization criterion.

This identifies potentially interesting drug classes but does not establish therapeutic benefit.

### Colocalization

High PP.H4 support strengthens the argument that the metabolite and CAD signals share a locus-level causal variant.

Colocalization alone does not identify the causal gene or prove a drug mechanism.

### VEP

VEP provides functional consequence and gene annotation.

A VEP-linked gene is not automatically:

```text
a causal CAD gene
a metabolite-regulating gene
a pharmacological target
or
a clinically actionable gene
```

Functional annotation is therefore preserved as a separate evidence layer.

---

## Reproducibility

For publication-grade analyses, record:

- R version;
- package versions;
- TwoSampleMR version;
- MRPRESSO version;
- PLINK version;
- Ensembl VEP version;
- VEP cache release;
- GRCh37 reference FASTA;
- CAD GWAS accession;
- mGWAS accessions;
- LD reference panel;
- colocalization source and threshold;
- Nature Medicine workbook version;
- exact configuration file;
- exact VEP command;
- generated QC outputs.

Large GWAS and mGWAS summary-statistic files are not intended to be redistributed through the repository.

---

## Scope

The current implementation focuses on CAD, but the architecture is intended to be reusable for other diseases by replacing the disease-specific GWAS and corresponding colocalization layer while retaining the same evidence-integration principles.

Potential future applications include other cardiovascular and complex diseases.

---

## Code availability

Repository:

```text
https://github.com/doganlab/Framework_v1
```

---

## Citation

If you use Framework_v1 or adapt the workflow, please cite the associated publication once the final bibliographic information is available.

---

## Status

Framework_v1 is a research workflow for evidence integration and candidate prioritization.

It is **not** a clinical decision-support system.

**Research use only.**
