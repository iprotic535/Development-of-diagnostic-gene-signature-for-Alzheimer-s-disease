# Development of Diagnostic Gene Signatures for Alzheimer's Disease (AD) and Mild Cognitive Impairment (MCI)

## **Background**

- The cause of Mild Cognitive Impairment (MCI) remains unknown.
- Individuals with MCI have a higher risk of developing Alzheimer's Disease (AD) or other neurological disorders.
- Some scientists believe MCI represents a transitional phase between normal cognition and AD.
- Not all individuals with MCI progress to AD; some may even recover.

## **Project Aims**

1. Develop a **blood-based gene signature** to differentiate **AD patients from controls**.
2. Develop a **blood-based gene signature** to differentiate **MCI patients from controls**.
3. Develop a **blood-based gene signature** to differentiate **patients with cognitive impairment (AD/MCI) from controls**.

## **Study Cohorts**

- **Cohort A**: 249 controls, 204 AD patients, 134 MCI patients.
  - Data files: `expr_cohort_A.txt`, `pheno_cohort_A.txt`
- **Cohort B**: 104 controls, 145 AD patients, 80 MCI patients.
  - Data files: `expr_cohort_B.txt`, `pheno_cohort_B.txt`
- **Cohort C**: 134 controls, 139 AD patients, 109 MCI patients (6 "OTHER" status individuals).
  - Data files: `expr_cohort_C.txt`, `pheno_cohort_C.txt`

## **Requirements**

1. **Select one cohort** as the **discovery cohort** to develop gene signatures.
2. **Validate** the gene signature(s) in the **other two cohorts**.
3. **Generate figures** showcasing the performance of gene signatures in both discovery and validation cohorts.
4. **Required Files** for this project can be accessed from this link: https://drive.google.com/file/d/1sb0ceVLUW69pceEBU8Qw3Wk3GwVJWp4T/view?usp=sharing

## **Installation & Dependencies**

Ensure you have R installed along with the following required packages:

```r
install.packages("survival")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("caret")
```

## **Running the Analysis**

1. **Load the provided R script** (`Alzheimers_GeneSignature.R`).
2. Set the working directory to the folder containing the expression and phenotype files.
3. Run the script to **develop the gene signature** in the discovery cohort.
4. Validate the performance in the other two cohorts.
5. Visualize the survival analysis results.

## **Outputs**

- Gene signature file (`gene_signature.txt`).
- Survival analysis plots (`KM_survival_plot.png`).
- Performance metrics (`performance_metrics.txt`).

## **Example Usage**

```r
source("Alzheimers_GeneSignature.R")
```

## **License**

This project is open-source and available under the MIT License.

## **Contact**

For any queries, reach out to [**protic535@outlook.com**](mailto\:protic535@outlook.com) or create an issue in the repository.

---

**Created for Alzheimer's Disease Gene Signature Development Project**

