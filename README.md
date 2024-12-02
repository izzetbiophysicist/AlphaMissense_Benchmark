# Alphamissense Evaluation on ClinVar Subsets

This repository provides a pipeline for evaluating Alphamissense performance on specific subsets of ClinVar data, focusing on genes or diseases of interest. The pipeline processes ClinVar data, matches mutations with Alphamissense predictions, and evaluates its performance metrics.

---

## Requirements

### Software and Libraries
- **Python 3.8+**
- **Pandas**
- **Matplotlib**
- **Seaborn**
- **Scikit-learn**
- **Tqdm**
- **Biopython**
- **Requests**

### External Resources
1. **ClinVar Data**: A `.tsv` file containing mutation annotations from ClinVar.
2. **Alphamissense File**: The Alphamissense predictions file, specifically formatted for hg38 (`AlphaMissense_hg38.tsv`).

---

## Pipeline Overview

This pipeline is composed of **three sequential steps** that must be executed in the following order:

---

### Step 1: Canonicalize ClinVar Data

**Script**: `step1_canonize_clinvar.py`

This script extracts canonical mutations from the ClinVar data by aligning them to Ensembl's canonical transcripts.

#### Input
- **ClinVar Data File**: A `.tsv` file containing ClinVar annotations (e.g., mutation name, reference sequences, gene symbols).

#### Processing
1. Fetches canonical transcripts for each mutation using Ensembl REST API.
2. Translates canonical transcripts into proteins and applies mutations.
3. Determines the resulting protein mutation.

#### Output
- **Canonized ClinVar File**: `canonized_clinvar.csv` with added canonical mutation column.

#### Usage
```bash
python step1_canonize_clinvar.py <clinvar_file>
```

---

### Step 2: Match Alphamissense Predictions

**Script**: `step2_match_AlphaMiss.py`

This script matches ClinVar mutations with Alphamissense predictions based on chromosome, position, and protein variant.

#### Input
- **Canonized ClinVar File**: `canonized_clinvar.csv` (output of Step 1).
- **Alphamissense Predictions File**: `AlphaMissense_hg38.tsv`.

#### Processing
1. Normalizes and aligns data columns for matching.
2. Iterates through Alphamissense predictions to find matches in the ClinVar data.
3. Adds Alphamissense pathogenicity scores and classifications to the matched ClinVar entries.

#### Output
- **Matched ClinVar File**: `updated_clinvar_alphamiss.csv`.

#### Usage
```bash
python step2_match_AlphaMiss.py
```

### Step 3: Performance Analysis

**Script**: `step3_analysis.py`

This script evaluates Alphamissense performance on the matched ClinVar data by generating various visualizations and performance metrics. It provides an in-depth analysis of how Alphamissense pathogenicity predictions align with ClinVar classifications.

#### Input
- **Matched ClinVar File**: `updated_clinvar_alphamiss.csv` (output of Step 2).

#### Specific Analyses
1. **Frequency Distribution**:
   - Counts the frequency of germline classifications and visualizes it using a log-transformed bar plot.

2. **Pathogenicity Distribution**:
   - Creates violin plots to display the distribution of Alphamissense pathogenicity scores for each germline classification.
   - Overlays strip plots to highlight individual data points.

3. **Classification Accuracy**:
   - Evaluates the accuracy of Alphamissense predictions (`AM Class`) for `Benign` and `Pathogenic` variants, grouped by `Germline Review Status`.
   - Displays results as a bar plot with accuracy percentages.

4. **Confusion Matrix**:
   - Constructs a confusion matrix comparing true classifications (`Benign` or `Pathogenic`) against Alphamissense predictions (`Likely Benign`, `Likely Pathogenic`, `Ambiguous`).
   - Visualizes the matrix as a heatmap.

5. **ROC Curve**:
   - Computes the Receiver Operating Characteristic (ROC) curve for Alphamissense pathogenicity scores.
   - Identifies the optimal threshold using the Youden Index.
   - Visualizes the curve and annotates the best threshold.

#### Output
- Various performance plots saved in the current directory:
  - **Log Frequency Bar Plot**: Highlights the frequency of germline classifications.
  - **Violin Plot**: Shows Alphamissense pathogenicity distribution for each germline classification.
  - **Classification Accuracy Bar Plot**: Displays accuracy rates for `Benign` and `Pathogenic` predictions.
  - **Confusion Matrix Heatmap**: Compares actual and predicted classifications.
  - **ROC Curve**: Displays performance of Alphamissense pathogenicity scoring.

#### Customization
The script is modular and can be modified for different types of analyses or visualizations. You can:
- Add or remove analyses as needed.
- Change plot aesthetics (e.g., colors, sizes, labels).
- Use other metrics to evaluate performance.
- Filter the data for specific genes, diseases, or classifications.

#### Usage
```bash
python step3_analysis.py

