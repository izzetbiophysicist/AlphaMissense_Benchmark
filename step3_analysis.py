import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.metrics import auc

# Load the dataset
data = pd.read_csv("/home/lucas/AlphaMissense_Benchmark/cataract_benchmark/updated_clinvar_alphamiss.csv")

# Select the required columns and create a new DataFrame
selected_columns = ['Gene(s)', 'Canonic mutation', 'Germline classification', 
                    'Germline review status', 'am_class', 'am_pathogenicity']
subset_df = data[selected_columns].copy()

# Rename the columns for clarity
subset_df.columns = ['Gene', 'Variant', 'Germline Classification', 
                     'Germline Review Status', 'AM Class', 'AM Pathogenicity']

# Clean the data: remove rows with NaN in key columns
subset_df = subset_df.dropna(subset=['Germline Classification', 'AM Pathogenicity', 'AM Class', 'Germline Review Status'])

# Ensure AM Pathogenicity is numeric
subset_df['AM Pathogenicity'] = pd.to_numeric(subset_df['AM Pathogenicity'], errors='coerce')

# Remove rows where AM Pathogenicity is NaN or outside the range [0, 1]
subset_df = subset_df[(subset_df['AM Pathogenicity'] >= 0) & (subset_df['AM Pathogenicity'] <= 1)]

# Define the custom order of classes
class_order = [
    'Uncertain significance',
    'Likely benign',
    'Benign/Likely benign',
    'Benign',
    'Likely pathogenic',
    'Pathogenic/Likely pathogenic',
    'Pathogenic',
    'Conflicting classifications of pathogenicity',
    'Uncertain significance/Uncertain risk allele',
    'not provided',
    'Pathogenic/Likely pathogenic/Pathogenic, low penetrance'
]

# Reorder the Germline Classification column
subset_df['Germline Classification'] = pd.Categorical(
    subset_df['Germline Classification'], 
    categories=class_order, 
    ordered=True
)

# Count frequencies of each Germline Classification, ensuring all categories are included
freq_data = subset_df['Germline Classification'].value_counts().reindex(class_order, fill_value=0)

# Replace NaN or inf values with zero in the log-transformed data
log_freq_data = np.log10(freq_data.replace(0, np.nan))  # Replace 0 with NaN for log calculation
log_freq_data = log_freq_data.replace([np.inf, -np.inf, np.nan], 0)  # Replace inf or NaN with 0

# Create the figure with two subplots
fig, axes = plt.subplots(2, 1, figsize=(18, 16), gridspec_kw={'height_ratios': [1, 3]})

# Bar plot for frequencies (log scale with zero for 0 frequencies)
axes[0].bar(range(len(class_order)), log_freq_data, color='skyblue', alpha=0.8)
axes[0].set_title('Log Frequency of Germline Classifications', fontsize=16)
axes[0].set_ylabel('Log10(Frequency)', fontsize=14)
axes[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)  # Hide x-axis labels
axes[0].set_xlim(-0.5, len(class_order) - 0.5)  # Align x-axis with violin plot

# Violin plot for AM Pathogenicity
sns.violinplot(
    data=subset_df, 
    x='Germline Classification', 
    y='AM Pathogenicity', 
    scale='width',
    cut=0,
    bw_adjust=0.5,
    palette="muted",
    order=class_order,
    inner=None,
    ax=axes[1]
)

# Overlay strip plot for individual points
sns.stripplot(
    data=subset_df, 
    x='Germline Classification', 
    y='AM Pathogenicity', 
    color='k', 
    size=2.5,
    alpha=0.5,
    dodge=True,
    order=class_order,
    ax=axes[1]
)

axes[1].set_title('Distribution of AM Pathogenicity by Germline Classification', fontsize=16)
axes[1].set_xlabel('Germline Classification', fontsize=14)
axes[1].set_ylabel('AM Pathogenicity', fontsize=14)
axes[1].set_ylim(-0.1, 1.1)

# Adjust x-axis labels to be angled
axes[1].tick_params(axis='x', labelsize=12)
for label in axes[1].get_xticklabels():
    label.set_rotation(45)
    label.set_horizontalalignment('right')

# Increase margins around the figure
plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.25)

# Save the plots
plt.savefig('germline_classification_violin_bar_plot_updated.png', dpi=300, bbox_inches='tight')
plt.show()

########################
## Confusion Matrix
########################

# Filter data for "reviewed by expert panel" and relevant classes
expert_panel = subset_df[
    (subset_df['Germline Review Status'] == 'reviewed by expert panel') &
    (subset_df['Germline Classification'].isin(['Benign', 'Pathogenic']))
].copy()

# Map classes to simplified labels
germline_label_mapping = {
    'Benign': 'Benign',
    'Pathogenic': 'Pathogenic'
}
am_class_label_mapping = {
    'likely_benign': 'Likely Benign',
    'likely_pathogenic': 'Likely Pathogenic',
    'ambiguous': 'Ambiguous'
}
expert_panel['Germline Classification'] = expert_panel['Germline Classification'].map(germline_label_mapping)
expert_panel['AM Class'] = expert_panel['AM Class'].map(am_class_label_mapping)

# Ground truth labels and predicted labels
true_labels = expert_panel['Germline Classification']
predicted_labels = expert_panel['AM Class']

# Define unique labels for rows (ground truth) and columns (predictions)
true_classes = ['Benign', 'Pathogenic']
predicted_classes = ['Likely Benign', 'Likely Pathogenic', 'Ambiguous']

# Re-index confusion matrix data to match the desired structure
conf_matrix = pd.crosstab(
    true_labels,
    predicted_labels,
    rownames=['True Labels'],
    colnames=['Predicted Labels'],
    dropna=False
).reindex(index=true_classes, columns=predicted_classes, fill_value=0)

# Plot the confusion matrix
plt.figure(figsize=(10, 8))
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', xticklabels=predicted_classes, yticklabels=true_classes)
plt.title('Confusion Matrix for "Reviewed by Expert Panel"', fontsize=16)
plt.xlabel('Predicted Labels', fontsize=14)
plt.ylabel('True Labels', fontsize=14)
plt.xticks(rotation=45, ha='right', fontsize=12)
plt.yticks(rotation=0, fontsize=12)

# Save the plot
plt.savefig('confusion_matrix_expert_panel_updated.png', dpi=300, bbox_inches='tight')
plt.show()

########################
## Discrepancy Analysis
########################

# Filter for discrepancies where AM Class is "Likely Benign" and Ground Truth is "Pathogenic"
discrepancies = expert_panel[
    (expert_panel['AM Class'] == 'Likely Benign') & (expert_panel['Germline Classification'] == 'Pathogenic')
]

# Select relevant columns for output
discrepancy_details = discrepancies[['Gene', 'Variant', 'AM Class', 'Germline Classification']]
print("Discrepancy Details:\n", discrepancy_details)

########################
## ROC Curve
########################

# Filter data for "Pathogenic" and "Benign" ground truth
roc_data = subset_df[subset_df['Germline Classification'].isin(['Benign', 'Pathogenic'])].copy()
roc_data['Binary Classification'] = roc_data['Germline Classification'].map({'Benign': 0, 'Pathogenic': 1})

# Extract true labels and prediction scores
true_labels = roc_data['Binary Classification']
predicted_scores = roc_data['AM Pathogenicity']

# Define custom thresholds at fixed intervals
custom_thresholds = np.arange(0.0, 1.01, 0.01)  # From 0.0 to 1.0 in steps of 0.01
tpr, fpr = [], []

for threshold in custom_thresholds:
    predicted_classes = (predicted_scores >= threshold).astype(int)
    tp = np.sum((true_labels == 1) & (predicted_classes == 1))
    fp = np.sum((true_labels == 0) & (predicted_classes == 1))
    fn = np.sum((true_labels == 1) & (predicted_classes == 0))
    tn = np.sum((true_labels == 0) & (predicted_classes == 0))
    tpr.append(tp / (tp + fn) if (tp + fn) > 0 else 0.0)
    fpr.append(fp / (fp + tn) if (fp + tn) > 0 else 0.0)

# Calculate AUC
roc_auc = auc(fpr, tpr)
youden_index = np.array(tpr) - np.array(fpr)
best_index = np.argmax(youden_index)
best_threshold = custom_thresholds[best_index]

# Plot ROC curve
plt.figure(figsize=(10, 8))
plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
plt.scatter(fpr[best_index], tpr[best_index], color='red', label=f'Youden Index (Threshold = {best_threshold:.2f})')
plt.plot([0, 1], [0, 1], color='gray', linestyle='--', lw=1, label='Random Guess')
plt.title('ROC Curve for Pathogenicity Prediction', fontsize=16)
plt.xlabel('False Positive Rate (FPR)', fontsize=14)
plt.ylabel('True Positive Rate (TPR)', fontsize=14)
plt.legend(fontsize=12)
plt.grid(alpha=0.3)
plt.savefig('roc_curve_updated.png', dpi=300, bbox_inches='tight')
plt.show()

# Print Youden Index and Best Threshold
print(f"Best Youden Index: {youden_index[best_index]:.2f}")
print(f"Best Threshold: {best_threshold:.2f}")


### By year
#
## Cumulative Accuracy by Year of Deposition
########################

# Extract the year from the "Germline Date Last Evaluated" column
subset_df['Year of Deposition'] = pd.to_datetime(subset_df['Germline Date Last Evaluated'], errors='coerce').dt.year

# Filter data for only Benign and Pathogenic ground truths, and drop rows with NaN in the AM Class column
yearly_data = subset_df[
    subset_df['Germline Classification'].isin(['Benign', 'Pathogenic']) & 
    subset_df['AM Class'].notna()
].copy()

# Calculate Correct Predictions
yearly_data['Correct Prediction'] = (
    ((yearly_data['Germline Classification'] == 'Benign') & (yearly_data['AM Class'] == 'likely_benign')) |
    ((yearly_data['Germline Classification'] == 'Pathogenic') & (yearly_data['AM Class'] == 'likely_pathogenic'))
)

# Initialize cumulative DataFrame
cumulative_data = pd.DataFrame()

# Process each year to include cumulative data
for year in sorted(yearly_data['Year of Deposition'].dropna().unique()):
    # Filter data up to the current year
    data_up_to_year = yearly_data[yearly_data['Year of Deposition'] <= year]
    
    # Calculate total and correct predictions for the cumulative data
    total_variants = len(data_up_to_year)
    correct_variants = data_up_to_year['Correct Prediction'].sum()
    
    # Calculate accuracy
    accuracy = (correct_variants / total_variants) * 100 if total_variants > 0 else 0
    
    # Append to cumulative DataFrame
    cumulative_data = pd.concat([cumulative_data, pd.DataFrame({
        'Year of Deposition': [year],
        'Total Variants': [total_variants],
        'Correct Variants': [correct_variants],
        'Accuracy': [accuracy]
    })])

# Reset index of the cumulative data
cumulative_data.reset_index(drop=True, inplace=True)

# Plot cumulative accuracy as a function of year
plt.figure(figsize=(12, 6))
plt.plot(cumulative_data['Year of Deposition'], cumulative_data['Accuracy'], marker='o', linestyle='-', color='skyblue')

# Add titles and labels
plt.title('Cumulative Accuracy by Year of Deposition', fontsize=16)
plt.xlabel('Cutoff Date of Last Evaluation', fontsize=14)
plt.ylabel('Cumulative Accuracy (%)', fontsize=14)
plt.grid(alpha=0.3)

# Add annotations for total variants per year, ensuring spacing and including first and last labels
spacing = max(len(cumulative_data) // 10, 1)  # Adjust label spacing dynamically
for index, row in cumulative_data.iterrows():
    # Adjust placement for first label to ensure it doesn’t go above the plot
    if index == 0:
        y_offset = -2  # Move the first label down
    elif index == len(cumulative_data) - 1:
        y_offset = 4  # Move the last label slightly higher
    else:
        y_offset = 2  # Default offset for other labels

    # Ensure first, last, and spaced labels are shown
    if index == 0 or index == len(cumulative_data) - 1 or index % spacing == 0:
        plt.text(row['Year of Deposition'], row['Accuracy'] + y_offset, 
                 f'N={int(row["Total Variants"])}', ha='center', fontsize=10, color='black')

# Save the plot
plt.savefig('cumulative_accuracy_by_year_of_deposition.png', dpi=300, bbox_inches='tight')

# Display the plot
plt.show()
