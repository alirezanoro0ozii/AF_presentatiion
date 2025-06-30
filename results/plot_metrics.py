import pandas as pd
import matplotlib.pyplot as plt

# Read the mean metrics CSV files
af2_df = pd.read_csv('results/AF2_mean_metrics.csv')
mp5_df = pd.read_csv('results/MP5_mean_metrics.csv')

# Exclude non-numeric metrics (like 'sequence' if present)
def filter_numeric(df):
    return df[pd.to_numeric(df['mean_value'], errors='coerce').notnull()]

af2_df = filter_numeric(af2_df)
mp5_df = filter_numeric(mp5_df)

# Merge the two dataframes on the 'metric' column
merged_df = pd.merge(af2_df, mp5_df, on='metric', suffixes=('_AF2', '_MP5'))

# Set metric as index for easier plotting
merged_df.set_index('metric', inplace=True)

# Create separate plots for each metric
for metric in ['tm_score', 'rmsd_value', 'seq_id', 'CA_rmsd', 'lddt', 'gdt']:
    if metric in merged_df.index:
        plt.figure(figsize=(8, 6))
        
        metric_data = merged_df.loc[metric, ['mean_value_AF2', 'mean_value_MP5']]
        bars = plt.bar(['AF2', 'MP5'], metric_data.values, color=['skyblue', 'lightcoral'])
        
        # Set title based on metric
        if metric == 'rmsd_value':
            plt.title(f'{metric.upper()} Comparison (Lower is Better)')
        elif metric in ['tm_score', 'seq_id', 'lddt', 'gdt']:
            plt.title(f'{metric.upper()} Comparison (Higher is Better)')
        else:
            plt.title(f'{metric.upper()} Comparison')
        
        plt.ylabel(metric.upper())
        plt.ylim(0, max(metric_data.values) * 1.2)
        
        # Add value labels on bars
        for bar, value in zip(bars, metric_data.values):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
                    f'{value:.3f}', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(f'results/{metric}_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
