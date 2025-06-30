import os
import csv
import pandas as pd
from collections import defaultdict

def get_mean_values_from_csv_folder(folder_path):
    """
    Iterate over CSV files in the given folder and calculate mean values for each metric.
    
    Args:
        folder_path (str): Path to the folder containing CSV files
    
    Returns:
        dict: Dictionary with metric names as keys and mean values as values
    """
    metric_values = defaultdict(list)
    
    # Check if folder exists
    if not os.path.exists(folder_path):
        print(f"Folder {folder_path} does not exist")
        return {}
    
    # Get all CSV files in the folder
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]
    
    if not csv_files:
        print(f"No CSV files found in {folder_path}")
        return {}
    
    # Read each CSV file and collect values
    for csv_file in csv_files:
        file_path = os.path.join(folder_path, csv_file)
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    metric = row['metric']
                    try:
                        value = float(row['value'])
                        metric_values[metric].append(value)
                    except (ValueError, KeyError):
                        # Skip non-numeric values or missing keys
                        continue
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue
    
    # Calculate mean for each metric
    mean_values = {}
    for metric, values in metric_values.items():
        if values:  # Only calculate mean if we have values
            mean_values[metric] = sum(values) / len(values)
    
    return mean_values

def print_mean_results(folder_path):
    """
    Print the mean results in a formatted way
    """
    mean_values = get_mean_values_from_csv_folder(folder_path)
    
    if not mean_values:
        print("No valid data found")
        return
    
    # Save mean values to CSV
    output_csv_path = f'results/{folder_path.split("/")[-1]}_mean_metrics.csv'
    with open(output_csv_path, 'w', newline='') as csvfile:
        fieldnames = ['metric', 'mean_value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for metric, mean_val in mean_values.items():
            writer.writerow({'metric': metric, 'mean_value': mean_val})
    
    print(f"Mean values saved to {output_csv_path}")
    
    return mean_values

# Example usage:
print_mean_results('results/MP5')
print_mean_results('results/AF2')

