import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import minimize
import random


residue_names = ["LYS", "GLU", "PHE", "TRP"]

def save_pdb(points, filename):
    """Save points to PDB file"""
    with open(filename, 'w') as f:
        for i, point in enumerate(points):
            f.write(f"ATOM      {i+1}  CA   {residue_names[i]}  {i+1}  {point[0]:.3f} {point[1]:.3f} {point[2]:.3f}  1.00  0.00           C\n")

def generate_random_points(n_points=4, bounds=(-10, 10)):
    """Generate n random points in 3D space within given bounds"""
    return np.random.uniform(bounds[0], bounds[1], (n_points, 3))

def compute_pairwise_distances(points):
    """Compute pairwise distances between all points"""
    return pdist(points)

def compute_distance_matrix(points):
    """Compute full distance matrix"""
    return squareform(pdist(points))

def objective_function(new_points, target_distances):
    """Objective function to minimize: difference between target and current distances"""
    new_points = new_points.reshape(-1, 3)
    current_distances = pdist(new_points)
    return np.sum((current_distances - target_distances) ** 2)

def reproduce_points_from_distances(target_distances, n_points=4, max_attempts=100):
    """Reproduce 4 points that have the exact pairwise distances"""
    best_points = None
    best_error = float('inf')
    
    for attempt in range(max_attempts):
        # Generate random initial guess
        initial_points = generate_random_points(n_points)
        
        # Flatten for optimization
        x0 = initial_points.flatten()
        
        # Optimize
        result = minimize(
            objective_function, 
            x0, 
            args=(target_distances),
            method='L-BFGS-B',
            bounds=[(-10, 10)] * (n_points * 3)
        )
        
        if result.success and result.fun < best_error:
            best_error = result.fun
            best_points = result.x.reshape(-1, 3)
    
    return best_points, best_error

# Generate original 4 points
print("Generating original 4 points...")
original_points = generate_random_points(4)
save_pdb(original_points, "original_points.pdb")
print("Original points:")
print(original_points)

# Compute pairwise distances
original_distances = compute_pairwise_distances(original_points)
print(f"\nPairwise distances: {original_distances}")

# Compute distance matrix
distance_matrix = compute_distance_matrix(original_points)
print(f"\nDistance matrix:")
print(distance_matrix)

# Reproduce points from distances
print("\nReproducing points from distances...")
reproduced_points, error = reproduce_points_from_distances(original_distances)
save_pdb(reproduced_points, "reproduced_points.pdb")

print(f"\nReproduced points:")
print(reproduced_points)
print(f"Error: {error}")

# Verify the reproduction
reproduced_distances = compute_pairwise_distances(reproduced_points)
print(f"\nReproduced pairwise distances: {reproduced_distances}")
print(f"Distance difference: {np.abs(original_distances - reproduced_distances)}")

# Test multiple reproductions
print("\n" + "="*50)
print("Testing multiple reproductions...")

# for i in range(5):
#     print(f"\nReproduction {i+1}:")
#     new_points, new_error = reproduce_points_from_distances(original_distances)
#     new_distances = compute_pairwise_distances(new_points)
#     print(f"Error: {new_error}")
#     print(f"Max distance difference: {np.max(np.abs(original_distances - new_distances))}")
