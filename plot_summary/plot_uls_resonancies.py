import numpy as np
import matplotlib.pyplot as plt

# Define the boundaries of the first box
out_from1, out_to1 = 0, 0.46
side_from1, side_to1 = 0.36, 0.62
long_from1, long_to1 = 0, 0.34

# Define the ranges for the boundaries of the second box
out_ranges = [(0.14, 0.18), (0.46, 0.54), (0.14, 0.54)]
side_ranges = [(0.02, 0.02), (0.1, 0.14), (0.02, 0.14)]
long_ranges = [(0.02, 0.02), (0.06, 0.1), (0.02, 0.1)]

# Define the number of points per unit size
points_per_unit = 100

# Calculate the size of the region for the first box
out_size1 = out_to1 - out_from1
side_size1 = side_to1 - side_from1
long_size1 = long_to1 - long_from1

# Calculate the number of points in each dimension for the first box
num_points_out1 = int(out_size1 * points_per_unit)
num_points_side1 = int(side_size1 * points_per_unit)
num_points_long1 = int(long_size1 * points_per_unit)

# Create a 3D grid for the first box
out1 = np.linspace(out_from1, out_to1, num_points_out1)
side1 = np.linspace(side_from1, side_to1, num_points_side1)
long1 = np.linspace(long_from1, long_to1, num_points_long1)
grid_out1, grid_side1, grid_long1 = np.meshgrid(out1, side1, long1, indexing='ij')

# Calculate the energy for each point in the first grid
energies1 = np.sqrt(grid_out1**2 + grid_side1**2 + grid_long1**2)
plt.hist(energies1.flatten(), bins=100, alpha=0.5, label=f'Box 1')

# Iterate over the different sizes for the second box
for i in range(len(out_ranges)):
    # Calculate the size of the region in each dimension for the second box
    out_size2 = out_ranges[i][1] - out_ranges[i][0]
    side_size2 = side_ranges[i][1] - side_ranges[i][0]
    long_size2 = long_ranges[i][1] - long_ranges[i][0]

    # Calculate the number of points in each dimension for the second box
    num_points_out2 = int(out_size2 * points_per_unit)
    num_points_side2 = int(side_size2 * points_per_unit)
    num_points_long2 = int(long_size2 * points_per_unit)

    # Create a 3D grid for the second box
    out2 = np.linspace(out_ranges[i][0], out_ranges[i][1], num_points_out2)
    side2 = np.linspace(side_ranges[i][0], side_ranges[i][1], num_points_side2)
    long2 = np.linspace(long_ranges[i][0], long_ranges[i][1], num_points_long2)
    grid_out2, grid_side2, grid_long2 = np.meshgrid(out2, side2, long2, indexing='ij')

    # Calculate the energy for each point in the second grid
    energies2 = np.sqrt(grid_out2**2 + grid_side2**2 + grid_long2**2)

    # Flatten the energies arrays and create the histograms
    plt.hist(energies2.flatten(), bins=100, alpha=0.5, label=f'Box 2 - Size {i+1}')
    total_energies = np.concatenate((energies1.flatten(), energies2.flatten()))

    # Create the histogram for the total energies
    plt.hist(total_energies, bins=100, alpha=0.5, label=f'Total {i+1}')

# Show the legend
plt.legend()

# Show the plot
plt.show()