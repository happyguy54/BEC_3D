import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Number of samples
N = 100000
sigma = 1  # Standard deviation of the Gaussian distribution

# Generate random points in 3D following a Gaussian distribution
x = np.random.normal(0, sigma, N)
y = np.random.normal(0, sigma, N)
z = np.random.normal(0, sigma, N)

# Calculate the radial distance r for each point
r = np.sqrt(x**2 + y**2 + z**2)

# Plot the histogram
plt.hist(r, bins=100, density=True, alpha=0.6, color='g', edgecolor='black')

# Plot the theoretical distribution for comparison
r_values = np.linspace(0, 4*sigma, 1000)
P_r = (r_values**2 / sigma**3) * (2 / np.sqrt(2*np.pi)) * np.exp(-r_values**2 / (2*sigma**2))
plt.plot(r_values, P_r, 'r-', lw=2)

plt.title('Histogram of Radial Distances with Gaussian Density Distribution')
plt.xlabel('Radial distance r')
plt.ylabel('Probability density')
plt.show()

# Calculate the bounds of the 68% central interval around the most probable value
# Since we are dealing with a 3D Gaussian, we find the 16th and 84th percentiles of the chi distribution with 3 degrees of freedom
from scipy.stats import chi
lower_bound = chi.ppf(0.16, df=3, scale=sigma)
upper_bound = chi.ppf(0.84, df=3, scale=sigma)

# Filter points to include only those in the positive octant
positive_octant = (x >= 0) & (y >= 0) & (z >= 0)
x = x[positive_octant]
y = y[positive_octant]
z = z[positive_octant]
r = r[positive_octant]

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot all points in the positive octant
ax.scatter(x, y, z, c='blue', s=0.1, alpha=0.3)

# Highlight points where the radial distance r is within the central 68% interval
highlight = (r >= lower_bound) & (r <= upper_bound)
ax.scatter(x[highlight], y[highlight], z[highlight], c='red', s=0.1, alpha=0.6)

# Draw the boundary of the sphere with radius 3 sigma for reference
u = np.linspace(0, np.pi / 2, 100)
v = np.linspace(0, np.pi / 2, 100)
u, v = np.meshgrid(u, v)
X = 3 * sigma * np.sin(u) * np.cos(v)
Y = 3 * sigma * np.sin(u) * np.sin(v)
Z = 3 * sigma * np.cos(u)
ax.plot_wireframe(X, Y, Z, color='k', alpha=0.2)

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('One-Eighth of a Sphere with Highlighted Region Within 68% Interval Around Most Probable r')

plt.show()