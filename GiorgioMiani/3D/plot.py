import numpy as np
import matplotlib.pyplot as plt

# Read the matrix_sparsity.dat file
data = np.loadtxt('matrix_sparsity.dat')

# Extract the row and column indices of non-zero elements
row_indices = data[:, 0]
col_indices = data[:, 1]

# Determine the size of the matrix
num_rows = int(np.max(row_indices)) + 1
num_cols = int(np.max(col_indices)) + 1

# Create a sparse matrix using a scatter plot
plt.scatter(col_indices, num_rows - row_indices, marker='.', color='black')

# Set the plot limits and labels
plt.xlim(0, num_cols)
plt.ylim(0, num_rows)
plt.xlabel('Column Index')
plt.ylabel('Row Index')
plt.title('Sparsity Structure of Matrix')

# Show the plot
plt.show()
