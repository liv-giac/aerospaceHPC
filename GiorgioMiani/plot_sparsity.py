import pandas as pd
import matplotlib.pyplot as plt

# Read the matrix from CSV file
matrix = pd.read_csv('matrix.csv')

# Convert matrix values to absolute values
matrix = matrix.abs()

# Plot sparsity
plt.spy(matrix)
plt.title('Matrix')
plt.xlabel('Columns')
plt.ylabel('Rows')
plt.xlim(0, matrix.shape[1])
plt.show()
