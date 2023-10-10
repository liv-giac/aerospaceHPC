import matplotlib.pyplot as plt
import pandas as pd

# Load the data from the CSV file
data = pd.read_csv("memory_used.csv")

# Extract the first two columns
x = data.iloc[:, 0]
y = data.iloc[:, 1]

# Plot the data
plt.plot(x, y)
plt.plot(x, x**3)
plt.yscale('log')
plt.xscale('log')
plt.xlabel("Number of elements")
plt.ylabel("Memory used (MB)")
plt.legend(["Memory used","x^3"])
plt.savefig('Memory_used.png')
