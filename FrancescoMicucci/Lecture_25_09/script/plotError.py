import matplotlib.pyplot as plt
import pandas as pd

inputData = pd.read_csv('../data.csv', sep = ",")

plt.plot(inputData.deltaX, inputData.err, 'b--o', label="error")
plt.plot(inputData.deltaX, inputData.deltaX**1, 'k--', label="deltaX")
plt.plot(inputData.deltaX, inputData.deltaX**2, 'k-', label="deltaX^2")
plt.xlabel('delta X')
plt.ylabel('Error')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()