import matplotlib.pyplot as plt
import pandas as pd

inputData = pd.read_csv('data.csv')

plt.plot(inputData['n'].tolist(), inputData['SecondsPerPoint'].tolist(), 'r--o')
plt.xlabel('Number of points')
plt.ylabel('Seconds per point')
plt.xscale('log')
plt.show()