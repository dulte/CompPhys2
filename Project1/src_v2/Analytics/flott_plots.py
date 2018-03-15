import numpy as np
import matplotlib.pyplot as plt

alpha = np.fromfile("../output/alphas.bin",sep=" ")
data = np.fromfile("../output/data.bin",sep=" ")
plt.plot(alpha,data,".")

plt.show()

