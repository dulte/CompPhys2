import numpy as np
import matplotlib.pyplot as plt

gradient = np.fromfile("../output/gradient_data.bin",sep=" ")
alphas = np.fromfile("../output/gradient_stamp.bin",sep=" ")

plt.plot(alpha,gradient,".")
plt.title(r"$\frac{\partial \rangle E_L \langle}{\partial \alpha}$ for differnt $\alpha$",fontsize=15)
plt.xlabel(r"$\alpha$",fontsize=20)
plt.ylabel(r"$\frac{\partial \rangle E_L \langle}{\partial \alpha}$",fontsize=20)
plt.show()
