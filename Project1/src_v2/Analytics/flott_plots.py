import numpy as np
import matplotlib.pyplot as plt

alpha = np.fromfile("../output/alphas.bin",sep=" ")
data = np.fromfile("../output/data.bin",sep=" ")

print("min alpha: ",alpha[np.argmin(data)])
print("min data: ",np.min(data))

plt.plot(alpha,data,".")

plt.show()


r = np.fromfile("../output/r_positions.bin",sep=" ")
rho = np.fromfile("../output/density.bin",sep=" ")



rs = np.linspace(r[0],r[1],int(r[2]),endpoint=True)

print(len(rs),len(rho[:len(rs)]))

plt.plot(rs,rho[:len(rs)])
plt.show()
