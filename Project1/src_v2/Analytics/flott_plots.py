import numpy as np
import matplotlib.pyplot as plt

alpha = np.fromfile("../output/alphas.bin",sep=" ")
data = np.fromfile("../output/data.bin",sep=" ")

print("min alpha: ",alpha[np.argmin(data)])
print("min data: ",np.min(data))

plt.plot(alpha,data,".")

plt.show()
exit()

r = np.fromfile("../output/r_positions.bin",sep=" ")
psi2 = np.fromfile("../output/psi_squared.bin",sep=" ")



psi2 = psi2[r.argsort()]
r = np.sort(r)

N  = 0

for i in range(1,len(r)):
    N += psi2[i]*(r[i]-r[i-1])

print(N)

plt.plot(r,psi2)

plt.show()
