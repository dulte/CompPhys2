import numpy as np
import matplotlib.pyplot as plt

alpha = np.fromfile("../output/stamp0.bin",sep=" ")
data = np.fromfile("../output/data.bin",sep=" ")

print("min alpha: ",alpha[np.argmin(data)])
print("min data: ",np.min(data))

plt.plot(alpha,data,".")

plt.show()


exit()

r = np.fromfile("../output/r_positions.bin",sep=" ")
rho = np.fromfile("../output/density.bin",sep=" ")
rho_non = np.fromfile("../output/density_non.bin",sep=" ")



rs = np.linspace(r[0],r[1],int(r[2]),endpoint=True)

print(len(rs),len(rho[:len(rs)]))

step = r[3]

N = np.sum(rho*step)

N_non = np.sum(rho_non*step)

print("For interacting N = ",N)
print("For noninteracting N = ",N_non)

print(10*rho[0]*step)

plt.plot(rs,10*rho[:len(rs)]*step,label="Interacting")
plt.plot(rs,10*rho_non[:len(rs)]*step,label="Noninteracting")
plt.legend()
plt.show()
