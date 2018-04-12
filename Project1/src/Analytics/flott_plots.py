import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

# alpha = np.fromfile("../output/stamp0.bin",sep=" ")
# data = np.fromfile("../output/data.bin",sep=" ")
#
# print("min alpha: ",alpha[np.argmin(data)])
# print("min data: ",np.min(data))
#
# plt.plot(alpha,data,".")
#
# plt.show()
#
#
# exit()

r = np.fromfile("../output/r_positions.bin",sep=" ")
rho = np.fromfile("../output/density.bin",sep=" ")
#rho_non = np.fromfile("../output/density_non.bin",sep=" ")
volume_factor = np.fromfile("../output/volume.bin",sep=" ")


rs = np.linspace(r[0],r[1],int(r[2]),endpoint=True)



step = r[3]

N = np.sum(rho[:len(rs)]*volume_factor[:len(rs)]*rs)

#N_non = np.sum(rho_non*step)

print("For interacting N = ",N)
#print("For noninteracting N = ",N_non)

print(10*rho[0]*step)
print(np.mean(10*rho))
print(np.sum(rho[:len(rs)]*rs))

plt.plot(rs,10*rho[:len(rs)],label="Interacting")
#plt.plot(rs,10*rho_non[:len(rs)],label="Noninteracting")
plt.title(r"Onebody density for N=10",fontsize=15)
plt.xlabel("$r$",fontsize=20)
plt.ylabel(r"$\rho(r)$",fontsize=20)
plt.xlim(0,4)
plt.legend()
plt.show()
