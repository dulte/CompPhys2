import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


sns.set_style("darkgrid")
plt.rcParams.update({'font.size':12})
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.family']='STIXGeneral'



rates = [0.500000,0.100000,0.050000,0.010000]
Ns = [1,2,3,4,5]

colors = ["r","g","b","y"]
ticks = ["*","s",".","8","D"]

for rate in range(len(rates)):
    for N in range(len(Ns)):
        filename = "../output/gradient_data_%s_%.6f" %(Ns[N],rates[rate])
        grad = np.fromfile(filename,sep=" ")
        iterations = np.arange(1,len(grad)+1)
        plt.plot(iterations,grad,color=colors[rate],marker=ticks[N],label="N=%s;Rate=%s" %(Ns[N],rates[rate]))

plt.legend()
plt.title("Convergence for Different Ns and Learning Rates with Importance Sampling",fontsize=40)
plt.xlabel("Iterations",fontsize=40)
plt.ylabel(r"$|\nabla E_L|$",fontsize=40)
plt.ylim(0,10)
plt.show()
