import numpy as np
import matplotlib.pyplot as plt


rates = [0.500000,0.100000,0.050000,0.010000]
Ns = [1,2,3,4,5]

colors = ["r","g","b","y"]
ticks = ["*","s",".","8","D"]

for rate in range(len(rates)):
    for N in range(len(Ns)):
        filename = "../output/Gradient Data/gradient_data_%s_%.6f" %(Ns[N],rates[rate])
        grad = np.fromfile(filename,sep=" ")
        iterations = np.arange(1,len(grad)+1)
        plt.plot(iterations,grad,color=colors[rate],marker=ticks[N],label="N=%s;Rate=%s" %(Ns[N],rates[rate]))

plt.legend()
plt.show()
