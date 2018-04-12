#from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
from time import time
import seaborn as sns
import pandas as pd


def read_parameters(folder):
    parameters = {}
    with open(folder+"metadata.txt") as f:
        for line in f:
            words = line.split()
            if len(words) == 0:
                continue
            if words[0][0] == "#":
                continue
            parameters[words[0]] = eval(words[1])

    return parameters

class Analytics:
    def __init__(self,folder,num_proc=1,flatten=False):
        self.folder = folder
        self.num_proc = num_proc
        self.parameters = read_parameters(folder)

        self.read_stamp_data()
        self.read_data(flatten)

        self.do_statistics(flatten)




    def read_stamp_data(self):
        #Most of this is reduntant, since all the simulations use the same alpahs
        #But it is added just in case
        self.num_alphas = len(np.fromfile((self.folder + "stamp0.bin"),sep=" "))

        if self.num_proc != 1:
            self.stamps = np.zeros((self.num_alphas,self.num_proc))
            for i in range(self.num_proc):
                self.stamps[:,i] = np.fromfile(self.folder + "stamp%i.bin"%i,sep=" ")
        else:
            self.stamps = np.fromfile((self.folder + "stamp0.bin"),sep=" ")

    def read_data(self,flatten=False):
        if flatten:
            if self.num_proc != 1:
                self.data = np.zeros((self.num_alphas,self.parameters["MC_cycle"]*self.num_proc))
                for i in range(self.num_proc):
                    data_flat = np.fromfile(self.folder + "data%i.bin"%i,sep=" ")
                    print(len(data_flat)/self.num_alphas)
                    self.data[:,int(i*len(data_flat)/float(self.num_alphas)):(i+1)*int(len(data_flat)/float(self.num_alphas))] = np.reshape(data_flat,(self.num_alphas,self.parameters["MC_cycle"]))
            else:
                data_flat = np.fromfile(self.folder + "data0.bin",sep=" ")
                self.data = np.reshape(data_flat,(self.num_alphas,self.parameters["MC_cycle"]))
        else:
            if self.num_proc != 1:
                self.data = np.zeros((self.num_alphas,self.parameters["MC_cycle"],self.num_proc))
                for i in range(self.num_proc):
                    data_flat = np.fromfile(self.folder + "data%i.bin"%i,sep=" ")

                    self.data[:,:,i] = np.reshape(data_flat,(self.num_alphas,self.parameters["MC_cycle"]))
            else:
                data_flat = np.fromfile(self.folder + "data0.bin",sep=" ")
                self.data = np.reshape(data_flat,(self.num_alphas,self.parameters["MC_cycle"]))


    def do_statistics(self,flatten=False):
        if flatten:
            self.meta_average = np.zeros(self.num_alphas)
            self.meta_error = np.zeros(self.num_alphas)
            if self.num_proc != 1:
                for index,a in enumerate(self.stamps[:,0]):
                    self.meta_average[index] = np.mean(self.data[index,:])
                    self.meta_error[index] = block(self.data[index,:])
            else:
                for index,a in enumerate(self.stamps):
                    self.meta_average[index] = np.mean(self.data[index,:])
                    self.meta_error[index] = block(self.data[index,:])

        else:
            self.meta_average = np.zeros(self.num_alphas)
            self.meta_error = np.zeros(self.num_alphas)
            if self.num_proc != 1:
                self.averages = np.zeros((self.num_alphas,self.num_proc))
                self.errors = np.zeros((self.num_alphas,self.num_proc))
                for i in range(self.num_proc):
                    for index,a in enumerate(self.stamps[:,i]):
                        self.averages[index,i] = np.mean(self.data[index,:,i])
                        self.errors[index,i] = block(self.data[index,:,i])
                self.meta_average = np.mean(self.averages, axis=1)
                self.meta_error = np.mean(self.errors, axis=1)
            else:
                self.averages = np.zeros(self.num_alphas)
                self.errors = np.zeros(self.num_alphas)
                for index,a in enumerate(self.stamps):
                    self.averages[index] = np.mean(self.data[index,:])
                    self.errors[index] = block(self.data[index,:])
                self.meta_average = np.copy(self.averages)
                self.meta_error = np.copy(self.errors)

    def plot_average(self):
        sns.set_style("darkgrid")
        if self.num_proc != 1:
            alphas = self.stamps[:,0]
        else:
            alphas = self.stamps

        title_text = r"Average $E_L$ for Various $\alpha$ for N = {}".format(self.parameters["N"])
        if(self.parameters["D"] != 0):
            title_text += " with Importance Sampling."

        plt.errorbar(alphas, self.meta_average,fmt='b.',yerr=np.sqrt(self.meta_error),ecolor='r')
        plt.title(title_text,fontsize=15)
        plt.xlabel(r"$\alpha$",fontsize=20)
        plt.ylabel(r"$E_L [\hbar \omega]$",fontsize=20)
        plt.show()


class SingleAlphaAnalytics:
    def __init__(self,folder,num_proc=1):
        self.folder = folder
        self.num_proc = num_proc
        self.parameters = read_parameters(folder)

        self.read_stamp_data()
        self.read_data()

        self.do_statistics()

    def read_data(self):
        if self.num_proc != 1:
            self.data = np.zeros(self.parameters["MC_cycle"]*self.num_proc)
            for i in range(self.num_proc):
                data_flat = np.fromfile(self.folder + "data%i.bin"%i,sep=" ")

                self.data[int(i*len(data_flat)):(i+1)*int(len(data_flat))] = data_flat
        else:
            self.data = np.fromfile(self.folder + "data0.bin",sep=" ")


    def read_stamp_data(self):
        self.alpha = np.fromfile((self.folder + "stamp0.bin"),sep=" ")[0]

    def do_statistics(self):

        print("Average: ", np.mean(self.data))
        print("Error: ", np.sqrt(block(self.data)))

        with open(self.folder + "benchmark.txt","w+") as f:
            f.write("Average: "+ str(np.mean(self.data)) + "\n")
            f.write("Error: "+ str(np.sqrt(block(self.data))) + "\n")



class getVariance:
    def __init__(self,folder):
        self.folder = folder
        self.parameters = read_parameters(folder)
        self.list_of_dx = [0.0005,0.001,0.005,0.01,0.05,0.1]
        self.file_names = ["data00005.bin","data0001.bin","data0005.bin","data001.bin",\
                "data005.bin","data01.bin"]
        #self.list_of_dx = [0.05,0.1,0.5,1]
        #self.file_names = ["data005.bin","data01.bin","data05.bin","data1.bin"]

        self.data = self.data = np.zeros((self.parameters["MC_cycle"],len(self.list_of_dx)))
        self.read_data()
        self.plot_variance()

    def read_data(self):
        for index,f in enumerate(self.file_names):
            self.data[:,index] = np.fromfile(self.folder + f,sep=" ")


    def plot_variance(self):
        sns.set_style("darkgrid")

        for i,dx in enumerate(self.list_of_dx):
            length = int(self.data.shape[0]/2**10)
            self.variance = np.zeros(length-1)
            iterations = np.linspace(0,self.data.shape[0],length-1,endpoint=True)
            for j in range(length-1):
                print(i,j)
                self.variance[j] = np.std(self.data[0:(j+1)*2**10,i])
            plt.plot(iterations,self.variance,label="dx = %g"%dx)
        plt.title("Evolving Error with Importance Sampling",fontsize=15)
        plt.xlabel(r"Iterations",fontsize=20)
        plt.ylabel(r"Error",fontsize=20)
        plt.legend(loc=4)
        plt.show()




def block(x):
    # preliminaries
    n = len(x);
    if abs(np.log2(n) - int(np.log2(n))) > 1e-5:
        print("Number of MC_cycles are not on the form 2^n. Change that " \
                            +"and come back an other time")

        print(np.log2(n))
        exit()
    d = int(np.log2(n)); s, gamma = np.zeros(d), np.zeros(d);
    mu = np.mean(x); t0 = time()

    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in np.arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*np.sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = np.var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])

    # generate the test observator M_k from the theorem
    M = (np.cumsum( ((gamma/s)**2*2**np.arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =np.array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in np.arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    ans = s[k]/2**(d-k)
    # print("Runtime: %g sec" % (time()-t0)); print("Blocking Statistics :")
    # print("average            iterations      std. error")
    # print("%8g %20g %15g" % (mu, k, ans**.5))
    return ans

if __name__ == '__main__':
    an = Analytics("../output/",1,flatten=True)
    an.plot_average()

    #an = SingleAlphaAnalytics("../output/",4)

    #an = getVariance("../output/varIS/")
