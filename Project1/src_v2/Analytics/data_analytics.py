#from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, floor
import numpy as np
from numpy.linalg import inv
from time import time

class Analytics:
    def __init__(self,folder,num_proc=1):
        self.folder = folder
        self.num_proc = num_proc
        self.read_parameters(folder)

        self.read_stamp_data()
        self.read_data()

        self.do_statistics()


    def read_parameters(self,folder):
        parameters = {}
        with open(folder+"metadata.txt") as f:
            for line in f:
                words = line.split()
                if len(words) == 0:
                    continue
                if words[0][0] == "#":
                    continue
                parameters[words[0]] = eval(words[1])

        self.parameters = parameters

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

    def read_data(self):

        if self.num_proc != 1:
            self.data = np.zeros((self.num_alphas,self.parameters["MC_cycle"],self.num_proc))
            for i in range(self.num_proc):
                data_flat = np.fromfile(self.folder + "data%i.bin"%i,sep=" ")

                self.data[:,:,i] = np.reshape(data_flat,(self.num_alphas,self.parameters["MC_cycle"]))
        else:
            data_flat = np.fromfile(self.folder + "data0.bin",sep=" ")
            self.data = np.reshape(data_flat,(self.num_alphas,self.parameters["MC_cycle"]))


    def do_statistics(self):
        self.meta_average = np.zeros(self.num_alphas)
        self.meta_error = np.zeros(self.num_alphas)
        if self.num_proc != 1:
            self.averages = np.zeros((self.num_alphas,self.num_proc))
            self.errors = np.zeros((self.num_alphas,self.num_proc))
            for i in range(self.num_proc):
                for index,a in enumerate(self.stamps[:,i]):
                    self.averages[index,i] = np.mean(self.data[index,:,i])
                    self.errors[index,i] = self.block(self.data[index,:,i])
            self.meta_average = np.mean(self.averages, axis=1)
            self.meta_error = np.mean(self.errors, axis=1)
        else:
            self.averages = np.zeros(self.num_alphas)
            self.errors = np.zeros(self.num_alphas)
            for index,a in enumerate(self.stamps):
                self.averages[index] = np.mean(self.data[index,:])
                self.errors[index] = self.block(self.data[index,:])
            self.meta_average = np.copy(self.averages)
            self.meta_error = np.copy(self.errors)



    def block(self,x):
        # preliminaries
        n = len(x);
        if abs(np.log2(n) - int(np.log2(n))) > 1e-5:
            print("Number of MC_cycles are not on the form 2^n. Change that " \
                                +"and come back an other time")
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
    an = Analytics("../output/",10)
    #print(an.data[:,:,0])
    print(an.meta_average)
    print(an.meta_error)
