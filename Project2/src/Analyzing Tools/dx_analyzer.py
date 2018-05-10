import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def read_parameters(folder):
    """
    Take the folder where the parameter file is found, and returns
    a dictionary with the name of the variables as keys,
    and the value of parameter as the value.
    """
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

class getVariance:
    """
    Is is an ad hoc made class to get the Evolving variance of different data,
    with different dx(dt). This is in no way standarized, so to get it to work
    the data has to be in files with the same names as in the init.
    """
    def __init__(self,folder):
        self.folder = folder
        self.parameters = read_parameters(folder)
        self.list_of_dx = [1,0.5,0.1,0.05,0.01]
        self.file_names = []
        for dx in self.list_of_dx:
            self.file_names.append("data_%.6f.bin" %(dx))
        self.data = np.zeros((self.parameters["MC_cycle"],len(self.list_of_dx)))
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
        plt.title("Evolving Error",fontsize=15)
        #plt.title("Evolving Error with Importance Sampling",fontsize=15)
        plt.xlabel(r"Iterations",fontsize=20)
        plt.ylabel(r"Error",fontsize=20)
        plt.legend(loc=4)
        plt.show()


if __name__ == '__main__':
    getVar = getVariance("../output/")
