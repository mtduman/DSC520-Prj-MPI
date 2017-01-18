
# coding: utf-8

# In[4]:

import numpy as np
from matplotlib import pyplot as plt

def Homework5_p2c(path, data_filenames):
    ### import data ###
    data = np.loadtxt(path+data_filenames+'.txt', delimiter=',', usecols=range(8))

    ### Spedup Calculation #############
    # For spedup I used  http://stackoverflow.com/questions/7925510/parallel-speedup-with-openmp
    # sequential time measurement: seq_time = endTime-startTime;
    # Parallel time measurement : paralleltime = endTime-startTime; 
    # speedup = seq_time/paralleltime;
    seq_time = data[0,4] 
    new_col = np.zeros((len(data), 1))
    for i,j in  enumerate(data):
        new_col[i] = seq_time / data[i,4]
    data = np.c_[data, new_col]

    print(data)

    ### plot results ###
    fig = plt.figure()
    plt.loglog(data[:,0]*data[:,1],data[:,8],'r*',markersize=20,label='Weak Scaling Test')

    plt.legend(loc='upper right')
    plt.xlabel('Total Random Samples (N * cores)')
    plt.ylabel('Speedup')
    plt.title('MPI Scaling Test')
    plt.savefig(path+data_filenames+'_1.png' , bbox_inches='tight')

    fig = plt.figure()
    plt.loglog(data[:,0]*data[:,1], data[:,5],'r*',markersize=14,label='Ring Passing Time')
    plt.loglog(data[:,0]*data[:,1], data[:,6],'bo',markersize=14,label='Reduction Time')

    plt.legend(loc='upper left')
    plt.xlabel('Total Random Samples (N * cores)')
    plt.ylabel('Time')
    plt.title('MPI: Ring Passing Time & Reduction vs Total Random Samples', size=12)
    plt.savefig(path+data_filenames+'_2.png' , bbox_inches='tight')

    plt.close('all')
    return

def run():    
    path = "/Users/ekinezgi/Documents/UmassD/DSC520/umassd-hpc-mehmetduman/HW5/"
    fnames_p2c = "hw5_p2c"

    Homework5_p2c(path,fnames_p2c) 
    return
    
if __name__ == '__main__':
    run()
