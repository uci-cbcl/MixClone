#!/usr/bin/env python

import sys

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from mixclone.model.utils import *
from mixclone import constants

def main():
    mu_N = 0.5
    c_N = 2
    phi = np.linspace(0, 1, 100)
    
    fig = plt.figure()
    #ax = fig.gca(projection='3d')
    
    mu_G = 0.5
    c_H = 0
    mu_E = get_mu_E(mu_N, mu_G, c_N, c_H, phi)
    c_E = get_c_E(c_N, c_H, phi)
    #ax.plot(phi, mu_E, c_E, label='NULL', color='r')
    plt.plot(mu_E, c_E, label='NULL', color='r')
    
    mu_G = 0
    c_H = 1
    mu_E = get_mu_E(mu_N, mu_G, c_N, c_H, phi)
    c_E = get_c_E(c_N, c_H, phi)
    #ax.plot(phi, mu_E, c_E, label='A', color='b')
    plt.plot(mu_E, c_E, label='A', color='b')
    
    mu_G = 0.5
    c_H = 2
    mu_E = get_mu_E(mu_N, mu_G, c_N, c_H, phi)
    c_E = get_c_E(c_N, c_H, phi)
    #ax.plot(phi, mu_E, c_E, label='AB', color='k')
    plt.plot(mu_E, c_E, label='AB', color='k')
    
    mu_G = 0
    c_H = 2
    mu_E = get_mu_E(mu_N, mu_G, c_N, c_H, phi)
    c_E = get_c_E(c_N, c_H, phi)
    #ax.plot(phi, mu_E, c_E, label='AA', color='y')
    plt.plot(mu_E, c_E, label='AA', color='y')

    mu_G = 1./4
    c_H = 4
    mu_E = get_mu_E(mu_N, mu_G, c_N, c_H, phi)
    c_E = get_c_E(c_N, c_H, phi)
    #ax.plot(phi, mu_E, c_E, label='AAAB', color='c')
    plt.plot(mu_E, c_E, label='AAAB', color='m')
    
    mu_G = 1./3
    c_H = 3
    mu_E = get_mu_E(mu_N, mu_G, c_N, c_H, phi)
    c_E = get_c_E(c_N, c_H, phi)
    #ax.plot(phi, mu_E, c_E, label='AAB', color='g')
    plt.plot(mu_E, c_E, label='AAB', color='g')
    

    
    
    
    
    
    #ax.legend()
    #ax.set_xlabel('cellular prevalence')
    #ax.set_xlim(0, 1)
    #ax.set_ylabel('BAF')
    #ax.set_ylim(0, 0.5)
    #ax.set_zlabel('copy number')
    #ax.set_zlim(0, 4)
    plt.legend()
    plt.xlabel('BAF')
    plt.xlim(0, 0.5)
    plt.ylabel('copy number')
    plt.ylim(0, 4)
    
    plt.show()
    
    
    
    
    
    
if __name__ == '__main__':
    main()