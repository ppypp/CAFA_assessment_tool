
import math
import pickle as cp
import os
import argparse
import sys
import yaml
from multiprocessing import Process
import ICHelper
import numpy as np
import scipy
def createDAG():
    return

def createA():
    return

def main():
    '''
    Get a GAF, Convert, Calculate, and output. 
    '''
    # Read Config
    
    # Create DAG : m-by-m adjancency matrix
    # DAG(i, j) means i has relationship to j
    DAG = createDAG()
    # Create A   : n-by-m, Ontology Annotation matrix
    # A(i,j) = True, protein i  has annotation j
    A = createA()
    # Validate DAG is square
    if DAG.shape[0] != DAG.shape[1]:
        #Die
        raise ValueError() 
    # Get size of DAG
    m = DAG.shape[0] # 1st dimension of DAG
    
    # Validate A has m columns
    
    
   ############ make subDAGs    
    
    subDAG = np.zeros(m)
    # For each ontology
    for ontology in ['BPO','CCO','MFO']:
        k = subDAG.shape[0]
        # Create a zeroed 1-by-k
        subIC = np.zeros(k)
        # Add fake protein to SubA to prevent infinite IC
        # subA = [subA, np.ones(k)] # This is commented out in matlab code
        # Loop through all k terms in SubDAG        
        i = 0
        while i < k:
            i += 1
            p = subDAG[i,:]
            support = np.all(subA,1)
            
            num = np.sum(support)
            dom = np.sum(support and subA,1)
        
        
        IC = np.zeros(m)

    
if __name__ == '__main__':
   main()
   
