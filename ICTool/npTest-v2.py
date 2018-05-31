# -*- coding: utf-8 -*-
"""
Created on Wed May 23 15:04:36 2018

@author: mcgerten
"""
import numpy as np
import math
import scipy as sc
from scipy import sparse

''' MATLAB equivelent
  for i = 1 : k
      p        = subDAG(i, :); % parent term(s)
      support  = all(subA(:, p), 2);
      S        = sum(support);
      subia(i) = sum(support & subA(:, i)) / S;
  end
'''
def main():
    ic = {}

    DAG = np.matrix([[0,0,0,0,0],[1,0,0,0,0],[1,0,0,0,0],[0,1,1,0,0],[0,0,1,0,0]])  
    A = np.matrix([[1,1,1,1,0],[1,0,1,0,1],[1,1,1,0,0],[1,0,1,0,0]])
    
    print("This it TERM to TERM")
    print(DAG)
    print(DAG.shape[0])
    print(DAG.shape[1])
    #
    print ("This is protein to Term")
    print(A)
    print(A.shape[0])
    print(A.shape[1])
    # 0 is columns
    # 1 is rows
    
    #Convert to Sparse    
    DAG = sparse.csr_matrix(DAG)
    A   = sparse.csc_matrix(A)  
    
    for i in range(0, DAG.shape[0]):
        # Make holder matrix
        #B = np.zeros((A.shape[0], 1))
        B = sparse.csc_matrix((A.shape[0], 1))
        D = sc.sparse.dok_matrix((A.shape[0], 1))
        #print("Matrix B when created")
        #print(B)
        #print(B.shape[0])
        #print(B.shape[1])
        
        parents = DAG.getrow(i)
        #print (" These are the parents of the {} term".format(i+1))
        #print (parents)
        #print(parents.shape)
        s = parents.shape[1]
        #I know know who the direct ancestor of a Term is/are
        #For those parents, see how often they occur in A
        # ie. how many proteins are annotated for a term
        # Sum for all terms 
        
        
        #print("parents size should be 5")
        ##### BUILD B
        # Following matlab to a T -> 
        '''
        B = np.all(A.getcol(parents),axis = 1)
        '''
        # Modify

        # DOesn't work, What i really want is for all to be on all parents
        parents = parents.todok()
        #OLD WAY REALLY SLOW
        for j in range(0, s):
            # If it is a parent
            #print("This is the value of the {}th item".format(j+1))
            #print(parents.getcol(j))
            p = parents.get((0,j))
            #print("LOOK HERE")
            #print(np.all(A.getcol(j),0))
            if p != 0:
                # Count how many times that parent occurs in proteins
                C = A.getcol(j)#.todok()
                # Adding two CSC matrices
                B = B + C
                del C
            else:
                pass
        # Now i have the indexes of all parents
        # Also B will only have values where valid 
        # may be large values for overlap, but should work
        
        # For everything in B, I find largest m
        # FIND MAX
        m = B.max()
        #print("MAX")
        #print(m)
        B = B.todok()
        # Change m values to 1, all others to 0
        for k in range(0, B.shape[0]):
            # For the sparcity warning
            if B.get((k,0)) != m:
                pass
            else:
                D[k,0] = 1        
        # Now B is correct for math
        # Once have done all parents
        #print("Matrix B")                
        #print (B)
        #print(B.getformat())
        a = A.getcol(i).todok()
        b = D
        #print ("A column I")
        #print (a)
        #print ("B")
        #print (b)
        #print ("Multiply")
        c = a.multiply(b)
        #print (c)
        # Get Num/Dom
        numerator   = np.sum(c)
        denominator = np.sum(b)
        
        print("NUMERATOR")
        print (numerator)
        print("DENOMINATOR")
        print (denominator)
        if numerator > denominator:
            print("Weird")
            ic[i] = 0
        else:
            # Calculate IC
            print("IC Value")
            if denominator != 0 and numerator != 0:
                ic[i] = numerator / denominator
                try:
                    print(-math.log(ic[i], 2))
                    ic[i] = -math.log(ic[i], 2)
                except ValueError: 
                    print("Log Error")
            else:
                print ("0")
                ic[i] = 0
    print (ic)
    return ic
        
    
    
if __name__ == '__main__':
   main()
   