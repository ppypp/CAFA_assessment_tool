# -*- coding: utf-8 -*-
"""
Created on Wed May 23 15:04:36 2018

@author: mcgerten
"""
import numpy as np
import math
    

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
    #
    print ("This is protein to Term")
    print(A)
    # 0 is columns
    # 1 is rows
        
    for i in range(0, DAG.shape[0]):
        # Make holder matrix
        B = np.zeros((A.shape[0],1))
        parents = DAG[i, :]
        print (" These are the parents of the {} term".format(i+1))
        print (parents)
        #I know know who the direct ancestor of a Term is/are
        #For those parents, see how often they occur in A
        # ie. how many proteins are annotated for a term
        # Sum for all terms 
        t = 0
        #print("parents size should be 5")
        print(parents.shape)
        for j in range(0, parents.shape[1]):
            # If it is a parent
            #print("I printed the loop index")
            #print(i)
            print("This is the value of the {}th item".format(j+1))
            print(parents.item((0,j)))
            if parents.item((0,j)) == 1:
                # Count how many times that parent occurs in proteins
                B = B + A[:, j]
                print("Matrix B")                
                print (B)
                ####t += np.sum(A[:, j]) #################
                # HERE IS THE PROBLEM
                # I am double counting 
                # need to combine first, letting duplicates get merged
                # then do the count
                
                #print("T is {} currently".format(t))
                
            else:
                pass
        # Now i have the indexes of all parents
        # Also B will only have values where valid 
        # may be large values for overlap, but should work
        
        # For everything in B, ifind largest m
        m = 0        
        for k in range(0, B.shape[0]):
            print("K is {}".format(k))
            if B.item((k,0)) > 0:
                if B.item((k,0)) > m:
                    m = B.item((k,0))
                print("Value is currently {}".format(B.item((k,0))))
                #B[k][0] = 1
            else:
                B[k][0] = 0
        # Change m values to 1, all others to 0
        for k in range(0, B.shape[0]):
            if B.item((k,0)) == m:
                B[k][0] = 1
            else:
                B[k][0] = 0        
        # Now B is correct for math
        # Once have done all parents
        print("Matrix B")                
        print (B)
        denominator = np.sum(B)
        
        # Get self
        Self = np.sum(A[:, i])
        #print("SELF VALUE")
        #print (Self)
        # Get numerator
        numerator = Self 
        print("NUMERATOR")
        print (numerator)
        print("DENOMINATOR")
        print (denominator)
        # Calculate IC
        print("IC Value")
        if denominator != 0:
            ic[i] = numerator / denominator
            try:
                print(-math.log(ic[i], 2))
                ic[i] = -math.log(ic[i], 2)
            except ValueError: 
                print("Log Error")
        else:
            print ("0")
            ic[i] = 0
    return ic
        
    
    
if __name__ == '__main__':
   main()
   