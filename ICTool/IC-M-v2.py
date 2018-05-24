import math
import pickle as cp
import os
import argparse
import sys
import yaml
import ICHelper
import numpy as np
import scipy as sc

def main():
    # Read OBO
    # Build DAG
    DAG, LIST = makeDAG()
    #Read GAF-equivelent 
    # Build A
    A = makeA(LIST)
    IC = calculateIC(DAG, A, LIST)
    
    
def makeDAG():
    '''
    Read in OBO and make a list of terms for each ontology
    '''
    # Parse OBO
    obo_file = open("/home/mcgerten/cafa3/CAFA3_test_data/gene_ontology_edit.obo.2016-06-01","r")
    allGOterms = obo_file.read().split("[Term]")
    print("Parsed OBO")
    # Intialize variables
    GO_term   = ""
    namespace = ""
    # Lists for each namespace
    BPO  = []
    CCO  = []
    MFO  = []
    # 
    LIST = {}
    
    # Add terms to respective namspace lists
    for term in allGOterms[1:]:
        split_term = term.split("\n")
        for line in split_term:
            if "id:" in line and "GO:" in line and "alt_id" not in line:
                GO_term = "GO:"+line.split("GO:")[-1].strip()
            if   "namespace: biological_process" in line:
                namespace = "BPO"
                BPO.append(GO_term)
            elif "namespace: cellular_component" in line:
                namespace = "CCO"   
                CCO.append(GO_term)
            elif "namespace: molecular_function" in line:
                namespace = "MFO"
                MFO.append(GO_term)
    LIST["BPO"] = BPO
    LIST["CCO"] = CCO
    LIST["MFO"] = MFO
    print("Done making lists")
    ######### DONE MAKING LIST ###########################    
    alt_id_to_id = dict()
    DAG = {}    
    # 3 1D arrays to build COO sparse matrix
        
    BPO = LIST["BPO"] # makes sure same index is used for a term thoughout
    ROW_INDEX_BPO = np.empty(len(BPO), dtype = int)
    COL_INDEX_BPO = np.empty(len(BPO), dtype = int)
    DATA_BPO      = np.empty(len(BPO), dtype = int)
    
    CCO = LIST["CCO"]
    ROW_INDEX_CCO = np.empty(len(CCO), dtype = int)
    COL_INDEX_CCO = np.empty(len(CCO), dtype = int)
    DATA_CCO      = np.empty(len(CCO), dtype = int)
  
    MFO = LIST["MFO"]    
    ROW_INDEX_MFO = np.empty(len(MFO), dtype = int)
    COL_INDEX_MFO = np.empty(len(MFO), dtype = int)
    DATA_MFO      = np.empty(len(MFO), dtype = int)

    # Make list of all GO:terms that occur  (split by ontology)  
    for term in allGOterms[1:]:
        split_term = term.split("\n")
        # Read term info in by line
        for line in split_term:
            # Get the Go Term
            if "id:" in line and "GO:" in line and "alt_id" not in line:
                GO_term = "GO:" + line.split("GO:")[-1].strip()
            # Add to ALT_ID if needed
            if "GO:" in line and "alt_id" in line:
                alt_id = "GO:" + line.split("GO:")[-1].strip()
                alt_id_to_id[alt_id] = GO_term
            # Get the Namespace for this entry
            if   "namespace: biological_process" in line:
                namespace = "BPO"
            elif "namespace: cellular_component" in line:
                namespace = "CCO"   
            elif "namespace: molecular_function" in line:
                namespace = "MFO"
            # Add edges based on the relationships in the entry
            if(("is_a:" in line or "relationship: part_of" in line) and "GO" in line):
                
                if "is_a" in line:
                    relation = "is_a"
                else:
                    relation = "part_of"
                try:
                    # Parent GO Term is listed if in a relationship
                    parent_GO_term = "GO:" + line.split("GO:")[1][:7]
                except IndexError:
                    print(line)
                    exit()
                
                # For every is_A | part_of
                # make matric DAG[i,j] = 1 
                # where i is Term, j is Term that is in relationship 
                try:
                    if   namespace == "BPO" and parent_GO_term in BPO:
                        np.append(ROW_INDEX_BPO, BPO.index(GO_term))
                        np.append(COL_INDEX_BPO, BPO.index(parent_GO_term))
                        np.append(DATA_BPO, 1)
                        
                    elif namespace == "CCO" and parent_GO_term in CCO:
                        np.append(ROW_INDEX_CCO, CCO.index(GO_term))
                        np.append(COL_INDEX_CCO, CCO.index(parent_GO_term))
                        np.append(DATA_CCO, 1)
                        
                    elif namespace == "MFO" and parent_GO_term in MFO:
                        np.append(ROW_INDEX_MFO, MFO.index(GO_term))
                        np.append(COL_INDEX_MFO, MFO.index(parent_GO_term))
                        np.append(DATA_MFO, 1)
                        
                except(ValueError):
                    # Term crosses ontolgies or is invalid in some other way
                    print(GO_term)
    # Make COO sparse matrix                
    BPO_DAG = sc.sparse.coo_matrix((DATA_BPO, (ROW_INDEX_BPO, COL_INDEX_BPO)))                
    CCO_DAG = sc.sparse.coo_matrix((DATA_CCO, (ROW_INDEX_CCO, COL_INDEX_CCO)))               
    MFO_DAG = sc.sparse.coo_matrix((DATA_MFO, (ROW_INDEX_MFO, COL_INDEX_MFO)))               
                    
    # Store subDAGS into dictionary
    DAG["BPO"] = BPO_DAG
    DAG["CCO"] = CCO_DAG
    DAG["MFO"] = MFO_DAG
    print("DAG is complete")
    return DAG, LIST

    
def makeA(LIST_Terms):
    '''
    
    
    LIST is needed to verify lengths
    '''
    
    BPO = []
    CCO = []
    MFO = []
    LIST_Protein = {}
    # Load list to memory
    reader = open("/home/mcgerten/cafa3/CAFA3_test_data/uniprot_exp_20170118_withTAS.txt","r")
    # For every annotation, create an id and add to data
    for line in reader:
        Protein, GOTerm, Aspect = line.split("\t")
        Aspect, Empty =  Aspect.split("\n")
        if   Aspect == "P":
            BPO.append([Protein, GOTerm])
        elif Aspect == "C":
            CCO.append([Protein, GOTerm])
        elif Aspect == "F":
            MFO.append([Protein, GOTerm])
        else:
            print("Not in valid Ontology")
    LIST_Protein["BPO"] = BPO
    LIST_Protein["CCO"] = CCO
    LIST_Protein["MFO"] = MFO     
    
    ################ DONE READING LIST (GAF)
    # Set and print Lengths of Protien/Term Lists
    print   (len(LIST_Protein["BPO"]))
    LP_BPO = len(LIST_Protein["BPO"])
    print   (len(LIST_Protein["CCO"]))
    LP_CCO = len(LIST_Protein["CCO"])
    print   (len(LIST_Protein["MFO"]))
    LP_MFO = len(LIST_Protein["MFO"])
    # Set and print Lengths of Term/Term Lists
    print   (len(LIST_Terms["BPO"]))
    LT_BPO = len(LIST_Terms["BPO"])
    print   (len(LIST_Terms["CCO"]))
    LT_CCO = len(LIST_Terms["CCO"])
    print   (len(LIST_Terms["MFO"]))
    LT_MFO = len(LIST_Terms["MFO"])
    # Intilize arrays
    ROW_INDEX_BPO = np.empty(LP_BPO, dtype = int)
    COL_INDEX_BPO = np.empty(LP_BPO, dtype = int)
    DATA_BPO      = np.empty(LP_BPO, dtype = int)
    
    ROW_INDEX_CCO = np.empty(LP_CCO, dtype = int)
    COL_INDEX_CCO = np.empty(LP_CCO, dtype = int)
    DATA_CCO      = np.empty(LP_CCO, dtype = int)
    
    ROW_INDEX_MFO = np.empty(LP_MFO, dtype = int)
    COL_INDEX_MFO = np.empty(LP_MFO, dtype = int)
    DATA_MFO      = np.empty(LP_MFO, dtype = int)
    
    A = {}
    # {i, j} i is protein, j is GO:term
    '''
    # BPO 
    i = 0
    for item in BPO:
        # Find the index that the term has in the DAG
        try:
            j = LIST_Terms["BPO"].index(item[1])
            np.append(ROW_INDEX_BPO, i)
            np.append(COL_INDEX_BPO, j)
            np.append(DATA_BPO, 1)
        except ValueError:
            print(item[1])
        # Increment for each new protein
        i += 1
    '''
    # CCO    
    i = 0
    for item in CCO:
        # Find the index that the term has in the DAG
        print (item)
        try:
            j = LIST_Terms["CCO"].index(item[1])
            print (j)
            np.append(ROW_INDEX_CCO, i)
            np.append(COL_INDEX_CCO, j)
            np.append(DATA_CCO, 1)
        except ValueError:
            print(item[1])
        # Increment for each new protein
        i += 1
    '''
    # MFO    
    i = 0
    for item in MFO:
        # Find the index that the term has in the DAG
        try:
            j = LIST_Terms["MFO"].index(item[1])
            np.append(ROW_INDEX_MFO, i)
            np.append(COL_INDEX_MFO, j)
            np.append(DATA_MFO, 1)
        except ValueError:
            print(item[1])
        # Increment for each new protein
        i += 1
    '''
        
    # Make A matrices    
    #A_BPO = sc.sparse.coo_matrix((DATA_BPO, (ROW_INDEX_BPO, COL_INDEX_BPO)))                
    A_CCO = sc.sparse.coo_matrix((DATA_CCO, (ROW_INDEX_CCO, COL_INDEX_CCO)))               
    #A_MFO = sc.sparse.coo_matrix((DATA_MFO, (ROW_INDEX_MFO, COL_INDEX_MFO)))               
    # Store in a dictionary    
    #A["BPO"] = A_BPO
    A["CCO"] = A_CCO
    #A["MFO"] = A_MFO
    print("A is complete")
    return A
    
def calculateIC(DAGs, As, LIST):
    IC = {}
    #for ontology in ['BPO','CCO','MFO']:
    for ontology in ['CCO']:
        I = {}
        DAG  = DAGs[ontology]
        A    = As[ontology]
        
        ic = ICcalculator(DAG, A)
        for term in LIST[ontology]:
            prob = ic.index(term)
            I[term] = [prob, -math.log( prob, 2 )]
        # Store in IC dictionary 
        IC[ontology] = I
        cp.dump(IC[ontology], open("ICdata/ia_{}.map".format(ontology), "wb"))
        convertToReadable(IC[ontology], "{}".format(ontology))
    return IC
    
    
''' MATLAB equivelent
  for i = 1 : k
      p        = subDAG(i, :); % parent term(s)
      support  = all(subA(:, p), 2);
      S        = sum(support);
      subia(i) = sum(support & subA(:, i)) / S;
  end
'''
def ICcalculator(DAG, A):
    ic = {}

    #DAG = np.matrix([[0,0,0,0,0],[1,0,0,0,0],[1,0,0,0,0],[0,1,1,0,0],[0,0,1,0,0]])  
    #A = np.matrix([[1,1,1,1,0],[1,0,1,0,1],[1,1,1,0,0],[1,0,1,0,0]])
     
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
        
    for i in range(0, DAG.shape[0]):
        # Make holder matrix
        B = np.empty((A.shape[0],1))
        print(B)
        print(B.shape[0])
        print(B.shape[1])
        parents = DAG[i, :]
        print (" These are the parents of the {} term".format(i+1))
        print (parents)
        print(parents.shape)
        print(parents.shape[0])
        # Throw error becuase [x,] instead of [1,x] notation
        #print(parents.shape[1])
        s = parents.shape[0]
        # I know know who the direct ancestor of a Term is/are
        # For those parents, see how often they occur in A
        # ie. how many proteins are annotated for a term
        # Sum for all terms 
        t = 0
        #print("parents size should be 5")
        for j in range(0, s):
            # If it is a parent
            #print("I printed the loop index")
            #print(i)
            print("This is the value of the {}th item".format(j+1))
            print(parents.item((j)))
            if parents.item((j)) == 1:
                # Count how many times that parent occurs in proteins
                C = A[:, j]                 
                B = B + C
                del C
                #print("Matrix B")                
                #print (B)

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

    
def convertToReadable(IC, ontology):
    '''
    Convert Map to human readable Values to manually check accuracy
    '''
    data  = open("ICdata/IAmap_{}.txt".format(ontology), 'w')
    for term in IC:
        data.write('{}\t {}\t {}\n'.format(term, IC[term][0], IC[term][1]))
    data.close


if __name__ == '__main__':
   main()
   