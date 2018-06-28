import math
import pickle as cp
import os
import argparse
import sys
import yaml
import ICHelper
import numpy as np
import scipy as sc
import time

def main():
    '''
    Comment/ uncomment to complete steps
    Much faster to make DAG, then A, then IC
    than to try and complete all at once
    '''
    # Read OBO
    # Build DAG
    '''
    DAG, LIST = makeDAG()
    cp.dump(DAG,  open("ICdata/DAG.map",  "wb"))
    cp.dump(LIST, open("ICdata/LIST.map", "wb"))
    '''
    #Read GAF-equivelent 
    # Build A
    
    DAG  = cp.load(open("ICdata/DAG.map", "rb"))
    LIST = cp.load(open("ICdata/LIST.map","rb"))
    print("DAG / LIST loaded")
    '''
    A, LIST_Protein = makeA(LIST)
    cp.dump(A,            open("ICdata/A.map",            "wb"))
    cp.dump(LIST_Protein, open("ICdata/LIST_Protein.map", "wb"))
    '''
    A            = cp.load(open("ICdata/A.map",           "rb"))
    LIST_Protein = cp.load(open("ICdata/LIST_Protein.map","rb"))
    print(" A is loaded")
    
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
    print(len(BPO))
    ROW_INDEX_BPO = []
    COL_INDEX_BPO = []
    DATA_BPO      = []
    
    CCO = LIST["CCO"]
    print(len(CCO))
    ROW_INDEX_CCO = []
    COL_INDEX_CCO = []
    DATA_CCO      = []
  
    MFO = LIST["MFO"]   
    print(len(MFO))
    ROW_INDEX_MFO = []
    COL_INDEX_MFO = []
    DATA_MFO      = []
    
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
                        r = BPO.index(GO_term)
                        c = BPO.index(parent_GO_term)
                        ROW_INDEX_BPO.append(r)
                        COL_INDEX_BPO.append(c)
                        DATA_BPO.append(1)
                        
                    elif namespace == "CCO" and parent_GO_term in CCO:
                        r = CCO.index(GO_term)
                        c = CCO.index(parent_GO_term)
                        ROW_INDEX_CCO.append(r)
                        COL_INDEX_CCO.append(c)
                        DATA_CCO.append(1)
                    elif namespace == "MFO" and parent_GO_term in MFO:
                        r = MFO.index(GO_term)
                        c = MFO.index(parent_GO_term)
                        ROW_INDEX_MFO.append(r)
                        COL_INDEX_MFO.append(c)
                        DATA_MFO.append(1)
                        
                        #print(len(ROW_INDEX_MFO))
                        
                except(ValueError):
                    # Term crosses ontolgies or is invalid in some other way
                    #print(GO_term)
                    pass
    print("Done with relationships")
    # Make COO sparse matrix                
    BPO_DAG = sc.sparse.coo_matrix((DATA_BPO, (ROW_INDEX_BPO, COL_INDEX_BPO)))                
    CCO_DAG = sc.sparse.coo_matrix((DATA_CCO, (ROW_INDEX_CCO, COL_INDEX_CCO)))               
    MFO_DAG = sc.sparse.coo_matrix((DATA_MFO, (ROW_INDEX_MFO, COL_INDEX_MFO)))               
                    
    # Store subDAGS into dictionary
    DAG["BPO"] = BPO_DAG
    DAG["CCO"] = CCO_DAG
    DAG["MFO"] = MFO_DAG
    print("DAG is complete")
    print(BPO_DAG.get_shape)
    print(CCO_DAG.get_shape)
    print(MFO_DAG.get_shape)
    
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
    print("DONE reading in proteins to list")
    # Intilize arrays
    ROW_INDEX_BPO = []
    COL_INDEX_BPO = []
    DATA_BPO      = []
    
    ROW_INDEX_CCO = []
    COL_INDEX_CCO = []
    DATA_CCO      = []
  
    ROW_INDEX_MFO = []
    COL_INDEX_MFO = []
    DATA_MFO      = [] 
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
    
    A = {}
    # {i, j} i is protein, j is GO:term
    
    # BPO
    i = 0
    for item in BPO:
        # Find the index that the term has in the DAG
        try:
            r = i
            c = LIST_Terms["BPO"].index(item[1])
            ROW_INDEX_BPO.append(r)
            COL_INDEX_BPO.append(c)
            DATA_BPO.append(1)
        except ValueError:   
            print(item[1])
        i+=1
        if (i%10000 == 0):
            print("I have done {} items in BPO".format(i))
    
    # CCO    
    i = 0
    for item in LIST_Protein["CCO"]:
        # Find the index that the term has in the DAG
        #print (item)
        try:
            
            r = i
            c = LIST_Terms["CCO"].index(item[1])
            #print(r)
            #print(c)
            
            ROW_INDEX_CCO.append(r)
            COL_INDEX_CCO.append(c)
            DATA_CCO.append(1)
        except ValueError:
            print(item[1])
        i+=1
        if (i%10000 == 0):
            print("I have done {} items in CCO".format(i))
    
    # MFO
    i = 0
    for item in LIST_Protein["MFO"]:
        # Find the index that the term has in the DAG
        #print (item)
        try:
            r = i
            c = LIST_Terms["MFO"].index(item[1])
            ROW_INDEX_MFO.append(r)
            COL_INDEX_MFO.append(c)
            DATA_MFO.append(1)
        except ValueError:
            print(item[1])
        i+=1
        if (i%10000 == 0):
            print("I have done {} items in MFO".format(i))
    
        
    # Make A matrices    
    A_BPO = sc.sparse.coo_matrix((DATA_BPO, (ROW_INDEX_BPO, COL_INDEX_BPO)))                
    A_CCO = sc.sparse.coo_matrix((DATA_CCO, (ROW_INDEX_CCO, COL_INDEX_CCO)))               
    A_MFO = sc.sparse.coo_matrix((DATA_MFO, (ROW_INDEX_MFO, COL_INDEX_MFO)))               
    # Store in a dictionary    
    A["BPO"] = A_BPO
    A["CCO"] = A_CCO
    A["MFO"] = A_MFO
    print(A_BPO.get_shape)
    print(A_CCO.get_shape)
    print(A_MFO.get_shape)
    print("A is complete")
    return A, LIST_Protein
    
    
def calculateIC(DAGs, As, LISTs):
    '''
    
    '''
    IC = {}
    for ontology in ['BPO','CCO','MFO']:
    #for ontology in ['BPO']:
        print(ontology)
        I = {}
        DAG  = DAGs[ontology]
        A    = As[ontology]
        LIST = LISTs[ontology]
        #ic = ICcalculator(DAG, A)
        #cp.dump(ic, open("ICdata/ic_temp_{}.map".format(ontology), "wb"))
        print("TIME TO CONVERT TO OLD STYLE")
        ic = cp.load(open("ICdata/ic_temp_{}.map".format(ontology),"rb"))
        for term in LIST:
            try:
                # Get probality
                prob = ic[LIST.index(term)]
                '''
                #Convert to base E           
                print(prob)
                N = math.log(prob)
                D = math.log(2)
                print(N)
                print(D)
                prob  = N /D
                print(prob)
                '''
                I[term] = [prob, prob]
            except:
                print("ERROR")
                I[term] = [0, 0]
        # Store in IC dictionary 
        IC[ontology] = I
        # Store dictionary where assessment will use it
        cp.dump(IC[ontology], open("ICdata/ia_{}.map".format(ontology), "wb"))
        # Save in human readable format
        convertToReadable(IC[ontology], "{}".format(ontology))
        print("Done with {}".format(ontology))
    return IC
    
    
''' MATLAB equivelent
  for i = 1 : k
M1      p        = subDAG(i, :); % parent term(s)
M2      support  = all(subA(:, p), 2);
M3      S        = sum(support);
M4      subia(i) = sum(support & subA(:, i)) / S;
  end
'''


def ICcalculator(DAG, A):
    ic = {}
    start_time = time.time()
    #Convert to CSR/CSC
    DAG = DAG.tocsr()
    A   = A.tocsc()    
    #paralellalize here? Create a pool?) 
    print("Converted to CSR/CSC")
    print(time.time() - start_time)
    # For all terms in ontology
    for i in range(0, DAG.shape[0]): # 30000 Times at most 
        # Make holder matrix
        B = sc.sparse.csc_matrix((A.shape[0], 1))
        D = sc.sparse.dok_matrix((A.shape[0], 1))
        # Get parents for term i (csr format)
        parents = DAG.getrow(i)
        s = parents.shape[1]
        # Convert to DOK for lookup
        parents = parents.todok()
        # For all Parent terms
        for j in range(0, s): # Takes 1/10th sec on a 4000 / 120000 set)
            p = parents.get((0, j))
            # If is a parent, p = 1. Use != for sparse matrix speed
            if p != 0:
                # Get the column in A for parent term
                C = A.getcol(j)
                # Adding two CSC matrices (fast)
                # Add parent values to B
                B = B + C
                del C
            else:
                pass 
        # Find the largest value in B
        m = B.max() 
        # Convert to DOK to enable fast O(1) lookup
        B = B.todok()
        # For all values in B, create D where a max value in B is 1 in D
        for k in range(0, B.shape[0]):
            # For the sparcity warning
            if B.get((k,0)) != m:
                # Do nothing to go fast               
                pass
            else:
                # Update new matrix with only values that matter
                D[k,0] = 1    
        # Convert to DOK        
        a = A.getcol(i).todok()
        b = D
        # Piece-wise multiplication, will zero out terms not in both
        c = a.multiply(b)
        # Get Num/Dom
        numerator   = np.sum(c)
        denominator = np.sum(b)
        if numerator > denominator:
            print("Weird")
            ic[i] = 0
        else:
            # Calculate IC
            #print("IC Value")
            if denominator != 0 and numerator != 0:
                ic[i] = numerator / denominator
                try:
                    
                    ic[i] = -math.log(ic[i], 2)
                    #print(ic[i])
                except ValueError: 
                    print("Log Error")
                    ic[i] = 0
            else:
                #print ("0")
                ic[i] = 0
        if (i%100 == 0):
            print("I have done {} items".format(i))
            print("TIME TOTAL ")
            print(time.time() - start_time)
            
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
   