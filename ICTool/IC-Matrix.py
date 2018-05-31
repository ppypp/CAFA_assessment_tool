
import math
import pickle as cp
import os
import argparse
import sys
import yaml
import ICHelper
import numpy as np
from scipy.sparse import coo_matrix
    
def createDAGs():
    '''
    
    '''
    # Parse OBO
    obo_file = open("/home/mcgerten/cafa3/CAFA3_test_data/gene_ontology_edit.obo.2016-06-01","r")
    allGOterms = obo_file.read().split("[Term]")
    # Intialize variables
    GO_term   = ""
    namespace = ""
    relation  = ""
    # Alt name dictionary 
    alt_id_to_id = dict()
    # Lists for each namespace
    BPO  = []
    CCO  = []
    MFO  = []
    # 
    LIST = {}
    DAG  = {}
    
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
    # Make Matrices         
    print(len(BPO))
    print(len(CCO))
    print(len(MFO))
    BPO_DAG = np.zeros((len(BPO), len(BPO)))
    CCO_DAG = np.zeros((len(CCO), len(CCO)))
    MFO_DAG = np.zeros((len(MFO), len(MFO)))
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
                        BPO_DAG[BPO.index(GO_term)][BPO.index(parent_GO_term)] = 1
                    elif namespace == "CCO" and parent_GO_term in CCO:
                        CCO_DAG[CCO.index(GO_term)][CCO.index(parent_GO_term)] = 1
                    elif namespace == "MFO" and parent_GO_term in MFO:
                        MFO_DAG[MFO.index(GO_term)][MFO.index(parent_GO_term)] = 1
                except(ValueError):
                    print(GO_term)
    # Store subDAGS into dictionary
    DAG["BPO"] = BPO_DAG
    DAG["CCO"] = CCO_DAG
    DAG["MFO"] = MFO_DAG
    print("DAG is complete")
    return LIST, DAG


def createA(LIST):
    '''
    
    '''
    counter = 0
    BPO = []
    CCO = []
    MFO = []
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
            # Die
            exit
        if (counter % 1000 == 0):
            print ("I have done {} annotations".format(counter))
        counter += 1
 

    print ("I have done {} annotations".format(counter))
    # Now have a List inported
    print(len(BPO))
    print(len(CCO))
    print(len(MFO))
    print(len(LIST["BPO"]))
    print(len(LIST["CCO"]))
    print(len(LIST["MFO"]))
    # Intilize variables
    #BPO_A = np.zeros((len(BPO), len(LIST["BPO"])))
    CCO_A = np.zeros((len(CCO), len(LIST["CCO"])))
    MFO_A = np.zeros((len(MFO), len(LIST["MFO"])))
    A = {}
    # {i, j} i is protein, j is GO:term
    i = 0
    '''
    for item in BPO:
        # Find the index that the term has in the DAG
        j = LIST["BPO"].index(item[1])
        # Increment for each new protein
        BPO_A[i][j] = 1
        i += 1
    '''
    i = 0
    for item in CCO:
        # Find the index that the term has in the DAG
        try:
            j = LIST["CCO"].index(item[1])
        except ValueError:
            print(item[1])
        # Increment for each new protein
        CCO_A[i][j] = 1
        i += 1
        
    i = 0
    for item in MFO:
        # Find the index that the term has in the DAG
        try:
            j = LIST["MFO"].index(item[1])
        except ValueError:
            print(item[1])
        # Increment for each new protein
        MFO_A[i][j] = 1
        i += 1
    
    # Store in a dictionary    
    #A["BPO"] = BPO_A
    A["CCO"] = CCO_A
    A["MFO"] = MFO_A
    print("A is complete")
    return A

def main():
    '''
    Get a GAF, Convert, Calculate, and output. 
    '''
    # Read Config
    
    # Create DAG : m-by-m adjancency matrix
    # DAG(i, j) means i has relationship to j
    LIST, DAGs = createDAGs()
    # Create A   : n-by-m, Ontology Annotation matrix
    # A(i,j) = True, protein i  has annotation j
    As = createA(LIST)
    IC = {}
    BPO_IC = {}
    CCO_IC = {}
    MFO_IC = {}
    # For each ontology
    #for ontology in ['BPO','CCO','MFO']:
    for ontology in ['CCO', 'MFO']:
        DAG  = DAGs[ontology]
        A    = As[ontology]
        # Validate DAG is square
        if DAG.shape[0] != DAG.shape[1]:
            #Die
            raise ValueError() 
        
        ic = ICcalculator(DAG, A)
        # Store in IC dictionary 
        IC[ontology] = ic
    # Print / Write IC   
    #Convert to old format
    #OLD:         Key: Term, Value: [prob, log]  
    for term in LIST["BPO"]:
        prob = ic.index(term)
        BPO_IC[term] = [prob, -math.log( prob, 2 )]
    for term in LIST["CCO"]:
        prob = IC["CCO"][LIST["CCO"].index(term)]
        #CCO_IC[term] = [prob, -math.log( prob, 2 )]
        CCO_IC[term] = [prob, prob]
    for term in LIST["MFO"]:
        prob = IC["MFO"][LIST["MFO"].index(term)]
        MFO_IC[term] = [prob, prob]
            
        
    # Save IA Map
    cp.dump(BPO_IC, open("ICdata/ia_BPO.map","wb"))
    convertToReadable(BPO_IC, "BPO")
        
    cp.dump(CCO_IC, open("ICdata/ia_CCO.map","wb"))
    convertToReadable(CCO_IC, "CCO")
    
    cp.dump(MFO_IC, open("ICdata/ia_MFO.map","wb"))
    convertToReadable(MFO_IC, "MFO") 
    
    print("I'm done")
    
def convertToReadable(IC, ontology):
    '''
    Convert Map to human readable Values to manually check accuracy
    '''
    data  = open("ICdata/IAmap_{}.txt".format(ontology), 'w')
    for term in IC:
        data.write('{}\t {}\t {}\n'.format(term, IC[term][0], IC[term][1]))


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
        B = np.zeros((A.shape[0],1))
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
   
