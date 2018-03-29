# -*- coding: utf-8 -*-

GAF21FIELDS = [
                'DB',
                'DB_Object_ID',
                'DB_Object_Symbol',
                'Qualifier',
                'GO_ID',
                'DB:Reference',
                'Evidence',
                'With',
                'Aspect',
                'DB_Object_Name',
                'Synonym',
                'DB_Object_Type',
                'Taxon_ID',
                'Date',
                'Assigned_By',
                'Annotation_Extension',
                'Gene_Product_Form_ID'
                ]

from Bio.UniProt  import GOA
import collections
import math
import pickle as cp
import gc
import os
import argparse
import sys
import yaml

import ICHelper

  
    
    
def LISTtoDICT(list_path):
    '''
    Converts an imported List into a dictionary
    
    Input: 
    gaf : GAF 2.1 object generator
    
    Ouput: 
    [0] : Dictionary KEY: AnnotationID, VALUE: Annotation 
    '''
    
    counter = 0
    data = dict()
    # Load list to memory
    reader = open(list_path,"r")
    # For every annotation, create an id and add to data
    for line in reader:
        Protein, GOTerm, Aspect = line.split("\t")
        #print Protein
        #print GOTerm
        Aspect, Empty =  Aspect.split("\n")
        #print Aspect
        
        data[counter] = [Protein, GOTerm, Aspect]
        #print data[counter]
        if (counter % 1000 == 0):
            print "I have done {} annotations".format(counter)
        counter += 1
 

    print "I have done {} annotations".format(counter)
    return data
    
   
def ProteinToGO(data):
    '''
    This function creates a dictionary where key is a protein. 
    Each protein refers to a list of GO_TERMS.
    
    Input: 
    data  : Dictionary KEY: AnnotationID, VALUE: 
    
    Output: 
    [0]   : Dictionary KEY: Proteins,     VALUE: GO Terms
    [1]   : Set        all GO terms in the data
    '''
    
     # Load Alternate map
    alt_id_to_id = cp.load( open("ALTERNATE_ID_TO_ID.map", "rb" ) )
    protein_to_go = dict()
    all_GO = []
    
    for annotationID in data:
        annotation = data[annotationID]
        #print annotation
        protein = annotation[0] 
        GO_term = annotation[1]
        aspect  = annotation[2]
        #print protein
        #print GO_term
        #print Aspect
        # If in alternate map, replace
        if GO_term in alt_id_to_id:
            GO_term = alt_id_to_id[GO_term]
        # Add GO to total list
        all_GO.append( GO_term )
        # Check if protein has been added
        if protein not in protein_to_go:
            #print "added protein"
            protein_to_go[protein] = []
        # Check if aspect has been added    
        if [GO_term, aspect] not in protein_to_go[protein]:
                protein_to_go[protein].append( [GO_term, aspect] )
    return protein_to_go, list( set( all_GO ) )


def propagateOntologies(Protein_to_GO):
    '''
    Takes each annotation and constructs the ancestors of that term from their aspect
    
    Input:
    Protein_to_GO : Dictionary KEY: Proteins, VALUE: GO Terms
    
    Output:    
    [0]           : Dictionary KEY: Proteins, VALUE: GO Terms, including ancestors 
    '''
    
    Prot_to_GO_new = dict()
    
    mf_ancestors = cp.load(open("MFO_ANCESTORS.graph","rb"))
    bp_ancestors = cp.load(open("BPO_ANCESTORS.graph","rb"))
    cc_ancestors = cp.load(open("CCO_ANCESTORS.graph","rb"))
    
    for protein in Protein_to_GO:
        ancestors = []
        annotations = Protein_to_GO[protein]
        for annotation in annotations:
            aspect = annotation[1]
            GO_term = annotation[0]
            try:
                if aspect == 'F':
                    ancestors.extend(mf_ancestors[GO_term])
                if aspect == 'P':
                    ancestors.extend(bp_ancestors[GO_term])
                if aspect == 'C':
                    ancestors.extend(cc_ancestors[GO_term])
            except KeyError: # Key doesn't exist
                pass
        ancestors = list( set( ancestors ) )
        Prot_to_GO_new[protein] = ancestors
        
    return Prot_to_GO_new
    

def assignProbabilitiesToOntologyTree(graph, Protein_to_GO, all_GO_Terms, ontology_to_ia, aspect ):
    '''
    Calculates probalities for a given aspect graph
    
    Input:
    graph           : 
    Protein_to_GO   : Dictionary KEY: Proteins,    VALUE: GO Terms
    all_GO_Terms    : Set        all GO terms in the data
    ontology_to_ia  : Dictionary KEY: Node (Term), VALUE: List[Probability, Log]
    aspect          : String     (MFO, BPO, CCO)
    
    Output:
    [0]             : Dictionary KEY: Node (Term), VALUE: List[Probability, Log]
    '''
    for node_num, node in enumerate( graph.nodes() ):
        if( node not in all_GO_Terms ):
            #Node :[Probability, Log]
            ontology_to_ia[node] = [0, 0]
            continue
        if node_num % 100 == 0:
            print( node_num , " proteins processed for ", aspect )
            
        predecessor = graph.successors( node )
        predecessor_with_node = []
        predecessor_with_node.extend( predecessor )
        predecessor_with_node.append( node )
        denominator = findFrequency( predecessor, Protein_to_GO )
        numerator   = findFrequency( predecessor_with_node, Protein_to_GO )
        
        if( denominator == 0.0 ):
            prob = 0.0
        else:
            prob = float(numerator) / float(denominator)
        if (prob != 0.0):
            ontology_to_ia[node] = [prob, -math.log( prob, 2 )]
        else:
            ontology_to_ia[node] = [prob, 0]
    #Needed?
    return ontology_to_ia

def assignProbabilitiesToOntologyGraphs( Protein_to_GO, all_GO_Terms):
    '''
    Calculates probalities for all aspects
    
    Input:
    Protein_to_GO : Dictionary KEY: Proteins,    VALUE: GO Terms
    all_GO_Terms  : Set        all GO terms in the data
    
    Output:
    [0]           : Dictionary KEY: Node (Term), VALUE: List[Probability, Log]
    '''
    mf_g = cp.load( open( "MFO.graph", "rb" ) )
    bp_g = cp.load( open( "BPO.graph", "rb" ) )
    cc_g = cp.load( open( "CCO.graph", "rb" ) )
    
    ontology_to_ia = dict()
    
    assignProbabilitiesToOntologyTree( mf_g, Protein_to_GO, all_GO_Terms, ontology_to_ia, 'MFO' )
    assignProbabilitiesToOntologyTree( bp_g, Protein_to_GO, all_GO_Terms, ontology_to_ia, 'BPO' )
    assignProbabilitiesToOntologyTree( cc_g, Protein_to_GO, all_GO_Terms, ontology_to_ia, 'CCO' )
    
    return ontology_to_ia
    
   
    
def findFrequency(annotations, Protein_to_GO ):
    '''
    Find the frequency of annotations in the protein data
    '''
    
    count = 0.0
    # Check for no annotations
    if annotations == None:
        return 0
        
    for protein in Protein_to_GO:
        if set( annotations ).issubset( set( Protein_to_GO[protein] ) ):
            count += 1.0
    return count


def countProteins(data):
    '''
    Find the total number of proteins
    
    Input:
    data  : Dictionary KEY: AnnotationID, VALUE: Annotation (Use GAFtoDICT to get)
    
    Output: 
    [0]   : Integer    # of proteins
    '''
    
    allproteins = []
    
    for annotation in data:
        allproteins.append(data[annotation]['DB'] + "_" + data[annotation]['DB_Object_ID'])
        
    return len(set(allproteins))



def freqGO_TERM(data):
    '''
    Find the frequency of every GO Term
    
    Input:
    data  : Dictionary KEY: AnnotationID, VALUE: Annotation (Use GAFtoDICT to get)
    
    Output:
    [0]   : Dictionary KEY: GO_term_ID,   VALUE: # of occurances
    '''
    
    go_freq = dict()
    
    for annotation in data:
        if data[annotation]['GO_ID'] in go_freq:
            go_freq[data[annotation]['GO_ID']] += 1
        else:
            go_freq[data[annotation]['GO_ID']] = 1
    return go_freq


def PhillipLordIC(data):
    '''
    Calculate Phillip Lord's Information Content
    
    Input: 
    data           : Dictionary KEY: AnnotationID, VALUE: Annotation (Use GAFtoDICT to get)
    
    Output:
    [0]            : Dictionary KEY: Node (Term), VALUE: IC]
    '''
    
    go_terms = []
    # Create a list of all terms
    for annotationID in data:
        #print(data[annotationID]["GO_ID"])
        go_terms.append(data[annotationID]["GO_ID"])
    # Put it in a collection    
    PL_info = collections.Counter(go_terms)
    
    ic = dict()
    
    for term in PL_info:
        #print(term , PL_info[term], float(PL_info[term]) / len(go_terms)) #DEBUG STATEMENT
        #print(len(go_terms))
        # Calculate PL IC and add to dictionary
        ic[term] = -math.log(float(PL_info[term]) / len(go_terms), 2)
    
    return ic	
    
    
def WyattClarkIC(data):
    '''
    Calculate Wyatt Clark Information Content 
    
    Input:
    data            : Dictionary KEY: AnnotationID, VALUE: Annotation (Use GAFtoDICT to get)
    
    Output:
    [0]             : Dictionary KEY: Node (Term), VALUE: List[Probability, Log]
    '''
    
    # Get Protein to GO term Dictionary
    Protein_to_GO, all_GO_Terms_in_corpus = ProteinToGO(data)
    # Propagate ancestors 
    Protein_to_GO_propagated = propagateOntologies(Protein_to_GO)
    # Calculate IA
    ontology_to_ia = assignProbabilitiesToOntologyGraphs(Protein_to_GO_propagated, all_GO_Terms_in_corpus)
    
    # Save IA Map
    cp.dump(ontology_to_ia, open("ia.map","wb"))
    
    #protInfoAccretion = calculateInformationAccretion( Prot_to_GO_Map_propagated, ontology_to_ia_map )

    return ontology_to_ia
    

def extant_file(x):
    '''
    Description - extant?
    
    Input:
    x   : 
    
    Output:
    [0] :
    '''
    
    if not os.path.isfile(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    else:
        return(open(x,'r'))
        
        
def read_config():
    '''
    Read in the configuration file
    
    Output:
    [0] : String   OBO         file path
    [1] : String   GAF         file path
    '''
    
    parser = argparse.ArgumentParser(description='Precision- Recall assessment for CAFA predictions.', )
    
    parser.add_argument('config_stream', type = extant_file, help = 'Configuration file')            
    args = parser.parse_args()
    # Load config file to dictionary
    try:
        config_dict = yaml.load(args.config_stream)['assessment']
    except yaml.YAMLError as exc:
        print(exc)
        sys.exit()
    # Store config into itermediate variables    
    obo_path = config_dict['obo_path']
    list_path = config_dict['gaf_path']
     
    return(obo_path, list_path)
    

def main():
    '''
    Get a GAF, Convert, Calculate, and output. 
    '''
    
    # Use the assessment configuration file to grab the OBO file.
    obo_path, list_path = read_config()
    print "Read the config"
    ICHelper.setupGraphs(obo_path)
    print "Graphs are done"
    data = LISTtoDICT(list_path)
    print "I have made the dictionary"
    WyattClarkIC(data)
    print("IC values created")
    # WyattClarkIC stores to disk, we are done!
    
    
if __name__ == '__main__':
   main()