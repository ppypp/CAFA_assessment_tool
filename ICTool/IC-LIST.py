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

import math
import pickle as cp
import os
import argparse
import sys
import yaml

import ICHelper

  
    
    
def LISTtoDICT(list_path):
    '''
    Converts an imported List into a dictionary
    
    Input: 
    gaf : filename of List
    
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
        Aspect, Empty =  Aspect.split("\n")
        
        data[counter] = [Protein, GOTerm, Aspect]
        #print data[counter]
        if (counter % 1000 == 0):
            print ("I have done {} annotations".format(counter))
        counter += 1
 

    print ("I have done {} annotations".format(counter))
    return data
    
   
def ProteinToGO(data):
    '''
    This function creates a dictionary where key is a protein. 
    Each protein refers to a list of GO_TERMS.
    
    Input: 
    data  : Dictionary KEY: AnnotationID, VALUE: {Protein, Go_Term, aspect}
    
    Output: 
    [0]   : Dictionary KEY: Proteins,     VALUE: GO Terms
    [1]   : Set        all GO terms in the data
    '''
    
     # Load Alternate map
    alt_id_to_id  = cp.load( open("ICdata/ALTERNATE_ID_TO_ID.map", "rb" ) )
    # Intialize variables
    protein_to_go = dict()
    all_GO        = []
    # For all annotations
    for annotationID in data:
        annotation = data[annotationID]
        protein = annotation[0] 
        GO_term = annotation[1]
        aspect  = annotation[2]
        # If in alternate map, replace 
        if GO_term in alt_id_to_id:
            GO_term = alt_id_to_id[GO_term]
        # Add GO to total list
        all_GO.append( GO_term )
        # Check if protein has been added
        if protein not in protein_to_go:
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
    
    mf_ancestors = cp.load(open("ICdata/MFO_ANCESTORS.graph","rb"))
    bp_ancestors = cp.load(open("ICdata/BPO_ANCESTORS.graph","rb"))
    cc_ancestors = cp.load(open("ICdata/CCO_ANCESTORS.graph","rb"))
    
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
    

def assignProbabilitiesToOntologyTree(graph, Protein_to_GO, all_GO_Terms, IC, aspect ):
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
    for count, term in enumerate( graph.nodes() ):
        if( term not in all_GO_Terms ):
            #Node :[Probability, Log]
            IC[term] = [0, 0]
            continue
        if count % 100 == 0:
            print( count , " proteins processed for ", aspect )
            
        predecessor = graph.successors( term )
        predecessor_with_term = []
        # Get predecessor with the current term
        predecessor_with_term.extend( predecessor )
        predecessor_with_term.append( term )
        denominator = findFrequency( predecessor          , Protein_to_GO )
        numerator   = findFrequency( predecessor_with_term, Protein_to_GO )
        
        # Catching Div by 0 Error
        if( denominator == 0.0 ):
            prob = None
        else:
            prob = float(numerator) / float(denominator)
        # Catching Log(0) Error    
        if   (prob != 0.0 and prob is not None):
            IC[term] = [prob, -math.log( prob, 2 )]
        elif (prob == 0.0):
            IC[term] = [prob, 0]
        else: # prob is None
            IC[term] = [None, None]
    return IC
    

def assignProbabilities( Protein_to_GO, all_GO_Terms):
    '''
    Calculates probalities for all aspects
    
    Input:
    Protein_to_GO : Dictionary KEY: Proteins,    VALUE: GO Terms
    all_GO_Terms  : Set        all GO terms in the data
    
    Output:
    [0]           : Dictionary KEY: Node (Term), VALUE: List[Probability, Log]
    '''
        
    IC = dict()
    
    mf_g = cp.load( open( "ICdata/MFO.graph", "rb" ) )
    IC = assignProbabilitiesToOntologyTree( mf_g, Protein_to_GO, all_GO_Terms, IC, 'MFO' )
    del mf_g
    
    bp_g = cp.load( open( "ICdata/BPO.graph", "rb" ) )
    IC = assignProbabilitiesToOntologyTree( bp_g, Protein_to_GO, all_GO_Terms, IC, 'BPO' )
    del bp_g
    
    cc_g = cp.load( open( "ICdata/CCO.graph", "rb" ) )
    IC = assignProbabilitiesToOntologyTree( cc_g, Protein_to_GO, all_GO_Terms, IC, 'CCO' )
    del cc_g
        
    return IC
   
    
def findFrequency(terms, Protein_to_GO ):
    '''
    Find the frequency of terms in the protein data
    '''
    
    count = 0.0
    # Check for no annotations
    if terms == None:
        return 0
    # For all proteins    
    for protein in Protein_to_GO:
        if set( terms ).issubset( set( Protein_to_GO[protein] ) ):
            count += 1.0
    return count


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
    IC = assignProbabilities(Protein_to_GO_propagated, all_GO_Terms_in_corpus)
    
    # Save IA Map
    cp.dump(IC, open("ICdata/ia.map","wb"))
    convertToReadable(IC)
 
 
def convertToReadable(ontology_to_ia):
    '''
    Convert Map to human readable Values to manually check accuracy
    '''
    data  = open(".ICdata/IAmap.txt", 'w')
    for term in ontology_to_ia:
        data.write('{}\t {}\t {}\n'.format(term, ontology_to_ia[term][0], ontology_to_ia[term][1]))


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
    obo_path  = config_dict['obo_path']
    list_path = config_dict['gaf_path']
     
    return(obo_path, list_path)
    

def main():
    '''
    Get a GAF, Convert, Calculate, and output. 
    '''
    
    # Use the assessment configuration file to grab the OBO file.
    obo_path, list_path = read_config()
    print ("Read the config")
    ICHelper.setupGraphs(obo_path)
    print ("Graphs are done")
    data = LISTtoDICT(list_path)
    print ("I have made the dictionary")
    WyattClarkIC(data)
    print("IC values created")
    # WyattClarkIC stores to disk, we are done!
    
    
if __name__ == '__main__':
   main()