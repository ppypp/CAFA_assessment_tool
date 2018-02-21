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
from dateutil     import parser
import collections
import math
import pickle as cp
import ICHelper


def GAFtoDICT(gaf):
    '''
    Converts an imported GAF 2.1 file into a dictionary
    
    Input: 
    gaf : GAF 2.1 object
    
    Ouput: 
    [0] : Dictionary KEY: AnnotationID, VALUE: Annotation 
    '''
    
    # Load Alternate map
    alt_id_to_id = cp.load( open("ALTERNATE_ID_TO_ID.map", "rb" ) )
    counter = 1
    data = dict()
    # For every annotation, create an id and add to data
    for annotation in gaf:
        id = "annotation" + str( counter )
        # If in alternate map, replace
        if annotation['GO_ID'] in alt_id_to_id:
            annotation['GO_ID'] = alt_id_to_id[annotation['GO_ID']]
        #print(id) # DEBUG
        annotation['DB:Reference'] = annotation['DB:Reference'][0]
        annotation['Date'] = parser.parse(annotation['Date']).date()
        data[id] = annotation
        counter += 1
        
    return data
    

def ProteinToGO(data):
    '''
    This function creates a dictionary where key is a protein. 
    Each protein refers to a list of GO_TERMS.
    
    Input: 
    data  : Dictionary KEY: AnnotationID, VALUE: Annotation (Use GAFtoDICT to get)
    
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
        proteinID = annotation['DB'] + '_' + annotation['DB_Object_ID']
        GO_term = annotation['GO_ID']
        # If in alternate map, replace
        if GO_term in alt_id_to_id:
            GO_term = alt_id_to_id[GO_term]
        # Add GO to total list
        all_GO.append( GO_term )
        # Check if protein has been added
        if proteinID not in protein_to_go:
            protein_to_go[proteinID] = []
        # Check if aspect has been added    
        if [GO_term, annotation['Aspect']] not in protein_to_go[proteinID]:
                protein_to_go[proteinID].append( [GO_term, annotation['Aspect']] )
            
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
            if aspect == 'F':
                ancestors.extend(mf_ancestors[GO_term])
            if aspect == 'P':
                ancestors.extend(bp_ancestors[GO_term])
            if aspect == 'C':
                ancestors.extend(cc_ancestors[GO_term])
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
    

def InformationAccretionProtein(Protein_to_GO, ontology_to_ia ): #NEVER USED###############################################
    '''
    Calculate the IA for each protein
    
    Input:
    Protein_to_GO   : Dictionary KEY: Proteins,     VALUE: GO Terms
    ontology_to_ia  : Dictionary KEY: Node (Term),  VALUE: List[Probability, Log]
    
    Output:
    [0]             : Dictionary KEY: Proteins,     VALUE: IA value
    '''
    
    # Load alternate map
    alt_id_to_id = cp.load( open("ALTERNATE_ID_TO_ID.map", "rb" ) )
    infoAccr = dict()
    
    for protein in Protein_to_GO:
        annotations = Protein_to_GO[protein]
        ia = 0
        for annotation in annotations:
            if len( annotation ) == 2:
                GO_term = annotation[0]
            else:
                GO_term = annotation
            # If in alternate, replace    
            if GO_term not in ontology_to_ia:
                GO_term = alt_id_to_id[GO_term]
                
            ia += ontology_to_ia[GO_term][1]
        infoAccr[protein] = ia
    return infoAccr


def writeToFile(data, filename):
    '''
    This function will write the content of the data structure 'data' to the output file. 
    
    Input:
    data     : Dictionary KEY: AnnotationID, VALUE: Annotation (Use GAFtoDICT to get)
    filename : String     Name of file to save data to
    Output: 
    [0]      : NONE
    '''
        
    f = open( filename + ".out", "w" )
    f.write("This is a fake header")
    for key in data:
        per_annotation = data[key]
        per_annotation['Qualifier']='|'.join(per_annotation['Qualifier'])
        per_annotation['With']='|'.join(per_annotation['With'])
        per_annotation['Synonym']='|'.join(per_annotation['Synonym'])
        per_annotation['Taxon_ID']='|'.join(per_annotation['Taxon_ID'])
        per_annotation['Date']=''.join(str(per_annotation['Date']).split("-"))
        #print(per_annotation)
        string = ""
        for field in GAF21FIELDS:
            try:
                string += per_annotation[field] + "\t"
            except TypeError:
                print("Exception has occurred in function writeToFile")
                print(per_annotation)
                print(field)
                print(per_annotation[field])
                exit()
        string += '\n'
        f.write( string )
    f.close()
    
    
def printDetailsAboutData(data):
    '''
    Helper Function that outputs info about data set
    '''    
    
    print("Total number of annotations in the provided Database ",len(data))
    prots=[]
    ref=[]
    for attnid in data:
        annotation=data[attnid]
        prots.append(annotation['DB']+"_"+annotation['DB_Object_ID'])
        ref.append(annotation['DB:Reference'])
    print("Total number of unique proteins in the provided Database ",len(set(prots)))
    print("Total number of unique references in the provided Database ",len(set(ref)))


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
    [0]            : 
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
    Calculate Wyatt Clark Information Content NEEDS WORK GONNA SKIP FOR A BIT
    
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
    

    
    

def main():
    '''
    Get a GAF, Convert, Calculate, and output. 
    '''
    ICHelper.setupGraphs()
    gaf = GOA._gaf20iterator(open("goa_human.gaf"))
    data = GAFtoDICT(gaf)
    ic = WyattClarkIC(data)
    print ic
    # I now have the IC values by term in PL_info
    #PL_info = PhillipLordIC(data, 10, .70)
    #f = open( "Output.out", "w" )
    #f.write("This is a fake header")
    #for item in PL_info:
    #    f.write(item)
    # Can use this to weigh terms
    #writeToFile(data, "Ouput")
    
    
if __name__ == '__main__':
   main()