# -*- coding: utf-8 -*-
"""
IC from scratch

Created on Tue Feb  6 13:48:00 2018

@author: mcgerten
"""

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


# SO I want to import the GAF file 
from Bio.UniProt import GOA
from dateutil import parser
import collections
import math
import os
import sys

def GAFtoDICT(gaf):
    """
    Description Here    
    
    Input: GAF 2.1 object
    Ouput: Dictionary
    """
    
    counter = 1
    data = dict()
    # For every annotation, create an id and add to data
    for annotation in gaf:
        id = "annotation" + str( counter )
        #print(id) # DEBUG
        annotation['DB:Reference'] = annotation['DB:Reference'][0]
        annotation['Date'] = parser.parse(annotation['Date']).date()
        data[id] = annotation
        counter += 1
    return data


def PhillipLordIC(data, crisp, percentile_val):
    '''
    Calculate Phillip Lord's Information Content
    
    Input:
    Output:
    '''
    
    go_terms=[]
    # Create a list of all terms
    for annotationID in data:
        #print(data[annotationID]["GO_ID"])
        go_terms.append(data[annotationID]["GO_ID"])
    # Put it in a collection    
    PL_info=collections.Counter(go_terms)
    
    ic=[]
    
    for term in PL_info:
        print(term,PL_info[term],float(PL_info[term])/len(go_terms)) #DEBUG STATEMENT
        print(len(go_terms))
        PL_info[term]=-math.log(float(PL_info[term])/len(go_terms), 2)
        ic.append(PL_info[term])
        print PL_info[term]
    if crisp==None:
        threshold=((max(ic)-min(ic))*float(percentile_val)/100)+min(ic)
    else:
        threshold=float(crisp)
    # Status Printing
    print("The maximum value of information content is ", max(ic))
    print("The minimum value of information content is ", min(ic))
    print("The chosen threshold is ",threshold)
    
    new_data=dict()
    # Only takes annotations above the threshold
    for annotationID in data:
        annotation = data[annotationID]
        if PL_info[annotation["GO_ID"]] >= threshold:
            new_data[annotationID]=data[annotationID]
            
    return PL_info
    #print(collections.Counter(go_terms)) #DEBUG Statement

def WyattClarkIC(data, recal, crisp, percentile_val, aspects, outputfiles, input_num):
    """
    Calculate Wyatt Clark Information Content NEEDS WORK GONNA SKIP FOR A BIT
    
    Input:
    Output:
    """
    #vprint(outputfiles[0].split("_"))
    ontology_to_ia_map_filename="ontology_to_ia_map_"+"_".join(outputfiles[input_num].split("/")[-1].split("_")[:-2])+".txt"
    #vprint(ontology_to_ia_map_filename)
    #exit()
    Prot_to_GO_Map, all_GO_Terms_in_corpus = createProteinToGOMapping( data )
    Prot_to_GO_Map_propagated = propagateOntologies( Prot_to_GO_Map )
    if(recal==1):
        print("Recalculating Information Accretion for Wyatt Clark Information Content. This may take a long time depending on the size of input")
        ontology_to_ia_map=assignProbabilitiesToOntologyGraphs(Prot_to_GO_Map_propagated,all_GO_Terms_in_corpus,aspects)
        if os.path.isdir("data/temp/")==False:
            os.makedirs("data/temp/")
        cp.dump(ontology_to_ia_map,open("data/temp/"+ontology_to_ia_map_filename,"wb"))
    else:
        print("Skipping recalculation of Information Accretion for Wyatt Clark")
    try:
        ontology_to_ia_map = cp.load( open( "data/temp/"+ontology_to_ia_map_filename, "rb" ) )
    except IOError as e:
        print("File for GO_Term to ia NOT FOUND. Please rerun the program with the argument -recal 1")
        exit()
    #protInfoAccretion = calculateInformationAccretion( Prot_to_GO_Map_propagated, ontology_to_ia_map )
    ia=[]
    for mapping in ontology_to_ia_map:
        if ontology_to_ia_map[mapping][0]!=0:
            ia.append(ontology_to_ia_map[mapping][1])
    #vprint(sorted(ia))
    vprint(len(ia))
    # Doing Some statistical analysis with the distribution of information content
   
    
    if crisp==None:
        threshold=(max(ia)-min(ia))*float(percentile_val)/100+min(ia)
    else:
        threshold=float(crisp)
    #print("Wyatt Clark Threshold",threshold,min(ia),max(ia))
    new_data=dict()
    
            
    for attnid in data:
        annotation=data[attnid]
        #print(ontology_to_ia_map[annotation["GO_ID"]])
        if ontology_to_ia_map[annotation["GO_ID"]][1]>=threshold:
            new_data[attnid]=data[attnid]
    #print(threshold)
    return new_data
    
    
def writeToFile(data, filename):
    """
    This function will write the content of the data structure 'data' to the output file. 
    
    It requires the input file to read the header. Inclusion of the header is mandatory.
    Input:
    Output: The data on disk
    """
    
    
        
    
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
    
    
def main():
    '''
    Get a GAF, Convert, Calculate, and output. 
    '''
    gaf = GOA._gaf20iterator(open("goa_human.gaf"))
    data = GAFtoDICT(gaf)
    # I now have the IC values by term in PL_info
    PL_info = PhillipLordIC(data, 10, .70)
    f = open( "Output.out", "w" )
    f.write("This is a fake header")
    for item in PL_info:
        f.write(item)
    # Can use this to weigh terms
    #writeToFile(data, "Ouput")
    
    
if __name__ == '__main__':
   main()