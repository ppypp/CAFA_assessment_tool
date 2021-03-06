# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
 CAFA module
 
"""

import sys
from collections import defaultdict
import os
from Ontology.IO import OboIO
import helper
import assessment.FMAX  as F
import assessment.WFMAX as W
import assessment.SMIN  as S
import assessment.NSMIN as N
import assessment.AUC   as A

legal_species = [
"all",
"eukarya",
"prokarya",
 'HELPY',
 'ECOLI',
 'RAT',
 'DANRE',
 'SULSO',
 'DROME',
 'PSEPK',
 'STRPN',
 'PSEAE',
 'BACSU',
 'MYCGE',
 'HUMAN',
 'METJA',
 'DICDI',
 'YEAST',
 'SCHPO',
 'ARATH',
 'XENLA',
 'MOUSE',
 'PSESM',
 'SALTY',
 'CANAX',
 'SALCH']

root_terms      = ['GO:0008150', 'GO:0005575', 'GO:0003674']
legal_ontology  = ["bpo", "mfo", "cco", "hpo"] 
legal_types     = ["type1", "type2", "typex"]
legal_subtypes  = ["easy", "hard"]
legal_modes     = ["full", "partial"]

######################################BENCHMARK START#################################
def go_ontology_split(ontology):
    """
    Split an GO obo file into three ontologies
    
    Input:
    ontology : Ontology File generate by OBOReader
    
    Output:
    [0]      : Set   MFO Terms
    [1]      : Set   BPO Terms
    [2]      : Set   CCO Terms
    """
    
    mfo_terms = set({})
    bpo_terms = set({})
    cco_terms = set({})
    
    for term in ontology.get_ids(): # loop over node IDs and alt_id's
        if ontology.namespace[term]   == "biological_process": # P
            bpo_terms.add(term)
        elif ontology.namespace[term] == "cellular_component": # C
            cco_terms.add(term)
        elif ontology.namespace[term] == "molecular_function": # F
            mfo_terms.add(term)
        else:
            raise(ValueError,"{} has no namespace".format(term))
    return (mfo_terms, bpo_terms, cco_terms)
    
    
def go_ontology_ancestors_split_write(obo_path):
    """
    Takes an OBO file and writes 3 files with ancestors corresponding to the namespace
    
    Input:
    obo_path : filepath of the obo file
    
    Output:
    [0]      : List [Length of BPO, Length of CCO, Length of MFO]
    """
    
    go = OboIO.OboReader(open(obo_path)).read()
    # Split the OBO into the three namespaces 
    mfo_terms, bpo_terms, cco_terms = go_ontology_split(go)
    # BPO
    obo_bpo_out = open("{}_ancestors_bpo.txt".format(os.path.splitext(obo_path)[0]), "w")
    for term in bpo_terms:
        ancestors = go.get_ancestors(term)
        obo_bpo_out.write("%s\t%s\n" % (term, ",".join(ancestors)))
    obo_bpo_out.close()
    # CCO
    obo_cco_out = open("{}_ancestors_cco.txt".format(os.path.splitext(obo_path)[0]), "w")
    for term in cco_terms:
        ancestors = go.get_ancestors(term)
        obo_cco_out.write("%s\t%s\n" % (term, ",".join(ancestors)))
    obo_cco_out.close()    
    # MFO
    obo_mfo_out = open("{}_ancestors_mfo.txt".format(os.path.splitext(obo_path)[0]), "w") 
    for term in mfo_terms:
        ancestors = go.get_ancestors(term)
        obo_mfo_out.write("%s\t%s\n" % (term, ",".join(ancestors)))    
    obo_mfo_out.close()
    
    return([len(bpo_terms), len(cco_terms), len(mfo_terms)])
    
    
def read_benchmark(ontology, species, types, benchmark_directory, obopath):
    '''
    Read Benchmark.
    
    Input:
    ontology            : String        {bpo, cco, mfo, hpo}
    species             : String        List above of accepted options
    types               : String        {full, partial}
    fullbenchmarkfolder : String        Directory
    obopath             : String        File 
    
    Output:
    [0]                 : Benchmark
    [1]                 : obocountDict 
    '''
    
    # Error Checking
    if ontology not in legal_ontology:
        sys.stderr.write("Namespace not accepted, choose from 'bpo', 'cco', 'mfo' and 'hpo'\n")
    elif (species not in legal_species) and (species not in legal_subtypes):
        sys.stderr.write('Species not accepted')
    elif types not in legal_types:
        sys.stderr.write('Type not accepted, choose from "type1","type2" and "typex"\n')
    else:
        # Create name for finding benchmark file
        matchname = ontology + '_' + species + '_' + types + '.txt'
    # Generate ancestor files
    # Store size of each ontology
    obocounts = go_ontology_ancestors_split_write(obopath)
    obocountDict = {'bpo':obocounts[0],'cco':obocounts[1],'mfo':obocounts[2]}
    # Ontology-specific calculations
    if   ontology == 'bpo':
        full_benchmark_path = benchmark_directory + '/groundtruth/' + 'leafonly_BPO.txt'
        ancestor_path = os.path.splitext(obopath)[0] + "_ancestors_bpo.txt"
    elif ontology == 'cco':
        full_benchmark_path = benchmark_directory + '/groundtruth/' + 'leafonly_CCO.txt'
        ancestor_path = os.path.splitext(obopath)[0] + "_ancestors_cco.txt"
    elif ontology == 'mfo':
        full_benchmark_path = benchmark_directory + '/groundtruth/' + 'leafonly_MFO.txt'
        ancestor_path = os.path.splitext(obopath)[0] + "_ancestors_mfo.txt"
        
    benchmarkListPath = benchmark_directory + '/lists/' + matchname
    # If the List file exists and has data
    if os.path.isfile(benchmarkListPath) and os.path.getsize(benchmarkListPath)>0:
        handle    = open(benchmark_directory + '/lists/' + matchname, 'r')
        proteins  = set()
        for line in handle:
            proteins.add(line.strip())
        handle.close()
        tempfilename = 'temp_%s_%s_%s.txt' % (ontology, species, types)
        tempfile = open(benchmark_directory + '/' + tempfilename , 'w')
        for line in open(full_benchmark_path,'r'):
            protein = line.split('\t')[0]
            if protein in proteins:
                tempfile.write(line)
        tempfile.close()
        bench = benchmark(ancestor_path, tempfile.name)
        bench.propagate()
        #os.remove(tempfile.name)
    else:
        print('Benchmark set is empty.\n')
        bench = None
    return(bench, obocountDict)   
    
########################END OF FUNCTIONS OUTSIDE OF A CLASS ##############################
    
class benchmark:
    '''
    Defines a benchmark for comparing the predictions to
    '''
    
    
    def __init__(self, ancestor_path, benchmark_path):
        '''
        Initialize the benchmark.
        
        Input: 
        benchmark_path : String    ontology specific file location
        ancestor_path  : String    ontology specific file location
        '''
        
        # Key: protein
        # Value: set of benchmark leaf terms
        self.ancestors = defaultdict(set)
        # Read GO ancestors file generated with go_ontology_ancestors_split_write()
        # File format: 
        # go_term <tab> ancestor_1,ancestor_2,..,ancestor_n
        with open(ancestor_path) as ancestors_input:
            for inline in ancestors_input:
                inrec = inline.strip().split('\t')
                term = inrec[0]
                if len(inrec) == 1:
                    self.ancestors[term] = set({})
                else:
                    term_ancestors = inrec[1]
                    self.ancestors[term] = set(term_ancestors.split(','))
                
        self.true_base_terms = defaultdict(set)
        with open(benchmark_path) as benchmark_input:
            for inline in benchmark_input:
                protein, term = inline.strip().split('\t')
                self.true_base_terms[protein].add(term)
                
                
    def propagate(self):
        '''
        Progate Benchmark terms.
        '''
        
        self.true_terms = defaultdict(set)
        # Key: protein
        # Value: set of benchmark propagated terms
        
        for protein in self.true_base_terms:
            for term in self.true_base_terms[protein]:
                try:
                    
                    ancestors = self.ancestors[term].difference(root_terms)    
                    
                except KeyError:
                    sys.stderr.write("%s not found \n" % term) 
                # Add the term    
                self.true_terms[protein].add(term)
                # Set Union 
                self.true_terms[protein] |= ancestors
        

#####################################BENCHMARK END ############################

###############################################################################
#################################### START OF INFO ############################
###############################################################################

class Info:
    ''' Holds all the stuff we need, allows for just importing once '''
    
    def writeBenchwithIC(self, ontology, Type):
        '''
        Write to file the benchmark IC values for manual validation of methods
        '''
        
        data  = open(self.path + '/benchmarkIC_{}_{}.txt'.format(ontology, Type), 'w')
        
        for protein in self.true_terms:
            for term in self.true_terms[protein]:
                try:
                    data.write('{}\t {}\t {}\n'.format( protein, term, self.ic[term][1])) 
                except KeyError:
                    data.write('{}\t {}\t {}\n'.format( protein, term, "Key Error"))
        
        data.close() 
        
        #data  = open(self.path + '/test.txt', 'w')
        
        #data.write('{}\n'.format(self.ancestors)) 
                       
        #data.close()
        
        
    def __init__(self, benchmark, os_prediction_path, obocounts, ic, results_path):
        '''
        Initalize.
        
        Input:
        benchmark          : instance of the benchmark class
        os_prediction_path : !ontology-specific! prediction file 
                              separated in GOPred (without headers)
        obocounts          : total number of terms in the ontology
        ic                 : Dictionary KEY: Term VALUE: IC value
        '''
        
        # Initialize variables
        self.path                             = results_path
        self.exist                            = True
        self.ancestors                        = benchmark.ancestors
        self.true_terms                       = benchmark.true_terms
        self.obocount                         = obocounts
        # Set of all obsolete terms found
        self.obsolete                         = set()
        # count_above_threshold is the number of proteins 
        # with at least one term above threshold
        self.count_above_threshold            = defaultdict()
        # count_predictions_in_benchmark is the number of predicted proteins 
        # in this file that are in the benchmark file (for coverage)
        self.count_predictions_in_benchmark   = 0
        self.count_true_terms                 = len(self.true_terms)
        # self.data is the dictionary for all predicted terms
        self.data                             = defaultdict(list)
        # self.predicted_bench is the dictionary 
        # for all predicted terms that are benchmarks
        self.predicted_bench                  = defaultdict(defaultdict)
        # key:protein
        # value: list of dictionaries
        # key: GO term
        # value: tuple(confidence, Boolean) Boolean is whether in self.true_terms     
        self.ic                               = ic   
        # key: term
        # value: IC
        
        # Read in prediction file        
        if os.path.getsize(os_prediction_path) > 0:
            for inrec in open(os_prediction_path,'r'):
                fields = [i.strip() for i in inrec.split()]
                self.data[fields[0]].append({'term': fields[1], 'confidence': float(fields[2])})

            
            # Propagated prediction
            for protein in self.data:
                # Take only protein that are in benchmark
                if self.true_terms[protein]:
                    '''
                    benchmark.true_terms[protein] not an empty set
                    The protein is in the benchmark file
                    i.e. gained experimental annota
                    '''
                    self.count_predictions_in_benchmark += 1
                    for tc in self.data[protein]:
                        try:
                            ancterms = self.ancestors[tc['term']].difference(root_terms)
                        except KeyError:
                            # Add unknown term to the obsolete set
                            self.obsolete.add(tc['term'])
                            continue
                        # For each term
                        if tc['term'] in self.predicted_bench[protein]:
                            # Term already exists, update confidence
                            self.update_confidence(protein,tc)
                        else:
                            # Add term to self.predicted_bench
                            # Add confidence and compare with self.true_terms
                            # Regardless of comparision, propagate
                            if tc['term'] not in root_terms:
                                self.predicted_bench[protein][tc['term']] = self.compare(protein, tc)
                                for ancterm in ancterms:
                                    # Make a new TC with ancestors
                                    newtc = {'term':ancterm,'confidence':tc['confidence']}
                                    if ancterm in self.predicted_bench[protein]:
                                        # Term already exists, update confidence
                                        self.update_confidence(protein, newtc)
                                    else:
                                        # Add term to self.predicted_bench
                                        self.predicted_bench[protein][ancterm] = self.compare(protein, newtc)
                                        
            if self.count_predictions_in_benchmark == 0:
                self.exist = False
                print("No protein in this predicted set became a benchmark\n")
        else:
            # File does not exist
            self.exist = False
            print('No prediction made in this ontology.\n')
            
            
    def coverage(self):
        ''' 
        Determine the coverage. 
        
        Output:
        [0] : Float       coverage of the prediction
        '''
        
        return float(self.count_predictions_in_benchmark)/self.count_true_terms  
       
     
    def update_confidence(self, protein, tc):
        '''
        Update Confidence for given protein and propagate.
        
        This function compares the confidence value in tc to the confidence in self.predicted
        If tc is larger, than it overwrites the confidence in self.predicted
        
        Input:
        protein : chosen protein 
        tc      : Dictionary KEY: Confidence VALUE: GO_Term
        '''
        
        # Defined for readablity
        confidence = tc['confidence'] 
        term = tc['term']
        
        if confidence > self.predicted_bench[protein][term][0]:
            # Update the confidence
            self.predicted_bench[protein][term][0] = confidence
            # Propagate changes if necessary
            for ancterm in self.ancestors[term].difference(root_terms):
                if confidence > self.predicted_bench[protein][ancterm][0]:
                    # Update the confidence
                    self.predicted_bench[protein][ancterm][0] = confidence
                    
                    
    def compare(self, protein, tc):
        '''
        Check if tc['term'] is a True term.
        
        This function compares if tc['term'] is in self.true_terms
        
        Input:
        protein  : Protein for comparision
        tc       : Dictionary KEY: Confidence VALUE: GO_Term
        
        Output:
        [0]      : List [Confidence value, Boolean]
        '''
        
        if tc['term'] in self.true_terms[protein]:
            return [tc['confidence'], True]
        else:
            return [tc['confidence'], False]
            
        
    def getObsolete(self):
        ''' 
        Get obsolete terms used by the prediction team. 
        
        Output:
        [0] : Set      GOTerms that are obsolute
        '''
        
        return(self.obsolete)
        
        
    def check(self, tool, ontology, Type, mode):
        '''
        Calls subroutines for specified tool
        
        Inputs:
        tool     : {FMAX, WFMAX, SMIN, NSMIN, AUC}
        ontology : {BPO, CCO, MFO} #HPO?
        Type     : {NK, LK}{Type1, Type2}
        mode     : {full, partial}
        
        Output:
        [0]   : Float         Output Value
        [1]   : List[Float]   First Values
        [2]   : List[Float]   Second values
        [3]   : Float         Threshold
        '''
        
        info = self
        
        if  (tool == "FMAX"):
            return F.output(info, ontology, Type, mode)
        elif(tool == "WFMAX"):
            return W.output(info, ontology, Type, mode)
        elif(tool == "SMIN"):
            return S.output(info, ontology, Type, mode)
        elif(tool == "NSMIN"):
            return N.output(info, ontology, Type, mode)
        elif(tool == "AUC"):
            return A.output(info, ontology, Type, mode)
        else:
            raise ValueError("Not a valid tool")