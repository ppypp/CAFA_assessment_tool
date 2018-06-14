from collections import defaultdict
import numpy
import pickle as cp
import sys
import os
import helper
from Ontology.IO import OboIO

legal_species = [
 "all", "eukarya", "prokarya",
 'HELPY', 'ECOLI', 'RAT',
 'DANRE', 'SULSO', 'DROME',
 'PSEPK', 'STRPN', 'PSEAE',
 'BACSU', 'MYCGE', 'HUMAN',
 'METJA', 'DICDI', 'YEAST',
 'SCHPO', 'ARATH', 'XENLA',
 'MOUSE', 'PSESM', 'SALTY',
 'CANAX', 'SALCH']

root_terms      = ['GO:0008150', 'GO:0005575', 'GO:0003674']
legal_ontology  = ["BPO", "MFO", "CCO", "HPO"] 
legal_types     = ["type1", "type2", "typex"]
legal_subtypes  = ["easy", "hard"]
legal_modes     = ["full", "partial"]




class Info:
    def __init__(self, ontology, Type):
        # KEY: Protein              Value: GoTerms
        self.terms                  = defaultdict(set)
        # KEY: Protein              Value: GoTerms
        self.truth                  = defaultdict(set)
        # KEY: Protein              Value: GoTerms
        # KEY: GoTerm               Value: Confidence
        self.prediction             = defaultdict(defaultdict)
        # KEY: Protein              Value: MAX Confidence
        self.maxConfidence          = defaultdict()
        # KEY: GoTerm               Value: IC Value
        self.IC                     = defaultdict()
        # KEY: Protein              Value: Thresholds  
        # KEY: Threshold            Value: {FP, TN, TP, FN, POS, NEG, TRUE, FALSE, TOT}
        self.data                   = defaultdict(defaultdict)
        # KEY: Protein              Value: Thresholds   
        # KEY: Threshold            Value: {FP, TN, TP, FN, POS, NEG, TRUE, FALSE}
        self.data_unweighted        = defaultdict(defaultdict)
        # KEY: threshold            Value: Count of Proteins 
        self.ProteinInPrediction    = defaultdict()
        self.ProteinInBenchmark     = 0
        # FOR AUTHOR / FILE INFO
        self.author                 = ""
        self.model                  = ""
        self.keywords               = ""
        self.taxon                  = ""
    
        # FOR Choosing metric / data
        self.ontology               = ontology
        self.Type                   = Type
        self.mode                   = ""
        
        
        # FOR Benchmark inporting
        self.Benchmark              = None
        self.OBOCount               = 0
        # KEY: GoTerm                  Value: GoTerm
        self.TermAncestors          = defaultdict()
        # KEY: Protein                 Value: GoTerm
        self.ProteinTrueTerms       = defaultdict()
        self.TrueTermCount          = 0
        
        # For filepaths
        self.ResultPath             = ""
        self.PredictionPath         = ""
        self.OBOPath                = ""
        self.Path                   = ""
        
    def setMode(self, mode):
        self.mode = mode
        
        
    def setTeamInfo(self, author, model, keywords, taxon):
            '''
            
            '''
            
            vprint("Saving Team Info", 9)
            self.author   = author
            vprint("Author: {}".format(self.author), 9)
            self.model    = model
            vprint("Model: {}".format(self.model), 9)
            self.keywords = keywords
            vprint("Keywords: {}".format(self.keywords), 9)
            self.taxon    = taxon
            vprint("Taxon: {}".format(self.taxon), 9)
                    
                    
    def setBenchmark(self, OBOCounts, benchmark_directory):
        '''
        
        '''
        
        Benchmark = read_benchmark(self.ontology, helper.taxon_name_converter(self.taxon), self.Type, benchmark_directory, self.OBOPath)
        self.Benchmark = Benchmark
        vprint("Benchmark Done", 10)
        self.OBOCount = OBOCounts[self.ontology]
        vprint("OBO Count for {} : {}".format(self.ontology, self.OBOCount), 10)
        self.TermAncestors = Benchmark.ancestors
        self.ProteinTrueTerms = Benchmark.ProteinTrueTerms

        
    def setIC(self, ic_path):
        '''
        
        '''
        
        IC_old = cp.load(open("{}ia_{}.map".format(ic_path, self.ontology), "rb"))
        IC = {}
        for term in IC_old:
            vprint("Term : {}".format(term), 10 )
            vprint("IC   : {}".format(IC_old[term][1]), 10)
            IC[term] = IC_old[term][1]
        self.IC = IC    
        vprint(self.IC, 15)
        
        
    def setPaths(self, results_path, prediction_path, obo_path, path):  
        '''
        
        '''
        
        self.ResultPath = results_path
        vprint("Results Path: {}".format(self.ResultPath), 9)
        # Make directory 
        helper.mkdir_results(self.ResultPath)
        vprint("Result directory structure complete", 9)
        self.PredictionPath = prediction_path
        vprint("Prediction Path: {}".format(self.PredictionPath), 9)
        self.OBOPath = obo_path
        vprint("OBO Path: {}".format(self.OBOPath), 9)
        self.Path = path
        vprint("Path: {}".format(self.Path), 9)
        
    def setPrediction(self):
        self.ProteinInBenchmark = 0
        PredictionProteinTemp = defaultdict(list)     
        if os.path.getsize(self.Path) > 0:
            # Read in prediction file   
            for inrec in open(self.Path,'r'):
                fields     = [i.strip() for i in inrec.split()]
                protein    = fields[0]
                vprint("Protein: {}".format(protein), 15)
                term       = fields[1]
                vprint("Term: {}".format(term), 15)
                confidence = float(fields[2])
                vprint("Confidence: {}".format(confidence), 15)
                # data[Protein] = {Term, Confidence}
                PredictionProteinTemp[protein].append(
                        {'term': term, 'confidence': confidence})
            # Propagated prediction
            for protein in PredictionProteinTemp:
                vprint("Protein: {}".format(protein), 15)
                # Take only Protein that has terms in the benchmark
                if self.ProteinTrueTerms[protein]:
                    self.ProteinInBenchmark += 1
                    vprint(PredictionProteinTemp[protein], 15)
                    for value in PredictionProteinTemp[protein]:
                        term       = value['term']
                        confidence = value['confidence']
                        try:
                            ancterms = self.TermAncestors[term].difference(root_terms)
                            vprint(ancterms,15)
                        except KeyError:
                            # Add unknown term to the obsolete set
                            self.obsolete.add(term)
                            continue
                        # For each term
                        if term in self.prediction[protein]:
                            # Term already exists, update confidence
                            self.update_confidence(protein, term, confidence)
                        else:
                            # Add term to self.predicted_bench
                            if term not in root_terms:
                                self.prediction[protein][term] = confidence
                                # Propagate for ancestors
                                for ancterm in ancterms:
                                    if ancterm in self.prediction[protein]:
                                        # Term already exists, update confidence
                                        self.update_confidence(protein, term, confidence)
                                    else:
                                        # Add term to self.predicted_bench
                                        self.prediction[protein][ancterm] = confidence            
            if self.ProteinInBenchmark == 0:
                vprint("No protein in this predicted set became a benchmark\n", 1)
            else:
                vprint("Protein in Benchmark: {}".format(self.ProteinInBenchmark), 10)
        else:
            vprint('No prediction made in this ontology.\n', 1)
           
            
    def coverage(self):
        ''' 
        Determine the coverage. 
        
        Output:
        [0] : Float       coverage of the prediction
        '''
        
        return float(self.ProteinInPrediction[0.00])/self.ProteinInBenchmark  
        
        
    def update_confidence(self, protein, term, confidence):
        '''
        Update Confidence for given protein and propagate.
        
        This function compares the confidence value in tc to the confidence in self.predicted
        If tc is larger, than it overwrites the confidence in self.predicted
        
        Input:
        protein : choosen protein 
        '''
        
        if confidence > self.prediction[protein][term]:
            # Update the confidence
            self.prediction[protein][term] = confidence
            # Propagate changes if necessary
            for ancterm in self.TermAncestors[term].difference(root_terms):
                if confidence > self.prediction[protein][ancterm]:
                    # Update the confidence
                    self.prediction[protein][ancterm] = confidence
                            
            
        
    def getObsolete(self):
        ''' 
        Get obsolete terms used by the prediction team. 
        
        Output:
        [0] : Set      GOTerms that are obsolute
        '''
        
        return(self.obsolete)   
    
    
    
    
    def setInfo(self):
        # Find Max confidence for each protein
        for protein in self.prediction:
            self.maxConfidence[protein] = 0
            for term in self.prediction[protein]:
                if self.prediction[protein][term] > self.maxConfidence[protein]:
                    self.maxConfidence[protein] = self.prediction[protein][term]
        vprint("Max Confidence : {} ".format(self.maxConfidence), 15)    
        # Generate protein in prediction count for each threshold
        
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            threshold = numpy.around(threshold, decimals = 2)
            total = 0
            for protein in self.prediction:
                if self.maxConfidence[protein] >= threshold:
                    total += 1
            self.ProteinInPrediction[threshold] = total
            vprint("Proteins at {} : {}".format(threshold, self.ProteinInPrediction[threshold]), 10)
        
        # Generate Terms and truths
        terms = set()
        truths = set()
        for protein in self.prediction:
            for term in self.prediction[protein]:
                terms.add(term)
                if term in self.ProteinTrueTerms[protein]:
                    truths.add(term)
            self.terms[protein] = terms
            self.truth[protein] = truths
            vprint(protein, 15)
            vprint("Terms: {}".format(terms), 15)
            vprint("Truth: {}".format(truths), 15)
        
        print("Info has been set")
        self.generateIntermediateData()
        # Ready to run metrics now
        return True 
    
    ################### Actual work for metrics ##################   
    def generateIntermediateData(self):
        for protein in self.prediction:
            terms       = self.terms[protein]
            truth       = self.truth[protein]
            prediction  = self.prediction[protein]
            ic          = self.IC
            #vprint(ic,5)
            terms_with_IC = terms & ic.keys()
            #vprint(list(ic.keys()), 5)            
            
            #vprint(terms_with_IC, 5)
            badTerms = terms - terms_with_IC
            vprint("Bad terms Set: {}".format(badTerms), 15)
            # Bad terms are those without IC
            terms = terms - badTerms
            vprint("Terms Set: {}".format(terms), 15)
            countPositive = 0
            countNegative = 0
            positive = truth
            vprint("Positive Set: {}".format(positive), 15)
            # Set difference
            negative = terms - truth
            vprint("Negative Set: {}".format(negative), 15)
            for term in terms:
                if term in positive:
                    countPositive += 1
                if term in negative:
                    countNegative += 1  
            for threshold in numpy.arange(0.00, 1.01, 0.01, float):
                threshold = numpy.around(threshold, decimals = 2)
                FP  = 0
                TN  = 0
                TP  = 0
                FN  = 0 
                   
                for term in prediction:
                    if (prediction[term] > threshold 
                        and term in negative):
                        FP += 1
                    if (prediction[term] > threshold 
                        and term in positive):
                        TP += 1       
                            
                FP = FP
                TN = countNegative - FP
                TP = TP
                FN = countPositive - TP
                POS = TP + FP
                NEG = FN + TN
                TRUE = TP + FN
                FALSE = FP + TN
                #vprint(POS,1)
                self.data_unweighted[protein][threshold] = {
                 'FP':FP, 'TN':TN, 'TP':TP, 'FN':FN, 
                 'POS':POS, 'NEG':NEG, 'TRUE':TRUE,'FALSE':FALSE}
                
            #COUNT FOR WEIGHTED METRICS
            
            weightedPositive = 0
            weightedNegative = 0
            # term is the index 
            for term in terms:
                
                if term in positive:
                    try:
                        weightedPositive += ic[term]
                    except KeyError:
                        vprint("{} has no IC value".format(term), 5) 
                if term in negative:
                    try:
                        weightedNegative += ic[term]   
                    except KeyError:
                        vprint("{} has no IC value".format(term), 5)  
            for threshold in numpy.arange(0.00, 1.01, 0.01, float):
                threshold = numpy.around(threshold, decimals = 2)
                FP = 0
                TN = 0
                TP = 0
                FN = 0
                for term in prediction:
                   if (prediction[term] > threshold 
                       and term in negative):
                           
                       try:
                           FP += ic[term]
                       except KeyError:
                           vprint("{} has no IC value".format(term), 5) 
                   if (prediction[term] > threshold 
                       and term in positive):
                           
                       try:        
                           TP += ic[term]       
                       except KeyError:
                           vprint("{} has no IC value".format(term), 5) 
                FP = FP
                TN = weightedNegative - FP
                TP = TP
                FN = weightedPositive - TP
                POS = TP + FP
                NEG = FN + TN
                TRUE = TP + FN
                FALSE = FP + TN
                TOT = Info.tot(ic, truth, prediction)    
                self.data[protein][threshold] = {
                  'FP':FP, 'TN':TN, 'TP':TP, 'FN':FN, 
                  'POS':POS, 'NEG':NEG, 'TRUE':TRUE,'FALSE':FALSE, 'TOT':TOT}
        
        print ("Done with Intermediate value creation")      
        return True
    
    def tot(IC, T, P):
        '''
        Sum the IC of all terms
        
        Input:
        info       : Object
        T   : Set [ Truth      ]
        P   : Set [ Prediction ]
        
        Output:
        [0] : Float
        '''
        
        #vprint("Finding TOT", 5) 
        # Sum of IC values of terms meeting criteria   
        total = 0.0
        # Sum the IC values of every element in T or P
        # Add all truths
        for term in T:
            try:
                total += IC[term]
            except KeyError:
                vprint("{} has no IC value".format(term), 5)  
        # Add predictions that are not truths
        for term in P:
            if term not in T:
                try:
                    total += IC[term]
                except KeyError:
                    vprint("{} has no IC value".format(term), 5)  
        return total

###################### HELPER FUNCTIONS ###################
def vprint(message, priority):
    # Change number for different verbosity
    if priority < 11:
        print(message)
def clear(location):
    # Delete old results/ files
    os.remove(location)
    vprint("Old Results have been cleared", 1)
    
def vwrite(message, location, priority):
    if priority < 5:
        # Append to file
        data = open(location, 'a+')
        data.write(message)
        data.close()


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
        # Build Truth        
        self.ProteinTrueTermsBase = defaultdict(set)
        with open(benchmark_path) as benchmark_input:
            for inline in benchmark_input:
                protein, term = inline.strip().split('\t')
                self.ProteinTrueTermsBase[protein].add(term)
                
                
    def propagateTrueTerms(self):
        '''
        Progate Benchmark terms.
        '''
        
        self.ProteinTrueTerms = defaultdict(set)
        # Key: protein
        # Value: set of benchmark propagated terms
        
        for protein in self.ProteinTrueTermsBase:
            for term in self.ProteinTrueTermsBase[protein]:
                try:
                    ancestors = self.ancestors[term].difference(root_terms)    
                except KeyError:
                    sys.stderr.write("{} not found \n".format(term)) 
                # Add the term    
                self.ProteinTrueTerms[protein].add(term)
                # Set Union 
                self.ProteinTrueTerms[protein] |= ancestors


def read_benchmark(ontology, species, types, benchmark_directory, obopath):
    '''
    Read Benchmark.
    
    Input:
    ontology            : String        {BPO, CCO, MFO, HPO}
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
        sys.stderr.write("Namespace not accepted, choose from 'BPO', 'CCO', 'MFO' and 'HPO'\n")
    elif (species not in legal_species) and (species not in legal_subtypes):
        sys.stderr.write('Species not accepted')
    elif types not in legal_types:
        sys.stderr.write('Type not accepted, choose from "type1","type2" and "typex"\n')
    else:
        # Create name for finding benchmark file
        matchname = ontology.lower() + '_' + species + '_' + types + '.txt'
        vprint("Matchname: {}".format(matchname), 9)
    # Ontology-specific calculations
    full_benchmark_path = benchmark_directory + '/groundtruth/' + 'leafonly_{}.txt'.format(ontology)
    vprint("full_benchmark_path: {}".format(full_benchmark_path), 9)
    ancestor_path       = os.path.splitext(obopath)[0] + "_ancestors_{}.txt".format(ontology)
    vprint("ancestor_path: {}".format(ancestor_path), 9)
    benchmarkListPath   = benchmark_directory + '/lists/' + matchname
    vprint("benchmarkListPath: {}".format(benchmarkListPath), 9)
    # If the List file exists and has data
    if os.path.isfile(benchmarkListPath) and os.path.getsize(benchmarkListPath) > 0:
        handle    = open(benchmark_directory + '/lists/' + matchname, 'r')
        proteins  = set()
        for line in handle:
            proteins.add(line.strip())
        handle.close()
        tempfilename = 'temp_{}_{}_{}.txt'.format(ontology, species, types)
        tempfile = open(benchmark_directory + '/' + tempfilename , 'w')
        for line in open(full_benchmark_path,'r'):
            protein = line.split('\t')[0]
            if protein in proteins:
                tempfile.write(line)
        tempfile.close()
        bench = benchmark(ancestor_path, tempfile.name)
        bench.propagateTrueTerms()
        #os.remove(tempfile.name)
    else:
        print('Benchmark set is empty.\n')
        bench = None
    return bench 


############### Read OBO ################  
def go_ontology_split(OBO):
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
    
    for term in OBO.get_ids(): # loop over node IDs and alt_id's
        if OBO.namespace[term]   == "biological_process": # P
            bpo_terms.add(term)
        elif OBO.namespace[term] == "cellular_component": # C
            cco_terms.add(term)
        elif OBO.namespace[term] == "molecular_function": # F
            mfo_terms.add(term)
        else:
            raise(ValueError,"{} has no namespace".format(term))
    return (mfo_terms, bpo_terms, cco_terms)
 
   
def readOBO(obo_path):
    """
    Takes an OBO file and writes 3 files with ancestors corresponding to the namespace
    
    Input:
    obo_path : filepath of the obo file
    
    Output:
    [0]      : List [Length of BPO, Length of CCO, Length of MFO]
    """
    
    OBO = OboIO.OboReader(open(obo_path)).read()
    # Split the OBO into the three namespaces 
    MFO_terms, BPO_terms, CCO_terms = go_ontology_split(OBO)
    # BPO
    obo_BPO_out = open("{}_ancestors_BPO.txt".format(os.path.splitext(obo_path)[0]), "w")
    for term in BPO_terms:
        ancestors = OBO.get_ancestors(term)
        obo_BPO_out.write("%s\t%s\n" % (term, ",".join(ancestors)))
    obo_BPO_out.close()
    # CCO
    obo_CCO_out = open("{}_ancestors_CCO.txt".format(os.path.splitext(obo_path)[0]), "w")
    for term in CCO_terms:
        ancestors = OBO.get_ancestors(term)
        obo_CCO_out.write("%s\t%s\n" % (term, ",".join(ancestors)))
    obo_CCO_out.close()    
    # MFO
    obo_MFO_out = open("{}_ancestors_MFO.txt".format(os.path.splitext(obo_path)[0]), "w") 
    for term in MFO_terms:
        ancestors = OBO.get_ancestors(term)
        obo_MFO_out.write("%s\t%s\n" % (term, ",".join(ancestors)))    
    obo_MFO_out.close()
    
    OBOcountDict = {'BPO':len(BPO_terms),'CCO':len(CCO_terms),'MFO':len(MFO_terms)}
    return OBOcountDict