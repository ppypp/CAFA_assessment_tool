import numpy
import pickle as cp
import helper
import sys
from collections import defaultdict
import os
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

class Data:
    '''
    New info object
    
    '''
    def __init__(self, ontology, Type):
        '''
        
        '''
        
        self.Benchmark       = None
        self.OBOCount        = 0
        self.obsolete        = None
        self.ready           = False
        
        
        # Dictionary
        self.IC_dict              = defaultdict() 
        # KEY: Protein   Value: [Confidence values in same order as GoTerms]
        self.Prediction            = defaultdict(defaultdict)
        
        # List of GoTerms
        self.Terms                 = []
        self.IC                    = []
        self.PredictionConfidence  = defaultdict(list) 
        self.Truth                 = defaultdict(list) 
        # Determines 
        self.ontology        = ontology
        self.Type            = Type
        vprint ("Running {}, {}".format(ontology, Type), 1)
        self.Info            = defaultdict(list)
        

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
                
                
    def setBenchmark(self, taxon, OBOCounts, benchmark_directory):
        '''
        
        '''
        
        Benchmark = read_benchmark(self.ontology, helper.taxon_name_converter(taxon), self.Type, benchmark_directory, self.OBOPath)
        self.Benchmark = Benchmark
        vprint("Benchmark Done", 9)
        self.OBOCount = OBOCounts[self.ontology]
        vprint("OBO Count for {} : {}".format(self.ontology, self.OBOCount), 9)
        self.TermAncestors = Benchmark.ancestors
        vprint("Term Ancestors", 10)def total(info, T, P):
    '''
    Sum the IC of all terms
    
    Input:
    info       : Object
    T   : Set [ Truth      ]
    P   : Set [ Prediction ]
    
    Output:
    [0] : Float
    '''
    
    # Sum of IC values of terms meeting criteria   
    total = 0.0
    # Sum the IC values of every element in T or P
    # Add all truths
    for term in T:
        try:
            total += info.ic[term][1]
 'w')
        except KeyError:
            pass
    # Add predictions that are not truths
    for term in P:
        if term not in T:
            try:
                total += info.ic[term][1]
            except KeyError:
                pass
    return total
    
        self.ProteinTrueTerms = Benchmark.ProteinTrueTerms
        vprint("Protein True Terms", 10)
        self.TrueTermsCount = len(Benchmark.ProteinTrueTerms)
        
        
    def setIC(self, ic_path):
        '''
        
        '''
        
        IC_old = cp.load(open("{}ia_{}.map".format(ic_path, self.ontology), "rb"))
        IC = {}
        for term in IC_old:
            IC[term] = IC_old[term][1]
        self.IC_dict = IC    
        vprint(self.IC_dict, 15)
        
        
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
        '''
        
        '''
        
        self.ProteinInBenchmark = 0
        self.PredictionProtein = defaultdict(list)     
        if os.path.getsize(self.Path) > 0:
            # Read in prediction file   
            for inrec in open(self.Path,'r'):
                fields     = [i.strip() for i in inrec.split()]
                protein    = fields[0]
                vprint("Protein: {}".format(protein), 10)
                term       = fields[1]
                vprint("Term: {}".format(term), 10)
                confidence = float(fields[2])
                vprint("Confidence: {}".format(confidence), 10)
                # data[Protein] = {Term, Confidence}
                self.PredictionProtein[protein].append({'term': term, 'confidence': confidence})
            # Propagated prediction
            for protein in self.PredictionProtein:
                # Take only Protein that has terms in the benchmark
                if self.ProteinTrueTerms[protein]:
                    self.ProteinInBenchmark += 1
                    for value in self.PredictionProtein[protein]:
                        term       = value['term']
                        confidence = value['confidence']
                        try:
                            ancterms = self.TermAncestors[term].difference(root_terms)
                        except KeyError:
                            # Add unknown term to the obsolete set
                            self.obsolete.add(term)
                            continue
                        # For each term
                        if term in self.Prediction[protein]:
                            # Term already exists, update confidence
                            self.update_confidence(protein, value)
                        else:
                            # Add term to self.predicted_bench
                            # Add confidence and compare with self.true_terms
                            # Regardless of comparision, propagate
                            if term not in root_terms:
                                self.Prediction[protein][term] = self.compare(protein, value)
                                for ancterm in ancterms:
                                    # Make a new TC with ancestors
                                    newValue = {'term' : ancterm, 'confidence' : confidence}
                                    if ancterm in self.Prediction[protein]:
                                        # Term already exists, update confidence
                                        self.update_confidence(protein, newValue)
                                    else:
                                        # Add term to self.predicted_bench
                                        self.Prediction[protein][ancterm] = self.compare(protein, newValue)             
            if self.ProteinInBenchmark == 0:
                print("No protein in this predicted set became a benchmark\n")
        else:
            print('No prediction made in this ontology.\n')

'''            
    def buildLists(self):
        '''

        '''      
        
        for term in self.IC_dict:
            self.Terms.append(term)
            self.IC.append(self.IC_dict[term])
        for protein in self.Prediction:   
            vprint("Protein: {}".format(protein), 3)
            PL    = [0 for i in range(len(self.Terms))]
            Truth = [False for i in range(len(self.Terms))]
            for term in self.Prediction[protein]:
                try:
                    PL[self.Terms.index(term)] = self.Prediction[protein][term]
                    if term in self.ProteinTrueTerms[protein]:
                        Truth[self.Terms.index(term)] = True
                        vprint("True added", 15)
                except ValueError:
                    vprint("Term not in Terms: {}".format(term),3)
            self.PredictionConfidence[protein] = PL
            self.Truth[protein]                = Truth
'''            
            
    def setInfo(self):
        vprint("Set Info", 1)
        for protein in self.Prediction: 
            for term in self.Prediction[protein]:
                self.Info[protein].append(
                { 'Term'      : term, 
                  'IC'        : self.IC_dict[term], 
                  'Confidence': self.Prediction[protein][term][0], 
                  'Truth'     : self.Prediction[protein][term][1]
                  })
                vprint(self.Info, 1)
                
        
    def coverage(self):
        ''' 
        Determine the coverage. 
        
        Output:
        [0] : Float       coverage of the prediction
        '''
        
        return float(self.count_predictions_in_benchmark)/self.count_true_terms  
        
        
    def update_confidence(self, protein, value):
        '''
        Update Confidence for given protein and propagate.
        
        This function compares the confidence value in tc to the confidence in self.predicted
        If tc is larger, than it overwrites the confidence in self.predicted
        
        Input:
        protein : choosen protein 
        '''
        
        # Defined for readablity
        confidence = value['confidence'] 
        term       = value['term']
        
        if confidence > self.Prediction[protein][term][0]:
            # Update the confidence
            self.Prediction[protein][term][0] = confidence
            # Propagate changes if necessary
            for ancterm in self.TermAncestors[term].difference(root_terms):
                if confidence > self.Prediction[protein][ancterm][0]:
                    # Update the confidence
                    self.Prediction[protein][ancterm][0] = confidence
                    
                    
    def compare(self, protein, value):
        '''
        Check if tc['term'] is a True term.
        
        This function compares if tc['term'] is in self.true_terms
        
        Input:
        protein  : Protein for comparision
        tc       : Dictionary KEY: Confidence VALUE: GO_Term
        
        Output:
        [0]      : List [Confidence value, Boolean]
        '''
        
        if value['term'] in self.ProteinTrueTerms[protein]:
            return [value['confidence'], True]
        else:
            return [value['confidence'], False]
            
        
    def getObsolete(self):
        ''' 
        Get obsolete terms used by the prediction team. 
        
        Output:
        [0] : Set      GOTerms that are obsolute
        '''
        
        return(self.obsolete)   
         

########################## END OF DATA #######################################        
############### START OF EVAUATION CODE ################################
def assessMetrics(data, mode):
    Results = {}
    # Call assessAllProteins
    # Evaluate metrics (For this particular Ontology, Type, and both Modes)
    assessAllProteins(data)
    # for each metric
    # Make result object
    # Run over all threshold values from 0 to 1, two signifigant digits
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        threshold = numpy.around(threshold, decimals = 2)
        
        
        
        
    #Do work
    #Store in Results
    
    return Results


def assessAllProteins(data):
    
    
    
    
    #TermsDict = {} # KEY: Term, VALUE: [IC, Truth, Prediction]
    
    
        
    
    

        
        
        
        #Values, WeightedValues = assessProtein[data.Terms, data.IC, Truth, data.Prediction[protein]]
    print("Have made FP, TN, TP, FN for all proteins")
        
        
        
        
        
        
       
        
def assessProtein(Terms, IC, Truth, Prediction):
    # I want IC as a n by 1 vector (ie. a list)
    '''
    # List of Terms 
    Terms = []
    # IC values of terms in respective index
    IC = []
    # Truth values of terms in respective index
    Truth = []
    # Prediction confidence values of terms in respective index
    Prediction = []
    
    
    '''
    Values = {}
    WeightedValues = {}
    # Intailize 
    Positive = Truth
    Negative = ~Truth
    
    # COUNT FOR NON-Weighted
    countPositive = 0
    countNegative = 0
    
    for i in Truth:
        if Positive[i] is True:
            countPositive += 1
        if Negative[i] is True:
            countNegative += 1    
    
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
       FP  = 0
       TN  = 0
       TP  = 0
       FN  = 0 
       
       for term in Prediction:
           if (Prediction[term] > threshold and Negative[Terms.index(term)] is True):
               FP += 1
           if (Prediction[term] > threshold and Positive[Terms.index(term)] is True):
               TP += 1       
           
       FP = FP
       TN = countNegative - FP
       TP = TP
       FN = countPositive - TP
           
       Values[threshold] = [FP, TN, TP, FN] 
        
    
    #COUNT FOR WEIGHTED METRICS
    
    weightedPositive = 0
    weightedNegative = 0
    # term is the index 
    for term in Truth:
        if Positive[term] is True:
            weightedPositive += IC[term]
        if Negative[term] is True:
            weightedNegative += IC[term]   
    
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
       FP = 0
       TN = 0
       TP = 0
       FN = 0
       # term is the index  
       for term in Prediction:
           if (Prediction[term] > threshold and Negative[term] is True):
               FP += IC[term]
           if (Prediction[term] > threshold and Positive[term] is True):
               TP += IC[term]       
           
       FP = FP
       TN = weightedNegative - FP
       TP = TP
       FN = weightedPositive - TP
           
       WeightedValues[threshold] = [FP, TN, TP, FN] 
       
    return Values, WeightedValues
  
########### END of EVALUATION CODE #################
########### START OF RESULT CLASS ###############  
class result:
    def __init__(self, results_path):
        ''' 
        State all variables needed 
        
        Input:
        results_path : String       The directory to store results in
        '''
        self.path           = results_path        

        self.Opt            = 0.0
        self.Value1         = []
        self.Value2         = []
        self.Value3         = []
        self.OptThreshold   = 0.0
        
        self.coverage       = 0.0
        
        
    def update(self, subval1, subval2, subval3, value, threshold):
        '''
        Stores the outpu of each metric in a single object
        
        Input:
        subval1   : List[Float]
        subval2   : List[Float]
        subval3   : List[Float]
        value     : Float
        threshold : Float
        '''

        self.Opt            = float(value)
        self.Value1         = subval1
        self.Value2         = subval2
        self.Value3         = subval3
        self.OptThreshold   = float(threshold)  
        
        
    def teamInfo(self, author, model, keywords, taxon, predictionFile):
        '''
        Store team info
        
        Input:
        author           : String
        model            : String
        keywords         : List[String]
        taxon            : String
        predictionFile   : String
        '''        
        self.author         = author
        self.model          = model
        self.keywords       = keywords
        self.taxon          = taxon
        self.predictionFile = predictionFile    
        
    
    def toFile(self, Ontology, Type, Mode, tool):
        '''
        Method to dump data to file for evaluation
        
        Input:
        Ontology   : String
        Type       : String
        Mode       : String
        tool       : String
        '''
        # Set Type to correct format
        Type = helper.typeConverter(Type) 
        
        p = self.path
        handle  = open(p + "/{}_{}_{}_{}_{}{}_{}_results.txt".format(Ontology, self.taxon, Type, Mode, self.author, self.model, tool), 'w')
        handle.write('!AUTHOR:      \t {}  \n'.format(self.author))
        handle.write('!MODEL:       \t {}  \n'.format(self.model))
        handle.write('!KEYWORDS:    \t {}  \n'.format(self.keywords))
        handle.write('!SPECIES:     \t {}  \n'.format(self.taxon))
        handle.write('!ONTOLOGY:    \t {}  \n'.format(Ontology))
        handle.write('!TYPE:        \t {}  \n'.format(Type))
        handle.write('!MODE:        \t {}  \n'.format(Mode))
        handle.write('<OPTIMAL:\t {:.6f}   \n'.format(self.Opt))
        handle.write('<THRESHOLD:\t {:.2f} \n'.format(self.OptThreshold))
        handle.write('<COVERAGE:\t {:.2f}  \n'.format(self.coverage))

        # Write data
        handle.write('{}\t {}\t {}\t {}\n'.format('Threshold', 'PR', 'RC', 'F')) 
        index = 0
        for threshold in numpy.arange(0.00, 1.01, 0.01, float):
            
            threshold = numpy.around(threshold, decimals = 2)
            try:
                handle.write('>{:.2f}\t {:.6f}\t {:.6f}\t {:.6f}\n'.format(threshold, self.Value1[index], self.Value2[index], self.Value3[index]))
            except (IndexError, TypeError):
                pass
            index += 1
        
        handle.close()       
        
#### END OF RESULT ####
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
        
        # Key: term
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
                
                
    def propagate(self):
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
        bench.propagate()
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
    
def vprint(message, priority):
    # Change number for different verbosity
    if priority < 5:
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