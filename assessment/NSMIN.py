'''
Class to calculate the normalized semantic distance mininum
'''

import numpy
import assessment.SMIN as S


def output(info, ontology, Type, mode):
    '''
    Main method
    
    k = 2 by convention
    
    Input:
    info     : Object 
    ontology : String     {bpo, cco, mfo}
    Type     : String     {type1, type2}
    mode     : String     {partial, full}
    
    Output:
    [0]      : List[RU, MI, NS, nsmin, nsmin_threshold]
    
    '''
    
    k = 2.0
    # Intialize Variables
    nsmin = float("inf")
    nsmin_threshold = 0.0
    RU   = []
    MI   = []
    NS   = []
    # Run over all threshold values from 0 to 1, two signifigant digits
    for threshold in numpy.arange(0.00, 1.01, 0.01, float):
        
        threshold = numpy.around(threshold, decimals = 2)
        # Run RUMI on given threshold
        ru, mi = rumi_average(info, k, threshold, ontology, Type, mode)
        if ru is None:
            # No prediction above this threshold 
            sval = None
            break
        else:
            RU.append(ru)
            MI.append(mi)
            # Find the S-value for this particular threshold
            sval = S.s(k, ru, mi)
        NS.append(sval)
        if (sval is not None and sval <= nsmin): 
            nsmin           = sval
            nsmin_threshold = threshold
    # Have found the NSmin at this point       
    return ([RU, MI, NS, nsmin, nsmin_threshold])


def norm(info, T, P):
    '''
    Normalize the values
    '''
    # Sum of IC values of terms meeting criteria   
    total = 0.0
    # Sum the IC values of every element in T or P
    # Add all truths
    for term in T:
        try:
            total += info.ic[term][1]
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
    
    

def rumi_average(info, k, threshold, ontology, Type, mode):
    '''
    Remaing Uncertainity, Misinformation   
    
    At a particular threshold, averaged over all proteins
    
    Input:
    info       : 
    k          : Integer  
    threshold  : Float      
    ontology   : String     {bpo, cco, mfo}
    Type       : String     {type1, type2}
    mode       : String     {partial, full}
    
    Output:
    [0]        : Float    Remaining Uncertainity
    [1]        : Float    Misinformation
    '''
    
    RU = 0.0  
    MI = 0.0
    count = 0
    
    info.count_above_threshold[threshold] = 0 # Ne in paper 
    
    for protein in info.predicted_bench:
        data  = open(info.path + "/NSMIN/{}/{}/{}/{}/{}.txt".format(ontology, Type, mode, threshold, protein), 'w')
        
        T = set()
        P = set()
        # Find T
        T = info.true_terms[protein]
        # Find P (within threshold)
        for term in info.predicted_bench[protein]:
            # If it is above the threshold, add to the P set
            if info.predicted_bench[protein][term][0] >= threshold:
                try:                
                    data.write('{}\t {}\t {}\n'.format(term, info.predicted_bench[protein][term][0], info.ic[term][1]))
                except KeyError:
                    data.write('{}\t {}\t {}\n'.format(term, "ERROR", "ERROR"))   
                P.add(term) 
        data.close()
        # Calculate ru & mi     
        r = S.ru(info, T, P)
        m = S.mi(info, T, P)
        n = norm(info, T, P)
        if r is not None and m is not None:
            try:
                rn = r/n
                mn = m/n
            except ZeroDivisionError:
                rn = None
                mn = None
        # Check both to ensure calculation worked
        if rn is not None and mn is not None:
            RU += rn
            MI += mn 
            
        if mode == 'partial':
            count = info.count_predictions_in_benchmark
        elif mode == 'full':
            count = info.count_true_terms
             
    try:
        remain  = RU / count
        misinfo = MI / count
        
    except ZeroDivisionError:
        remain  = None
        misinfo = None
        print("No prediction is made above the %.2f threshold\n" % threshold)
            
    return remain, misinfo        