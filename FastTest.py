from assessment.Fmetric import FMAX
from assessment.WeightedFmetric import WFMAX
from assessment.Smetric import SMIN
from assessment.NormalizedSmetric import NSMIN
from assessment.Tools import vprint
import pickle as cp

if __name__=='__main__':
    '''
    Once you have read in the prediction and benchmark for a particular file, 
    you can use the temp files saved to skip that part of the process
    Only work for 1 Ontology/Type at a time,
    but signifgantly speeds up the iterative cycle
    '''
    
    info = cp.load(open("Info.info", "rb"))
    #for mode in ['partial']:
    for mode in ['partial', 'full']:
        info.setMode(mode)
        vprint("Mode: {}".format(mode), 9)
        # RUN METRICS
        
        print("FMAX")
        print(FMAX(info))
        print("WFMAX")
        print(WFMAX(info))
        print("SMIN")
        print(SMIN(info))
        print("NSMIN")
        print(NSMIN(info))
        vprint("Done", 8)