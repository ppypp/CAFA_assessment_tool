from assessment.Fmetric import FMAX
from assessment.WeightedFmetric import WFMAX
from assessment.Smetric import SMIN
from assessment.NormalizedSmetric import NSMIN
from assessment.Tools import vprint
import pickle as cp

if __name__=='__main__':
    info = cp.load(open("Info.info", "rb"))
    for mode in ['partial', 'full']:
        info.setMode(mode)
        vprint("Mode: {}".format(mode), 9)
        # RUN METRICS
        #print(info.data_unweighted)
        print("FMAX")
        print(FMAX(info))
        #print("WFMAX")
        #print(WFMAX(info))
        #print("SMIN")
        #print(SMIN(info))
        #print("NSMIN")
        #print(NSMIN(info))
        vprint("Done", 8)