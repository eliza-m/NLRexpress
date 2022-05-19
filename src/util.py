import sys
import numpy as np
from .FeaturesData import *


def printResultsOutput(inputRoot: str, outdir:Path, inputData: FeaturesData, results:dict, outformat="all", cutoff=0.2) :

    countpos = {motif: 0 for motif in results};

    # try:

    if outformat in ["all", "long"]:
        longFile = open( str(outdir) + "/" + str(inputRoot) + '.long.output.txt', 'w' )
        #  Print header
        print('{0:<50} {1:>6} {2:>3}'.format("#ProtName", "resid", "aa"), sep='\t', end='\t', file=longFile)
        for motif in results:
            print('{0:>8}'.format(motif), sep='\t', end='\t', file=longFile)
        print(file=longFile)

    if outformat in ['all', "short"]:
        shortFile = open( str(outdir) + "/" + str(inputRoot) + '.short.output.txt', 'w' )
        #  Print header
        print('{0:<50} {1:>6} {2:>15} {3:>8}   | {4:>5} | {5:>15} | {6:5>} | '.format("#ProtName", "ResId", "Motif", "Proba", "-5pos", "MotifSeq", "+5pos"), sep='\t', file=shortFile)


    for prot in inputData.seqData:

            seq = inputData.seqData[prot]
            seqLength = len(seq)
            for i, aa in enumerate( seq ):

                if outformat in ["all", "long"]:
                    print('{0:<50} {1:>6} {2:>3}'.format(prot, i + 1, aa), sep='\t', end='\t', file=longFile)

                for motif in results:
                    # anno = annotations[prot][i] if len(annotations) != 0 else '-'
                    if i >= allMotifs[motif]["windLeft"] and i < seqLength - (allMotifs[motif]["motifSpan"] + allMotifs[motif]["windRight"]):

                        if outformat in ["all", "long"]:
                            print( '{0:8.2f}'.format(100 * results[motif][countpos[motif]][1]), end='\t', file=longFile)

                        if outformat in ["all", "short"]:

                            if round(results[motif][countpos[motif]][1], 4) >= cutoff :
                                print('{0:<50} {1:>6} {2:>15} {3:8.2f}   | {4:>5} | {5:>15} | {6:>5} | '.format(
                                  prot, i + 1, motif, 100 * results[motif][countpos[motif]][1],
                                  seq[ i-5 : i ],
                                  seq[ i : i + allMotifs[motif]["motifSpan"] ],
                                  seq[ i + allMotifs[motif]["motifSpan"] : i + allMotifs[motif]["motifSpan"] + 5]),
                                  sep='\t', file=shortFile)
                        countpos[motif] += 1
                    elif outformat in ["all", "long"]:
                        print('{0:>8}'.format('-'), end='\t', file=longFile)

                print(file=longFile)

    # except :
    #     print("Something went wronfg ", sys.exc_info()[0])
    #     # raise
    #



