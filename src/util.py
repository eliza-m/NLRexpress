import sys
import numpy as np
from .FeaturesData import *


def printResultsOutput(inputRoot: str, outdir:Path, inputData: FeaturesData, results:dict, outformat="all", cutoff=0.2) :

    countpos = {motif: 0 for motif in results};

    # try:

    if outformat in ["all", "long"]:
        longFile = open( str(outdir) + "/" + str(inputRoot) + '.long.output.txt', 'w' )
        #  Print header
        print("#ProtName", "ResId", "ResName", sep='\t', end='\t', file=longFile)
        for motif in results:
            print(motif, sep='\t', end='\t', file=longFile)
        print(file=longFile)

    if outformat in ['all', "short"]:
        shortFile = open( str(outdir) + "/" + str(inputRoot) + '.short.output.txt', 'w' )
        #  Print header
        print("#ProtName", "ResID_start", "Motif", "Proba", " | 5 positions upstream | Motif | 5 positions downstream | ", sep='\t', file=shortFile)


    for prot in inputData.hmmData:

            seqLength = len(inputData.hmmData[prot]["seq"])
            for i, aa in enumerate( inputData.hmmData[prot]["seq"] ):

                if outformat in ["all", "long"]:
                    print(prot, i + 1, aa, sep='\t', end='\t', file=longFile)

                for motif in results:
                    # anno = annotations[prot][i] if len(annotations) != 0 else '-'
                    if i >= allMotifs[motif]["windLeft"] and i < seqLength - (allMotifs[motif]["motifSpan"] + allMotifs[motif]["windRight"]):

                        if outformat in ["all", "long"]:
                            print( round(results[motif][countpos[motif]][1], 4), end='\t', file=longFile)

                        if outformat in ["all", "short"]:

                            if round(results[motif][countpos[motif]][1], 4) >= cutoff :
                                print(prot, i + 1, aa, motif, round(results[motif][countpos[motif]][1], 4),
                                  inputData.hmmData[prot]["seq"][ i-5 : i ],
                                  inputData.hmmData[prot]["seq"][ i : i + allMotifs[motif]["motifSpan"] + 1 ],
                                  inputData.hmmData[prot]["seq"][ i + allMotifs[motif]["motifSpan"] + 1 : i + allMotifs[motif]["motifSpan"] + 6],
                                  sep='\t', file=shortFile)
                        countpos[motif] += 1
                    elif outformat in ["all", "long"]:
                        print('-', end='\t', file=longFile)

                print(file=longFile)

    # except :
    #     print("Something went wronfg ", sys.exc_info()[0])
    #     # raise
    #



