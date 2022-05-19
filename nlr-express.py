import os
import click
import subprocess
from bioservices.apps import FASTA
from pathlib import Path

from src.ModelData import *
from src.ModuleData import *
from src.FeaturesData import *
from src.util import *

sys.path.insert(0, os.path.abspath("."))

@click.group(chain=True, invoke_without_command=True)
def cli():
    pass


@click.command()
@click.option('--input', required=True, help='Input FASTA file ')
@click.option('--outdir', required=True, help='Output folder')
@click.option('--module', required=True, help="""\b
                        Predifined prediction modules :
                          - cc  :  CCexpress contains motif predictors for the CC domain:
                                                    extendedEDVID: rdhhhdhEDVID
                          - nbs :  NBSexpress contains motif predictors for the NBS/NBARC domain:
                                                    VG:            bbGRx
                                                    P-loop:        GbGGbGKTT
                                                    RNSB-A:        FDbrhWhshs
                                                    Walker-B:      KRbbbbDD
                                                    RNSB-B:        KbbbTTR
                                                    RNSB-C:        LseeeSWeLF
                                                    RNSB-D:        CFLYCSLFP
                                                    GLPL:          GLPLA
                                                    MHD:           bHD
                          - lrr  :  LRRexpress contains motif predictors for the NBS/NBARC domain:
                                                    LRR pattern:   LxxLxL
                          - all  :  All modules above: CC, NBS and LRR espress """)
@click.option('--outformat', required=False, default="all", help="""\b
                        Output format layout :
                          - short          : Only the hits with more than 20% probability will be printed ( one line per motif )
                          - long           : All predicted residues will be printed ( one line per residue )
                          - all [default]  : Both short and long output formats are generated  """)



def predict( input: Path, outdir: Path, module: str, outformat:str, skipJhmmer=True ):
    """Predict NLR-related motifs"""

    if module == 'cc':
        motifs = {key: allMotifs[key] for key in ['extendedEDVID']}
    elif module == 'nbs':
        motifs = { key: allMotifs[key] for key in
                                     ["VG", "P-loop", "Walker-B", "RNSB-A", "RNSB-B", "RNSB-C", "RNSB-D", "GLPL", "MHD" ] }
    elif module == 'lrr':
        motifs = { key: allMotifs[key] for key in ['LxxLxL'] }
    elif module == 'all':
        motifs = allMotifs
    else:
        print("No such module present. Please use one of the followin: cc, nbs, lrr, all")
        raise()

    inputData = generateFeatures( inputFasta=Path(input), outdir=Path(outdir), motifs=motifs, skipJhmmer=skipJhmmer, annotations={})

    results = {}

    if module in ("cc", "all") :
        CCexpress = ModuleData.loadModels(
            modelsPath = { 'extendedEDVID': 'models/MLP_CC_extendedEDVID.pkl' })
        for p in CCexpress.predictors:
            results[p] = CCexpress.predictors[p].model.predict_proba( inputData.X[ p ] )

    if module in ("nbs", "all") :
        NBSexpress = ModuleData.loadModels(
            modelsPath={
                          'VG': 'models/MLP_NBS_VG.pkl',
                         'P-loop': 'models/MLP_NBS_P-loop.pkl',
                         'RNSB-A': 'models/MLP_NBS_RNSB-A.pkl',
                         'RNSB-B': 'models/MLP_NBS_RNSB-B.pkl',
                         'RNSB-C': 'models/MLP_NBS_RNSB-C.pkl',
                         'RNSB-D': 'models/MLP_NBS_RNSB-D.pkl',
                         'Walker-B': 'models/MLP_NBS_Walker-B.pkl',
                         'GLPL': 'models/MLP_NBS_GLPL.pkl',
                         'MHD': 'models/MLP_NBS_MHD.pkl'
                         })
        for p in NBSexpress.predictors:
            results[p] = NBSexpress.predictors[p].model.predict_proba( inputData.X[ p ] )


    if module in ("lrr", "all") :
        LRRexpress = ModuleData.loadModels(
            modelsPath={'LxxLxL': 'models/MLP_LRR_LxxLxL.pkl'})
        for p in LRRexpress.predictors:
            results[p] = LRRexpress.predictors[p].model.predict_proba( inputData.X[ p ] )

    printResultsOutput(Path(input).stem, outdir, inputData, results, outformat)


if __name__ == '__main__':
    predict()














