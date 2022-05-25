import click
from src.ModuleData import *
from src.util import *
import logging
from datetime import datetime
import os


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
                                                    extEDVID: rdhhhdhEDVID
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

@click.option('--cpunum', required=False, default=4,  help="""\b
                        Set the number of CPU threads to be used. """)

@click.option('--writeinputfile', required=False, default=True,  help="""\b
                        Write input file to be used in case of job intreruption. """)

@click.option('--skipjhmmer', required=False, default="False", hidden=True)


def predict( input: Path, outdir: Path, module: str, outformat:str, cpunum:int, skipjhmmer:bool, writeinputfile:bool ):
    """Predict NLR-related motifs"""

    setupLogger(input, outdir)
    logger = logging.getLogger("")
    logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + '############ NLRexpress started ############ ')
    logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Input FASTA: ' + str(input) )

    scriptDir = Path(__file__).resolve().parent

    if module == 'cc':
        motifs = {key: allMotifs[key] for key in ['extEDVID']}
    elif module == 'nbs':
        motifs = { key: allMotifs[key] for key in
                                     ["VG", "P-loop", "Walker-B", "RNSB-A", "RNSB-B", "RNSB-C", "RNSB-D", "GLPL", "MHD" ] }
    elif module == 'lrr':
        motifs = { key: allMotifs[key] for key in ['LxxLxL'] }
    elif module == 'all':
        motifs = allMotifs
    else:
        print("No such module present. Please use one of the following: cc, nbs, lrr, all")
        raise()



    inputData = generateFeatures( inputFasta=Path(input), outdir=Path(outdir), motifs=motifs, skipJhmmer=skipjhmmer, cpuNum=cpunum, annotations={}, writeInputFile=writeinputfile)
    results = {}

    if module in ("cc", "all") :
        logger.info(
            datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running CCexpress : started')
        CCexpress = ModuleData.loadModels(
            modelsPath = { 'extEDVID': str(scriptDir) + '/models/MLP_CC_extEDVID.pkl' })
        for p in CCexpress.predictors:
            X = generateXmat( inputData, p)
            results[p] = CCexpress.predictors[p].model.predict_proba( X )
            logger.info(
                datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running CCexpress : Running ' + p + ' predictor: ...done')
        logger.info(
            datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running CCexpress : done')

    if module in ("nbs", "all") :
        logger.info(
            datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running NBSexpress : started')
        NBSexpress = ModuleData.loadModels(
            modelsPath={
                          'VG': str(scriptDir) + '/models/MLP_NBS_VG.pkl',
                         'P-loop': str(scriptDir) + '/models/MLP_NBS_P-loop.pkl',
                         'RNSB-A': str(scriptDir) + '/models/MLP_NBS_RNSB-A.pkl',
                         'RNSB-B': str(scriptDir) + '/models/MLP_NBS_RNSB-B.pkl',
                         'RNSB-C': str(scriptDir) + '/models/MLP_NBS_RNSB-C.pkl',
                         'RNSB-D': str(scriptDir) + '/models/MLP_NBS_RNSB-D.pkl',
                         'Walker-B': str(scriptDir) + '/models/MLP_NBS_Walker-B.pkl',
                         'GLPL': str(scriptDir) + '/models/MLP_NBS_GLPL.pkl',
                         'MHD': str(scriptDir) + '/models/MLP_NBS_MHD.pkl'
                         })
        for p in NBSexpress.predictors:
            X = generateXmat(inputData, p)
            results[p] = NBSexpress.predictors[p].model.predict_proba( X )
            logger.info(
                datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running NBSexpress : Running ' + p + ' predictor: ...done')
        logger.info(
            datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running NBSexpress ...done')


    if module in ("lrr", "all") :
        logger.info(
            datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running LRRexpress ...started')
        LRRexpress = ModuleData.loadModels(
            modelsPath={'LxxLxL': str(scriptDir) + '/models/MLP_LRR_LxxLxL.pkl'})
        for p in LRRexpress.predictors:
            X = generateXmat(inputData, p)
            results[p] = LRRexpress.predictors[p].model.predict_proba( X )
            logger.info(
                datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running LRRexpress : Running ' + p + ' predictor: ...done')
        logger.info(
            datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running LRRexpress ...done')

    logger.info(
        datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Printing final results...started')


    printResultsOutput(Path(input).stem, outdir, inputData, results, outformat)

    logger.info(
        datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Printing final results...done')

    logger.info( datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + '############ NLRexpress finished ############ ')






def setupLogger(input: Path, outdir: Path) :
    logging.basicConfig(filename= str(outdir) + '/' + Path(input).stem + '.log', level=logging.DEBUG)
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    logger = logging.getLogger('').addHandler(console)
    return logger



if __name__ == '__main__':
    predict()










