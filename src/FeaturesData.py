from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import logging
from datetime import datetime
import sys
import subprocess
import re

# allMotifs= {
#     "rdhhhdhEDVID" : { "windLeft": 5, "windRight": 5, "motifSpan": 12},
#     "hhGRE"        : { "windLeft": 5, "windRight": 5, "motifSpan":  5},
#     "GmGGvGKTT"    : { "windLeft": 5, "windRight": 5, "motifSpan":  9},
#     "FDhrhWhshs"   : { "windLeft": 5, "windRight": 5, "motifSpan": 10},
#     "KRhhhhDD"     : { "windLeft": 5, "windRight": 5, "motifSpan":  8},
#     "KhhhTTR"      : { "windLeft": 5, "windRight": 5, "motifSpan":  7},
#     "LseeeSWeLF"   : { "windLeft": 5, "windRight": 5, "motifSpan": 10},
#     "CFLYCSLFP"    : { "windLeft": 5, "windRight": 5, "motifSpan":  9},
#     "GLPLA"        : { "windLeft": 5, "windRight": 5, "motifSpan":  5},
#     "MHD"          : { "windLeft": 5, "windRight": 5, "motifSpan":  3},
#     "LxxLxL"          : { "windLeft": 5, "windRight": 5, "motifSpan":  6}
# }

allMotifs = {
    "extEDVID": {"windLeft": 5, "windRight": 5, "motifSpan": 12},

    "bA": {"windLeft": 5, "windRight": 5, "motifSpan": 10},
    "aA": {"windLeft": 5, "windRight": 5, "motifSpan": 7},
    "bC": {"windLeft": 5, "windRight": 5, "motifSpan": 8},
    "aC": {"windLeft": 5, "windRight": 5, "motifSpan": 6},
    "bDaD1": {"windLeft": 5, "windRight": 5, "motifSpan": 16},
    #"aD1": {"windLeft": 5, "windRight": 5, "motifSpan": 5},
    "aD3": {"windLeft": 5, "windRight": 5, "motifSpan": 13},

    "VG": {"windLeft": 5, "windRight": 5, "motifSpan": 5},
    "P-loop": {"windLeft": 5, "windRight": 5, "motifSpan": 9},
    "RNSB-A": {"windLeft": 5, "windRight": 5, "motifSpan": 10},
    "Walker-B": {"windLeft": 5, "windRight": 5, "motifSpan": 8},
    "RNSB-B": {"windLeft": 5, "windRight": 5, "motifSpan": 7},
    "RNSB-C": {"windLeft": 5, "windRight": 5, "motifSpan": 10},
    "RNSB-D": {"windLeft": 5, "windRight": 5, "motifSpan": 9},
    "GLPL": {"windLeft": 5, "windRight": 5, "motifSpan": 5},
    "MHD": {"windLeft": 5, "windRight": 5, "motifSpan": 3},
    "LxxLxL": {"windLeft": 5, "windRight": 5, "motifSpan": 6}
}



@dataclass
class FeaturesData:
    """

    """
    seqData: dict
    hmmData: dict


def generateFeatures( inputFasta:Path, outdir:Path, motifs:dict, cpuNum:int, skipJhmmer:bool, writeInputFile:bool, annotations={}) -> FeaturesData :
    """

    :param inputFasta:
    :param outdir:
    :param motifs:
    :param skipJhmmer:
    :return:
    """

    logger = logging.getLogger("")

    # STEP 1: Check FASTA file

    logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Checking FASTA file - started')

    processesInputFasta = Path(str(outdir) + "/" + str(inputFasta.stem) + '.fasta_proc')
    seqData = processFastaFile( inputFasta, processesInputFasta )

    logger.info( datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Checking FASTA file - done')


    # STEP 2: Run Jhmmer if needed

    if skipJhmmer == 'False':

        # not all Jhmmer params were added, as the training was done with specific params and therefore for
        # new data the same params should be used. For now only cpuNum and targetDB path are customizable
        # from here
        params = { 'cpuNum':cpuNum,
                    'targetDB': "hmmer_db/targetDB.fasta" }

        logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running JackHMMER - started')
        jhmmerLog = runJhmmer(processesInputFasta, outdir, params)
        logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Running JackHMMER - done')


    # STEP 3: Parse and organize HMM data; generate input file for later use;

    logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Preparing features: Parsing HMM profile - started')

    try:
        hmmFile1 = str(outdir) + "/" + str(processesInputFasta.stem) + '-1.hmm'
        hmm_it1 = parse_hmm_multiprot( hmmFile1 )
    except FileNotFoundError:
        raise FileNotFoundError('Preparing features: HMM profile iteration 1 was not found at. Execution stopped ')

    try:
        hmmFile2 = str(outdir) + "/" + str(processesInputFasta.stem) + '-2.hmm'
        hmm_it2 = parse_hmm_multiprot( hmmFile2 )

    except FileNotFoundError:
        logger.warning(
            datetime.now().strftime(
                "%d/%m/%Y %H:%M:%S") + ':\t' + 'Preparing features: HMM profile iteration 2 was not found at: '
                + hmmFile2 + '. The first iteration HMM profile will be used alone.')
        hmm_it2 = parse_hmm_multiprot( hmmFile1 )


    hmmData = generateInputFile( seqData, hmm_it1, hmm_it2, processesInputFasta, outdir, annotations, writeInputFile )

    return FeaturesData( seqData=seqData, hmmData=hmmData)




def generateXmat (FeaturesData:dict, motif:str):
    ###############################################################
    # STEP 4: Create the X matrices for the requested motifs specs

    X = []

    logger = logging.getLogger("")
    logger.info( datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Preparing features: NN input for motif ' + motif + ' started')

    for prot in FeaturesData.hmmData:
        seqLength = len( FeaturesData.seqData[prot] )

        for i in range(seqLength):

                windLeft = allMotifs[motif]["windLeft"]
                windRight = allMotifs[motif]["windRight"]
                motifSpan = allMotifs[motif]["motifSpan"]

                if i >= windLeft and i < seqLength - ( motifSpan + windRight) :
                    features = FeaturesData.hmmData[prot]

                    X.append([])
                    for w in range(windLeft * (-1), motifSpan + windRight + 1):
                        X[-1] += features[i+w]

    logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' + 'Preparing features: NN input for motif ' + motif + ' done')

    return X

def generateInputFile( seqData:dict, hmm_it1:dict, hmm_it2:dict, inputFasta:Path, outdir:Path, annotations:dict, writeInputFile:bool  ) -> dict:
    """
    :param seqData:
    :param hmm_it1:
    :param hmm_it3:
    :param outdir:
    :return:
    """

    data = {}
    # try:

    if writeInputFile == 'True':
        outfile = open( str(outdir) + "/" + str(inputFasta.stem) + '.input', 'w' )



    for name in seqData:
        seq = seqData[name]
        data[name] = []

        for i, aa in enumerate(seq):
            if writeInputFile == 'True':
                anno = annotations[name][i] if len(annotations) != 0 else '-'
                print( name, i+1, aa, anno, sep='\t', end='\t', file=outfile )
            data[name].append([])

            for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                if writeInputFile == 'True':
                    print(val, end='\t', file=outfile)
                data[name][-1].append(val)

            if name in hmm_it2:
                for k, (key, val) in enumerate( hmm_it2[name][i].items() ):
                    if writeInputFile == 'True':
                        print(val, end='\t', file=outfile)
                    data[name][-1].append(val)
            else:
                for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                    if writeInputFile == 'True':
                        print(val, end='\t', file=outfile)
                    data[name][-1].append(val)

            if writeInputFile == 'True': print(file=outfile)
    # except :
    #     print("Something went wronfg ", sys.exc_info()[0])
    #     raise

    return data



def processFastaFile( input:Path, output:Path) -> dict :

    # TODO implement a better FASTA file check

    seqData = {}

    try:
        with open( input, 'r' ) as inputFile :
            lines = inputFile.readlines()
            for i, line in enumerate(lines):
                if line[0] == ">":
                    if i>0:
                        seqData[ name ] = seq
                    name = line.split()[0][1:]
                    seq = ''
                else:
                    seq += line[:-1]
            seqData[name] = seq

        if len(seqData) > 1001:
            raise("The input FASTA file contains more than the maximum allowed of 100 sequences. Please use splitFasta.py to split your input file.")

        with open(output, 'w') as outputFile:
            for name in seqData:
                print('>', name, sep='', file=outputFile)
                print(''.join(seqData[name].split()), file=outputFile)

    except OSError as e:
        print("FASTA file error:", sys.exc_info()[0])
        raise

    return seqData







def runJhmmer( inputFasta:Path, outdir:Path, params:dict ) -> str :
    """
    Runs JackHMMER for the input sequence and checks the log.
    :param inputFasta: Path to the FASTA file (multi sequence FASTA allowed)
    :param outdir: Path to the output directory where the HMM profiles will be saved (two files for iterations 1 and 2)
    :param params: Dictorary with custom parameters for JackHMMER (more will be added)
    :return: the stdout and stderr of the JackHMMER output
    """
    # TODO add more customizable params
    # not all Jhmmer params were added, as the training was done with specific params and therefore for
    # new data the same params should be used. For now only cpuNum and targetDB path are customizable
    # from here

    scriptDir = Path(__file__).resolve().parents[1]
    #jhhmerLog = subprocess.run(["/usr/local/bin/jackhmmer",
    jhhmerLog = subprocess.run(["jackhmmer",
                                "--cpu", str(params["cpuNum"]),
                                "-o", "/dev/null",
                                "-N", "2",
                                "-E", "1e-5",
                                "--domE", "1e-5",
                                "--noali",
                                "--chkhmm", str(outdir) + "/" + str(inputFasta.stem),
                                inputFasta,
                                str(scriptDir) + '/' + params["targetDB"],
                                ],
                   stdout=subprocess.PIPE)

    valideteJhmmerLog(jhhmerLog)

    return jhhmerLog




def valideteJhmmerLog(jhhmerLog:str):
    # not implemented
    # TODO check the Jhmmer log for red flags

    return 1;




def parse_hmm_multiprot( hmmFile:Path ) -> dict:
    """

    :param hmmFile:
    :return:
    """

    # TODO Better try catch system to check for inconsistancies

    header1 = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    header2 = ["m->m", "m->i", "m->d", "i->m", "i->i", "d->m", "d->d"]

    hmm = {}

    logger = logging.getLogger("")

    try:
        with open(hmmFile, 'r') as f:
            lines = f.readlines()

            start = 0
            hascomp = 0
            length = 0

            for i, line in enumerate(lines):
                if line[0:4] == "NAME":
                    name = line.split()[1]
                    if name[-3:]=="-i1": name = name[:-3]

                    hmm[ name ] = []
                    start = 0
                    hascomp = 0

                elif line[0:4] == "LENG":
                    length = int( line.split()[1] )

                elif line[0:3] == 'HMM':
                    start = i + 5

                elif "COMPO" in line:
                    hascomp = 1

                elif start > 0 and i >= start and i <= 3*(start + length) :
                    if not hascomp:
                        logger.error(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' +
                                     "Problem parsing HMM file: no COMP line was found for prot " + name + line )
                        raise Exception( "Problem parsing HMM file: no COMP line was found for prot ", name, line )

                    elif length == 0:
                        logger.error(datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' +
                                     "Problem parsing HMM file: no LENGTH line was found for prot " + name)
                        raise Exception( "Problem parsing HMM file: no LENGTH line was found for prot ", name )

                    else:
                        l = line.split()
                        if len(l) > 3:

                            if (i - start) % 3 == 0:
                                hmm[name].append({})
                                for i, val in enumerate(l[1:21]):
                                    if val == '*': val = 'inf'
                                    hmm[name][-1][header1[i]] = float(val)

        return hmm

    except FileNotFoundError:
        logger.error(
            datetime.now().strftime(
                "%d/%m/%Y %H:%M:%S") + ':\t' + 'Preparing features: HMM profile not found at: ' + str(hmmFile))
        raise FileNotFoundError('Preparing features: HMM profile not found at: ' + str(hmmFile))

    except:
       logger.error( datetime.now().strftime("%d/%m/%Y %H:%M:%S") + ':\t' +
                     "Something went wrong when parsing the HMM file " + str(hmmFile) )
       raise Exception( "Something went wrong when parsing the HMM file ", str(hmmFile) )







