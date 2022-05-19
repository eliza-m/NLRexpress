from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from .MultiFASTA import *
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
    "extendedEDVID": {"windLeft": 5, "windRight": 5, "motifSpan": 12},
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
    seqData: MultiFASTA
    hmmData: dict
    X: dict


def generateFeatures( inputFasta:Path, outdir:Path, motifs:dict, skipJhmmer=False, annotations={} ) -> FeaturesData :
    """

    :param inputFasta:
    :param outdir:
    :param motifs:
    :param skipJhmmer:
    :return:
    """

    # STEP 1: Check FASTA file
    seqData = parseFastaFile( inputFasta )

    # STEP 2: Run Jhmmer if needed
    if not skipJhmmer:

        # not all Jhmmer params were added, as the training was done with specific params and therefore for
        # new data the same params should be used. For now only cpuNum and targetDB path are customizable
        # from here
        params = { 'cpuNum':4,
                    'targetDB': "../hmmer_db/targetDB.fasta" }

        runJhmmer(inputFasta, outdir, params)


    # STEP 3: Parse and organize HMM data; generate input file for later use;

    hmm_it1 = parse_hmm_multiprot( str(outdir) + "/" + str(inputFasta.stem) + '-1.hmm' )
    hmm_it2 = parse_hmm_multiprot( str(outdir) + "/" + str(inputFasta.stem) + '-2.hmm' )

    hmmData = generateInputFile( seqData, hmm_it1, hmm_it2, inputFasta, outdir, annotations )


    ###############################################################
    # STEP 4: Create the X matrices for the requested motifs specs

    X = { motif: [] for motif in motifs }
    # data[name] = {'seq': seq,
    #               'features': []}

    for prot in hmmData:
        seqLength = len( hmmData[prot]['seq'] )

        for i in range(seqLength):
            for motif in motifs:
                windLeft = motifs[motif]["windLeft"]
                windRight = motifs[motif]["windRight"]
                motifSpan = motifs[motif]["motifSpan"]

                if i >= windLeft and i < seqLength - ( motifSpan + windRight) :
                    features = hmmData[prot]['features']

                    X[motif].append([])
                    for w in range(windLeft * (-1), motifSpan + windRight + 1):
                        X[motif][-1] += features[i+w]

    return FeaturesData( seqData=seqData, hmmData=hmmData, X=X)




def generateInputFile( seqData:MultiFASTA, hmm_it1:dict, hmm_it2:dict, inputFasta:Path, outdir:Path, annotations:dict  ) -> dict:
    """
    :param seqData:
    :param hmm_it1:
    :param hmm_it3:
    :param outdir:
    :return:
    """

    data = {}
    # try:
    with open( str(outdir) + "/" + str(inputFasta.stem) + '.input', 'w' ) as outfile :

            for n in range(len(seqData.df.Accession)):
                name = seqData.df.Accession[n]
                seq = re.sub( r"[\n\t\s]*", "", seqData.df.Sequence[n] )
                data[name] = {'seq':seq,
                              'features':[]}

                for i, aa in enumerate(seq):
                    anno = annotations[name][i] if len(annotations) != 0 else '-'
                    print( name, i+1, aa, anno, sep='\t', end='\t', file=outfile )
                    data[name]['features'].append([])

                    for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                        print(val, end='\t', file=outfile)
                        data[name]['features'][-1].append(val)

                    for k, (key, val) in enumerate( hmm_it2[name][i].items() ):
                        print(val, end='\t', file=outfile)
                        data[name]['features'][-1].append(val)

                    print(file=outfile)
    # except :
    #     print("Something went wronfg ", sys.exc_info()[0])
    #     raise

    return data



def parseFastaFile( input:Path) -> MultiFASTA :

    # TODO implement a better FASTA file check

    seqData = MultiFASTA()

    try:
        seqData.read_fasta(input)

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

    jhhmerLog = subprocess.run(["jackhmmer",
                                "--cpu", str(params["CpuNum"]),
                                "-o", "/dev/null",
                                "-N", "2",
                                "-E", "1e-5",
                                "--domE", "1e-5",
                                "--noali",
                                "--chkhmm", str(outdir) + "/" + str(inputFasta.stem),
                                inputFasta,
                                params["targetDB"],
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
                        raise Exception( "Problem parsing HMM file: no COMP line was found for prot ", name, line )
                    elif length == 0:
                        raise Exception( "Problem parsing HMM file: no LENGTH line was found for prot ", name )

                    else:
                        l = line.split()
                        if len(l) > 3:

                            if (i - start) % 3 == 0:
                                hmm[name].append({})
                                for i, val in enumerate(l[1:21]):
                                    if val == '*': val = 'inf'
                                    hmm[name][-1][header1[i]] = float(val)

    except:
       raise Exception( "Something went wrong when parsing the HMM file ", hmmFile )

    return hmm





