from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from .MultiFASTA import *
import sys
import subprocess
import re

@dataclass
class FeaturesData:
    """
        deals with generating the features
    """
    seqData: MultiFASTA
    X: list



def generate_features( input:Path, outdir:Path, skipJhmmer=False ) -> FeaturesData :



    #################################
    # STEP 1: Check FASTA file
    # TODO implement a proter FASTA file check

    seqData = MultiFASTA()
    X = []

    try:
        seqData.read_fasta(input)


    except OSError as e:
        print("FASTA file error:", sys.exc_info()[0])
        raise





    #################################
    # STEP 2: Run Jhmmer
    # TODO check the Jhmmer log for red flags

    cpuNum = 4
    targetDB = "../hmmer_db/targetDB.fasta"

    if not skipJhmmer:

        jhhmerLog = subprocess.run(["jackhmmer",
                                "--cpu", str(cpuNum),
                                "-o", "/dev/null",
                                "-N", "2",
                                "-E", "1e-5",
                                "--domE", "1e-5",
                                "--noali",
                                "--chkhmm", str(outdir) + "/" + str(input.stem),
                                input,
                                targetDB,
                                ],
                   stdout=subprocess.PIPE)





    ####################################
    # STEP 3: Parse and organize data
    # TODO check for inconsistancies

    hmm_it1 = parse_hmm_multiprot( str(outdir) + "/" + str(input.stem) + '-1.hmm' )
    hmm_it2 = parse_hmm_multiprot( str(outdir) + "/" + str(input.stem) + '-2.hmm')


    # try:

    with open( str(outdir) + "/" + str(input.stem) + '.input', 'w' ) as outfile :

            for n in range(len(seqData.df.Accession)):
                name = seqData.df.Accession[n]
                seq = re.sub( r"[\n\t\s]*", "", seqData.df.Sequence[n] )

                for i, aa in enumerate(seq):

                    print( name, i+1, aa, '-', sep='\t', end='\t', file=outfile )
                    X.append([])

                    for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                        print(val, end='\t', file=outfile)
                        X[-1].append(val)

                    for k, (key, val) in enumerate( hmm_it2[name][i].items() ):
                        print(val, end='\t', file=outfile)
                        X[-1].append(val)

                    print(file=outfile)
    # except :
    #     print("Something went wronfg ", sys.exc_info()[0])
    #     raise

    return FeaturesData( seqData=seqData, X=X)








def parse_hmm_multiprot( hmmFile:Path ) -> dict:
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



