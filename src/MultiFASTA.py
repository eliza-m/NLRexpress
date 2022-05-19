from bioservices.apps import FASTA, MultiFASTA
import logging
from datetime import datetime
import re


class MultiFASTA(MultiFASTA):
    """Inherits bioservices.apps.MultiFASTA and introduces additional features."""

    def read_fasta(self, filename):

        """
        Load several FASTA from a filename.
        Overrides base method. Temporary fix for bioservices MultiFASTA read_fasta funtion
        that recognises only particular headers such as Swissprot(sp) but not Trembl (tr) part of Uniprot...
        """
        fh = open(filename, "r")
        data = fh.read()
        fh.close()

        # we split according to ">2 character
        for thisfasta in data.split(">")[1:]:
            f = FASTA()
            f._fasta = f._interpret(thisfasta)

            # temporary fix for header recognition
            if f.accession == None :
                if thisfasta[0:3] == 'tr|':
                    thisfasta = 'sp|' + thisfasta[3:]
                if '|' not in thisfasta:
                    temp = thisfasta.split('\n')
                    newheader = 'sp|' + temp[0].split()[0] + '|blabla\n'
                    seq = ''
                    for l in range( 1, len(temp) ):
                        seq = seq + temp[l]
                    thisfasta = newheader + seq
                f._fasta = f._interpret(thisfasta)


            if f.accession != None and f.accession not in self.ids:
                self._fasta[f.accession] = f
            else:
                print("Accession %s is already in the ids list or could not be interpreted. skipped" % str(
                    f.accession))

    def write_fasta(self, filename):
        with open(filename, 'w') as outFastaFile:
            for n in range(len(self.df.Accession)):
                name = self.df.Accession[n]
                seq = re.sub(r"[\n\t\s]*", "", self.df.Sequence[n])
                print('>', name, sep='', file=outFastaFile)
                print(seq, file=outFastaFile)


