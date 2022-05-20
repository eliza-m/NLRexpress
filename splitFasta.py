from src.FeaturesData import *
import click

@click.group(chain=True, invoke_without_command=True)
def cli():
    pass


@click.command()
@click.option('--input', required=True, help='Input FASTA file ')
@click.option('--outdir', required=True, help='Output folder')
@click.option('--batchsize', required=False, default=1000, help='Number of sequences to be split per file. [default: 1000]')

def splitFastaFile( input:Path, outdir:Path, batchsize:int)  :

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
    except :
        print("Input FASTA file error")
        raise


    count = 1
    part = 1

    outputFile = open( str(outdir) + '/' + Path(input).stem + '_part_' +  str(part) + '.fasta', 'w')

    for name in seqData:
        if count <= batchsize:
            print('>', name, sep='', file=outputFile)
            print(seqData[name], file=outputFile)
            count += 1
        else:
            outputFile.close()
            part += 1
            count = 1
            outputFile = open(str(outdir) + '/' + Path(input).stem + '_part_' + str(part) + '.fasta', 'w')

            print('>', name, sep='', file=outputFile)
            print(seqData[name], file=outputFile)



if __name__ == '__main__':

    splitFastaFile()
