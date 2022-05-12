import os
import click
import subprocess
from bioservices.apps import FASTA
from pathlib import Path

from .src.ModelData import *
from .src.ModuleData import *
from .src.FeaturesData import *

# @click.group(chain=True, invoke_without_command=True)
# def cli():
#     pass


@click.command()
@click.option('--input', required=True, help='Input FASTA file ')
@click.option('--outdir', required=True, help='Output folder')
@click.option('--module', required=True, help="""\b
                        Predifined prediction modules : 
                          - cc             :  CC motifs - EDVID
                          - nbs            :  NBS motifs - .....
                          - lrr            :  LRR motifs - LxxLxL
                          - all            :  All modules """)


def predict( input: Path, outdir: Path, module: str ):
    """Print FILENAME."""

    inputData = FeaturesData().generate_features(input, outdir)

    results = {}

    if module in ("cc", "all") :
        CCexpress = ModuleData().loadModels(
            modelsPath={ 'extendedEDVID': 'models/' })

        for model in CCexpress.models:
            results[model.name] = model.predict_proba( inputData.X )

    if module in ("nbs", "all") :
        NBSexpress = ModuleData().loadModels(
            modelsPath={'extendedEDVID': 'models/'})
        for model in NBSexpress.models:
            results[model.name] = model.predict_proba( inputData.X )


    if module in ("lrr", "all") :
        LRRexpress = ModuleData().loadModels(
            modelsPath={'LxxLxL': 'models/'})
        for model in LRRexpress.models:
            results[model.name] = model.predict_proba( inputData.X )


if __name__ == '__main__':
    predict()



