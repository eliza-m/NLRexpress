import sys, os

sys.path.insert(0, os.path.abspath(".."))
from src.FeaturesData import *


# data = generateFeatures( Path("input_data/2prot.fasta"), Path("output_data"),
#                          motifs={"extendedEDVID": {"windLeft": 5, "windRight": 5, "motifSpan": 12}},
#                          skipJhmmer=True )



print( data.X )
