import sys, os

sys.path.insert(0, os.path.abspath(".."))
from src.FeaturesData import *

data = generate_features( Path("input_data/2prot.fasta"), Path("output_data"), skipJhmmer=True )

print( data.X )
