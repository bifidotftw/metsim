from modules.metsim_core import *

from modules.global_vars import *



def pyruvate_to_oxaloacetate(number_of_molecules, source_pyruvate):
    if verbose:
        print('running pyruvate_to_oxaloacetate')

    # Randomly choose rows from serine.tmp and remove them from serine.tmp
    tmp = source_pyruvate.tmp.sample(n = number_of_molecules)
    source_pyruvate.tmp = source_pyruvate.tmp.loc[~source_pyruvate.tmp.index.isin(tmp.index)]

    pyruvate_to_oxaloacetate.tmp = tmp.assign(C4 = number_of_molecules*0)

    if verbose:
        print(pyruvate_to_oxaloacetate.tmp)
        print('')
