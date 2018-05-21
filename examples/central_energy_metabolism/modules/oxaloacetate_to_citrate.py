from modules.metsim_core import *

from modules.global_vars import *

def oxaloacetate_to_citrate(number_of_molecules, source_oxaloacetate, source_pyruvate):
    if verbose:
        print('running oxaloacetate_to_citrate')

    tmp_oxaloacetate = source_oxaloacetate.tmp.sample(n = number_of_molecules)
    source_oxaloacetate.tmp = source_oxaloacetate.tmp.loc[~source_oxaloacetate.tmp.index.isin(tmp_oxaloacetate.index)] # TODO: check if working
    tmp_oxaloacetate = tmp_oxaloacetate.rename(columns={'C1': 'C6', 'C2': 'C3', 'C3': 'C4', 'C4': 'C5'}) # Rename carbons according to new position in citrate, side chain is assigned C6
    tmp_oxaloacetate = tmp_oxaloacetate.reset_index(drop = True)

    tmp_pyruvate = source_pyruvate.tmp.sample(n = number_of_molecules)
    source_pyruvate.tmp = source_pyruvate.tmp.loc[~source_pyruvate.tmp.index.isin(tmp_pyruvate.index)] # TODO: check if working
    tmp_pyruvate = tmp_pyruvate.drop(['C1'], axis = 1)
    tmp_pyruvate = tmp_pyruvate.rename(columns = {'C2': 'C1', 'C3': 'C2'})
    tmp_pyruvate = tmp_pyruvate.reset_index(drop = True)

    oxaloacetate_to_citrate.tmp = tmp_oxaloacetate.join(tmp_pyruvate)

    if verbose:
        print(oxaloacetate_to_citrate.tmp)
        print('')
