from modules.metsim_core import *

from modules.global_vars import *



def citrate_to_alphaKG(number_of_molecules, source_citrate):
    if verbose:
        print('running citrate_to_alphaKG')

    tmp = source_citrate.tmp.sample(n = number_of_molecules)
    source_citrate.tmp = source_citrate.tmp.loc[~source_citrate.tmp.index.isin(tmp.index)]
    citrate_to_alphaKG.tmp = tmp.drop(['C6'], axis = 1)

    if verbose:
        print(citrate_to_alphaKG.tmp)
        print('')
