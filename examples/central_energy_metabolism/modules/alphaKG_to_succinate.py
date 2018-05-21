from modules.metsim_core import *

from modules.global_vars import *



def alphaKG_to_succinate(number_of_molecules, source_alphaKG):
    if verbose:
        print('running alphaKG_to_succinate')

    tmp = source_alphaKG.tmp.sample(n = number_of_molecules)
    source_alphaKG.tmp = source_alphaKG.tmp.loc[~source_alphaKG.tmp.index.isin(tmp.index)]

    alphaKG_to_succinate.tmp = tmp.drop(['C5'], axis = 1)

    if verbose:
        print(alphaKG_to_succinate.tmp)
        print('')
