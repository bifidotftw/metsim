from modules.metsim_core import *

from modules.global_vars import *
#import modules.global_vars

def serine_to_glycine(number_of_molecules, source_serine):
    if verbose:
        print('running serine_to_glycine')

    # Randomly choose rows from serine.tmp and remove them from serine.tmp
    serine_to_glycine.tmp_to_split = source_serine.tmp.sample(n = number_of_molecules)
    source_serine.tmp = source_serine.tmp.loc[~source_serine.tmp.index.isin(serine_to_glycine.tmp_to_split.index)]

    serine_to_glycine.tmp = serine_to_glycine.tmp_to_split.drop(['C3'], axis = 1)

    if verbose:
        print(serine_to_glycine.tmp)
        print('')

def serine_to_meTHF():
    if verbose:
        print('running serine_to_meTHF')

    serine_to_meTHF.tmp = serine_to_glycine.tmp_to_split.drop(['C1', 'C2'], axis = 1)
    serine_to_meTHF.tmp = serine_to_meTHF.tmp.rename(columns = {'C3': 'C1'})

    if verbose:
        print(serine_to_meTHF.tmp)
        print('')
