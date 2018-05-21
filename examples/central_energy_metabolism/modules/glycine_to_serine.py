from modules.metsim_core import *

from modules.global_vars import *



def glycine_to_serine(number_of_molecules, source_glycine, source_meTHF):
    if verbose:
        print('running glycine_to_serine')
    
    # Randomly choose rows from glycine.tmp and remove them from glycine.tmp
    tmp1 = source_glycine.tmp.sample(n = number_of_molecules)
    source_glycine.tmp = source_glycine.tmp.loc[~source_glycine.tmp.index.isin(tmp1.index)]

    # Randomly choose rows from meTHF and remove them from meTHF.tmp
    tmp2 = source_meTHF.tmp.sample(n = number_of_molecules)
    source_meTHF.tmp = source_meTHF.tmp.loc[~source_meTHF.tmp.index.isin(tmp2.index)]

    # meTHF C1 -> C3
    tmp2 = tmp2.rename(columns = {'C1': 'C3'})

    # Join
    tmp1 = tmp1.reset_index(drop = True)
    tmp2 = tmp2.reset_index(drop = True)
    glycine_to_serine.tmp = tmp1.join(tmp2)

    if verbose:
        print(glycine_to_serine.tmp)
        print('')
