from modules.metsim_core import *

from modules.global_vars import *



def glycine_to_serine(number_of_molecules, source_glycine, source_meTHF):
    if verbose:
        print('running glycine_to_serine')
    
    # Randomly choose rows from glycine.tmp and remove them from glycine.tmp
    tmp_glycine = source_glycine.tmp.sample(n = number_of_molecules)
    source_glycine.tmp = source_glycine.tmp.loc[~source_glycine.tmp.index.isin(tmp_glycine.index)]

    # Randomly choose rows from meTHF and remove them from meTHF.tmp
    tmp_meTHF = source_meTHF.tmp.sample(n = number_of_molecules)
    source_meTHF.tmp = source_meTHF.tmp.loc[~source_meTHF.tmp.index.isin(tmp_meTHF.index)]

    # meTHF C1 -> C3
    tmp_meTHF = tmp_meTHF.rename(columns = {'C1': 'C3'})

    # Join
    tmp_glycine = tmp_glycine.reset_index(drop = True)
    tmp_meTHF = tmp_meTHF.reset_index(drop = True)
    glycine_to_serine.tmp = tmp_glycine.join(tmp_meTHF)

    if verbose:
        print(glycine_to_serine.tmp)
        print('')
