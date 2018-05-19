from metsim_core import *

verbosity(True)

# Initialize pools

pyruvate = metabolite_pool('pyruvate', 3, 5) # will be replenished infinitely with labeled carbons
pyruvate.initialize_pool()
#
#citrate = metabolite_pool('citrate', 6, 50)
#citrate.initialize_pool()
#
#glutamate = metabolite_pool('glutamate', 5, 50)
#glutamate.initialize_pool()
#
#succinate = metabolite_pool('succinate', 4, 50)
#succinate.initialize_pool()
#
#oxaloacetate = metabolite_pool('oxaloacetate', 4, 50)
#oxaloacetate.initialize_pool()
#
#
## Sequence
#
## Start with labeled pyruvate pool
#pyruvate.to_tmp(5)
#pyruvate.introduce_molecules(5, {'C1': 1, 'C2': 1, 'C3': 1})
#
#for i in range(100):
#    # Move to tmp
#    pyruvate.to_tmp(5)
#    citrate.to_tmp(5)
#    glutamate.to_tmp(5)
#    succinate.to_tmp(5)
#    oxaloacetate.to_tmp(5)
#
#    # Move from tmp
#    pyruvate.introduce_molecules(5, {'C1': 1, 'C2': 1, 'C3': 1})
#    oxaloacetate_to_citrate(5)
#    citrate_to_glutamate(5)
#    glutamate_to_succinate(5)
#    oxaloacetate.from_tmp(succinate.tmp)
#
#    # Calculate enrichment
#    citrate.calculate_enrichment()
#    glutamate.calculate_enrichment()
#    oxaloacetate.calculate_enrichment()
