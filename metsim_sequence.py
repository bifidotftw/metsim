from modules.metsim_core import *
from modules.serine_to_glycine_meTHF import *
from modules.glycine_to_serine import *
from modules.pyruvate_to_oxaloacetate import *
from modules.oxaloacetate_to_citrate import *
from modules.citrate_to_alphaKG import *
from modules.alphaKG_to_succinate import *


scale = 10

# Initialize pools

PG3 = metabolite_pool('3PG', 3, 10*scale) # 3-phosphoglycerate
PG3.initialize_pool()
PG3.to_tmp(40)
PG3.introduce_molecules(40, '111')


serine = metabolite_pool('serine', 3, 400*scale)
serine.initialize_pool()

meTHF = metabolite_pool('meTHF', 1, 200*scale)
meTHF.initialize_pool()

glycine = metabolite_pool('glycine', 2, 1000*scale)
glycine.initialize_pool()

pyruvate = metabolite_pool('pyruvate', 3, 10*scale)
pyruvate.initialize_pool()

lactate = metabolite_pool('lactate', 3, 500*scale)
lactate.initialize_pool()

alanine = metabolite_pool('alanine', 3, 200*scale)
alanine.initialize_pool()

citrate = metabolite_pool('citrate', 6, 10*scale)
citrate.initialize_pool()

alphaKG = metabolite_pool('alphaKG', 5, 10*scale)
alphaKG.initialize_pool()

glutamate = metabolite_pool('glutamate', 5, 1000*scale)
glutamate.initialize_pool()

glutamine = metabolite_pool('glutamine', 5, 14) # will be replenished infinitely with unlabeled carbons
glutamine.initialize_pool()

succinate = metabolite_pool('succinate', 4, 50*scale)
succinate.initialize_pool()

fumarate = metabolite_pool('fumarate', 4, 50*scale)
fumarate.initialize_pool()

malate = metabolite_pool('malate', 4, 250*scale)
malate.initialize_pool()

oxaloacetate = metabolite_pool('oxaloacetate', 4, 10*scale)
oxaloacetate.initialize_pool()

aspartate = metabolite_pool('aspartate', 4, 750*scale)
aspartate.initialize_pool()

steps = 1
for i in range(steps):
    # Move to tmp
    PG3.to_tmp(40)
    serine.to_tmp(15)
    glycine.to_tmp(10)
    meTHF.to_tmp(10)
    pyruvate.to_tmp(45)
    lactate.to_tmp(30)
    alanine.to_tmp(5)
    citrate.to_tmp(9)
    alphaKG.to_tmp(22)
    glutamate.to_tmp(23)
    glutamine.to_tmp(14)
    succinate.to_tmp(13)
    fumarate.to_tmp(26)
    malate.to_tmp(39)
    oxaloacetate.to_tmp(32)
    aspartate.to_tmp(10)


    # Move from tmp
    PG3.introduce_molecules(40, '111')
    PG3.check_pool_size()


    serine.from_tmp(10, PG3)
    glycine_to_serine(5, glycine, meTHF)
    serine.from_tmp(5, glycine_to_serine)
    serine.check_pool_size()

    serine_to_glycine(10, serine)
    glycine.from_tmp(10, serine_to_glycine)
    glycine.check_pool_size()

    serine_to_meTHF()
    meTHF.from_tmp(10, serine_to_meTHF)
    meTHF.check_pool_size()

    pyruvate.from_tmp(30, PG3)
    pyruvate.from_tmp(15, lactate)
    pyruvate.check_pool_size()

    lactate.from_tmp(30, pyruvate)
    lactate.check_pool_size

    alanine.from_tmp(5, pyruvate)
    alanine.check_pool_size()

    oxaloacetate_to_citrate(9, oxaloacetate, pyruvate)
    citrate.from_tmp(9, oxaloacetate_to_citrate)
    citrate.check_pool_size()

    citrate_to_alphaKG(4, citrate)
    alphaKG.from_tmp(4, citrate_to_alphaKG)
    alphaKG.from_tmp(18, glutamate)
    alphaKG.check_pool_size()

    glutamate.from_tmp(14, glutamine)
    glutamate.from_tmp(9, alphaKG)
    glutamate.check_pool_size()

    glutamine.introduce_molecules(14, '00000')
    glutamine.check_pool_size()

    alphaKG_to_succinate(13, alphaKG)
    succinate.from_tmp(13, alphaKG_to_succinate)
    succinate.check_pool_size()


    fumarate.from_tmp(13, succinate)
    fumarate.from_tmp(13, malate)
    fumarate.check_pool_size()

    malate.from_tmp(26, fumarate)
    malate.from_tmp(13, oxaloacetate)
    malate.check_pool_size()

    oxaloacetate.from_tmp(26, malate)
    oxaloacetate.from_tmp(5, aspartate)
    pyruvate_to_oxaloacetate(1, pyruvate)
    oxaloacetate.from_tmp(1, pyruvate_to_oxaloacetate)
    oxaloacetate.check_pool_size()

    aspartate.from_tmp(10, oxaloacetate)
    aspartate.check_pool_size()


    print('step ' + str(i+1) + '(' + str((i+1)/steps*100) + ' %)')
    serine.calculate_enrichment()
    glycine.calculate_enrichment()
    pyruvate.calculate_enrichment()
    lactate.calculate_enrichment()
    alanine.calculate_enrichment()
    citrate.calculate_enrichment()
    alphaKG.calculate_enrichment()
    glutamate.calculate_enrichment()
    succinate.calculate_enrichment()
    fumarate.calculate_enrichment()
    malate.calculate_enrichment()
    oxaloacetate.calculate_enrichment()
    aspartate.calculate_enrichment()


serine.export()
#glycine.export()
#pyruvate.export()
#lactate.export()
#alanine.export()
#citrate.export()
#alphaKG.export()
#glutamate.export()
#succinate.export()
#fumarate.export()
#malate.export()
#oxaloacetate.export()
#aspartate.export()



## Check tmp
#PG3.check_tmp_size()
#serine.check_tmp_size(5)
#glycine.check_tmp_size(5)
#meTHF.check_tmp_size(5)
#pyruvate.check_tmp_size(0)
#lactate.check_tmp_size(15)
#alanine.check_tmp_size(5)
#citrate.check_tmp_size(5)
#alphaKG.check_tmp_size(0)
#glutamate.check_tmp_size(5)
#glutamine.check_tmp_size(0)
#succinate.check_tmp_size(0)
#fumarate.check_tmp_size(0)
#malate.check_tmp_size(0)
#oxaloacetate.check_tmp_size(0)
#aspartate.check_tmp_size(5)
