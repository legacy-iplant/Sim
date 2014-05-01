import simuPOP as sim
import random
from simuPOP.utils import saveCSV

pop = sim.Population(size=[6000], loci=1000, infoFields=["qtrait"])


def qtrait(geno):
    trait = random.normalvariate(sum(geno)*5, random.uniform(0.0001, 3))
    if trait <= 0:
        trait = random.uniform(0.0001, 1)
    return trait

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(prop=[0.7, 0.3])
    ],
    matingScheme=sim.RandomMating(),
    postOps=[sim.PyQuanTrait(loci=(0, 1, 2, 3, 4, 5, 6, 10, 100), func=qtrait, infoFields=["qtrait"])],
    gen=10
)

geno = list()
for i in pop.individuals():
    geno.append(i.genotype())

pheno = list()
for i in pop.individuals():
    pheno.append(i.qtrait)

f = open("qtrait1.txt", "w")
f.write("\n".join(map(lambda x: str(x), pheno)))
f.close()

saveCSV(pop,"sample.csv")