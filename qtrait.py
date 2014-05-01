import simuPOP as sim
import random

pop = sim.Population(size=5000, loci=1, infoFields=['qtrait1','age'])
pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=[40]))

def qtrait(geno, age):
    'Return two traits that depends on genotype and age'
    return random.normalvariate(age * sum(geno), 10)

pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(prop=[0.8,0.2]),
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        # use random age for simplicity
        sim.InitInfo(lambda:random.randint(20, 75), infoFields='age'),
        sim.PyQuanTrait(loci=(1), func=qtrait, infoFields=['qtrait1']),
        sim.Stat(meanOfInfo=['qtrait1'], subPops=[(0, sim.ALL_AVAIL)],
            vars='meanOfInfo_sp'),
        sim.PyEval(r"'Mean of trait1: %.3f (age < 40), %.3f (age >=40)\n' % "
            "(subPop[(0,0)]['meanOfInfo']['qtrait1'], subPop[(0,1)]['meanOfInfo']['qtrait1'])"),
    ],
    gen = 100
)

qtrait1_ls = list()
for ind in pop.individuals():
    qtrait1_ls.append(ind.qtrait1)

f = open("qtrait1.txt", "w")
f.write("\n".join(map(lambda x: str(x), qtrait1_ls)))
f.close()

saveCSV(pop,"sample.csv")