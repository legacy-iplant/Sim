import getopt
import random
import sys
import simuPOP as sim
from simuPOP.utils import saveCSV

def usage():
	print "\n"
	print "-h or --help for help"
	print "-v or --verbose for verbose"
	print "-s or --size to specify population size"
	print "-n or --number to specify number of loci"
	print "-l or --loci to specify loci with effects (separated by commas)"
	print "-e or --effects to specify corresponding loci effects (separated by commas)"
	print "-f or --filename for naming output in CSV format"


def trait(geno):
	my_sum = 0
	my_total_sum = 0
	true_count = 0
	snp_count = 0
	for each in geno:
		my_sum += each * effects[snp_count]
		true_count += 1
		if true_count % 2 is 0:
			snp_count += 1
			my_total_sum += my_sum
			my_sum = 0
	my_trait = random.normalvariate(my_total_sum, 0.2)
	return my_trait


def main():
	## Check for arguments passed
	try:
		opts, args = getopt.getopt(sys.argv[1:], shortopts="vhs:n:l:e:f:", longopts=["verbose", "help", "size=",
		                                                                             "number=", "loci=", "effect=",
		                                                                             "filename="])
	except getopt.GetoptError as err:
		print(err)
		usage()
		sys.exit()

	verbose = False
	has_filename = False

	for o in opts:
		if o[0] in ("-v", "--verbose"):
			verbose = True
			print ("Verbose mode")
	for o in opts:
		if o[0] in ("-h", "--help"):
			usage()
			sys.exit()
		elif o[0] in ("-s", "--size"):
			individuals = o[1]
			if verbose:
				print "Population size is set at", individuals
		elif o[0] in ("n", "--number"):
			number = o[1]
			if verbose:
				print "Number of loci per individual is set at", number
		elif o[0] in ("-l", "--loci"):
			loci = o[1].split(",")
			loci = map(int, loci)
			if verbose:
				print "Loci positions per individual are:", loci
		elif o[0] in ("-e", "--effect"):
			global effects
			effects = o[1].split(",")
			effects = map(int, effects)
			if verbose:
				print "Effects for loci per individual are:", effects
		elif o[0] in ("-f", "--filename"):
			filename = o[1]
			has_filename = True
			if verbose:
				print "File will be saved as:", filename

	## Start quantitative trait simulation
	if verbose:
		print "Creating population..."

	pop = sim.Population(size=int(individuals), loci=int(number), infoFields=["qtrait"])

	if verbose:
		print "Evolving population..."

	pop.evolve(initOps=[sim.InitSex(), sim.InitGenotype(prop=[0.7, 0.3])], matingScheme=sim.RandomMating(),
	           postOps=[sim.PyQuanTrait(loci=loci, func=trait, infoFields=["qtrait"])],
	           gen=5)

	genotypes = list()
	for i in pop.individuals():
		genotypes.append(i.genotype())

	phenotypes = list()
	for i in pop.individuals():
		phenotypes.append(i.qtrait)

	if has_filename is False:
		filename = "my"

	f = open(filename + "_qtrait.txt", "w")
	f.write("\n".join(map(lambda x: str(x), phenotypes)))
	f.close()

	saveCSV(pop, filename + "_genomes.csv")


if __name__ == "__main__":
	main()