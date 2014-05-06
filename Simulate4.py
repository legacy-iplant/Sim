import getopt
import random
import sys
import math
from scipy.stats.stats import pearsonr
import numpy
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
	print "-i or --herit for specifying heritability on a scale of 0 to 1"


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
	my_trait = random.normalvariate(my_total_sum, estimated_sigma2)
	return my_trait


def exact_trait(genotype):
	true_count = 0
	snp_count = 0
	new = list()
	for each in genotype:
		if snp_count in loci:
			new.append(each * effects[loci.index(snp_count)])
		if snp_count not in loci:
			new.append(0.0)
		true_count += 1
		if true_count % 2 is 0:
			snp_count += 1
	return sum(new)


def prob_trait(genotype, sigma2):
	true_count = 0
	snp_count = 0
	new = list()
	for each in genotype:
		if snp_count in loci:
			new.append(each * effects[loci.index(snp_count)])
		if snp_count not in loci:
			new.append(0.0)
		true_count += 1
		if true_count % 2 is 0:
			snp_count += 1
	return random.normalvariate(sum(new), math.sqrt(sigma2))


def main():
	## Check for arguments passed
	try:
		opts, args = getopt.getopt(sys.argv[1:], shortopts="vhs:n:l:e:f:i:", longopts=["verbose", "help", "size=",
		                                                                             "number=", "loci=", "effect=",
		                                                                             "filename=", "herit="])
	except getopt.GetoptError as err:
		print(err)
		usage()
		sys.exit()

	verbose = False
	has_filename = False
	print "\n"

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
		elif o[0] in ("-n", "--number"):
			number = o[1]
			if verbose:
				print "Number of loci per individual is set at", number
		elif o[0] in ("-l", "--loci"):
			global loci
			loci = o[1].split(",")
			loci = map(int, loci)
			if verbose:
				print "Loci positions per individual are:", loci
		elif o[0] in ("-e", "--effect"):
			global effects
			effects = o[1].split(",")
			effects = map(float, effects)
			if verbose:
				print "Effects for loci per individual are:", effects
		elif o[0] in ("-f", "--filename"):
			filename = o[1]
			has_filename = True
			if verbose:
				print "File will be saved as:", filename
		elif o[0] in ("-i", "--herit"):
			heritability = float(o[1])
			has_heritability = True
			if verbose:
				print "Heritability for simulation specified as:", heritability

	## Start quantitative trait simulation
	if verbose:
		print "Creating population..."

	pop = sim.Population(size=int(individuals), loci=int(number), infoFields=["qtrait"])

	if verbose:
		print "Evolving population..."

	pop.evolve(initOps=[sim.InitSex(), sim.InitGenotype(prop=[0.7, 0.3])], matingScheme=sim.RandomMating(),
	           #postOps=[sim.PyQuanTrait(loci=loci, func=trait, infoFields=["qtrait"])],
	           gen=5)

	if verbose:
		print "Population evolved."

	genotypes = list()
	for i in pop.individuals():
		genotypes.append(i.genotype())
		#print i.genotype()

	def fun(sigma2, h):
		exact_traits = list()
		for i in genotypes:
			exact_traits.append(exact_trait(i))
		prob_traits = list()
		for i in genotypes:
			prob_traits.append(prob_trait(i, sigma2))
		r = pearsonr(exact_traits, prob_traits)[0]
		return r - math.sqrt(h)

	def newton(p):
		xn = 100
		p_d = p.deriv()
		count = 0
		while abs(p(xn)) > 0.01:
			if count > 1000:
				print "Unable to converge after 1000 iterations..."
				sys.exit()
			count += 1
			xn = xn - p(xn)/p_d(xn)
		if verbose:
			print "Estimated variance using Newton's method for solving roots is: ", xn
		return xn

	## Make sure to change to 100 points around the average "effects"
	my_points = list()
	for i in range(100):
		my_points.append(fun(i, heritability))
	z = numpy.polyfit(x=my_points, y=range(100), deg=2)
	#print z
	p = numpy.poly1d(z)
	if verbose:
		print "Polynomial fit for finding tuning parameter to match heritability: \n", p
	estimated_variance = newton(p)
	#print estimated_variance

	phenotypes = list()
	for i in pop.individuals():
		phenotypes.append(prob_trait(i.genotype(), estimated_variance))
	#print phenotypes

	if has_filename is False:
		filename = "my"

	f = open(filename + "_qtrait.txt", "w")
	f.write("\n".join(map(lambda x: str(x), phenotypes)))
	f.close()

	saveCSV(pop, filename + "_genomes.csv")
	print "\n\n"


if __name__ == "__main__":
	main()