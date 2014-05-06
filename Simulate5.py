import getopt
import random
import sys
import math
from scipy.stats.stats import pearsonr
import numpy
import simuPOP as sim
from simuPOP.utils import saveCSV


def drange(start, stop, step):
	r = start
	while r < stop:
		yield r
		r += step


def usage():
	#print "\n"
	print "-h or --help for help"
	print "-v or --verbose for verbose"
	print "-s or --size to specify population size"
	print "-n or --number to specify number of loci"
	print "-l or --loci to specify loci with effects (separated by commas)"
	print "-e or --effect to specify corresponding loci effects (separated by commas)"
	print "-f or --filename for naming output in CSV format"
	print "-i or --herit for specifying heritability (between 0 and 1)"


def additive_model(geno):
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
	my_trait = my_total_sum
	return my_trait


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
	has_heritability = False
	filename = "my"
	heritability = 0.2
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
	           postOps=[sim.PyQuanTrait(loci=loci, func=additive_model, infoFields=["qtrait"])],
	           gen=5)

	if verbose:
		print "Coalescent process complete. Population evolved."

	genotypes = list()
	for i in pop.individuals():
		genotypes.append(i.genotype())
		#print i.genotype()

	phenotypes = list()
	for i in pop.individuals():
		phenotypes.append(i.qtrait)
		#print i.qtrait

	def fun(sigma, h):
		x_exact = phenotypes
		x_random = list()
		for each in phenotypes:
			x_random.append(random.normalvariate(each, sigma))
		r = pearsonr(x_exact, x_random)[0]
		return r - math.sqrt(h)

	#print fun(2.25, 0.25)

	if verbose:
		print "Building polynomial model for variance tuning..."

	points = list()
	for i in drange(0, max(effects)*10, 0.001):
		points.append(i)
	y_points = list()
	for i in points:
		y_points.append(fun(i, heritability))
	z = numpy.polyfit(x=points, y=y_points, deg=3)
	p = numpy.poly1d(z)

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
			print "Estimated variance of phenotypes using Newton's method for specified heritability: ", xn
		return xn

	if verbose:
		print "Using Newton's method to find polynomial roots..."

	estimated_variance = newton(p)
	new_phenotypes = list()
	for each in phenotypes:
		new_phenotypes.append(random.normalvariate(each, estimated_variance))

	f = open(filename + "_qtrait.txt", "w")
	f.write("\n".join(map(lambda x: str(x), new_phenotypes)))
	f.close()

	saveCSV(pop, filename + "_genomes.csv")
	print "\n\n"


if __name__ == "__main__":
	main()