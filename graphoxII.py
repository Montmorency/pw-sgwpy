import numpy as np
import itertools 
import os
from scipy.misc import comb

#Class stuff(object):
#	def __init__(self, value):
#		self.label = 'c'
#		self.node = None

# USER ADJUST PARAMETERS
n = 6
m = 6
a0 = 1.42
num_ox = 2

def brute_force(ox_sites, num_ox):
	for config in itertools.combinations(ox_sites,num_ox):
		if test_valid(config):
			valid_config.append(list(config))
		else:
			'not a valid config.'
	return valid_config

# for picked oxygen atom prune the graphene flake
# recursively call the algorithm and place 
# the oxygen in a valid spot

def recursive_decorate(config, ox_sites, num_ox):
	configs=[]
	ox_sites_prune=[]
	if (num_ox > 1): 
		for ox_site in ox_sites:
			if (dist(config, ox_site) == False): 
				pass
			else:
				ox_sites_prune.append(ox_site)
		configs.extend(recursive_decorate(config, ox_sites_prune, num_ox-1))
	elif (num_ox ==1):
		return ox_sites
	return configs

def gen_carbon_lattice(f):
# lattice vectors
	a = np.array([a0*np.sqrt(3), 0])
	b = np.array([a0*np.sqrt(3)/2.0, a0*3.0/2.0])
# carbon unit cell
	c1 =  np.array( [0.0, 0.0] )
	c2 =  np.array( [a0*np.sqrt(3)/2.0, a0/2.0] )
# Generate carbon lattice:
	for i in range(n):
		for j in range(m):
			ca = c1 + float(i)*a + float(j)*b
			cb = c2 + float(i)*a + float(j)*b
			print >>f, 'C {0:2.5f} {1:2.5f} {2:2.5f} '.format(ca[0], ca[1], 0.0)
			print >>f, 'C {0:2.5f} {1:2.5f} {2:2.5f} '.format(cb[0], cb[1], 0.0)

def test_valid(ox_sites):
	boole = True
# Generates all combinations of oxygen sites
# and tests if they are far enough apart
	for pair in itertools.combinations(ox_sites,2):
		boole = boole and dist(pair[0], pair[1])
		if boole == False:
			break
	return boole

def dist(a,b):
	dist_param = 1.3
	if(np.sqrt((a-b).dot(a-b)) > dist_param):
		bool = True
	else:
		bool = False
	return bool

# lattice vectors
a = np.array([a0*np.sqrt(3), 0])
b = np.array([a0*np.sqrt(3)/2.0, a0*3.0/2.0])

# carbon unit cell
c1 =  np.array( [0.0, 0.0] )
c2 =  np.array( [a0*np.sqrt(3)/2.0, a0/2.0] )

#oxygen lattice vectors
o1 = np.array([ a0*np.sqrt(3.0)/4.0, a0/4.0]) 
o2 = np.array([-a0*np.sqrt(3.0)/4.0, a0/4.0]) 
o3 = np.array([0.0, -a0/2]) 

try:
	os.mkdir('structs')
except:
	print 'dir structs already exists'

#Generate oxygen interstitial lattice
ox_sites=[]
for i in range(0,n):
	for j in range(0,m):
		oa = o1 + float(i)*a + float(j)*b
		ob = o2 + float(i)*a + float(j)*b
		oc = o3 + float(i)*a + float(j)*b
		ox_sites.append(oa)
		ox_sites.append(ob)
		ox_sites.append(oc)

print '{0} possible choices of {1} O on {2} lattice '.format(comb(len(ox_sites),num_ox), num_ox, len(ox_sites))
print '{0} combinations of O-O pairs'.format(comb(num_ox, 2))

valid_configs = []
print len(ox_sites)
for config in ox_sites:
	print config

first_o = ox_sites
for i, config in enumerate(first_o):
	print 'site', i
	oxlist = ox_sites
	configs = recursive_decorate(config, oxlist, num_ox)
	valid_configs.append(configs)

print len(valid_configs)
print 'There are {0} valid configurations'.format(len(valid_configs))

for i, isite in enumerate(ox_sites):
	for jsite in valid_configs[i]:
		f = open('./structs/c{0}.xyz'.format(i),'w')
		print >> f, 2*n*m+num_ox
		print >> f, ''
		gen_carbon_lattice(f)
		print >>f, 'O {0:2.5} {1:2.5} {2:2.5f} '.format(isite[0], isite[1], 0.0)
		print >>f, 'O {0:2.5} {1:2.5} {2:2.5f} '.format(jsite[0], jsite[1], 0.0)
		f.close()
