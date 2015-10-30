import numpy as np
import itertools 
import os
import argparse
from scipy.misc import comb
from pprint import pprint

# This function should return a list of lists
# For the simplest case we only consider first option at
# each level of the tree and two oxygen atoms.
# We should find the function returns
# a list nested as such: 
# [site1, [valid_sites]]
# The next step would be three oxygen 
# atoms (again only considering first site at 
# each level of the tree)
# at this level we should find a nested list:
#   O1,     O2,       O3
# [site1, [site2, [valid_sites]] ]
# Now let's consider the case where we loop 
# over all the possible configurations at each level
# How would the list look then?
# For two oxygen's would look like this?
# [[site1, [valid_sites]], [site2, [valid_sites]] ]
# Now for the case of three oxygen atoms?
# [[site1, [site11, [valid_sites]]], [site1, [site2,[valid_sites]]], ...]

class GraphOx(object):
    def __init__(self, configs=[]):
        self.configs = []

    def backtrack_recurse(self, ox_sites, config, num_ox):
        n = 0
        if num_ox > 2:
            config = config + [ox_sites[0]]
            for site in ox_sites[1]:
                if len(site) > 2 :
                    config += [site[0]]
                    self.backtrack_recurse(site[1], config, num_ox-1)
                    config.pop()
                else : 
                    self.backtrack_recurse(site, config, num_ox-1)
        elif num_ox == 2:
            n = 0
            config = config + [ox_sites[0]]
            for site in ox_sites[1]:
                config = config + [site]
                #should sort list here:
                #print 'config', config
                a = sorted(config, key=lambda x : (x[0],x[1])) 
                a = (tuple(vec) for vec in a)
                #print 'a', a
                self.configs.append(a)
                config.pop()
        return

    def prune_configs(self): 
        #self.configs = set(self.configs)
        self.configs = set(map(tuple, self.configs))
        return

class node(object):
  def __init__(self, value, children=[]):
  	self.value = value
  	self.children = children

def brute_force(ox_sites, num_ox):
  for config in itertools.combinations(ox_sites,num_ox):
  	if test_valid(config):
  		valid_config.append(list(config))
  	else:
  		'not a valid config.'
  return valid_config

def expand_sites_mut(sites):
    for site in sites:
        while len(site[-1]) > 1:
            site.extend(site.pop())
        site.extend(site.pop())

def backtrack(ox_sites, num_ox):
    k = 0
    sites = ox_sites
    solt  = []
    while k < (num_ox-1):
        k = k+1
        solt.append(sites[0])
        sites = sites[1]
    for site in sites:
        yield solt + [site]

def iterate(ox_sites,num_ox):
    for ox_site in ox_sites:
        if (dist(config, ox_site) == False): 
            pass
        else:
            ox_sites_prune.append(ox_site)
        yield recursive_decorate(ox_sites_prune, num_ox-1)

def recursive_decorate(ox_sites, num_ox):
  try:
    config  = ox_sites[0]
  except IndexError:
      print 'Not sure that many Oxygens will fit.'

  configs = [] 
  ox_sites_prune = []
  if (num_ox > 1): 
  	for ox_site in ox_sites:
  		if (dist(config, ox_site) == False): 
  			pass
  		else:
  			ox_sites_prune.append(ox_site)
  	configs.extend([config, recursive_decorate(ox_sites_prune, num_ox-1)])
  elif (num_ox ==1):
  	return ox_sites
  return configs

def recursive_decorate_stack(ox_sites, num_ox):
  n = 0
  configs = [] 
  while (n < len(ox_sites)): 
        ox_sites_prune = []
        config  = ox_sites[n]
        if (num_ox > 1) :
  	    for ox_site in ox_sites:
  	    	if (dist(config, ox_site) == False): 
  	        	pass
  		else:
  			ox_sites_prune.append(ox_site)
            configs.append([config, recursive_decorate_stack(ox_sites_prune, num_ox-1)])
        elif (num_ox ==1):
  	    return ox_sites
        n += 1
  return configs

def backtrack_recurse(ox_sites, config, num_ox):
    n = 0
    if num_ox > 2:
        config = config + [ox_sites[0]]
        for site in ox_sites[1]:
            if len(site) > 2 :
                config += [site[0]]
                backtrack_recurse(site[1], config, num_ox-1)
                config.pop()
            else : 
                backtrack_recurse(site, config, num_ox-1)
    elif num_ox == 2:
        n = 0
        config = config + [ox_sites[0]]
        for site in ox_sites[1]:
            config = config + [site]
    #        print config
            config.pop()
    return

def expand_list(ox_list, n, num_ox):
    configs = ox_list[0]
    if n < (num_ox-1):
        configs.append(configs, expand_list(ox_list[n], n+1, num_ox))
    else:
        return list(itertools.product(*ox_list))

#from black cat
def expand_sites(ox_list):
    return [getsubsite_rec(site) for site in ox_list]
def getsubsite_rec(site):
    if len(site) == 1:
        return site
    else:
        return site[:1] + getsubsite_rec(site[1])
def getsubsite(site):
    site = site[:] #copy the list
    while len(site[-1]) > 1:
        site.extend(site.pop())
    site.extend(site.pop())
    return site

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


parser = argparse.ArgumentParser(description='Pattern Hexagonal Flake with Oxygen')
parser.add_argument('--cell', metavar='m,n', type=int, nargs='+', help='size of supercell')
parser.add_argument('--numox', metavar='numox', type=int, help='number of oxygen atoms')
args = parser.parse_args()
n = args.cell[0] 
m = args.cell[1]
a0 = 1.42
num_ox = args.numox

# lattice vectors
a = np.array([a0*np.sqrt(3), 0])
b = np.array([a0*np.sqrt(3)/2.0, a0*3.0/2.0])

# carbon unit cell
c1 =  np.array( [0.0, 0.0] )
c2 =  np.array( [a0*np.sqrt(3)/2.0, a0/2.0] )

#oxygen lattice vectors
#o1 = np.array([ a0*np.sqrt(3.0)/4.0, a0/4.0]) 
#o2 = np.array([-a0*np.sqrt(3.0)/4.0, a0/4.0]) 
#o3 = np.array([0.0, -a0/2]) 
#These ones should always fall inside flake
o1 = np.array([ a0*np.sqrt(3.0)/4.0, a0/4.0]) 
o2 = np.array([a0*np.sqrt(3.0)/2.0, a0]) 
o3 = np.array([3*a0*np.sqrt(3.0)/4.0,a0/4]) 

print ''
print ''
print ''
print '\t PATTERNING GRAPHENE FLAKE WITH OXYGEN '
if args.numox > 1:
    print '\t Number of oxygen atoms:', args.numox, ' Supercell:', m, n
elif args.numox == 1: 
    print '\t {0} oxygen atom', args.numox, 'Supercell', m, n
else:
    print "Number of oxygen atoms must be greater than or equal to 1."
print ''
print ''
print ''

try:
  os.mkdir('structs')
except:
  print '\t dir structs already exists'

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

print '\t Roughly {0} possible choices of {1} oxygens on a lattice with {2} sites'.format(comb(len(ox_sites),num_ox), num_ox, len(ox_sites))
print '\t {0} combinations of O-O pairs'.format(comb(num_ox, 2))

#pprint(ox_sites)
valid_configs = recursive_decorate_stack(ox_sites, num_ox)
#pprint(valid_configs[0])
configs =[]
#results = backtrack_recurse(valid_configs[0], configs, num_ox )
configs = []
flake = GraphOx()
for config in valid_configs:
    flake.backtrack_recurse(config, configs, num_ox )

print ''
print '\t There are {0} valid configurations'.format(len(flake.configs))
print ''
flake.prune_configs()
print ''
print '\t There are {0} unique configurations'.format(len(flake.configs))
print ''

for i, config in enumerate(flake.configs):
    f = open('./structs/c{0}.xyz'.format(str(i)), 'w')
    print >> f, 2*n*m+num_ox
    print >> f, ''
    gen_carbon_lattice(f)
    for osite in config:
        print >>f, 'O {0:2.5} {1:2.5} {2:2.5f} '.format(osite[0], osite[1], 0.0)
    f.close()
