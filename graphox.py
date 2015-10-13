import numpy as np
import itertools 

def grouper(n, iterable,fillvalue=None):
	args = [iter(iterable)]*n
	return itertools.izip_longest(fillvalue=fillvalue, *args)

def decorate_plane(ox_sites, oxys):
# Generates iterator of all possible
# combinations of lattice vector and
# Oxygen position.
  cartproduct = itertools.product(ox_sites,oxys)
  sites = [a+b for a,b in cartproduct]
#  itertools.permutations(ox_sites,num_ox)    
  print map(list,itertools.permutations(ox_sites,2))

def test_valid(ox_sites):
	boole = True
# Generates all combinations of oxygen sites
# and tests if they are far enough apart
	for pair in itertools.combinations(ox_sites,2):
		boole= boole and dist(pair[0], pair[1])
	if boole:
		print 'True'
	elif not boole:
		print 'False'

def dist(a,b):
	dist_param = 1.3
#	print a,b
#	print np.sqrt((a-b).dot(a-b))
	if(np.sqrt((a-b).dot(a-b)) > dist_param):
		bool = True
	else:
		bool = False
	return bool

def istrue(a,b):
	bool = True
	for ox1 in a:
		for ox2 in b:
			if(b and dist(ox1,ox2)): 
				bool = True
			else: 
				bool = False
	return b

# size of supercell
n = 2
m = 2
a0 = 1.42

# lattice vectors
a = np.array([a0*np.sqrt(3), 0])
b = np.array([a0*np.sqrt(3)/2.0, a0*3.0/2.0])

# carbon unit cell
c1 =  np.array( [0.0, 0.0] )
c2 =  np.array( [a0*np.sqrt(3)/2.0, a0/2.0] )

# Generate carbon lattice:
for i,j in itertools.product(range(n), range(m)):
	ca = c1 + float(i)*a + float(j)*b
	cb = c2 + float(i)*a + float(j)*b
	print 'C {0:2.5f} {1:2.5f} {2:2.5f} '.format(ca[0], ca[1], 0.0)
	print 'C {0:2.5f} {1:2.5f} {2:2.5f} '.format(cb[0], cb[1], 0.0)

#generate oxygen interstices first cartesian 
#product of all possible locations.
#Oxygen sites
O1 = np.array([ a0*np.sqrt(3.0)/4.0, a0/4.0]) 
O2 = np.array([-a0*np.sqrt(3.0)/4.0, a0/4.0]) 
O3 = np.array([0.0, -a0/2]) 

num_ox = 2
ox_sites = [O1, O2, O3]
cart_prod=list(itertools.product(range(n), range(m),repeat=1))
oxys = [np.array([0.0, 0.0]) for n in range(num_ox)]

for ncell in range(n):
	for nox in range(num_ox):
		for nox2 in range(nox):
			oxys[nox2] =  oxys[nox2]  
			for mcell in range(m):
				oxys[nox] =  oxys[nox]  
				for i,j in cart_prod:
					oxys[-1] = float(i)*a + float(j)*b
					decorate_plane(oxys, ox_sites)
