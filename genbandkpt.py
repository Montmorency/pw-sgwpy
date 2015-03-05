import numpy as np

def gen_kline(lin, npoints, offset=np.array([0.0,0.0,0.0])):
  for i in range(npoints+1):
    tmp = offset + (lin/float(npoints))*float(i)
    print '{0:1.5f}  {1:1.5f}  {2:1.5f}  1.000'.format(tmp[0], tmp[1], tmp[2])

#High symmetry points For MgB2
#A   = np.array([0.0, 0.0, 0.5])
#Gam = np.array([0.0, 0.0, 0.0])
#K   = np.array([0.33333, 0.33333, 0.00000])
#H   = np.array([0.33333, 0.33333, 0.50000])
#M   = np.array([0.5, 0.0000, 0.0000])
#L   = np.array([0.5, 0.000, 0.5])
#Z   = np.array([0.000, 0.000, 0.5])

#High symmetry points For CaC6
T     = np.array([0.5, 0.5, -0.5])
Gam   = np.array([0.0, 0.0, 0.0])
L     = np.array([0.0, 0.5, 0.0])
Lambda = np.array([0.0, 0.0, 3.0])

#Reciprocal Lattice for MgB2
#b1 = np.array([  1.000000,  0.577350,  0.000000 ])
#b2 = np.array([  0.000000,  1.154701,  0.000000 ])
#b3 = np.array([  0.000000,  0.000000,  0.875604 ])

#Reciprocal Lattice for CaC6
b1 = np.array( [ 1.178111, -0.680183,  0.382401] )
b2 = np.array( [ 0.000000,  1.360365,  0.382401] )
b3 = np.array( [-1.178111, -0.680183,  0.382401] )

crysttocart = np.array([b1, b2, b3])
crysttocart = crysttocart.T

lin= T-Gam
norm=np.linalg.norm(lin)
npoints = 30
lin = L
a   = np.linalg.norm(lin)
normpts = int(npoints*(a/norm))
gen_kline(lin, normpts)

lin = Lambda
normpts = int(npoints)
gen_kline(lin, normpts)

lin = np.dot(crysttocart,T)
normpts = int(npoints)
gen_kline(lin, normpts, np.dot(crysttocart,T))

