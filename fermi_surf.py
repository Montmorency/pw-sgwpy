import sys
import re
import numpy as np
#f is a string with the contents of the file
class FermiSurface(object):
  def __init__(self):
    self.nx = 16
    self.ny = 16
    self.nz = 16
    self.dimvec = np.array([float(self.nx), float(self.ny), float(self.nz)])
    self.fermixyz = {}
    self.muk = {}
    self.prefix = 'mgb2'
    self.nbndmin = 1
    self.nbndmax = 4

  def __repr__(self):
    return 'Fermi Surface/Mu Tensor Object'

#should have the in pwstruct module, and also read the crystal,
#and reciprocal space coordinates direct from the pw output.
  def cryst_to_cart(self, kvec):
#pb crystal axes
    at1 = np.array([-0.500000,   0.000000,  0.500000])
    at2 = np.array([ 0.000000,   0.500000,  0.500000])
    at3 = np.array([-0.500000,   0.500000,  0.000000])
#mgb2 crystal axes
    at1 = np.array([ 1.000000,   0.000000,   0.000000])                          
    at2 = np.array([-0.500000,   0.866025,   0.000000])                           
    at3 = np.array([ 0.000000,   0.000000,   1.142069])                             
    at  = np.array([at1, at2, at3])
    outvec = np.dot(at, kvec)
    return outvec

  def pull_fermi(self,f):
    #fermi_regex  = re.compile(r'k\s=\s([0-9\.\-]+)\s([0-9\.\-]+)\s([0-9\.\-]+).*?:\n\n\s+([0-9\.\-\s]+)')
    fermi_regex  = re.compile(r'k\s=\s?(\-?[0-9\.]+)\s?(\-?[0-9\.]+)\s?(\-?[0-9\.]+).*?:\n\n\s+([0-9\.\-\s]+)')
    print len(fermi_regex.findall(f))
    for a, b, c, d in fermi_regex.findall(f):
      a = float(a)
      b = float(b)
      c = float(c)
      kvec = np.array([a,b,c])
      d = map(float, d.split())
#turn kpoint coordinates into integer indices
      kvec = self.cryst_to_cart(kvec)
      for i, a in enumerate(kvec):
        if (a<0.0): kvec[i] = kvec[i] + 1.0
      index = [round(a) for a in np.multiply(kvec, self.dimvec)]
      for i, a in enumerate(index):
        if index[i] == 16.0: 
          index[i] = 0.0
    
      self.fermixyz[tuple(index)] = d
#    for a in sorted(self.fermixyz.keys()):
#      print a 

#returns dictionary keys:xyz coordinates, values:eigenvalues.
  def print_xsf(self, surf, title='colour', band1=2, band2=5):
    for ibnd in range(band1, band2):
      #f1 = open('{0}.col.band{1}.xsf'.format(self.prefix, ibnd), 'w')
      f1 = open('{0}.{1}.{2}.xsf'.format(self.prefix, title, ibnd), 'w')
      print >>f1, "BEGIN_BLOCK_DATAGRID_3D" 
      print >>f1, "{0}_band_{1}".format(self.prefix, ibnd)     
      print >>f1, " BEGIN_DATAGRID_3D_{0}".format(self.prefix) 
      print >>f1, " {0}  {1}  {2} ".format(self.nx, self.ny, self.nz)
      #CaC6:
      print >>f1, "0.000000  0.000000  0.000000"   
     # print >>f1, "-1.000000 -1.000000  1.000000"
     # print >>f1, " 1.000000  1.000000  1.000000"
     # print >>f1, "-1.000000  1.000000 -1.000000"
     # print >>f1, ""
    #MgB2:
      print >>f1, "1.000000  0.577350  0.000000"
      print >>f1, "0.000000  1.154701  0.000000"
      print >>f1, "0.000000  0.000000  0.875604"
      print >>f1, ""
    #Pb
      for z in range(self.nz):
        for y in range(self.ny):
          for x in range(self.nx):
            try:
              print>>f1, surf[x,y,z][ibnd], " ",
            except TypeError:
              print>>f1, surf[x,y,z], " ",
              #print>>f1, 1.0, " ",
            except KeyError:
              print 'Missing key' 
          print >> f1, ""
        print >> f1, ""
      print >>f1, "END_DATAGRID_3D"  
      print >>f1, "END_BLOCK_DATAGRID_3D"  
      f1.close()

  def pull_muk(self, f):
#division along x,y,z
    f2 = open('matel','w')
    f1 = open('{0}.muk.colmap'.format(self.prefix), 'w')
    for line in f.split('\n'):
      try:
        a,b,c,d = map(float, line.split())
       #might need to duplicate values so the colour map is right.
        if(a<0): a = a+1.0
        if(b<0): b = b+1.0
        if(c<0): c = c+1.0
        ai = a*self.nx
        bi = b*self.ny
        ci = c*self.nz
        ai = round(ai)
        bi = round(bi)
        ci = round(ci)
        self.muk[ai,bi,ci] = d
        print >> f2, ai,bi,ci,d
      except:
        print "Couldn't read the following line:"
        print line

    for z in range(self.nx):
      for y in range(self.ny):
        for x in range(self.nz):
          try:
            print>>f1, self.muk[x,y,z], " ",
          except KeyError:
            print 'key error'
            print x,y,z
            pass
        print >> f1, ""
      print >> f1, ""

if __name__=="__main__":
  f = open(sys.argv[1]).read()
  fs = FermiSurface()
#  fs.pull_fermi(f)
#  fs.print_xsf(fs.fermixyz)
  fs.pull_muk(f)
  fs.print_xsf(fs.muk, 'muk')
