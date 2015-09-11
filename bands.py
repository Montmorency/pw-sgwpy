import sys
import re

"""
pass Quantum Espresso output file and prefix = name for output files 
as arguments on file line e.g.

python qe_bands.py --prefix si ./si.bands.out

file parses QE output for kpoints and eigenvalues and prints
a gnuplot script for plotting. 
"""

def parse_args(args):
  extra = []
  vars  = []
  current_var = None

  for arg in args:
    if arg.startswith('--'):
      current_var = arg[2:]
    else:
      if current_var is not None:
        vars.append((current_var, arg))
        current_var = None
      else:
        extra.append(arg)
  return (extra, vars)


class BandStruct(object):
  def __init__(self,
    fermilevel =   0.00,
    ylower     = -15.0,
    yupper     =   6.0,
    prefix     = 'si',
    title      = 'si.bands.out',
    qe_output  = 'si.scf.out',
    kpoints    = None,
    es         = None):

    self.fermilevel = fermilevel
    self.ylower     = ylower
    self.yupper     = yupper
    self.prefix     = prefix[0]
    self.title      = title
    self.qe_output  = qe_output

    if kpoints is None:
      self.kpoints = [] 
    else:
      self.kpoints = kpoints

    if es is None:
      self.es  = []
    else:
      self.es  = es

#extracts all bands, kpoints, and the fermi energy.
  def extract_bands(self, f):
    f = f.read()
    g = open('{0}.bands'.format(self.prefix),'w')
#  f is a string with the contents of the file
#  re.S match new lines as white space. re.M multiline...
    kpoint_regex = re.compile(r'\s+k =[\s-]([0-9\.]+)[\s-]([0-9\.]+)[\s-]([0-9\.]+).*?bands \(ev\):\n(.*?)\n\n', re.S)
    fermi_regex =  re.compile(r'Fermi energy is\s+([0-9\.]+) ev')
    self.fermilevel = float(fermi_regex.findall(f)[0])
    for a,b,c,es in kpoint_regex.findall(f):
      a = float(a)
      b = float(b)
      c = float(c)
      es = map(float, es.split())
      self.kpoints.append((a,b,c))
      self.es.append(es)

      print >> g, " ".join(map(str,es))

    g.close()

  def gnu_plot_bands(self):
    f = open('{0}.bands.gnu'.format(self.prefix),'w')
#default plot full ranges of kpoints and eigenvalues
    ylower = min(self.es[0])-5.0
    yupper = max(self.es[0])+5.0
    nks    = len(self.kpoints)
    nbnd   = len(self.es[0])

    print >> f, "unset key"
#print "unset xtics"
    print >> f, "set xr [0:{0}]".format(nks)
    print >> f, "set yr [{0}:{1}]".format(ylower, yupper)
    print >> f, "plot \'{0}.bands\' u ($13 - {1}) w l lt -1".format(self.prefix, self.fermilevel)

    for i in range(1,nbnd):
      print >>f, "replot \'{0}.bands\' u (${1} - {2}) w l lt -1".format(self.prefix, i, self.fermilevel)

    f.close()

def any(xs):
  for x in xs:
    if x:
      return True
  else:
    return False

def gen_dict(vars):
  vars_values = []
  for var_name, values in vars:
    values = values.split(",")
    try:
      if any(['.' in value for value in values]):
        values = map(float, values)
      else:
        values = map(int, values)
    except ValueError:
      pass
    vars_values.append((var_name, values))

  return dict(vars_values)

if __name__=='__main__':
  extra, vars = parse_args(sys.argv[1:])
  vars = gen_dict(vars)
  qe_output =  open(extra[0], 'r')
  bands = BandStruct(**vars)
  bands.extract_bands(qe_output)
  bands.gnu_plot_bands()

  print 'number of bands {0}'.format(len(bands.es[0]))
  print 'number of kpoints {0}'.format(len(bands.kpoints))
  print 'Fermi Energy {0}'.format(bands.fermilevel)


