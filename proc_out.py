import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import gen_files as gen



class Eps(object):
  def __init__(self):
    self.w    = [(0.0, 0.0)]
    self.eps  = [(0.0, 0.0)]
    self.eels = [(0.0, 0.0)]
    self.qvec = (0.0, 0.0, 0.0)
    self.nws  = 0

  def __repr__(self):
      return 'epsilon object for qpoint {0}'.format(self.qvec)

  def gen_eels(self):
   #unit cell volume
    omega  = 716.8086
   #nkstot = 100
    nkstot = 1
    #const  = 64.0*np.pi/(omega*nkstot) 
    const  = np.pi
    self.eels = [(const*(epsi)/(np.square(epsr) + np.square(epsi))) for epsr, epsi in self.eps]
    

def pull_eps(f):
#get all the frequencies, epsre, and qvecs out of the output file.
  freq_re    = re.compile(r'Imag. Frequencies:\s\n(.*?)nwgreen:',re.S)
  qpoints_re = re.compile(r'(?<=Calculation of q =)(.*?)\n')
  eps_re     = re.compile(r'(?<=inveps_{GG}\(q,w\) =)(.*?)\n')

  freqs     = freq_re.findall(f)[0].split('\n')
  qpoints   = qpoints_re.findall(f)
  eps       = eps_re.findall(f)

  try:
    freqs     = [map(float, (x.split())) for x in freqs[:-1]]
    qpoints   = [map(float, (x.split())) for x in qpoints]
    eps       = [map(float, (x.split())) for x in eps]
  except:
    print 'oh crap'
    pass

  epsw = Eps()
  epsw.w    = [(w[1], w[2]) for w in freqs]
  epsw.eps  = eps
  epsw.qvec = qpoints[0]
  epsw.nws  = len(epsw.w)

  return epsw

def plot_eps(eps):
    epsre, epsim = zip(*eps.eps)
    eps.gen_eels()
    wre, wim = zip(*eps.w)
    plt.plot(wre[1:eps.nws-1], epsre[1:eps.nws-1], 'black')
    plt.plot(wre[1:eps.nws-1], epsim[1:eps.nws-1], 'blue')
    plt.plot(wre[1:eps.nws-1], eps.eels[1:eps.nws-1],'red')
  
#return struct with freqno, re(w), im(w), re(eps(w)), im(eps(w))
if __name__=='__main__':
  extra, vars = gen.parse_args(sys.argv[1:])
  output_files = []

  if extra == []: 
    print "python proc_out.py ./dir suffix"
    #sys.exit()
    #gen.find_all_ext('MoS2perp', 'out', output_files)
    #gen.find_all_ext('MoS2perp2', 'out', output_files)
    gen.find_all_ext('MoS2K', 'out', output_files)
    gen.find_all_ext('MoS2K2', 'out', output_files)
    gen.find_all_ext('MoS2K3', 'out', output_files)
  else:
    gen.find_all_ext(extra[0], extra[1], output_files)

  eps_objs = []
#parse multiple
  for file in output_files:
    try:
      eps_tmp = pull_eps(open(file, 'r').read())
    except:
      print "couldn't parse file {0}".format(file)
      eps_tmp = None

    eps_objs.append(eps_tmp)

  for eps in eps_objs:
    try:
      plot_eps(eps)
    except:
      print sys.exc_info()[0]

  plt.show()

