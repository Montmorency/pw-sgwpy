import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import gen_files as gen

class Eps(object):
  def __init__(self):
    self.w    = [(0.0, 0.0)]
    self.eps  = [(0.0, 0.0)]
    self.qvec = (0.0, 0.0, 0.0)
    self.nws  = 0

  def __repr__(self):
      return 'epsilon object for qpoint {0}'.format(self.qvec)

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
    wre, wim = zip(*eps.w)

    plt.plot(wre, epsre[:eps.nws],'r')
    plt.plot(wre, epsim[:eps.nws],'b')
  
#return struct with freqno, re(w), im(w), re(eps(w)), im(eps(w))
if __name__=='__main__':
  extra, vars = gen.parse_args(sys.argv[1:])
  print extra, vars

  output_files = []
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

  for eps in eps_objs[1:-1]:
    print eps
    plot_eps(eps)

  plt.show()

