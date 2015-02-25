import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import gen_files as gen
import scipy
#from scipy.fftpack import dct

#Container for dielectric functions \epsilon(\q,\omega):
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
  def plot_eps(eps):
    epsre, epsim = zip(*eps.eps)
    eps.gen_eels()
    wre, wim = zip(*eps.w)
    plt.plot(wre[1:eps.nws-1], epsre[1:eps.nws-1], 'black')
    plt.plot(wre[1:eps.nws-1], epsim[1:eps.nws-1], 'blue')
    plt.plot(wre[1:eps.nws-1], eps.eels[1:eps.nws-1],'red')

#Container for selfenergy matrix elements band i,j: \sigma_{ij}(\w).
class Sigma(object):
  def __init__(self):
    self.sigmare  = [0.00]
    self.sigmaim  = [0.00]
    self.aspec    = [0.00]
    self.nws      = 0
    self.w        = [0.0]
    self.kvec     = (0.0, 0.0, 0.0)
    self.lda_e     = [0.0]
    self.vxc       = [0.0]
    self.exx       = [0.0]
    self.qp_energy = [0.0]
    self.qp_renorm = [0.0]

  def __repr__(self):
      return 'Sigma object for kpoint {0}'.format(self.kvec)
  def plot_sig(self):
    res   = open('resigma.dat', 'w')
    ims   = open('imsigma.dat', 'w')
    aspec = open('aspec.dat',   'w')
    for line in self.sigmare:
      print >>res, line
    for line in self.sigmaim:
      print >>ims, line
    for line in self.aspec:
      print >>aspec, line

def pull_sig(f):
  sigma_regex     = re.compile(r'GW qp renorm.*?\n\n(.*?)\Z', re.M | re.S)
  sigmare_regex   = re.compile(r'REsigma\n(.*?)IMsigma', re.M | re.S)
  sigmaim_regex   = re.compile(r'IMsigma\n(.*?)ASpec', re.M | re.S)
  sigmaspec_regex = re.compile(r'ASpec\n(.*?)\n\s{0,}\n', re.M | re.S)
  sigmakpoint_regex = re.compile(r'Sigma_k\s+(.*?)\n', re.S)

#all regexs related to matrix elements and qp renormalization
  lda_e_regex  = re.compile(r'(?<=LDA eigenval \(eV\)).*$', re.M)
  vxc_regex    = re.compile(r'(?<=Vxc expt val \(eV\)).*$', re.M)
  qp_sigmexch_regex  = re.compile(r'(?<=Sigma_ex val \(eV\)).*$',   re.M)
  qp_e_regex   = re.compile(r'(?<=GW qp energy \(eV\)).*$',   re.M)
  qp_renorm_regex = re.compile(r'(?<=QP renorm).*$',   re.M)

  try:
    lda_e       = lda_e_regex.findall(f)
    vxc         = vxc_regex.findall(f)
    qp_sigmexch = qp_sigmexch_regex.findall(f)
    qp_energy   = qp_e_regex.findall(f)
    qp_renorm   = qp_renorm_regex.findall(f)
  except:
    print "couldn't pull some matels"

  sig  = Sigma()

  sig.lda_e       = map(float, lda_e[0].split())
  sig.vxc         = map(float, vxc[0].split())
  sig.exx         = map(float, qp_sigmexch[0].split())
  sig.qp_energy   = map(float, qp_energy[0].split())
  sig.qp_renorm   = map(float, qp_renorm[0].split())

  kvec        = sigmakpoint_regex.findall(f)
  sig.kvec = [map(float,x.split()) for x in kvec]

  block        = sigmare_regex.findall(f)
  sig.sigmare  = [x for x in block[0].split('\n')]
  block        = sigmaim_regex.findall(f)
  sig.sigmaim  = [x for x in block[0].split('\n')]
  block        = sigmaspec_regex.findall(f)
  sig.aspec    = [x for x in block[0].split('\n')]
  return sig

def pull_eps(f):
# Get all the frequencies, epsre, and 
# qvecs out of the output file.
# return a list of of eps objects.
  freq_re    = re.compile(r'Imag. Frequencies:\s\n(.*?)eta',re.S)
  qpoints_re = re.compile(r'Calculation of q =    (.*?)\n')
  eps_re     = re.compile(r'(?<=inveps_{GG}\(q,w\) =)(.*?)\n')
  ecut_re       = re.compile(r'(?<=Ecut, Ngmsco:\n)(.*?)\n',re.S)

  freqs         = freq_re.findall(f)[0].split('\n')
  qpoints       = qpoints_re.findall(f)
  eps           = eps_re.findall(f)
  ecut          = ecut_re.findall(f) 
  energy, ngcut = ecut[0].split()
  print energy, ngcut
  ngcut = int(ngcut)

  try:
#    freqs     = [map(float, (x.split())) for x in freqs[:-1]]
    qpoints   = [map(float, (x.split())) for x in qpoints]
    eps       = [map(float, (x.split())) for x in eps]
  except:
    print 'oh crap'
    pass
#the last entry in the qpoints list refers to Calculation of q
#for sigma_matel so we just pop that out of the stack:
  del(qpoints[-1])
  list_eps = []

  for q in range(len(qpoints)):
    epsw = Eps()
    epsw.w    = [(w[1], w[2]) for w in freqs]
    epsw.nws  = len(epsw.w)
    epsw.eps  = eps[q*ngcut:q*ngcut+ngcut]
    epsw.qvec = qpoints[q]
    list_eps.append(epsw)
  return list_eps

def gen_spec(sig_objs):
  f = open('autospec.dat', 'w')
  nws = len(sig_objs[0].aspec)
  nqs = len(sig_objs) 
  aspec = np.zeros([nws, nqs])
  for i, sig in enumerate(sig_objs):
    for iw in range(nws):
      resig = map(float, sig.sigmare[iw].split())
      imsig = map(float, sig.sigmaim[iw].split())
      w_ev = resig[0]
      for ibnd in range(1, 19):
        res = (w_ev - sig.lda_e[ibnd] - (resig[ibnd] + sig.exx[ibnd] - sig.vxc[ibnd]))
        ims = imsig[ibnd] 
       #if just broadening the LDA spec fxn.
       #ims = 0.2
        aspec[iw][i] = aspec[iw][i] + (1.0/np.pi)*(abs(ims)/(np.square(res) + np.square(ims)))

  for iq in range(len(sig_objs)):
    for iw in range(nws):
      print >>f, iq, (sig.aspec[iw].split())[0], aspec[iw][iq]
    print >>f, ''
  f.close()


def fourier_interp(sig_objs):
  nws = len(sig_objs[0].aspec)
  nqs = len(sig_objs) 
  pads = 8*nqs
  aspec = np.zeros([nws, nqs+pads])
  aux = np.zeros(nqs)

  f = open('aspecinterp.dat', 'w')
  for iw in range(nws):
    for i, sig in enumerate(sig_objs):
      tmp = map(float, sig.aspec[iw].split())
#trace over valence manifold:
      aux[i] = sum(tmp[2:17])
#  now have fourier coefficients interpolate back on to dense grid
#dct not present on magneto!
#    aux[:] = dct(aux[:], 2, norm='ortho')
#dct type 3 is the 
    auxd = np.pad(aux, (0, pads), 'constant')
 #   auxd = dct(auxd[:], 3, norm='ortho')
    for iq in range(len(auxd)):
      aspec[iw][iq] = auxd[iq]
  for iq in range(nqs+pads):
    for iw in range(nws):
      print >>f, iq, (sig.aspec[iw].split())[0], aspec[iw][iq]
    print >>f, ''
  f.close()

def plot_spec(sig_objs):
  f = open('fullspec.dat', 'w')  
  for ik, sig in enumerate(sig_objs):
    for iw in range(len(sig.aspec)):
      tmp = map(float, sig.aspec[iw].split())
      print >> f, ik, tmp[0], sum(tmp[1:17])
    print >> f 

#return struct with freqno, re(w), im(w), re(eps(w)), im(eps(w))
if __name__=='__main__':
  extra, vars = gen.parse_args(sys.argv[1:])
  output_files = []

  if extra == []: 
    print "python proc_out.py ./dir suffix"
    sys.exit()
  else:
    gen.find_all_ext(extra[0], extra[1], output_files)

  sig_objs = []
  print output_files
  for file in output_files:
    try:
      sig_tmp = pull_sig(open(file, 'r').read())
      sig_objs.append(sig_tmp)
    except:
      print "couldn't pull {0}".format(file)

  sig_objs.sort(key=lambda x: np.linalg.norm(x.kvec))
  print sig_objs

  plot_spec(sig_objs)
  fourier_interp(sig_objs)
  gen_spec(sig_objs)

