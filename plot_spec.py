import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import gen_files as gen
import scipy
from scipy.fftpack import dct

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

#finds most 'appropriate frequency to 
#interpolate imaginary selfenergy 
#onto frequency grid.

def find_best_sigim(w_evs, sig_obj):
  nwsig  = len(sig_obj.aspec)
  w_corr = np.zeros(len(w_evs), dtype = np.int)
  w_sig = np.zeros(nwsig)
#generate list of freqs produced in SGW output:
  for iw in range(nwsig):
    w_sig[iw] = sig_obj.aspec[iw].split()[0]
#find correspondence in fine mesh:
  for i, w in enumerate(w_evs):
    for j in range(nwsig):
      if w_sig[j]-0.2 < w < w_sig[j]+0.2:
        w_corr[i] = int(j)
  return w_corr

def gen_qpeig_spec(sig_objs):
  #f = open('lda_spec.dat', 'w')
  f = open('qp_spec.dat', 'w')
  w_evs = np.arange(-7.0, 0.1, 0.1)
  sig = sig_objs[:0:-1]
  [sig.append(x) for x in sig_objs]
  w_corr = find_best_sigim(w_evs, sig[0])
  for iq, sig in enumerate(sig):
    eigs = sig.qp_energy
    #eigs = sig.lda_e
    for i, w_ev in enumerate(w_evs):
      aspec = 0.0
      if w_corr[i] != 0:
        imsig = map(float, sig.sigmaim[w_corr[i]].split())
      elif w_corr[i] == 0 :
        imsig = [2.94 for x in range(13)]
      #print imsig
      for j, eig in enumerate(eigs[2:9]):
        res = (w_ev - eig)
        #ims = 0.2
# with ab initio broadening:
        ims = imsig[j+3]
        aspec = aspec + (1.0/np.pi)*(abs(ims)/(np.square(res) + np.square(ims)))
        print >>f, iq, w_ev, aspec
    print >>f, ''
  f.close()

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

def gen_spec(sig_objs):
  f = open('autospec.dat', 'w')
  nws = len(sig_objs[0].aspec)
  nqs = len(sig_objs)
  aspec = np.zeros([nws, nqs])
  nbnd = 9

  print>>f, "#eig value band", sig_objs[0].lda_e[nbnd-1]

  for i, sig in enumerate(sig_objs):
    for iw in range(nws):
      resig = map(float, sig.sigmare[iw].split())
      imsig = map(float, sig.sigmaim[iw].split())
      w_ev = resig[0]

      for ibnd in range(nbnd-1, nbnd):
        res = (w_ev - sig.lda_e[ibnd] - (resig[ibnd] + sig.exx[ibnd] - sig.vxc[ibnd]))
        ims = imsig[ibnd] 
        aspec[iw][i] = aspec[iw][i] + np.pi*(abs(ims)/(np.square(res) + np.square(ims)))
  for iq in range(len(sig_objs)):
    for iw in range(nws):
      print >>f, float(iq)/float(nqs), (sig.aspec[iw].split())[0], aspec[iw][iq]
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
      aux[i] = sum(tmp[2:9])
#  now have fourier coefficients interpolate back on to dense grid
    aux[:] = dct(aux[:], 2, norm='ortho')
#dct type 3 is the 
    auxd = np.pad(aux, (0, pads), 'constant')
    auxd = dct(auxd[:], 3, norm='ortho')
    for iq in range(len(auxd)):
      aspec[iw][iq] = auxd[iq]
  for iq in range(nqs+pads)[::-1]:
#pads with zeros
    for iw in range(nws):
      print >>f,  -iq, (sig.aspec[iw].split())[0], aspec[iw][iq]
    print >>f, ''
  for iq in range(nqs+pads):
#pads with zeros
    for iw in range(nws):
      print >>f, iq, (sig.aspec[iw].split())[0], aspec[iw][iq]
    print >>f, ''
  f.close()

def plot_spec(sig_objs):
  f = open('fullspec.dat', 'w')  
  for ik, sig in enumerate(sig_objs):
    for iw in range(len(sig.aspec)):
      tmp = map(float, sig.aspec[iw].split())
      print >> f, ik, tmp[0], sum(tmp[1:9])
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
  for file in output_files:
    try:
      sig_tmp = pull_sig(open(file, 'r').read())
      sig_objs.append(sig_tmp)
    except:
      print "couldn't pull {0}".format(file)

  sig_objs.sort(key=lambda x: np.linalg.norm(x.kvec))
  gen_spec(sig_objs)
#  plot_spec(sig_objs)

#  fourier_interp(sig_objs)
#  gen_qpeig_spec(sig_objs)

