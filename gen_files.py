import sys
import itertools
import os, os.path
import numpy as np

"""
Script for varying a bunch of variables in a bunch of files to generate a bunch of calculations.
Makes a directory for each, labelled by varied variables, i.e. ones with more than one value.
gen_files.py --prefix --wstart --wstop --dw etc. [target root dir] [template/file].
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

def freq_grid(wstart, wstop, dw):
  wfs = []
#first freq needs to be 0.0 0.0
  wfs.append('0.0d0 0.0d0')
  b = 0.4
  print (wstop-(wstart+dw))/dw 
  for a in np.arange(wstart, wstop,dw):
    wfs.append('{0}  {1}'.format(a,b))
  return(wfs)

def k_path():
  kpts = []
  for kxy in np.arange(0, 0.3333, 0.03333):
    kpts.append('{0}  {1}  {2}'.format(kxy, kxy, 0.0))
#  for kxy in np.arange(0, 1.000, 0.1):
#    kpts.append('{0}  {1}  {2}'.format(0.0, 0.0, kxy))
  return kpts  

if __name__=='__main__':
  wfs=freq_grid(15.0, 23.0, 0.2)
  kpts = k_path()

  extra, vars = parse_args(sys.argv[1:])
  vars = dict(vars)
  prefix = vars['prefix']

  template_file = '{0}'.format(extra[1])
  f = open(template_file).read()
  pbs = open('./template/mos2.pbs').read()
  target_dir = extra[0]

#need dict nws, ws, nks kpts
  template_context        = {}
  template_context['nws'] = len(wfs)
  template_context['ws']  = '\n'.join(wfs)
  template_context['nks'] = 1

  for i, kpt in enumerate(kpts):
    dir_path=os.path.join(target_dir, '{0}_{1}'.format(prefix, i))
    dir_path=os.path.abspath(dir_path)
#   make dir if it doesn't already exist
    if not os.path.isdir(dir_path):
      os.mkdir(dir_path)
    else:
      print 'directory already made.'

    job_name = '{0}.{1}.in'.format(prefix, i)
    source_target = os.path.join(dir_path, job_name)
    template_context['kpts'] = kpt

    print >> open(source_target,'w'), f.format(**template_context)

    source_name = '{0}.{1}.pbs'.format(prefix,i)
    source_target = os.path.join(dir_path, source_name)

    print >> open(source_target,'w'), pbs.format(jobin=job_name, jdir=dir_path, jobout='{0}.{1}.out'.format(prefix,i))

    if i==0:
      print "qsub ", '{0}/{1}.{2}.pbs'.format(dir_path,prefix,i)
    else:
      print "qsub -hold_jid", '{0}.{1}.pbs'.format(prefix, i-1), '{0}/{1}.{2}.pbs'.format(dir_path,prefix,i)

