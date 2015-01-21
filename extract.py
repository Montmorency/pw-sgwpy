import sys
import re
import itertools
import os, os.path
import numpy as np

#extract output imsigma and aspec from files, orders them
#copies to new directory...


def recursive_find(dir):
  files = []
  for f in os.listdir(dir):
    if os.path.isdir(os.path.join(dir, f)):
      files = files + recursive_find(os.path.join(dir, f))
    elif f.endswith('.out'):
      files.append(os.path.join(dir,f))
#    elif not f.endswith('.pyc'):
#      files.append(os.path.join(dir, f))
  return files

if __name__=='__main__': 
  output = recursive_find('./')
  target_dir='./SiDat'

  print output
  print output.sort()
#pull out integer tag and sort the list
  reg = re.compile("[0-9]+")
  output = sorted(output, key=lambda x: int(reg.findall(x)[0]))
  print output

  for i, f in enumerate(output):
    print i, f
    tmp_dir = os.path.join(target_dir,'K{0}'.format(i))

    if not os.path.isdir(tmp_dir):
      os.mkdir(tmp_dir)
    else:
      print 'directory already made.'

    f = open(f).read()
    ims   = open(os.path.join(tmp_dir,'imsigma.dat'), 'w')
    aspec = open(os.path.join(tmp_dir,'aspec.dat'),   'w')

    sigmaim_regex   = re.compile(r'IMsigma\n(.*?)ASpec', re.M | re.S)
    sigmaspec_regex = re.compile(r'ASpec\n(.*?)\n\s{0,}\n', re.M | re.S)

    block = sigmaim_regex.findall(f)
    for line in block[0].split('\n'):
      print >>ims, line 

    block = sigmaspec_regex.findall(f)
    for line in block[0].split('\n'):
      print >>aspec, line 

