import sys
import itertools
import os, os.path
import numpy as np
from gen_files import parse_args

if __name__=='__main__':

  extra, vars = parse_args(sys.argv[1:])
  vars_values = []

  for var_name, values in vars:
    values = values.split(",")
    vars_values.append((var_name, values[0]))

  vars = dict(vars_values)
  vars['jdir'] =  os.getcwd()

  try:
    pbs = open(extra[0]).read()
  except:
    print 'python job_sub.py --jobin foo -jobout bar ./template.pbs ./actual_job.pbs'
    sys.exit()

  print >> open(extra[1], 'w'), pbs.format(**vars)

