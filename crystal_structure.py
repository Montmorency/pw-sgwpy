import sys
import re

#Object container for basic simulation parameters
#job_prefix, real space unit cell, reciprocal space unit cell.
class CrystStruct(object):
  def __init__(self):
    self.prefix = 'Si'
    self.rvec = [(1,0,0), (0,1,0), (0,0,1)]
    self.gvec = [(1,0,0), (0,1,0), (0,0,1)]
