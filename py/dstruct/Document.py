#! /usr/bin/env python

"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Document class
"""

import math
import copy
import re
from helper.easierlife import *
from dstruct.Sentence import *
from dstruct.Word import *
from dstruct.Box import *

class Document(object):

  docid = None
  sents = None

  def __init__(self, _docid):
    self.docid = _docid
    self.sents = []
    self.sents.append(Sentence())

  def push_word(self, word):
    if self.sents[-1].push_word(word) == False:
      self.sents.append(Sentence())
      self.sents[-1].push_word(word)


