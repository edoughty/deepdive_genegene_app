#! /usr/bin/env python

"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Word class
"""

class Word(object):
  
  insent_id = None
  word = None
  pos = None
  ner = None
  lemma = None
  dep_label = None
  dep_par = None
  sentid = None
  box = None

  def __init__(self, _insent_id, _word, _pos, _ner, _lemma, _dep_label, _dep_par, _sentid, _box):
    self.insent_id = int(_insent_id) - 1
    (self.word, self.pos, self.ner, self.lemma, self.dep_label) = (_word, _pos, _ner, _lemma, _dep_label)
    self.lemma = self.lemma.replace('"', "''") # If do not do this, outputing an Array in the language will crash
    self.lemma = self.lemma.replace('\\', "_") # If do not do this, outputing an Array in the language will crash
    self.ner = self.ner
    self.dep_par = int(_dep_par) - 1
    self.sentid = int(_sentid.split('_')[-1]) - 1
    #self.box = _box

  def __repr__(self):
    return self.word

  def get_feature(self):
    if self.ner == 'O':
      return self.lemma
    else:
      return self.ner
