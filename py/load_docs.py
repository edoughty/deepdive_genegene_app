#! /usr/bin/env python

"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Loads documents into database. 
""" 

from helper.easierlife import *
import fileinput

from dstruct.Document import *
from dstruct.Word import *

for row in fileinput.input():

  (docid, folder) = row.rstrip('\n').split('\t')
  log(docid)
  doc = Document(docid)

  #try:
  for l in open(folder):
    ss = l.rstrip().split('\t')
    if len(ss) < 3: continue
    (insent_id, word, pos, ner, lemma, deppath, deppar, sentid, box) = ss
    doc.push_word(Word(insent_id, word, pos, ner, lemma, deppath, deppar, sentid, box))

  print "\t".join(["\\N", docid, serialize(doc)])










