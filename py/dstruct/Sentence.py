#! /usr/bin/env python

"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Sentence class
"""

class Sentence(object):

  _MAXLEN = 1000 # to avoid bad parse tree that have self-recursion
  sentid = None
  words = None

  def __init__(self):
    self.words = []
    self.sentid = None

  def __repr__(self):
    return " ".join([w.word for w in self.words])
  
  def push_word(self, word):
    if self.sentid == None:
      self.sentid = word.sentid
      self.words.append(word)
      return True
    else:
      if self.sentid == word.sentid:
        self.words.append(word)
        return True
      else:
        return False

  def get_path_till_root(self, word_index):
    """ Get all words in the dependency tree from word_index to the root """

    path = []
    c = word_index
    MAXLEN = self._MAXLEN
    while MAXLEN > 0:
      MAXLEN = MAXLEN -1
      try:
        if c == -1: break
        path.append(c)
        c = self.words[c].dep_par
      except:
        break
    return path

  def get_common_ancestor(self, path1, path2):
    """ Given two paths returned by get_path_till_root, find the least common ancestor """
    
    parent = None
    for i in range(0, max(len(path1), len(path2))):
      tovisit = 0 - i - 1
      if i >= len(path1) or i >= len(path2):
        break
      if path1[tovisit] != path2[tovisit]:
        break
      parent = path1[tovisit]
    return parent
  
  def get_direct_dependency_path_between_words(self, idx1, idx2):
    """
    Given two word idx1 and idx2, where idx2 is the parent of idx1, return the
    words on the dependency path
    """

    words_on_path = []
    c = idx1
    MAXLEN = self._MAXLEN
    while MAXLEN > 0:
      MAXLEN = MAXLEN - 1
      try:
        if c == -1: break
        if c == idx2: break
        if c == idx1: 
          words_on_path.append(str(self.words[c].dep_label)) # we do not include the word of idx1
        else:
          words_on_path.append(str(self.words[c].dep_label) + "|" + self.words[c].get_feature())
        c = self.words[c].dep_par
      except:
        break
    return words_on_path

  def get_word_dep_path(self, idx1, idx2):
    """ Given two word idx1 and idx2, return the dependency path feature between them """
    path1 = self.get_path_till_root(idx1)
    path2 = self.get_path_till_root(idx2)

    parent = self.get_common_ancestor(path1, path2)

    words_from_idx1_to_parents = self.get_direct_dependency_path_between_words(idx1, parent)
    words_from_idx2_to_parents = self.get_direct_dependency_path_between_words(idx2, parent)

    return "-".join(words_from_idx1_to_parents) + "@" + "-".join(words_from_idx2_to_parents)

  def get_prev_wordobject(self, mention):
    """ Given two word idx1 and idx2, return the dependency path feature between them """
    begin = mention.prov_words[0].insent_id
    if begin - 1 < 0: 
      return None
    else: 
      return self.words[begin - 1]

  def dep_parent(self, mention):
    begin = mention.prov_words[0].insent_id
    end = mention.prov_words[-1].insent_id

    paths = []
    for i in range(begin, end+1):
      for j in range(0, len(self.words)):
        if j >= begin and j <= end: continue

        path = self.get_word_dep_path(i, j)
        paths.append(path)

    return sorted(paths, key=len)[0:min(5,len(paths))]
        
        
  def dep_path(self, entity1, entity2):
    begin1 = entity1.prov_words[0].insent_id
    end1 = entity1.prov_words[-1].insent_id
    begin2 = entity2.prov_words[0].insent_id
    end2 = entity2.prov_words[-1].insent_id
    
    paths = []
    for idx1 in range(begin1, end1+1):
      for idx2 in range(begin2, end2+1):
        paths.append(self.get_word_dep_path(idx1, idx2))

    # we pick the one that is shortest
    path = ""
    ll = 100000000
    for p in paths:
      if len(p) < ll:
        path = p
        ll = len(p)
    return path


