#! /usr/bin/env python

"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Loads sentences into database. 
"""

from helper.easierlife import *
from extractor.EntityExtractor_Drug import *

for row in get_inputs():
	doc = deserialize(row["document"])
	log(doc.docid)

	for sent in doc.sents:
		sentence_text = " ".join([x.word for x in sent.words])
		print json.dumps({
			"docid": doc.docid, 
			"sentid":sent.sentid, 
			"sentence": serialize(sent),
			"text": sentence_text
    })

