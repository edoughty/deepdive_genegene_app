#! /usr/bin/env python

"""
Author: Emily K Mallory and Ce Zhang
Date: 08/24/15
Contact: emily.mallory@stanford.edu

Loads docids into database. 
"""

from helper.easierlife import *

#location of the text corpus.
INPUT_FOLDER = BASE_FOLDER + "/data/plos_full_conll"
INPUT_NEG_FOLDER = BASE_FOLDER + "/data/NEGARTICLE"
INPUT_FILE_GS_SKIP = BASE_FOLDER + "/data/plos_journals_dip_mint_pmids.txt"
INPUT_FILE_10K_SKIP = BASE_FOLDER + "/data/plos_docids_sample_10000.txt"

pmids_skip = set()
for x in [x.strip().split("\t")[0] for x in open(INPUT_FILE_GS_SKIP).readlines()]:
	pass
	#pmids_skip.add(x.rstrip().lower())

for x in [x.strip() for x in open(INPUT_FILE_10K_SKIP).readlines()]:
    pass
	#pmids_skip.add(x.rstrip().lower())

log("~~~~~~~~~~~~~~~~~" + INPUT_FOLDER)

for docid in os.listdir(INPUT_FOLDER):
	if docid.startswith('journal.pbio.') or docid.startswith('journal.pgen.') or docid.startswith('journal.pone.'): 
		folder = INPUT_FOLDER + "/" + docid 
		if docid.split(".pdf")[0] not in pmids_skip and docid.split(".nlp")[0] not in pmids_skip:
			print json.dumps({"docid": docid.split(".nlp")[0], "folder": folder})


log("~~~~~~~~~~~~~~~~~" + INPUT_NEG_FOLDER)

for docid in os.listdir(INPUT_NEG_FOLDER):
	if docid.startswith('.'): continue
	folder = INPUT_NEG_FOLDER + "/" + docid
	if docid.split(".")[0] not in pmids_skip and docid not in pmids_skip:
		print json.dumps({"docid": docid, "folder": folder})
