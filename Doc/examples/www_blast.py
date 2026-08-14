#!/usr/bin/env python
# Copyright 2000 Brad Chapman.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Example showing how to deal with internet BLAST from Biopython.

This code is described in great detail in the BLAST section of the Biopython
documentation.
"""

# standard library
from io import BytesIO

# biopython
from Bio import Blast
from Bio import SeqIO

# first get the sequence we want to parse from a FASTA file
f_record = next(SeqIO.parse("m_cold.fasta", "fasta"))

print("Doing the BLAST and retrieving the results...")
result_stream = Blast.qblast("blastn", "nr", format(f_record, "fasta"))

# save the results for later, in case we want to look at it
blast_results = result_stream.read()
with open("m_cold_blast.out", "wb") as save_file:
    save_file.write(blast_results)

print("Parsing the results and extracting info...")

# option 1 -- open the saved file to parse it
# option 2 -- create a stream from the bytes and parse it
blast_record = Blast.read(BytesIO(blast_results))

# now get the alignment info for all e values greater than some threshold
E_VALUE_THRESH = 0.1

for hit in blast_record:
    for hsp in hit:
        if hsp.annotations["evalue"] < E_VALUE_THRESH:
            print("****Alignment****")
            print("sequence: %s %s" % (hit.target.id, hit.target.description))
            print("length: %i" % len(hit.target))
            print("e value: %f" % hsp.annotations["evalue"])
            print(hsp[1][0:75] + "...")  # aligned query
            print(hsp.annotations["midline"][0:75] + "...")
            print(hsp[0][0:75] + "...")  # aligned target
