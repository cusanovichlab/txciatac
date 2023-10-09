from __future__ import print_function

import pysam
import sys

inbam = sys.argv[1]
outbam = sys.argv[2]

readsin = pysam.AlignmentFile(inbam, "rb")
readsout = pysam.AlignmentFile(outbam, "wb", template=readsin)

last_reference = None
observed_reads = set()
for read in readsin:

    reference_id = read.reference_id

    if last_reference != reference_id:
        print('Deduplicating %s...' % reference_id)
        last_reference = reference_id
        observed_reads = set()

    readname = read.qname.split(':')[0]

    # This part is a little tricky. For each read, I ask for the molecule start and end (keeping strand in mind).
    # If a molecule has been seen with the same start/stop on the same strand for that cell, I skip the current read.
    # Effectively, if there are PCR duplicates of a given fragment, I pick a random representative read for each end of the molecule.
    if read.tlen < 0:
        fragstart = read.mpos - read.tlen
        fragend = read.mpos
    else:
        fragstart = read.pos
        fragend = read.pos + read.tlen

    read_id = (readname, fragstart, fragend, read.is_read1)

    if read_id not in observed_reads:
        observed_reads.add(read_id)
        readsout.write(read)

readsin.close()
readsout.close()
