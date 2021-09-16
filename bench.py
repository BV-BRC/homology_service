from Bio import SeqIO
from lempel_ziv_complexity import lempel_ziv_complexity
import sys
import time
from subprocess import Popen, PIPE, STDOUT

input_file = '/home/olson/ecoli-protein.fa'

fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)

    start = time.time()
    p = Popen(["blastp",
               "-query", "-",
               "-db", "/disks/tmp/t9.features.faa",
#               "-db", input_file,
#               "-db", "/vol/bvbrc/blast/ref/bacteria-archaea.features.faa",
               "-outfmt", "6",
               "-num_threads", "2",
               # "-out", "/dev/null",
               ], stdin=PIPE,stdout=PIPE)
    (out, err) = p.communicate(f">{name}\n{sequence}\n".encode())
    # print(out)
    count = out.count(b"\n")
    p.wait()
    end = time.time()

    elap = end - start
    if p.returncode != 0:
        print(f"error from {name}\n")

    c = lempel_ziv_complexity(sequence)
    print(f"{name}\t{c}\t{len(sequence)}\t{elap}")

