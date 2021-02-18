#!/usr/bin/env python

import sys
from Bio import SeqIO

class Cleaner(object):
    infile = None
    outfile = None
    min_length = 0
    n_frac = 1.0

    def parseArgs(self, args):
        prev = ""
        for a in args:
            if prev == "-l":
                self.min_length = int(a)
                prev = ""
            elif prev == "-f":
                self.n_frac = float(a)/100.0
                prev = ""
            elif a in ["-l", "-f"]:
                prev = a
            elif self.infile is None:
                self.infile = a
            else:
                self.outfile = a

    def clean(self):
        nin = nout = 0
        with open(self.outfile, "w") as out:
            for seq_record in SeqIO.parse(self.infile, "fasta"):
                nin += 1
                seq = str(seq_record.seq).upper()
                nc = seq.count("N")
                sl = len(seq)
                if len(seq) >= self.min_length and 1.0*nc/sl <= self.n_frac:
                    SeqIO.write(seq_record, out, "fasta")
                    nout += 1
        sys.stderr.write("{} sequences read, {} written ({} filtered).\n".format(nin, nout, nin-nout))
            
if __name__ == "__main__":
    args = sys.argv[1:]
    C = Cleaner()
    C.parseArgs(args)
    C.clean()
