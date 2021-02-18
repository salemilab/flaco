#!/usr/bin/env python

import os
import sys
import csv
import time
import glob
from collections import defaultdict
from Bio import SeqIO

# Clean
def main_clean(fastafile, words):
    nin = 0
    nout = 0

    for seq in SeqIO.parse(fastafile, "fasta"):
        nin += 1
        header = seq.description.replace(" ", "_")
        skip = False
        for w in words:
            if w in header:
                skip = True
        if skip:
            continue
        sys.stdout.write(">" + header + "\n")
        sys.stdout.write(str(seq.seq).replace("-", ""))
        sys.stdout.write("\n")
        nout += 1
    sys.stderr.write("{} sequences in input files, {} after filtering.\n".format(nin, nout))


# Split
def parse_header(line):
    parts = line.rstrip("\r\n").split("|")
    if len(parts[2]) > 7:
        return parts[1], parts[2][:7]
    else:
        return parts[1], parts[2]

def main_split(fastafile, outdir):
    out = None
    for seq in SeqIO.parse(fastafile, "fasta"):
        (seqid, month) = parse_header(seq.name)
        destdir = outdir + "/" + month
        outfile = destdir + "/" + seqid + ".fa"
        try:
            os.makedirs(destdir)
        except:
            pass
        with open(outfile, "w") as out:
            SeqIO.write(seq, out, "fasta")

def previous_month(month):
    y, m = month.split("-")
    y = int(y)
    m = int(m)
    if m == 1 or m == 0:
        return "{}-{:02d}".format(y-1, 12)
    else:
        return "{}-{:02d}".format(y, m-1)

def next_month(month):
    y, m = month.split("-")
    y = int(y)
    m = int(m)
    if m == 12:
        return "{}-{:02d}".format(y+1, 1)
    else:
        return "{}-{:02d}".format(y, m+1)

def getStream(streams, month, outdir):
    if month in streams:
        return streams[month]
    else:
        destdir = outdir + "/" + month
        outfile = destdir + "/DB"
        try:
            os.makedirs(destdir)
        except:
            pass
        out = open(outfile, "w")
        streams[month] = out
        return out

def main_splitdb(fastafile, cleanfile, outdir, excludefile, words):
    nin = nout = nexcl = 0
    streams = {}
    with open(excludefile, "w") as eout, open(cleanfile, "w") as cout:
        for seq in SeqIO.parse(fastafile, "fasta"):
            nin += 1
            header = seq.description.replace(" ", "_")
            bases = str(seq.seq).replace("-", "")
            excl = False
            for w in words:
                if w in header:
                    excl = True
            if excl:
                eout.write(">" + header + "\n" + bases + "\n")
                nexcl += 1
            else:
                cout.write(">" + header + "\n" + bases + "\n")
                nout += 1
                (seqid, month) = parse_header(header)
                pm = previous_month(month)
                nm = next_month(month)
                out1 = getStream(streams, month, outdir)
                out2 = getStream(streams, pm, outdir)
                out3 = getStream(streams, nm, outdir)
                out1.write(">" + header + "\n" + bases + "\n")
                out2.write(">" + header + "\n" + bases + "\n")
                out3.write(">" + header + "\n" + bases + "\n")

            sys.stderr.write("\rSeq: {}, DB: {}, Excluded: {}".format(nin, nout, nexcl))
            sys.stderr.flush()
    for m in streams.keys():
        streams[m].close()
    sys.stderr.write("\nSequences read from {}: {}\n  Written to {}: {}\n  Filtered to {}: {}\n".format(fastafile, nin, cleanfile, nout, excludefile, nexcl))

# Parse

def getQuery(filename):
    with open(filename, "r") as f:
        l1 = f.readline().split("\t")
    return l1[0].split("|")

def monthBefore(date):
    try:
        year, month, day = [ int(p) for p in date.split("-") ]
    except:
        sys.stderr.write("1 Bad date - {}\n".format(date))
        sys.exit(1)
    if month == 1:
        year = year - 1
        month = 12
    else:
        month = month - 1
    return "{}-{:02}-{:-2}".format(year, month, day)

def monthAfter(date):
    try:
        year, month, day = [ int(p) for p in date.split("-") ]
    except:
        sys.stderr.write("2 Bad date - {}\n".format(date))
        sys.exit(1)
    if month == 12:
        year = year + 1
        month = 1
    else:
        month = month + 1
    return "{}-{:02}-{:02}".format(year, month, day)

def g(l, n):
    if l:
        return l[n]
    else:
        return "-"

def same(h1, h2):
    return (h1[2] == h2[2]) and (h1[3] == h2[3])

def getBestHits(filename, out, extra=False):
    try:
        (name, seqid, date, country) = getQuery(filename)
    except:
        sys.stderr.write("Cannot read header from " + filename + "\n")
        return []

    if len(date) == 10:
        pass
    elif len(date) == 7:
        sys.stderr.write("3 Bad date - " + date + " - correcting\n")
        date = date + "-15"
    else:
        sys.stderr.write("3 Bad date - " + date + "\n")
        return []

    bef = monthBefore(date)
    aft = monthAfter(date)

    needb = True
    needa = True

    bestBefore = None
    bestAfter = None

    ids = []

    with open(filename, "r") as f:
        c = csv.reader(f, delimiter='\t')
        for line in c:
            subject = line[1].split("|")
            subdate = subject[2]
            
            if needb:
                if bef <= subdate <= date: # still within range?
                    if bestBefore is None:
                        bestBefore = line
                        ids.append(line[1])
                    elif extra and same(line, bestBefore):
                        ids.append(line[1])
                    else:
                        needb = False

            if needa:
                if date <= subdate <= aft:
                    if bestAfter is None:
                        bestAfter = line
                        ids.append(line[1])
                    elif extra and same(line, bestAfter):
                        ids.append(line[1])
                    else:
                        needa = False

            if (not needa) and (not needb):
                break
    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(seqid, date, g(bestBefore, 1), g(bestBefore, 2), g(bestBefore, 3), g(bestAfter, 1), g(bestAfter, 2), g(bestAfter, 3)))
    return ids

def main_parse(blastdir, outfile, matches):
    filenames = glob.glob(blastdir + "/*/*.blast.csv")
    allids = defaultdict(int)
    with open(outfile, "w") as out:
        out.write("#SeqId\tDate\tBefore\tBefIdent\tBefLength\tAfter\tAftIdent\tAftLength\n")
        sys.stderr.write("Finding best hits in:\n")
        for filename in filenames:
            sys.stderr.write("\r  " + filename)
            sys.stderr.flush()
            ids = getBestHits(filename, out)
            for i in ids:
                allids[i] += 1
    sys.stderr.write("\n{} best match sequences found.\n".format(len(allids)))
    with open(matches, "w") as out:
        for k in allids.keys():
            r = k.split("|")
            out.write("{}\t{}\t{}\t{}\t{}\n".format(k, r[1], r[0].split("/")[1], r[2], allids[k]))

# Extract

def readWanted(filename):
    wanted = []
    with open(filename, "r") as f:
        for line in f:
            wanted.append(line.split("\t")[0])
    sys.stderr.write("{} sequences wanted.\n".format(len(wanted)))
    return wanted

def main_extract(fastafile, matches, outfile):
    wanted = readWanted(matches)
    nin = n = 0
    with open(outfile, "w") as out:
        for seq in SeqIO.parse(fastafile, "fasta"):
            nin += 1
            hdr = seq.description.replace(" ", "_") # in the blast db we removed spaces from sequence names
            #name = hdr[1:].rstrip("\r\n")
            if hdr in wanted:
                out.write(">" + hdr + "\n")
                out.write(str(seq.seq).replace("-", "") + "\n")
                n += 1
            sys.stderr.write("\r{}/{}".format(n, nin))
            sys.stderr.flush()
    sys.stderr.write("\n{} sequences extracted.\n".format(n))

# Main

if __name__ == "__main__":
    cmd = sys.argv[1]
    args = sys.argv[2:]

    if cmd == "clean":
        infile = args[0]
        words = args[1:]
        main_clean(infile, words)

    elif cmd == "splitdb":
        infile = args[0]
        cleanfile = args[1]
        outdir = args[2]
        excludefile = args[3]
        words = args[4:]
        main_splitdb(infile, cleanfile, outdir, excludefile, words)

    elif cmd == "split":
        infile = args[0]
        splitdir = args[1]
        main_split(infile, splitdir)

    elif cmd == "parse":
        main_parse(args[0], args[1], args[2])

    elif cmd == "extract":
        main_extract(args[0], args[1], args[2])
