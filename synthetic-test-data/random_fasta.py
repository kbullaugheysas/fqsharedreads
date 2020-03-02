#!/usr/bin/env python

import argparse
import sys
import random
import string

parser = argparse.ArgumentParser(description='Generate some random fasta data')
parser.add_argument('-reads', type=int, help='Number of reads', required=True)
parser.add_argument('-len', type=int, help='Read length', default=75)
args = parser.parse_args()

def log(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def randomReadName():
    return "".join([string.ascii_lowercase[random.randint(0,25)] for x in range(12)])

qual = "E" * args.len

bases = ["A", "T", "G", "C"]

for i in range(args.reads):
    name = f"@{randomReadName()}"
    seq = "".join([bases[random.randrange(4)] for j in range(args.len)])
    print(name)
    print(seq)
    print("+")
    print(qual)

# END
