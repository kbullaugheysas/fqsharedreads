# fqmultioverlap

## Overview

Ordinarily one doesn't expect to see many full exact paired-end read matches
between biologically distinct samples in NGS projects. This utility compares
one sample against an arbitrary number of other samples to identify which reads
are shared between the focal sample and one or more of the other samples.

## Where do shared reads come from?

There are three main possibilities for shared reads:

1. Low complexity sequences that one might expect to see due to the sequencing
   technologies. For example reads that are entirely composed of G nucleotides
   are typical on Illumina sequences when the spot the optical sensors fail to
   detect anything at some spots on the flowcell.
2. Contamination from some other amplified library. If a library with many
   copies of the same sequence contaminates multiple other samples then one
   would expect to see copies of the same molecules in multiple samples.
3. De novo origin of the same molecules in multiple samples. If the process
   under investigation tends to produce essentially identical molecules
   independently in multiple samples, then this would presumably show up in the
   sequencing data. For example 3' biased RNA-seq of highly expressed genes may
   very well result in the same reads appearing independently in multiple samples.

## Usage

You can see the usage help for the program as follows:

    fqmultioverlap -help

Which produces something like this:

    usage: fqmultioverlap [options]
      -batches int
            process files in batches to avoid open file limits (default 1)
      -continue string
            file with output from an existing run we'll add to
      -files string
            file that contains the list of fastq files (required)
      -limit int
            only consider the first LIMIT fastq records in each sample
      -progress string
            write data after each batch to this file
      -sample string
            sample ID for this sample (required)

## Relevant files

### The -files argument

The file given with the `-files` argument shoudl be a tab-separated file with
three columns:

    1. Sample ID
    2. Path to fastq file for mate 1
    3. Path to fastq file for mate 2

For example, it could look like this:

    sampleA     sampleA_1.fq    sampleA_2.fq
    sampleB     sampleB_1.fq    sampleB_2.fq

The reference sample must be included as an entry in this file.

### Output and -continue argument

The output file consists of tab-separated header lines that start with a `#`
character followed by a space and then a label. The possible header lines are:

    1. `sample` - the sample ID that was given with the -sample argument.
    0. `ref1` - the path to the fastq file that was used as the mate 1 sequence for this sample
    0. `ref2` - the path to the fastq file that was used as the mate 2 sequence for this sample
    0. `overlap` - these lines correspond to all of the samples that were
        compared with the reference for possible shared reads. The three fields
        after the label are the same as those given in the list of fastq files
        specified with the `-files` argument

When adding more samples during a subsequent run of `fqmultioverlap`, the
output from a previous run is passed via the `-continue FILE` argument. The
list of fastq files given with the `-files FILE` argument must include an entry
for the reference sample, as before, and the paths must match the `ref1` and
`ref2` headers in the previous runs output.

## Testing & Examples

There is an included subdirectory, `synthetic-test-data`, that contains some
tests that can be run as follows:

    cd synthetic-test-data
    ./run_tests.sh

These tests create random fastq data for 3 samples that vary in how many shared
reads they have among them. 

Once you have run the test script, you'll have a tab-separated file, `fqlist`,
with contents like this:

    sampleA     sampleA_1.fq    sampleA_2.fq
    sampleB     sampleB_1.fq    sampleB_2.fq
    sampleC     sampleC_1.fq    sampleC_2.fq

This file lists the samples among which we want to examine overlap.

If we want to identify shared reads from sampleC that are found in either of
the other two samples, we can run the following command, which is part of the
test script:

    fqmultioverlap -files fqlist -sample sampleC > sampleC.shared

This will produce a file something like this:

    # sample        sampleC
    # ref1          sampleC_1.fq
    # ref2          sampleC_2.fq
    # overlap       sampleA sampleA_1.fq    sampleA_2.fq
    # overlap       sampleB sampleB_1.fq    sampleB_2.fq
    lhytsyljjuad    TTTCAACGCATTGGAGGTGTGTACTTGACTCGGCAAGCGAGAAGAGCCTAAGTGTTAGTAATATTATTTCCGGTGA    AACGACCCCTTATAAGTGAGGTCGGTGAGAGGGTCCTCACGGCGGATTGTTTATTCACCATGAATCCCACCTTTTT    sampleB
    dfqzvjbnralh    TTTGCTGAATATCTAGGGCCTCACATAAGCTTTGCCACCCCGCGACTGCGCAACTCTAATCAGAACTCGATTGCTA    AGTGAACTCCTCAGGGGTCTCGTACCCACACGCCTCCATGTTAACTCTGCACATAATCCTTGGAATTGGTGTCGCC    sampleB
    nwubztucacpz    ACGTTTTTCAGATTTCAGTACGCCATCCTCCCACCACAACTATCACGAACGACGAAGATCCCGATATGGTTAACTA    GCTAGTATTCGTCAACCCCAGTAGAGCAACCCCACGACTAGGGAACCGACAACCTATGATTACGTCTGAGAAGTTA    sampleB
    uwmgtyxavcvg    CAATGGGTGGTTGACCGGAGATTGGGCCACCGGTCCCCTTATACACAAGATACGGATAGATAGCTCAGGTGGATGG    AGCCTATCCCCTGCTGAACCTTCTGCCCCCGTAGGGGTCCGCCTTGCTGATTCTGCGTGTCCGTCAGGTGTCGCAT    sampleA
    tccjqxastdug    ATGCTACGGAACCAGGTGTACCGGCATTTCGCCAAAACGTCGGTCGCCTGGAATCGGCTCGTAACCCGGTAATCCT    GTCATAGGTCGCCAGGCTGACTCTCTAACTACGCGCAAATGCTGACATATGCGCCTAAGCAGGAAACGACTGGTAG    sampleA
    ekrluqbjcppp    TCCTCGTCGAGTCCGTGGCACGGGGGCGCGGGGAGTATCCTATGACCAGGTACCAATCTGGAATGCAGGGTCACAT    CTTGGTTAATGCGTAGTCTTAATCCGCGGAAGGCTCCTTAAGCCCGGACATGTAAAAATTCAGCCAGACCATAGCA    sampleA
    nwajkdqvbjfw    TCAAGTTGTACTGCGATTCCGAGCTTGTACCGGTGTTTATACGGTTAGCCTACTCCTGTCACAGGATACTTCACTG    GTAAGGATAATGCTTTCCACGGCGTAGGTAAATAGGCGCATTCCTAACTCTTGACCTTCGTTCAGTAAGGGCCCCG    sampleB
    noecghnccjrj    TTCTGAGGATCGGTGCGGATTCACCCTATTGAACTTCTGTGCGGGGAAAGCGTCTCATTCCCGCCCGTAACAGCAC    TATCACTGAGCAGTACATCTTATGGAAATTCGTACGTTAGTGTCACTCCTTCAAACTATCGTGCGCTAGGGAACCG    sampleB,sampleA


