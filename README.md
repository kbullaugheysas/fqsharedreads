# fqmultioverlap

## Overview

Ordinarily one doesn't expect to see many full exact paired-end read matches between biologically distinct samples in NGS projects. This utility compares one sample against an arbitrary number of other samples to identify which reads are shared between the focal sample and one or more of the other samples.

## Where do shared reads come from?

There are three main possibilities for shared reads:

1. Low complexity sequences that one might expect to see due to the sequencing technologies. For example reads that are entirely composed of G nucleotides are typical on Illumina sequences when the spot the optical sensors fail to detect anything at some spots on the flowcell.
2. Contamination from some other amplified library. If a library with many copies of the same sequence contaminates multiple other samples then one would expect to see copies of the same molecules in multiple samples.
3. De novo origin of the same molecules in multiple samples. If the process under investigation tends to produce essentially identical molecules independently in multiple samples, then this would presumably show up in the sequencing data. For example 3' biased RNA-seq of highly expressed genes may very well result in the same reads appearing independently in multiple samples.

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
      -ref1 string
            fastq file for read 1 of the reference sample (required)
      -ref2 string
            fastq file for read 2 of the reference sample (required)
      -sample string
            sample ID for this sample (required)

## Output

## Examples

## Implementation
