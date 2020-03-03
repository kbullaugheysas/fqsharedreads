package main

import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
	"sync"
)

/* This program takes a reference sample, and a file listing other fastq files
 * and outputs the sequences and the samples containing them for any of the
 * reference sample's sequences that are found in the other samples.
 *
 * The file supplied with the -file argument should have three tab-separated
 * columns giving sampleId, fastq1, fastq2.
 */

type Args struct {
	Sample    string
	FastqList string
	Ref1      string
	Ref2      string
	Limit     int
	Batches   int
	Progress  string
	Continue  string
}

var args = Args{}

// We keep track of sequences we're looking for by storing the mate 1 read
// name. After the reference seqeunce is read in, this will not be written to
// again. And in particular, it's treated as read-only by the goroutines, and
// so simultaneous access is okay.
var sampleSequences map[string]string

// This should only be accessed when setting up the initial set of keys and by
// the reader goroutine.
var refseq map[string]map[string]int

/* Provide an ambidexterous interface to files to read that may be gzipped */
type AmbiReader struct {
	fp *os.File
	gz *gzip.Reader
	r  io.Reader
}

func (a AmbiReader) Read(b []byte) (n int, err error) {
	return a.r.Read(b)
}

func (a *AmbiReader) Open(fn string) error {
	if a.r != nil {
		return fmt.Errorf("AmbiReader already open")
	}
	var err error
	// If no filename is given, then read from stdin
	if fn == "" {
		a.r = os.Stdin
		return nil
	}
	a.fp, err = os.Open(fn)
	if err != nil {
		return err
	}
	if strings.HasSuffix(fn, ".gz") {
		a.gz, err = gzip.NewReader(a.fp)
		if err != nil {
			return err
		}
		a.r = a.gz
	} else {
		a.r = a.fp
	}
	return nil
}

func (a *AmbiReader) Close() error {
	if a.gz != nil {
		if err := a.gz.Close(); err != nil {
			return err
		}
	}
	if err := a.fp.Close(); err != nil {
		return err
	}
	return nil
}

type PairedEndReader struct {
	Records  int
	fn1      string
	fn2      string
	mate1    *AmbiReader
	mate2    *AmbiReader
	scanner1 *bufio.Scanner
	scanner2 *bufio.Scanner
	lineNum  int
}

func (r *PairedEndReader) Open(fn1, fn2 string) error {
	r.fn1 = fn1
	r.fn2 = fn2
	inputs := make([]AmbiReader, 2)
	r.mate1 = &inputs[0]
	r.mate2 = &inputs[1]
	if err := r.mate1.Open(fn1); err != nil {
		log.Fatalf("Failed to open mate1 %s: %v\n", fn1, err)
	}
	if err := r.mate2.Open(fn2); err != nil {
		log.Fatalf("Failed to open mate1 %s: %v\n", fn2, err)
	}
	r.scanner1 = bufio.NewScanner(inputs[0])
	r.scanner2 = bufio.NewScanner(inputs[1])
	return nil
}

func (r *PairedEndReader) Close() error {
	if err := r.mate1.Close(); err != nil {
		return err
	}
	if err := r.mate2.Close(); err != nil {
		return err
	}
	return nil
}

// Returns a slice of four strings (mate 1 name, seq 1, mate 2 name, seq 2) and
// no error if it reads a fastq entry. On EOF it returns an empty slice.
func (r *PairedEndReader) Read() ([]string, error) {
	record := make([]string, 4)
	var leftMate, rightMate string
	for i := 0; i < 4; i++ {
		if r.scanner1.Scan() {
			leftMate = r.scanner1.Text()
			if r.scanner2.Scan() {
				rightMate = r.scanner2.Text()
			} else {
				return nil, fmt.Errorf("file %s truncated at line %d", r.fn2, r.lineNum+1)
			}
			if i == 0 {
				if !strings.HasPrefix(leftMate, "@") {
					return nil, fmt.Errorf("expecting %s line %d to start with '@'", r.fn1, r.lineNum+1)
				}
				if !strings.HasPrefix(rightMate, "@") {
					return nil, fmt.Errorf("expecting %s line %d to start with '@'", r.fn2, r.lineNum+1)
				}
				// Store the read names, stripping off the '@' character
				record[0] = leftMate[1:]
				record[2] = rightMate[1:]
			}
			if i == 1 {
				record[1] = leftMate
				record[3] = rightMate
				r.Records++
			}
			if i == 2 {
				if !strings.HasPrefix(leftMate, "+") {
					return nil, fmt.Errorf("expecting %s line %d to start with '+'", r.fn1, r.lineNum+1)
				}
				if !strings.HasPrefix(rightMate, "+") {
					return nil, fmt.Errorf("expecting %s line %d to start with '+'", r.fn2, r.lineNum+1)
				}
			}
		} else {
			return nil, nil
		}
		r.lineNum++
	}
	return record, nil
}

func init() {
	log.SetFlags(0)
	flag.StringVar(&args.Sample, "sample", "", "sample ID for this sample (required)")
	flag.StringVar(&args.FastqList, "files", "", "file that contains the list of fastq files (required)")
	flag.IntVar(&args.Limit, "limit", 0, "only consider the first LIMIT fastq records in each sample")
	flag.IntVar(&args.Batches, "batches", 1, "process files in batches to avoid open file limits")
	flag.StringVar(&args.Progress, "progress", "", "write data after each batch to this file")
	flag.StringVar(&args.Continue, "continue", "", "file with output from an existing run we'll add to")

	flag.Usage = func() {
		log.Println("usage: fqmultioverlap [options]")
		flag.PrintDefaults()
	}
}

func scanFastQ(sampleId, fn1, fn2 string, ch chan string, wg *sync.WaitGroup) {
	defer wg.Done()

	sample := PairedEndReader{}
	err := sample.Open(fn1, fn2)
	if err != nil {
		log.Fatalf("Failed to open fastq files for sample %s: %v", sampleId, err)
	}
	defer sample.Close()

	for {
		record, err := sample.Read()
		if err != nil {
			log.Fatal("Failed reading from sample %s fastq at record %d: %v", sampleId, sample.Records, err)
		}
		if len(record) == 0 {
			break
		}
		if len(record) != 4 {
			log.Fatalf("record should have 4 fields, got %d", len(record))
		}
		key := record[1] + ":" + record[3]
		_, present := sampleSequences[key]
		if present {
			serialized := key + "@" + sampleId
			ch <- serialized
		}
		if args.Limit > 0 && sample.Records >= args.Limit {
			return
		}
	}
}

func recordSamples(hits chan string, done chan int) {
	// Loop until our channel is closed
	numHits := 0
	for serialized := range hits {
		tuple := strings.Split(serialized, "@")
		if len(tuple) != 2 {
			log.Fatalf("received tuple of length %d from channel", len(tuple))
		}
		key := tuple[0]
		sampleId := tuple[1]
		refseq[key][sampleId]++
		numHits++
	}
	done <- numHits
}

func writeOutput(fp io.Writer) int {
	sharedReads := 0
	for key, sampleSet := range refseq {
		if len(sampleSet) > 0 {
			pieces := strings.Split(key, ":")
			readName := sampleSequences[key]
			list := make([]string, 0)
			for sampleId, _ := range sampleSet {
				list = append(list, sampleId)
			}
			fmt.Fprintf(fp, "%s\t%s\t%s\t%s\n", readName, pieces[0], pieces[1], strings.Join(list, ","))
			sharedReads++
		}
	}
	return sharedReads
}

func writeOverlapHeaders(fastqFiles [][]string) {
	for _, fields := range fastqFiles {
		fmt.Printf("# overlap\t%s\t%s\t%s\n", fields[0], fields[1], fields[2])
	}
}

func logArgs() {
	log.Println("command:", strings.Join(os.Args, " "))

	bytes, err := json.Marshal(args)
	if err != nil {
		log.Fatal("can't serialize Args")
	}
	log.Println(string(bytes))
}

func main() {
	flag.Parse()

	if args.FastqList == "" || args.Sample == "" {
		log.Println("arguments -files, and -sample are required")
		flag.Usage()
		os.Exit(1)
	}

	logArgs()

	// Read through file given with the -files argument
	fp, err := os.Open(args.FastqList)
	if err != nil {
		log.Fatalf("Failed to open list of fastq files %s: %v", args.FastqList, err)
	}
	listScanner := bufio.NewScanner(fp)
	lineNum := 0
	foundSelf := false
	// In the fastqFiles slice we store the split lines of the contents of the -files file list.
	// We store a map of these so we can easily identify which samples have already been
	// processed earlier in the -continue file.
	var fastqFiles [][]string
	fastqFilesIndex := make(map[string]bool)
	for listScanner.Scan() {
		line := listScanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) != 3 {
			log.Fatalf("malformed line %d in %s: %s", lineNum+1, args.FastqList, line)
		}
		lineNum++
		if fields[0] == args.Sample {
			args.Ref1 = fields[1]
			args.Ref2 = fields[2]
			foundSelf = true
			continue
		}
		// We don't expect to get the same sample a second time
		_, ok := fastqFilesIndex[fields[0]]
		if ok {
			log.Fatalf("already saw sample %s in files list", fields[0])
		}
		fastqFilesIndex[fields[0]] = true
		fastqFiles = append(fastqFiles, fields)
	}
	fp.Close()
	originalFastqCount := len(fastqFiles)

	if foundSelf {
		log.Println("found self in file list and read files match ref1 and ref2")
	} else {
		log.Fatalf("failed to find self, %s, in list of fastq files", args.Sample)
	}

	fmt.Printf("# sample\t%s\n", args.Sample)
	fmt.Printf("# ref1\t%s\n", args.Ref1)
	fmt.Printf("# ref2\t%s\n", args.Ref2)

	// Create our global refseq data structure.
	refseq = make(map[string]map[string]int, 0)
	sampleSequences = make(map[string]string, 0)

	wroteOverlapHeaders := false
	if args.Continue != "" {
		log.Println("Continuing from", args.Continue)
		// This file should be uncompressed. If file is compressed, can come
		// from a pipe, as this only does one pass to read in the data.
		fp, err := os.Open(args.Continue)
		if err != nil {
			log.Fatalf("Failed to open continue file %s: %v", args.Continue, err)
		}
		contScanner := bufio.NewScanner(fp)
		// The first two lines should give the file names used for ref1 and
		// ref2 of the previous run. These must match exactly.
		lineNum = 0
		foundRef1 := false
		foundRef2 := false
		foundSample := false
		samplesSeen := make(map[string]bool)
		for contScanner.Scan() {
			line := contScanner.Text()
			if strings.HasPrefix(line, "# ") {
				if strings.HasPrefix(line, "# sample") {
					sample := line[9:]
					if sample != args.Sample {
						log.Fatalf("continue file is for sample %s instead of sample %s", sample, args.Sample)
					}
					foundSample = true
				} else if strings.HasPrefix(line, "# ref1") {
					fn := line[7:]
					if fn != args.Ref1 {
						log.Fatal("ref1 in continue file %s is %s, expecting %s", args.Continue, fn, args.Ref1)
					}
					foundRef1 = true
				} else if strings.HasPrefix(line, "# ref2") {
					fn := line[7:]
					if fn != args.Ref2 {
						log.Fatal("ref2 in continue file %s is %s, expecting %s", args.Continue, fn, args.Ref2)
					}
					foundRef2 = true
				} else if strings.HasPrefix(line, "# overlap") {
					rest := line[10:]
					fields := strings.Split(rest, "\t")
					if len(fields) != 3 {
						log.Fatalf("malformed '# overlap' line (%d): %s", lineNum+1, line)
					}
					samplesSeen[fields[0]] = true
					// Ensure this sample isn't in the list of ones we plan to process, but give it an overlap header
					// indicating it's already been processed.
					delete(fastqFilesIndex, fields[0])
					fmt.Println(line)
				} else {
					log.Fatalf("Unrecognized comment line (%d) in continue file: %s", lineNum+1, line)
				}
				continue
			}
			if !wroteOverlapHeaders {
				// We create a new ordered subset of fastqFiles that have not yet been removed from `fastqFilesIndex`.
				// We've already removed all the ones we saw in the -continue file so the remaining ones are new ones
				// we wish to add.
				filteredFastqList := make([][]string, 0)
				for _, row := range fastqFiles {
					if _, present := fastqFilesIndex[row[0]]; present {
						// only include files we didn't see
						filteredFastqList = append(filteredFastqList, row)
					}
				}
				fastqFiles = filteredFastqList
				writeOverlapHeaders(fastqFiles)
				wroteOverlapHeaders = true
			}
			fields := strings.Split(line, "\t")

			if len(fields) != 4 {
				log.Fatalf("malformed line %d in %s: %s", lineNum, args.Continue, line)
			}
			// Sequences are keyed by the combined DNA sequences
			key := strings.Join(fields[1:3], ":")
			readName := fields[0]
			samples := strings.Split(fields[3], ",")
			// We should not see the same sequence twice in a continue file
			_, ok := refseq[key]
			if ok {
				log.Fatalf("Already saw key %s in continue file %s", key, args.Continue)
			}
			// Load in the sa samples associated with a key
			refseq[key] = make(map[string]int, 0)
			for _, sampleId := range samples {
				refseq[key][sampleId]++
			}
			sampleSequences[key] = readName
		}
		if !foundRef1 {
			log.Fatalf("Expecting continue file %s to have '# ref1' line", args.Continue)
		}
		if !foundRef2 {
			log.Fatalf("Expecting continue file %s to have '# ref2' line", args.Continue)
		}
		if !foundSample {
			log.Fatalf("Expecting continue file %s to have # sample' line", args.Continue)
		}
		fp.Close()
	} else {
		// If we haven't yet written the overlap headers, we do that now.
		writeOverlapHeaders(fastqFiles)
	}

	log.Println("Processing ref sequence")
	ref := PairedEndReader{}
	err = ref.Open(args.Ref1, args.Ref2)
	addedRef := 0
	foundRef := 0
	if err != nil {
		log.Fatalf("Failed to open reference fastq files: %v", err)
	}
	log.Println("Opened ref1", args.Ref1)
	log.Println("Opened ref2", args.Ref2)
	skipped := 0
	for {
		if (foundRef+addedRef)%1e06 == 0 && (foundRef+addedRef) > 0 {
			log.Println(foundRef, addedRef)
		}
		record, err := ref.Read()
		if err != nil {
			log.Fatal("Failed reading reference fastq: %v", err)
		}
		if len(record) == 0 {
			break
		}
		if len(record) != 4 {
			log.Fatalf("ref record should have 4 fields, got %d", len(record))
		}
		key := record[1] + ":" + record[3]
		// We may already have an entry from the continue file, so only create
		// a new map if this is a read that was not previously shared.
		_, ok := refseq[key]
		if ok {
			if sampleSequences[key] != record[0] {
				if skipped < 10 {
					log.Printf("existing entry for %s has name %s, which is different from %s\n",
						key, sampleSequences[key], record[0])
				}
				skipped++
				continue
			}
			foundRef++
		} else {
			refseq[key] = make(map[string]int, 0)
			sampleSequences[key] = record[0]
			addedRef++
		}
		if args.Limit > 0 && ref.Records >= args.Limit {
			log.Println("Warning: reached refseq limit")
			break
		}
	}
	ref.Close()
	log.Printf("Done processing ref sequence, %d cached from continue, added %d and skipped %d\n",
		foundRef, addedRef, skipped)

	// Make a channel that all our fastq goroutines will write to
	hits := make(chan string)

	// Make another channel that we'll wait on to detect that our reader has finished.
	done := make(chan int)

	// Lanuch a goroutine to read from the channel. This goroutine will be done
	// when the channel is closed.
	go recordSamples(hits, done)

	log.Printf("Will read %d files of %d listed in %s in %d batches\n",
		len(fastqFiles), originalFastqCount, args.FastqList, args.Batches)
	for b := 0; b < args.Batches; b++ {
		thisBatch := 0

		// Wait until all our goroutines in this batch are done
		var wg sync.WaitGroup
		for i, tuple := range fastqFiles {
			if i%args.Batches == b {
				go scanFastQ(tuple[0], tuple[1], tuple[2], hits, &wg)
				wg.Add(1)
				thisBatch++
			}
		}
		// Wait for all the goroutines to be done after which we close the channel.
		log.Printf("Processing %d samples in batch %d\n", thisBatch, b)
		wg.Wait()
		// Write intermediate progress, unless we're on the last batch
		if args.Progress != "" && b != args.Batches-1 {
			fp, err := os.Create(args.Progress)
			if err == nil {
				log.Println("writing intermediate progress to", args.Progress)
				_ = writeOutput(fp)
				fp.Close()
			} else {
				log.Println("can't write to", args.Progress, "skipping")
			}
		}
	}

	close(hits)

	// Now wait for our reader to be done.
	numHits := <-done

	log.Println("Writing output")
	sharedReads := writeOutput(os.Stdout)

	log.Println("Got", sharedReads, "shared reads with", numHits, "sharing events in aggregate")
}

// END
