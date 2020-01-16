package main

import (
	"bufio"
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
	"sync"
)

/* This program takes a reference sample, and a file lising other fastq files
 * and outputs the sequences and the samples containing them for any of the
 * reference sample's sequences that are found in the other samples.
 *
 * The file supplied with the -file argument should have three tab-separated
 * columns giving sampleId, fastq1, fastq2.
 */

type Args struct {
	Ref1      string
	Ref2      string
	FastqList string
	Limit     int
	Progress  bool
}

var args = Args{}

// We keep track of sequences we're looking for. After the reference seqeunce
// is read in, this will not be written to again. And in particular, it's
// treated as read-only by the goroutines, and so simultaneous access is okay.
var sampleSequences map[string]bool

// This should only be accessed when setting up the initial set of keys and by
// the reader goroutine.
var refseq map[string]string

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

// Returns a slice of two strings and no error if it reads a fastq entry.
// On EOF it returns an empty slice.
func (r *PairedEndReader) Read() ([]string, error) {
	mates := make([]string, 2)
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
			}
			if i == 1 {
				mates[0] = leftMate
				mates[1] = rightMate
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
	return mates, nil
}

func init() {
	log.SetFlags(0)
	flag.StringVar(&args.Ref1, "ref1", "", "fastq file for read 1 of the reference sample")
	flag.StringVar(&args.Ref2, "ref2", "", "fastq file for read 2 of the reference sample")
	flag.StringVar(&args.FastqList, "files", "", "file that contains the list of other fastq files")
	flag.IntVar(&args.Limit, "limit", 0, "only consider the first LIMIT fastq records in each sample")
	flag.BoolVar(&args.Progress, "progress", false, "show progress")

	flag.Usage = func() {
		log.Println("usage: fqseq [options]")
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
		mates, err := sample.Read()
		if err != nil {
			log.Fatal("Failed reading from sample %s fastq at record %d: %v", sampleId, sample.Records, err)
		}
		if len(mates) == 0 {
			break
		}
		key := mates[0] + ":" + mates[1]
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
	numHits := 1
	for serialized := range hits {
		tuple := strings.Split(serialized, "@")
		if len(tuple) != 2 {
			log.Fatalf("received tuple of length %d from channel", len(tuple))
		}
		key := tuple[0]
		sampleId := tuple[1]
		sep := ""
		if len(refseq[key]) > 0 {
			sep = ","
		}
		refseq[key] = refseq[key] + sep + sampleId
		numHits++
	}
	done <- numHits
}

func main() {
	flag.Parse()

	if args.Ref1 == "" || args.Ref2 == "" || args.FastqList == "" {
		log.Println("must give -ref1, -ref2, and -files arguments")
		flag.Usage()
		os.Exit(1)
	}

	fmt.Println("# ref1", args.Ref1)
	fmt.Println("# ref2", args.Ref2)

	// Read through file given with the -files argument
	fp, err := os.Open(args.FastqList)
	if err != nil {
		log.Fatalf("Failed to open list of fastq files %s: %v", args.FastqList, err)
	}
	listScanner := bufio.NewScanner(fp)
	lineNum := 0
	var fastqFiles [][]string
	for listScanner.Scan() {
		line := listScanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) != 3 {
			log.Fatalf("malformed line %d in %s: %s", lineNum+1, args.FastqList, line)
		}
		lineNum++
		fastqFiles = append(fastqFiles, fields)
	}
	fp.Close()

	log.Println("Processing ref sequence")
	ref := PairedEndReader{}
	err = ref.Open(args.Ref1, args.Ref2)
	if err != nil {
		log.Fatalf("Failed to open reference fastq files: %v", err)
	}
	refseq = make(map[string]string, 0)
	sampleSequences = make(map[string]bool, 0)
	for {
		mates, err := ref.Read()
		if err != nil {
			log.Fatal("Failed reading reference fastq: %v", err)
		}
		if len(mates) == 0 {
			break
		}
		key := mates[0] + ":" + mates[1]
		refseq[key] = ""
		sampleSequences[key] = true
		if args.Limit > 0 && ref.Records >= args.Limit {
			log.Println("Warning: reached refseq limit")
			break
		}
	}
	ref.Close()

	// Wait until all our goroutines are done
	var wg sync.WaitGroup
	wg.Add(len(fastqFiles))

	// Make a channel that all our fastq goroutines will write to
	hits := make(chan string)

	// Make another channel that we'll wait on to detect that our reader has finished.
	done := make(chan int)

	for _, tuple := range fastqFiles {
		go scanFastQ(tuple[0], tuple[1], tuple[2], hits, &wg)
	}

	// Lanuch a goroutine to read from the channel. This goroutine will be done
	// when the channel is closed.
	go recordSamples(hits, done)

	// Wait for all the goroutines to be done after which we close the channel.
	wg.Wait()
	close(hits)

	// Now wait for our reader to be done.
	numHits := <-done
	log.Println("Got", numHits, "hits")

	log.Println("Writing output")
	for key, list := range refseq {
		if list != "" {
			pieces := strings.Split(key, ":")
			fmt.Printf("%s\t%s\t%s\n", pieces[0], pieces[1], list)
		}
	}
}

// END
