# YAMS16
Simple 16S pipeline running DADA2

YAMS16 requires R to be installed and will install the libraries that
it needs if they aren't already. The db directory containing the SILVA
16S db should be kept in the same directory as the YAMS16.R script,
which can be a different directory from where you are running it from.

Running YAMS16.R by itself gives the following info.

Usage: ./YAMS16.R [options]


Options:
	-a ADAPTER, --adapter=ADAPTER
		file of Illumina adapters (default: ./db/adapters/adapter.fasta)

	-i INPUT, --input=INPUT
		path of directory where fastq files are

	-p PREFIX, --prefix=PREFIX
		filename prefix name for project

	-o OUTPUT, --output=OUTPUT
		path where to output files (default .)

	-t THREADS, --threads=THREADS
		number of threads to use (default 1)

	-x MAXLENGTH, --maxlength=MAXLENGTH
		maximum amplicon length (default 240)

	-n MINLENGTH, --minlength=MINLENGTH
		minimum amplicon length (default 180)

	-h, --help
		Show this help message and exit

A typical usage would be something like

YAMS16.R -i fastq -p Project1 -t 16 

This would run YAMS16 on the fastqs in directory fastq, set the
project name to Project1, and use 16 threads.

When it finishes, it should generate several files

unit_table_LKT.tsv -- this is a table of LKT counts with columns being samples, rows being taxa
taxa_16S_cons_LKT.tsv -- this is a table of taxonomic info with rows being taxa, columns being rank
rep_set_LKT.fa -- this is a fasta file of the predicted amplicons from DADA2


