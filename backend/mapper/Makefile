# Makefile for genome mapper
.PHONY: build clean

# Default test files
# REF = io/inputs/ref/short_reads_reference_genome.fasta
REF = io/inputs/ref/short_reads_ref_genome.fasta
FINALREF = io/inputs/final/short_reads_reference_genome.fasta
# READS1 = io/inputs/reads/challenging_dataset_1.fastq
# READS2 = io/inputs/reads/challenging_dataset_2.fastq
# TRUTH = io/inputs/readGroundTruth/challenging_ground_truth.txt

READS1 = io/inputs/reads/short_reads_1_2000_subset.fastq
READS2 = io/inputs/reads/short_reads_2_2000_subset.fastq
TRUTH = io/inputs/readGroundTruth/short_reads_2000_subset_ground_truth_1_base.txt
OUTPUT = io/outputs/SAMOutputFile.SAM

FINALREADS1 = io/inputs/final/challenging_dataset_1.fastq
FINALREADS2 = io/inputs/final/challenging_dataset_2.fastq
FINALTRUTH = io/inputs/final/challenging_ground_truth.txt
FINALOUTPUT = io/outputs/Final_SAMOutputFile.SAM

# Default run without memory tracking
run:
	touch $(OUTPUT)
	python3 main.py -r $(REF) -1 $(READS1) -2 $(READS2) --truth $(TRUTH) -o $(OUTPUT)

run-final:
	touch $(FINALOUTPUT)
	python3 main.py -r $(FINALREF) -1 $(FINALREADS1) -2 $(FINALREADS2) --truth $(FINALTRUTH) -o $(FINALOUTPUT)

# Run with memory tracking
run-mem:
	touch $(OUTPUT)
	python3 main.py -r $(REF) -1 $(READS1) -2 $(READS2) --truth $(TRUTH) -o $(OUTPUT) -m

# Run benchmark mode (no SAM output, with timing)
bench:
	python3 main.py -r $(REF) -1 $(READS1) -2 $(READS2) --truth $(TRUTH) --no-sam

# Run benchmark with memory tracking
bench-mem:
	python3 main.py -r $(REF) -1 $(READS1) -2 $(READS2) --truth $(TRUTH) --no-sam -m

# Run with custom parameters
run-custom:
	touch $(OUTPUT)
	python3 main.py -r $(REF) -1 $(READS1) -2 $(READS2) --truth $(TRUTH) -o $(OUTPUT) -k 15 -w 60

# Clean up
clean:
	find . -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name '*.pyc' -delete
	rm -rf io/outputs/*.SAM
	rm -f mapper

clean-cython:
	find . -type f -name '*.so' -delete
	find . -type f -name '*.c' -delete
	find . -type f -name '*.cpp' -delete
	find . -type f -name '*.o' -delete

# Run the Cython-compiled version
cython:
	@echo ">>> Building Cython extension..."
	python3 setup.py build_ext --inplace
# 	@echo ">>> Running Cythonized main..."
# 	touch $(OUTPUT)
# 	python3 main.py -r $(REF) -1 $(READS1) -2 $(READS2) --truth $(TRUTH) -o $(OUTPUT)

build:
	touch $(OUTPUT)
	@echo '#!/bin/bash' > mapper
	@echo 'python3 main.py "$$@"' >> mapper
	@chmod a+x mapper

	