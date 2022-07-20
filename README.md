# Shepherd

## Getting Started

Shepherd is a Python program for correcting substitution errors and single insertion and deletion errors in DNA barcode reads. These errors occur during PCR amplification and sequencing of the DNA barcodes. Shepherd is cross-platform and runs on any computer with Python 3.8 or later and the Scipy library. 

The program consists of two python scripts: shepherd_t0.py and shepherd_multi.py. These scripts are used via a command line interface described below.

### The shepherd_t0.py script

This script is designed to cluster the sequencing reads from a single time point to correct substitution errors and single insertion and deletion errors. 

### Inputs

#### Required Inputs

These inputs must be provided to run the script.

**-l:** (integer) The correct barcode length.

**-f:** (.txt file) The input file with a sequence and a sequence count in each row, separated by whitespace. Currently this is the only input file format supported by Shepherd. Only sequences in the file with lengths l (correct barcode length), l + 1 (single insertion errors) and l - 1 (single deletion errors) will be processed by Shepherd.

    Example file:   testdata_t0.txt
    
                    TCCCTTACTAATCGAAGAAG	5
                    ATAGTATGGATCTGGACCGC	10	
                    ATCCAGTGCTAGTTCAACTC	3
                    AATTTTGGAACAGGCCGTAG	200
    
#### Optional Inputs

These inputs are optional and we recommend using the default values determined by Shepherd.

**-e:** (float) An estimate of the substitution error rate of the sequencing protocol used to generate the input data. This is a floating point number, e.g. 0.01 if the estimated error rate is 1%. If not provided this parameter is automatically determined based on the input data.

**-eps:** (integer) The maximum Hamming distance considered for merging two sequencing into the same cluster. If not provided this parameter is automatically determined based on the input data.

**-k:** (integer) The substring length used to divide the sequences into partitions. If not provided this parameter is automatically determined based on the input data.

**-tau:** (integer) A distance threshold for frequency 1 sequences that determines if they should be merged with another sequence. If a frequency 1 sequence has Hamming distance less than or equal to this threshold to a candidate sequence it will be merged. If not provided this parameter is automatically determined based on the input data.

**-ft:** (integer) A frequency threshold for defining true barcodes. Any sequence with a frequency higher than this threshold is defined as a true barcode. If not provided this parameter is automatically determined based on the input data.

**-bft:** (float) The threshold for log Bayes factor. The default value is -4.

**-Nh:** (integer) Number of sequences used for estimation of the substitution error rate. The default value is min(number of sequences, 500).

### Outputs

**_seq_clust.csv:** A .csv file where the unique sequences are in the first column and the cluster labels are in the second column.

**_pb_freq.csv:** A .csv file where the putative barcodes are in the first column and the estimated counts are in the second column.

**_index:** The k-mer Index stored in the pickle format.
                                                                                             
**_params:** The parameters used to run the script stored in the pickle format.

### Usage

**Command line usage example:** <code>python3 shepherd_t0.py -f testdata_t0.txt -l 20 -e 0.01</code>

### The shepherd_multi.py script

This script is designed to use the the clustering from the first time point, i.e., the outputs of shepherd_t0.py, to estimate the counts of the putative barcodes at later time points, given the sequencing reads from each time point. If new barcodes that did not appear in the first time point emerge in later time points, the program is capable of identifying and tracking them. Note that the shepherd_t0.py script must be executed in the same folder prior to running shepherd_multi.py. 

### Inputs

**-f0:** (.txt file) The same input file used to run the shepherd_t0.py script containing the sequences and the sequence counts.                                             
**-fn:** (.txt files) Space separated list of .txt files containing the sequences and sequence counts for each time point. These files should have the same format as the input file to shepherd_t0.py (see testdata_t0.txt) and should be ordered by time point (see usage example below).\
**-o:** (string) The prefix of the final output file. By default set to 'multi_freqs' which produces an output file called 'multi_freqs.csv'.

### Outputs

**multi_freqs.csv:** A .csv file where each row is a putative barcode and the columns give the estimated counts for each time point.
**seq_clust.csv:** A .csv file for each time point where each row contains a sequence and the cluster ID it was assigned.

### Usage

**Command line usage example:**\
<code>python3 shepherd_multi.py -f0 testdata_t0.txt -fn testdata_t1.txt testdata_t2.txt</code>
