# k-seek
Version 4.
de novo identification and quantification of satellite DNA from illumina whole genome sequence data.

k.seek.pl identifies and quantifies satellite DNAs (tandemly repeating simple sequences) where the base motif is <= 20bp and generates two outputs .rep and .total. The former contains the the fastq reads from which the repeat was called with an additional 5th line documenting the the called repeat and the number of occurences. The latter contains the summation for all the called repeats.
Usage:
perl k.seek.pl input.fastq output_basename

k.compiler.pl will concatenate multiple .total files in a folder into a table where each row is a sample and each column a repeat. Different phases and reverse compliments of the same repeat (e.g AAC, ACA, CAA, TTG, etc) will be summed. The .toal files need to be put into a folder without other files. The output will be .rep.compiled
Usage:
perl k.compile.pl folder/ output_basename


For reference see:
https://www.ncbi.nlm.nih.gov/pubmed/25512552
