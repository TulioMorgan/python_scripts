Python Scripts

Repository containing most of my Python scripts for a range of Bioinformatic tasks.

  - Data science
      - multidimensional scaling and PCA
      
  - Biological sequence manipulation
  
      - Sequence alignments (align primers to target sequences, merge alignment files)
        * alignPrimers.py
        * concatAlignments.py -> 
        
      - Fasta files
        * dropSequences.py -> [drop nucleotide sequences that bear x % of undetermined or ambiguous characters (N, R, W, ...)]
        * extraiSeqFastaIDs.py -> report sequences according to a list of IDs (provided in a txt file)
        * extractseq.py -> report a sequence or a subsequence (search by sequence ID)
        * findStopCodons.py -> report sequences without stop codons in-frame (or report all sequences but only until reaches the first stop codon in each one)
        * formataFasta.py -> remove line breaks within each sequence in a multifasta file
        * seqsAleatorias.py -> retrieve n random sequences from a multifasta file
        
        
        
