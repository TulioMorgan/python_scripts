## Python Scripts

### Repository containing most of my Python scripts for a range of Bioinformatic tasks.

  - Data science
      - multidimensional scaling and PCA
      
  - Biological sequence manipulation
  
      - Sequence alignments (align primers to target sequences, merge alignment files)
        * __alignPrimers.py__ (performs alignments between primers/probe and target seqeunces. Uses Exonerate v.2.2.0)
        * __concatAlignments.py__ (concatenate multiple sequence alignment files (fasta) for multi-gene phylogeny)
        
      - Fasta files
        * __dropSequences.py__ (drop nucleotide sequences that bear x % of undetermined or ambiguous characters (N, R, W, ...))
        * __extraiSeqFastaIDs.py__ (report sequences according to a list of IDs (provided in a txt file))
        * __extractseq.py__ (report a sequence or a subsequence (search by sequence ID))
        * __findStopCodons.py__ (report sequences without stop codons in-frame (or report all sequences but only until reaches the first stop codon in each one))
        * __formataFasta.py__ (remove line breaks within each sequence in a multifasta file)
        * __getSeqFromCoord.py__ (extract fasta sequences from a GFF3 file or a generic coordiante file)
        * __seqsAleatorias.py__ (retrieve n random sequences from a multifasta file)
        
        
        
