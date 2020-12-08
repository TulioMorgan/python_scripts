## Python Scripts

### Repository containing some of my Python scripts for a range of Bioinformatic tasks.

  - Data science
      - multidimensional scaling and PCA
      
  - Biological sequence manipulation
  
      - Sequence alignments (align primers to target sequences, merge alignment files)
        * __alignPrimers.py__ (performs alignments between primers/probe and target seqeunces. Uses Exonerate v.2.2.0)
        * __concatAlignments.py__ (concatenate multiple sequence alignment files (fasta) for multi-gene phylogeny)
        * __runExonerate.v2.py__ (Script to run to run Exonerate in parallel. I use it mainly to predict genes in new genome assemblies. Uses Exonerate v.2.2.0)
        * __pairwiseGlobalAlignment_sequenceIdentities.py__ (Script to perform global alignments. The input consists of a multifasta file. Uses EMBOSS (needle))
        
      - Fasta files
        * __dropSequences.py__ (drop nucleotide sequences that bear x % of undetermined or ambiguous characters (N, R, W, ...))
        * __extraiSeqFastaIDs.py__ (report sequences according to a list of IDs (provided in a txt file))
        * __extractseq.py__ (report a sequence or a subsequence (search by sequence ID))
        * __findStopCodons.py__ (report sequences without stop codons in-frame (or report all sequences but only until reaches the first stop codon in each one))
        * __formataFasta.py__ (remove line breaks within each sequence in a multifasta file)
        * __getSeqFromCoord.py__ (extract fasta sequences from a GFF3 file or a generic coordinate file)
        * __orderSequencesBySize.py__ (Script to order sequences in a multifasta file according to their sizes, from highest to lowest (default))
        * __removeSmallSequences.py__ (Script to remove sequences with length less than a user-defined value (--comprimento integer))
        * __seqsAleatorias.py__ (retrieve n random sequences from a multifasta file)
        * __splitFasta.py__ (randomly divides a multifasta file into independent fasta files)
        
  - Genomics
      - genomic features
        * __genomeMetricsExonsIntrons.py__ (Analyzes the distribution of size of introns/exons and the number of introns/exons in a GFF3 file.)
      - protein function prediction
        * __blastpAnnotation.py__ (analyzes all BLASTp/DIAMOND subject descriptions and retrieve the most likely annotation for each query (majority consensus))
        * __interpro2Annotation.py__ (analyzes and filter the output of InterProScan. Print the protein identifier followed by all domain descriptions and coordinates)
