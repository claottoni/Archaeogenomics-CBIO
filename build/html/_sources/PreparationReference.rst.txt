#####################################
Preparation of the reference sequence
#####################################

After quality-control and pre-processing, we will align the merged and trimmed reads to the *Felis silvestris* mtDNA reference sequence, available in the RefSeq NCBI database (https://www.ncbi.nlm.nih.gov/). The first thing to do is the prepare the reference sequence for the alignment procedure. 


*************************************
Index the reference sequence with bwa
*************************************

For the alignment of reads to the reference sequence we will use BWA, in particular the BWA-aln algorithm. BWA first needs to construct the FM-index for the reference genome, with the command *bwa index*. FM-indexing in Burrows-Wheeler transform is used to efficiently find the number of occurrences of a pattern within a compressed text, as well as locate the position of each occurrence. It is in essential step for querying the DNA reads to the reference sequence. This command generates five files with different extensions.
::
     
  bwa index -a is reference.fasta
     
.. warning::
  
  The option -a indicates the algorithm to use for constructing the index. For genomes smaller than < 2 Gb use the **is** algorithm. For larger genomes (>2 Gb), use the **bwtsw** algorithm. 	

*****************************
Create a reference dictionary
*****************************

You need to do this command line in order to run later in the pipline the **gatk RealignerTargetCreator**. A sequence dictionary contains the sequence name, sequence length, genome assembly identifier, and other information about sequences.
::

  picard CreateSequenceDictionary R= referece.fasta O= reference.dict
 
.. note::

  In some environments we can call picard just by typing the program name. In other environments (including your laptop) you may have to call picard by providing the full path to the java file (.jar) of picard:
  ::
  
    java -jar path_to_picard.jar CreateSequenceDictionary R= referece.fasta O= ref.dict

******************************************
Index the reference sequence with samtools
******************************************

You need to do this command line in order to run later in the pipline the **gatk IndelRealigner**. We use **samtools faidx**, which enables efficient access to arbitrary regions within the reference sequence. The index file typically has the same filename as the corresponding FASTA file, with .fai appended.
::

  samtools faidx reference.fasta

