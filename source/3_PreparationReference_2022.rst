#####################################
Preparation of the reference sequence
#####################################

After quality-control and pre-processing of the reads (trimming of adapter sequences and merging of paired reads), we will align the merged and trimmed reads to the reference sequence of interest, available in the `RefSeq NCBI`_ database. The first thing to do is the prepare the reference sequence for the alignment procedure. 

  .. _RefSeq NCBI: https://www.ncbi.nlm.nih.gov/refseq/

*************************************
Index the reference sequence with BWA
*************************************

To align the reads to the reference sequence we will use the program `BWA`_, in particular the ``BWA aln`` algorithm. BWA first needs to construct the **FM-index** for the reference genome, with the command ``BWA index``. FM-indexing in Burrows-Wheeler transform is used to efficiently find the number of occurrences of a pattern within a compressed text, as well as locate the position of each occurrence. It is an essential step for querying the DNA reads to the reference sequence. This command generates five files with different extensions: ``amb``, ``ann``, ``bwt``, ``pac``, ``sa``.

  .. _BWA: https://github.com/lh3/bwa

::
     
  bwa index -a is reference.fasta
     
.. note::
  
  The option ``-a`` indicates the algorithm to use for constructing the index. For genomes smaller than < 2 Gb use the ``is`` algorithm. For larger genomes (>2 Gb), use the ``bwtsw`` algorithm. 	

*****************************
Create a reference dictionary
*****************************

A dictionary file (``dict``) is necessary to run later in the pipeline ``GATK RealignerTargetCreator``. A sequence dictionary contains the sequence name, sequence length, genome assembly identifier, and other information about the sequences. To create the ``dict`` file we use `Picard`_. 

  .. _Picard: https://broadinstitute.github.io/picard/

::

  picard CreateSequenceDictionary R= referece.fasta O= reference.dict
 
.. note:: In our server environment we can call Picard just by typing the program name. In other environments (including your laptop) you may have to call Picard by providing the full path to the java file ``jar`` of the program:
   ::
  
     java -jar /path/to/picard.jar CreateSequenceDictionary R= referece.fasta O= ref.dict

******************************************
Index the reference sequence with Samtools
******************************************

The reference sequence has to be indexed in order to run later in the pipeline ``GATK IndelRealigner``. To do that, we will use `Samtools`_, in particular the tool ``samtools faidx``, which enables efficient access to arbitrary regions within the reference sequence. The index file typically has the same filename as the corresponding reference sequece, with the extension ``fai`` appended.

  .. _Samtools: http://www.htslib.org/

::

  samtools faidx reference.fasta
