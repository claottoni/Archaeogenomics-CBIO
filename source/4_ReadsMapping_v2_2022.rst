################################################
Alignment of the reads to the reference sequence
################################################

*********************************************************************
Alignment of pre-processed reads to the reference genome with BWA aln
*********************************************************************

To align the reads to the reference genome we will use ``BWA aln``, which supports an end-to-end alignment of reads to the reference sequence. The alternative algorithm, ``BWA mem`` supports also local (portion of the reads) and chimeric alignments (resulting in a larger number of mapped reads than ``BWA aln``). ``BWA aln`` is more suitable for aliging short reads, like expected for ancient DNA samples. The following comand will generate a ``sai`` file.
::

  bwa aln reference.fasta filename.fastq.gz -n 0.01 -l 1000 -o 2 > filename.sai

Some of the options available in ``BWA aln``: 

================ ========
Option           Function
================ ========
**-n** *number*  Maximum edit distance if the value is an *integer*. If the value is *float* the edit distance is automatically chosen for different read lengths (default=0.04)
**-l** *integer* Seed length. If the value is larger than the query sequence, seeding will be disabled. 
**-o** *integer* Maximum number of gap opens. For aDNA, tolerating more gaps helps mapping more reads (default=1).
================ ========

.. note::

  - Due to the particular damaged nature of ancient DNA molecules, carrying deaminations at the molecules ends, we deactivate the ``seed-length`` option by giving it a high value (e.g. ``-l 1000``). 
  - To get more stringent mapping conditions (for examples when mapping the reads to a bacterial reference genome), we can increase the float ``edit-distance`` value to ``-n 0.1`` 
  
  
Once obtained the ``sai`` file, we align the reads (``fastq`` file) to the reference (``fasta`` file) using ``BWA samse``, to generate the alignment file ``sam``.
::

  bwa samse reference.fasta filename.sai filename.fastq.gz -f filename.sam


*******************************
Converting sam file to bam file
*******************************

For the downstream analyses we will work with the binary (more compact) version of the ``sam`` file, called ``bam``. To convert the ``sam`` file in ``bam`` we will use ``Samtools view``. 
::

  samtools view -Sb filename.sam > filename.bam

.. note::

  The conversion from ``sam`` to ``bam`` can be piped (``|``) in one command right after the alignment step:
  ::

    bwa samse reference.fasta filename.sai filename.fastq.gz | samtools view -Sb - > filename.bam

To view the content of a ``sam`` file we can just use standard commands like ``head``, ``tail``, ``less``, while to view the content of a ``bam`` file (binary format of ``sam``) we have to use ``Samtools view``:
::

  samtools view filename.bam
  
You may want to display on the screen one read/line (scrolling with the spacebar):
::

  samtools view filename.bam | less -S

while to display just the header of the ``bam`` file:
:: 

  samtools view -H filename.bam

*********************************
Sorting and indexing the bam file
*********************************

To go on with the analysis, we have to sort the reads aligned in the ``bam`` file by leftmost coordinates (or by read name when the option ``-n`` is used) with ``Samtools sort``. The option ``-o`` is used to provide an output file name:
::

  samtools sort filename.sam -o filename.sort.bam

The sorted bam files are then indexed with ``Samtools index``. Indexes allow other programs to retrieve specific parts of the ``bam`` file without reading through each sequence. The following command generates a ``bai`` file, a companion file of the ``bam`` which contains the indexes:
::

  samtools index filename.sort.bam

*********************************************
Adding Read Group tags and indexing bam files
*********************************************

A number of predefined tags may be appropriately assigned to specific set of reads in order to distinguish samples, libraries and other technical features. To do that we will use ``Picard``. You may want to use ``RGLB`` (library ID) and ``RGSM`` (sample ID) tags at your own convenience based on the experimental design. Remember to call ``Picard`` from the path of the ``jar`` file, here: ``/home/aurochs/Software/picard/picard.jar``
::

 picard AddOrReplaceReadGroups INPUT=filename.sort.bam OUTPUT=filename.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT

.. note::

  In some environments we can call ``Picard`` just by typing the program name. In other environments (including this server) you may have to call Picard by providing the full path to the java file (``jar``) of the program. Here, the path is: ``java -jar /home/aurochs/Software/picard/picard.jar``

  ::
   
    java -jar /home/aurochs/Software/picard/picard.jar CreateSequenceDictionary R= referece.fasta O= ref.dict

.. note::
  
  - In some instances ``Picard`` may stop running and return error messages due to conflicts with ``sam`` specifications produced by ``BWA`` (e.g. "MAPQ should be 0 for unmapped reads"). To suppress this error and allow ``Picard`` to continue, we pass the ``VALIDATION_STRINGENCY=LENIENT`` options (default is ``STRICT``).
  - Read Groups may be also added during the alignment with ``BWA`` using the option ``-R``. 


.. warning::
  
  If the analysis is stopped due to a "Java Runtime Environment" try to add to the command the options ``USE_JDK_DEFLATER=true USE_JDK_INFLATER=true``

Once added the Read Group tags, we index again the ``bam`` file:
:: 

  samtools index filename.RG.bam

*******************
Removing duplicates
*******************

Amplification through PCR of genomic libraries leads to duplication formation (reads originating from a single fragment of DNA). 
To remove duplicates we will use DeDup, which is specifically designed for ultra-short DNA (e.g. ancient DNA) paired-end and merged sequence data. 
The issue with ultra-short DNA is that during the typical sequencing chemistry cycles the entire molecule will be sequenced and therefore we will find both ends.
Compared to typical deduplication tools, that only look for reads with the same starting position at the 5' end of a read, DeDup will remove 'true' duplicates using BOTH (5' and 3') ends of the reads. This can help increase coverage in low-preservation samples such used in ancient DNA by being more exact as to what are duplicates or not.
DeDup generates a ``bam`` file with the suffix ``_rmdup.bam``
::

  dedup -i filename.RG.bam -o folder_name -m -u

=================================== ========
Option                              Function
=================================== ========
**-h**                  			Shows the help page
**-i** *string*                		Select your input, otherwise you may use pipes to pipe in your data
**-m**             					The input only contains merged reads - don't care about missing prefixes for merged/reverse/forward reads
**-o** *string*           			output folder
**-u**								Do not automatically sort the output
=================================== ========

.. note::

  The output can be put in a separate folder (``-o folder_name``), or more simply if you are working in the folder containing the input ``bam``, use ``-o .`` to have the output in that same location.  


.. note::

  Alternatively, you can use the ``MarkDuplicates`` tool of Picard, which marks the reads as duplicates when the 5'-end positions of both reads and read-pairs match. A metric file with various statistics is created, and reads are removed from the bam file by using the ``REMOVE_DUPLICATES=True`` option (the default option is ``False``, which simply 'marks' duplicate reads keep them in the ``bam`` file).
  :: 
  
    java -jar /home/aurochs/Software/picard/picard.jar MarkDuplicates I= filename.RG.bam O= filename.DR.bam M=output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT

Once removed the duplicates with DeDup, we sort the reads and index of the bam output file (``filename_rmdup.bam``):
::
  
  samtools sort filename_rmdup.bam -o filename_rmdup_sort.bam
  samtools index filename_rmdup_sort.bam


**************************
Local realignment of reads
**************************

The presence of insertions or deletions (indels) in the genome may be responsible of misalignments and bases mismatches that are easily mistaken as SNPs. For this reason, we locally realign reads to minimize the number of mispatches around the indels. The realignment process is done in two steps using two different tools of the program ``GATK``. These tools are called with the ``-T`` option. We first detect the intervals which need to be realigned with the ``RealignerTargetCreator``, and save the list of these intervals in a file that we name ``target.intervals``. Like Picard, we have to call ``GATK`` with the full path to the ``jar`` file:
::

  java -jar /home/aurochs/Software/gatkv3/GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I filename_rmdup_sort.bam -o target.intervals
 
.. warning::
  
  In  *version 4* of ``GATK`` the indel realigment tools have been retired from the best practices (they are unnecessary if you are using an assembly based caller like **Mutect2** or **HaplotypeCaller**). To use the indel realignment tools make sure to install *version 3* of ``GATK``.  

Then, we realign the reads over the intervals listed in the ``target.intervals`` file with the option ``-targetIntervals`` of the tool ``IndelRealigner`` in ``GATK``:
::

  java -jar /home/aurochs/Software/gatkv3/GenomeAnalysisTK.jar -T IndelRealigner -R reference.fasta -I filename_rmdup_sort.bam -targetIntervals target.intervals -o filename.final.bam --filter_bases_not_stored 

  
.. note::

  - If you want, you can redirect the standard output of the command into a ``log`` file by typing at the end of the command ``&> logFile.log`` 
  - The option ``--filter_bases_not_stored`` is used to filter out reads with no stored bases (i.e. with * where the sequence should be), instead of failing with an error

The final ``bam`` file has to be sorted and indexed as previously done:
::

  samtools sort filename.final.bam -o filename.final.sort.bam
  samtools index filename.final.sort.bam


**********************
Generate flagstat file
**********************

We can generate a file with useful information about our alignment with ``Samtools flagstat``. This file is a final summary report of the bitwise ``FLAG`` fields assigned to the reads in the ``sam`` file.
::

  samtools flagstat filename.final.sort.bam > flagstat_filename.txt

.. note::

  - You could generate a flagstat file for the two ``bam`` files before and after refinement and see the differences. 
  - You can decode each ``FLAG`` field assigned to a read on the `Broad Institute`_ website.
   
      .. _Broad Institute: https://broadinstitute.github.io/picard/explain-flags.html

********************************
Visualization of reads alignment
********************************

Once generated the final ``bam`` file,  you can compare the ``bam`` files before and after the refinement and polishing process (duplicates removal, realignment around indels and sorting). To do so, we will use the program ``IGV``, in which we will first load the reference fasta file from *Genomes --> Load genome from file* and then we will add one (or more) bam files with *File --> Load from file*:

.. image:: images/igv-bam_bam.png
