################################################
Alignment of the reads to the reference sequence
################################################

*********************************************
Alignment of pre-processed reads with BWA-aln
*********************************************

We are going to use the *BWA aln* algorithm. BWA-aln supports an end-to-end alignemnt of reads to the reference sequence, whereas the alternative algorithm, BWA-mem supports also local (portion of the reads) and chimeric alignments (then producing a larger number of reads mapped than aln). BWA-aln is more suitable for aliging short reads, like expected for ancient DNA samples. The following comand will generate a .sai file for downstream analysis and a .log file.
::

  bwa aln reference.fasta filename.fastq.gz -n 0.1 -l 1000 > filename.sai

Some of the options available in BWA-aln: 

============== ========
Option         Function
============== ========
*-n* <number>  maximum edit distance if the value is an *integer*. If the value is *float* the edit distance is automatically chosen for different read lengths. 
*-l* <integer> seed length. If the value is larger than the query sequence, seeding will be disabled. 
============== ========

.. note::

  - Due to the particular damaged nature of ancient DNA molecules, carrying deaminations and the molecules ends, we deactivate the **seed-length** option (-l) by giving it a high value (e.g. 1000). 
  - Here we are aligning reads to a bacterial reference genome. To reduce the impact of spurious alignemnts due to presence bacterial species closely related to the one that we are investigating, we will adopt stringent conditions by decreasing the **maximum edit distance** option (-n 0.1). For alignment of DNA reads to the human reference sequence, less stringent conditions can be used (-n 0.01). 

Once obtained the sai file, we align the reads (fastq file) to the reference (fasta file) using **bwa samse**, to generate the alignment file (sam) and a log file.
::

  bwa samse reference.fasta filename.sai filename.fastq.gz -f filename.sam


*******************************
Converting sam file to bam file
*******************************

For the downstream analyses we will work with the binary (more compact) version of the sam file. To convert the sam file in ban we will use **samtools view**. 
::

  samtools view -Sb filename.sam > filename.bam

.. tip::

  The conversion from sam to bam can be piped (|) in one command to the alignment step as follows:::

    bwa samse reference.fasta filename.sai filename.fastq.gz | samtools view -Sb - > filename.bam

To view the content of a sam file we can just use standart commands like **head**, **tail**, **less**. BUT, to view the content of a bam file (binary format of sam) we have to use **samtools view**. For example, to display on the screen one read/line (scrolling with the spacebar):
::

  samtools view filename.bam | less -S

while to display just the header of the bam file:
:: 

  samtools view -H filename.bam

*********************************
Sorting and indexing the bam file
*********************************

To go on with the analysis, we have to sort the reads aligned in the bam file by leftmost coordinates (or by read name when the option -n is used) with **samtools sort**. The option -o is used to provide an output file name:
::

  samtools sort filename.bam -o filename.sort.bam

The sorted bam files are then indexed with **samtools index**. Indexes allow other programs to retrieve specific parts of the bam file without reading through all of the sequences. The following command generates a **.bai** file, a companion file of the bam which contains the indexes:::

  samtools index filename.sort.bam

*********************************************
Adding Read Group tags and indexing bam files
*********************************************

A number of predefined tags may be appropriately assigned to specific set of reads in order to distinguish samples, libraries and other technical features. You may want to use RGLB (library ID) and RGSM (sample ID) tags at your own convenience.
::

  picard AddOrReplaceReadGroups INPUT=filename.sort.bam OUTPUT=filename.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT

.. note::
  
  - In some instances Picard may stop running and return error messages due to conflicts with SAM specifications produced by bwa (e.g. "MAPQ should be 0 for unmapped reads"). To suppress this error and allow the Picard to continue, we pass the VALIDATION_STRINGENCY=LENIENT options (default is STRICT).
  - Read Groups may be also added during the alignment with BWA using the option -R. 

Once added the Read Group tags, we index again the bam file:
:: 

  samtools index filename.RG.bam

*******************************
Marking and removing duplicates
*******************************

Amplification through PCR of genomic libraries leads to duplication formation, hence reads originating from a single fragment of DNA. The MarkDuplicates tool of Picard marks the reads as duplicates when the 5'-end positions of both reads and read-pairs match. A metric file with various statistics is created, and reads are removed from the bam file by using the REMOVE_DUPLICATES=True option (the default option is FALSE, which simply 'marks' duplicate reads keep them in the bam file).
:: 

  picard MarkDuplicates I=filename.RG.bam O=filename.DR.bam M=output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &> logFile.log

Once removed the duplicates, we index again the bam file:
::

  samtools index filename.DR.bam

**************************
Local realignment of reads
**************************

The presence of insertions or deletions (indels) in the genome may be responsible of misalignments and bases mismatches that are easily mistaken as SNPs. For this reason, we locally realign reads to minimize the number of mispatches around the indels. The realignment process is done in two steps using two different tools of GATK called with the -T option. We first detect the intervals which need to be realigned with the **GATK RealignerTargetCreator**, and save the list of these intevals in a file that we name targets.interval:
::

  gatk -T RealignerTargetCreator -R reference.fasta -I filename.DR.bam -o targets.intervals
 
.. note::

  like Picard, in our server environment we can call gatk just by typing the program name. In other environments (including your laptop) you may have to call gatk by providing the full path to the java file (.jar) of gatk:
  ::

    java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -h

Then, we realign the reads over the intervals listed in the targets.intervals file (*-targetIntervals* option) with **GATK IndelRealigner**:
::

  gatk -T IndelRealigner -R reference.fasta -I filename.RG.DR.bam -targetIntervals targets.intervals -o filename.final.bam --filter_bases_not_stored &> logFile.log

The final bam file has to be sorted and indexed as previously done:
::

  samtools sort filename.final.bam -o filename.final.sort.bam
  samtools index filename.final.sort.bam

**********************
Generate flagstat file
**********************

We can generate a file with useful information about our alignment with **samtools flagstat**. It represents a final summary report of the bitwise FLAG fields asigned to the reads in the sam file (to decode each FLAG field assigned to a read see https://broadinstitute.github.io/picard/explain-flags.html)
::

  samtools flagstat filename.final.sort.bam > flagstat_filename.txt

.. note::

  You could generate a flagstat file for the two bam files before and after refinement and see the differences. 


********************************
Visualization of reads alignment
********************************

Once generated the final bam file, let's compare the bam files before and after the refinement and polishing process (duplicates removal, realignment around indels and sorting). To do so, we will use the program IGV, in which we will first load the reference fasta file from *Genomes --> Load genome from file* and then we will add one (or more) bam files with *File --> Load from file*:

.. image:: images/igv-bam_bam.png
