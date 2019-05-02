#####################################################
Damage analysis and quality rescaling of the BAM file
#####################################################

To authenticate our analysis we will assess the post-mortem damage of the reads aligned to the reference sequence. We can track the post-portem damage accumulated by DNA molecules in the form of fragmentation due to depurination and cytosine deamiation, which generated the typical pattern of C-T and G-A variation ad the molecule ends. To assess the post-mortem damage patterns in our bam file we will use **MapDamage**, which analyses the size distribution of the reads and the base composition of the genomic regions located up- and downstream of each read, generating various plots and summary tables. To start the analysis we need the final bam and the reference sequence: 
::

  mapDamage -i filename.final.sort.bam -r reference.fasta

mapDamage creates a new folder where the output files are created. One of these files, is named _Fragmisincorporation\_plot.pdf_ which contains the following plots:

.. image:: images/damage.png

If DNA damage is detected, we can run mapDamage again using the `--rescale-only` option and providing the path to the MapDamage results folder. This command will downscale quality scores at positions likely affected by deamination according to their initial quality values, position in reads and damage patterns. 
A new rescaled bam file is generated. 
::

  mapDamage -i filename.final.sort.bam -r reference.fasta --rescale-only -d results_folder

You can rescale the bam file directly in the first MapDamage command: 
::

  mapDamage -i filename.final.sort.bam -r reference.fasta --rescale

.. note::

  Another useful tool for estimating post-mortem damage (PMD) is *PMDTools* (https://github.com/pontussk/PMDtools). This program uses a model incorporating PMD, base quality scores and biological polymorphism to assign a PMD score to the reads. PMD > 0 indicates support for the sequence being genuinely ancient. PMDTools filters the damaged reads (based on the selected score) in a separate bam file which can be used for downstream analyses (e.g. variant call).

The rescaled bam file has to be indexed, as usual.
::

  samtools index filename.final.sort.rescaled.bam