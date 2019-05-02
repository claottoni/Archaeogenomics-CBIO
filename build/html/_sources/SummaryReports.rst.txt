######################
Create Summary Reports
######################

We will use Qualimap to create summary report from the generated bam files. As mentioned in the Qualimap website, Qualimap examines sequencing alignment data in SAM/BAM files according to the features of the mapped reads and provides an overall view of the data that helps to the detect biases in the sequencing and/or mapping of the data and eases decision-making for further analysis. 
::

  qualimap bamqc -c -bam input.bam 

Here are some screenshots about Qualimap outputs:

.. image:: images/qualimap.png


At this stage we have created different type of summary report using fastqc and qualimap. In addition, other program that we used like atropos (Cutadapt) generates also summary information about reads filtring and adaptor clipping. A list of programs that generate output files recognized by MultiQC are availble in this link https://github.com/ewels/MultiQC To Create one unique summary that integrate and compare all the generated report, we will use **`MutiQC`**. If all generated reports are in the same directory and its sub-directories, you can run simply MultiQC as follow:
::

  multiqc . 


After compliting, multiqc will creates a summary report in html format that will let you compare all the summary reports for ech of your samples.

.. image:: images/multiqc.png

