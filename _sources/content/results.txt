=======
Results
=======

Sequence quality report
-----------------------
The tools from http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/ generate an excellent quality control report tool for high throughput sequence data. Here is an example plot:

.. image:: /images/per_base_quality.png

A zip archive named something like ``12_DW-9.fq_fastqc.zip`` containing all the plots will be created for every sample.

Digital profile plot
--------------------
This plot shows the result of the digital digestion.

.. image:: /images/profile.pdf
   :width: 1024 px

A PDF named something like ``12_DW-9_profile.pdf`` will be created for every sample.

Cross correlation plot
----------------------
In order to check that profiles are well aligned, the software computes the cross correlation between the digital profile and the wetlab profile. The digital profile is truncated for this plot.

.. image:: /images/correlation.pdf

A PDF named something like ``12_DW-9_profile.pdf`` will be created for every sample.
Plot only available when a wetlab profile is given.

Mirror plot
-----------
Displays the the digital profile and the wetlab profile, side by side. The digital profile is shifted for this plot.

.. image:: /images/mirror.pdf
   :width: 1024 px

A PDF named something like ``12_DW-9_mirror.pdf`` will be created for every sample.
Plot only available when a wetlab profile is given.

Annotation file
---------------
A comma separated file containing the result of the digital digestion is written.
A CSV file named something like ``12_DW-9_peaks.xls`` will be created for every sample.

Digital profile matrix
----------------------
A summary of the digital peak proportions for all samples in a comma separated format.
A single CSV file named something like ``HaeIII_peaks.xls`` will be created.

Digital profile matrix lagged and cut
-------------------------------------
Same thing but with the cutoff and lag applied.