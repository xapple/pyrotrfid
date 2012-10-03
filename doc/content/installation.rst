============
Installation
============

Dependencies
------------
The project has dependencies on the following third party python packages:

    * scipy
    * matplotlib
    * biopython
    * pysam

These can generally be obtained using the ``easy_install`` or ``pip`` commands on standard unix platforms.

In addition these two executables are required and should be accessible somewhere in your ``$PATH``:

    * bwa (http://bio-bwa.sourceforge.net)
    * fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

This public database of reference sequences needs to be downloaded, installed, and accessible (more details later):

    * Greengenes (http://greengenes.lbl.gov/cgi-bin/nph-index.cgi)

This library also needs to be installed if you want to use the ``qiime`` option when running the pipeline. It's not an easy library to install.

    * qiime (http://qiime.sourceforge.net)

Once you have done all this you are ready to download the latest version of pyrotrfid here: https://github.com/xapple/pyrotrfid/tags

and run the install script found in the pyrotrfid directory::

    $ python setup.py install

Making reference indexes
------------------------
Currently, Greengenes reference sequences used to map the dataset from sequencing experiments are provided by the Center for Environmental Biotechnology (Berkeley, USA) in the form of a FASTA file. You will need to obtains this file from http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/ and then modify it in a particular way::

    $ sed -i 's/ /_/g' current_GREENGENES_gg16S_unaligned.fasta

Next you need to prepare the file with the ``index`` command from ``BWA``. It takes typically a few minutes and is invoked as below::

    $ bwa index current_GREENGENES_gg16S_unaligned.fasta

Put the resulting files in the input directory with the rest of your samples.

Licensing
---------
The code is free and open source, licensed under the GNU General Public License version 3. You can download, copy or fork the source code on `github <https://github.com/xapple/pyrotrfid>`_. To download a copy of the latest version of the source code to your computer you would type::

    $ git clone git://github.com/xapple/pyrotrfid.git

Reporting bugs
--------------
As always with software that is released for the first time there will be bugs. We would be glad to hear from you if you have found any. You can post them in our `issue tracking system <https://github.com/xapple/pyrotrfid/issues>`_ that is found in the github repository. You will however need to create a github account if you don't already have one to open a new issue, sorry.
