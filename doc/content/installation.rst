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
Currently, Greengenes reference sequences used to map the dataset from sequencing experiments are provided by the Center for Environmental Biotechnology (Berkeley, USA) in the form of a FASTA file.

To make the plain text FASTA file into something usable, you would have to download appropriate files and prepare them with the ``index`` command from ``BWA``. It takes typically a few hours and is invoked as below::

    $ bwa index database.fa

Obviously, this step has to be done only if the index does not exists yet and has to be created or if there is a need for an update. Greengenes files are found here http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/. Put the resulting files in the input directory with the rest of your samples.

One could use other reference databases such as RDP, Silva, if ones adds the right parsing code in the python source files.

Graphical library
-----------------
The matplotlib library can use many different backends. If the default one does not work on your system, you can easy change it and select an other one (http://matplotlib.sourceforge.net/faq/usage_faq.html#what-is-a-backend). The safest bet is usually "AGG".

Licensing
---------
The code is free and open source, licensed under the GNU General Public License version 3. You can download, copy or fork the source code on `github <https://github.com/xapple/pyrotrfid>`_. To download a copy of the latest version of the source code to your computer you would type::

    $ git clone git://github.com/xapple/pyrotrfid.git

Reporting bugs
--------------
As always with software that is released for the first time there will be bugs. We would be glad to hear from you if you have found any. You can post them in our `issue tracking system <https://github.com/xapple/pyrotrfid/issues>`_ that is found in the github repository. You will however need to create a github account if you don't already have one to open a new issue, sorry.
