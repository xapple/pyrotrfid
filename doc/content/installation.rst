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

These can generally be obtained using the ``easy_install`` or ``pip`` command on standard unix platforms.

In addition these two executables are required and should be accessible somewhere in your ``$PATH``:

    * bwa (http://bio-bwa.sourceforge.net)
    * fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

This database needs to be installed and accessible (more on that later):

    * Green Genes (http://greengenes.lbl.gov/cgi-bin/nph-index.cgi)

This library also needs to be installed if you want to use the ``qiime`` option when running the pipeline. It's not an easy library to install.

     * qiime (http://qiime.sourceforge.net)

Once you have done all this you are ready to run the install script found in the pyrotrfid directory::

    $ python setup.py install

Making reference indexes
------------------------
Currently, reference sequences used to map the dataset from sequencing experiments are provided by the Center for Environmental Biotechnology (Berkeley, USA) in the form of a FASTA file. However, one can use other reference databases such as the RDP database, for instance.

To do this you would have to download appropriate files and prepare them as with the ``index`` command from ``BWA``. It takes typically a few hours and is invoked as below::

    $ bwa index database.fa

Obviously, this step has to be done only if the index does not exists yet and has to be created or if there is a need for an update. Green genes files are found here http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/. Put the resulting files in the input directory with the rest of your samples.

Development
-----------
Development is done using a git repository and can be cloned at::

    $ git clone git://github.com/xapple/pyrotrfid.git