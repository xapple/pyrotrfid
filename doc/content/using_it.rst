========
Using it
========

Input and output
----------------
The input directory has to be created by the user and filled with the appropriate files. The output directory will be created automatically.

The input directory can contain both FASTA or SFF files. If you give FASTA files, no sequence quality report can be made however.

The command line
----------------
You launch the analysis via the only command::

    $ pyrotrfid 112708

The number ``112708`` is the input directory. Change this to launch the analysis on different experiments.

There are extra optional parameters which will be detailed by calling the help as such::

    $ pyrotrfid --help

For instance, if you want to change the enzyme used, you could do the following::

    $ pyrotrfid 112708 --enzyme HaeIII

Other required files
--------------------
The input directory also needs to contain the original wetlab profiles in text format. These are expected to be contained in a file named "original_profile.csv". The first line indicate TRF lengths and the following lines indicate peak fractions with sample names. It should look something like this::

    ;123;126;128;153;177;180;185;186;189;190
    DW01;0;0;0;0.17300476;0;0;0.637599701;0;0.689988466;0
    DW02;0;0;0;0;0;0;0.592538187;0;0.660029974;0

If you wish to use qiime, the input directory also needs to contain the qiime mapping file in text format. This is expected to be contained in a file named "qiime_mapping.txt". And should look something like this::

    #SampleID       BarcodeSequence LinkerPrimerSequence    Description
    DW01    AGCAATGG        AGAGTTTGATCMTGGCTCAG    DW01
    DW02    AGCAAGCG        AGAGTTTGATCMTGGCTCAG    DW02

More information on the qiime website (http://qiime.org/documentation/file_formats.html#mapping-file-overview).