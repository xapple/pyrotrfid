==================
Executable options
==================

Usage: pyrotrfid SAMPLES_PATH [options]

A script to perform digital TRFLP. An example usage is the following::

   $ pyrotrfid 112708 --enzyme HaeIII

Options:
  -e ENZYME, --enzyme=ENZYME
                        The enzyme used for virtual digestion. Defaults to
                        'HaeIII'.
  -i REFERENCE, --reference=REFERENCE
                        The database to use to map the sequences. Defaults to
                        'greengenes'
  -x MIN_FRAG_SIZE, --min_frag_size=MIN_FRAG_SIZE
                        Minimum fragment length (inclusive). Defaults to 50.
  -y MAX_FRAG_SIZE, --max_frag_size=MAX_FRAG_SIZE
                        Maximum fragment length (inclusive). Defaults to 500.
  -p PRIMER_LENGTH, --primer_length=PRIMER_LENGTH
                        When SFF files are given, cut this much from every
                        sequence. Defaults to 20.
  -s SW_THRESHOLD, --sw_threshold=SW_THRESHOLD
                        Remove all short reads that mapped with a Smith-
                        Waterman score below this. Defaults to 150.
  -z HARD_LAG, --hard_lag=HARD_LAG
                        Manually specify the lag to apply between both
                        profiles. By default, lag calaucation is automatic.
  -f FILE_FORMAT, --file_format=FILE_FORMAT
                        The file format in which the graphs will be created.
                        Defaults to 'pdf'.
  -q, --qiime           Use QIIME to denoize the inputed SFF files. Defaults
                        to False.
  -l MIN_READ_LENGTH, --min_read_length=MIN_READ_LENGTH
                        Parameter used only in combination with the -q option.
                        Smaller reads are thrown away. Defaults to 300.
  -L MAX_READ_LENGTH, --max_read_length=MAX_READ_LENGTH
                        Parameter used only in combination with the -q option.
                        Larger reads are thrown away. Defaults to 500.
  --version             show program's version number and exit
  -h, --help            show this help message and exit
