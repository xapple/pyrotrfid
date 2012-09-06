"""
====================
Module **mapseq.py**
====================

Binding for some sequence mapping tools
"""

# Built-in modules #
import os, re

# Internal modules #
from pyrotrfid.cmd import command

################################################################################
@command
def sam_to_bam(sam_filename, bam_filename):
    """Convert *sam_filename* to a BAM file.

    *sam_filename* must obviously be the filename of a SAM file.
    Returns the filename of the created BAM file.

    Equivalent: ``samtools view -b -S -o ...``
    """
    return {"arguments": ["samtools","view","-b","-S","-o",bam_filename,sam_filename],
            "return_value": bam_filename}

################################################################################
@command
def bwa_sw(reads_path, reference_path, sam_path, z=7):
    """Calls the BWA-SW aligner.
    http://bio-bwa.sourceforge.net/

    :param reads_path: The path to the short reads.
    :type reads_path: str
    :param reference_path: The path to the reference to map to.
    :type reference_path: str
    :param sam_path: The place were the SAM file will be created.
    :type sam_path: str
    :param z: A balance between quality and speed. Small z is speed.
    :type z: int
    """
    return {'arguments': ["bwa", "bwasw", "-z", str(z), "-f", sam_path, reference_path, reads_path],
            'return_value': sam_path}

################################################################################
@command
def fastqc(fastqfile,outdir=None,options=None):
    """Binds ``fastqc`` (http://www.bioinformatics.bbsrc.ac.uk/) which
    generates a QC report of short reads present into the fastq file.
    """
    outfile = re.sub(".fastq","",os.path.basename(fastqfile))+'_fastqc.zip'
    if not(isinstance(options,list)): options = []
    if outdir and os.path.isdir(outdir):
        outfile = os.path.join(outdir,outfile)
        options += ["--outdir",outdir]
    return {'arguments': ["fastqc","--noextract"]+options+[fastqfile],'return_value': outfile}

################################################################################
@command
def sort_bam(bamfile, filename):
    """Sort a BAM file *bamfile* by chromosome coordinates.

    Returns the filename of the newly created, sorted BAM file.

    Equivalent: ``samtools sort ...``
    """
    return {'arguments': ['samtools','sort',bamfile,filename],
            'return_value': filename + '.bam'}

################################################################################
@command
def index_bam(bamfile):
    """Index a sorted BAM file.

    Returns the filename in *bamfile* with ``.bai`` appended, that is,
    the filename of the newly created index.  *bamfile* must be sorted
    for this to work.

    Equivalent: ``samtools index ...``
    """
    return {'arguments': ['samtools','index',bamfile],
            'return_value': bamfile + '.bai'}