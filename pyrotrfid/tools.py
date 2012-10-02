"""
====================
Module **tools.py**
====================

Binding for some sequence mapping tools amongst others.
"""

# Internal modules #
from pyrotrfid.cmd import command

################################################################################
@command
def sam_to_bam(sam_filename, bam_filename):
    """Convert *sam_filename* to a BAM file."""
    return {"arguments": ["samtools","view","-b","-S","-o",bam_filename,sam_filename]}

################################################################################
@command
def bwa_sw(reads_path, reference_path, sam_path, z=7):
    """Calls the BWA-SW aligner.
    http://bio-bwa.sourceforge.net/

    :param reads_path: The path to the short reads.
    :param reference_path: The path to the reference to map to.
    :param sam_path: The place were the SAM file will be created.
    :param z: A balance between quality and speed. Small z is speed.
    """
    return {'arguments': ["bwa", "bwasw", "-z", str(z), "-f", sam_path, reference_path, reads_path]}

################################################################################
@command
def fastqc(fastqfile, outdir="."):
    """Binds ``fastqc`` (http://www.bioinformatics.bbsrc.ac.uk/) which
    generates a QC report of short reads present into the fastq file"""
    return {'arguments': ["fastqc", "--noextract", "--outdir", outdir, fastqfile]}

################################################################################
@command
def sort_bam(bamfile, filename):
    """Sort a BAM file *bamfile* by chromosome coordinates."""
    return {'arguments': ['samtools', 'sort', bamfile, filename]}

################################################################################
@command
def index_bam(bamfile):
    """Index a sorted BAM file."""
    return {'arguments': ['samtools', 'index', bamfile]}