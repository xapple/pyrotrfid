"""
===================
Module **qiime.py**
===================

Binding for some of the tools found at:
http://qiime.sourceforge.net/scripts/

The tutorial for denoizing is here:
http://qiime.org/tutorials/denoising_454_data.html
"""

# Internal modules #
from pyrotrfid.cmd import command

################################################################################
@command
def check_id_map(mapping_file, output_dir):
    return {'arguments': ['check_id_map.py', '-m', mapping_file, '-o', output_dir]}

@command
def process_sff(input_dir):
    return {'arguments': ['process_sff.py', '-i', input_dir, '-f']}

@command
def split_libraries(fasta_path, qual_path, mapping_file, tag_length, output_dir, min_read_length, max_read_length):
    return {'arguments': ['split_libraries.py', '-f', fasta_path, '-q', qual_path, '-m', mapping_file, '-o', output_dir, '-b', str(tag_length), '-s', '20', '-l', min_read_length, '-L', max_read_length, '-n', '1000000']}

@command
def denoise_wrapper(sff_dump_path, fasta_path, mapping_file, output_dir):
    return {'arguments': ['denoise_wrapper.py', '-i', sff_dump_path, '-f', fasta_path, '-m', mapping_file, '-o', output_dir]}

@command
def inflate_denoiser_output(centroids, singletons, seqs, denoiser_mapping, denoised_seqs, output_path):
    return {'arguments': ['inflate_denoiser_output.py', '-c', centroids, '-s', singletons, '-f', seqs, '-d', denoiser_mapping, '-o', output_path]}