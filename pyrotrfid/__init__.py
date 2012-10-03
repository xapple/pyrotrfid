b'This module needs Python 2.6 or later.'

# Special variables #
__version__ = '0.9.1'

# Built-in modules #
import os, platform, glob, shutil, re, time, tempfile
from collections import defaultdict

# Plotting module #
import matplotlib
matplotlib.use('Agg', warn=False)
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator

# Third party modules #
import scipy, pysam, Bio.Seq, Bio.SeqIO, Bio.Restriction

# Internal modules #
from pyrotrfid.common import Color, property_cached, natural_sort
from pyrotrfid.common import wrap_string, andify_strings, shift_list
from pyrotrfid.cmd import pause_for_parralel_jobs
from pyrotrfid import qiime, tools

# Constants #
seperator = "\t"
delimiter = ";"

################################################################################
class Job(object):
    """Every time the ``pyrotrfid`` script is run, a job is created.

    :param samples_path: Path to sff files directory.
    :type  samples_path: str
    :param options: User provided options in a dictionary.
    :type  options: dict
    """

    def __init__(self, samples_path, options):
        # Copy variables #
        self.samples_path = samples_path
        self.__dict__.update(options)
        # Other useful stuff #
        self.version = __version__
        self.module = os.path.abspath(__file__)
        self.hostname = platform.node()
        # Set the plotting module format #
        matplotlib.use(self.file_format, warn=False)

    def run(self):
        """Run the job to completion."""
        self.print_info()
        #-------------------------------#
        print "Writing job parameters..."
        with open(self.output_dir + "job_parameters.txt", 'w') as file:
            lines = ((k + ":").ljust(32) + str(v) + "\n" for k,v in self.parameters.items())
            file.writelines(lines)
        #-------------------------------#
        print "Converting the samples..."
        sff_samples = [sample for sample in self.samples if sample.format == 'sff']
        for sample in sff_samples: sample.convert_to_fastq(self.output_dir, self.primer_length)
        #-------------------------------#
        print "Generating the quality reports..."
        for sample in sff_samples: sample.quality_report(self.output_dir)
        #-------------------------------#
        if self.qiime:
            print "Denoizing the samples with Qiime..."
            for sample in sff_samples: sample.qiime_split_sff_file(self.tmp_dir)
            for sample in sff_samples: sample.barcode = self.barcodes[sample.name]
            for sample in sff_samples: sample.add_barcode(self.tmp_dir)
            for sample in sff_samples: sample.qiime_filter(self.tmp_dir, self.qiime_mapping_file, self.min_read_length, self.max_read_length)
            futures = [sample.qiime_denoise_local(self.tmp_dir, self.qiime_mapping_file) for sample in self.samples]
            pause_for_parralel_jobs()
            for future in futures: future.wait()
            for sample in sff_samples: sample.qiime_inflate(self.tmp_dir)
        #-------------------------------#
        print "Mapping the reads..."
        futures = [sample.map_reads_local(self.reference_path, self.tmp_dir) for sample in self.samples]
        pause_for_parralel_jobs()
        for future in futures: future.wait()
        #for sample in self.samples: sample.sam_path = self.tmp_dir + sample.name + '.sam'
        #-------------------------------#
        print "Indexing the results..."
        for sample in self.samples: sample.create_bam_file(self.tmp_dir)
        #-------------------------------#
        print "Reporting non-mapping reads..."
        for sample in self.samples: sample.report_non_mapping_reads(self.output_dir)
        #-------------------------------#
        print "Digesting the reads..."
        for sample in self.samples: sample.digest(self.enzyme, self.primer_length, self.sw_threshold)
        #-------------------------------#
        print "Creating digital profile plots..."
        for sample in self.samples: sample.digital_profile_plot(self.output_dir, self.file_format)
        #-------------------------------#
        print "Loading wetlab profiles..."
        for sample in self.samples:
            try: sample.wetlab_profile = self.wetlab_profiles[sample.name]
            except KeyError: sample.wetlab_profile = None
            except TypeError: sample.wetlab_profile = False
        wetlab_available = [sample for sample in self.samples if sample.wetlab_profile]
        #-------------------------------#
        print "Creating correlation plots..."
        for sample in wetlab_available: sample.correlation_plot(self.output_dir, self.file_format)
        #-------------------------------#
        print "Creating mirror plots..."
        for sample in wetlab_available: sample.mirror_plot(self.output_dir, self.file_format)
        #-------------------------------#
        print "Exporting annotations..."
        for sample in self.samples: sample.create_annotation_file(self.output_dir)
        #-------------------------------#
        print "Exporting peak frequencies..."
        with open(self.output_dir + self.enzyme + "_peaks.xls", 'w') as file:
            file.writelines(";" + ";".join(map(str,range(self.max_frag_size+1))) + '\n')
            file.writelines(sample.text_digital_peaks + '\n' for sample in self.samples)
        #-------------------------------#
        print "Exporting lagged cut peak frequencies..."
        with open(self.output_dir + self.enzyme + "_peaks_cut_lagged.xls", 'w') as file:
            file.writelines(";" + ";".join(map(str,range(self.max_frag_size+1))) + '\n')
            file.writelines(sample.text_digital_peaks_cut_lagged + '\n' for sample in wetlab_available)
        #-------------------------------#
        print "Cleaning up..."
        shutil.rmtree(self.tmp_dir)
        #-------------------------------#
        print "Success. Results are in %s" % self.output_dir

    @property_cached
    def parameters(self):
        """All the job parameters."""
        return dict(((k,v) for k,v in sorted(self.__dict__.items()) if not k.startswith('__')))

    def print_info(self):
        """Display some info before the job starts"""
        print Color.f_pur + "Starting job on %s" % time.asctime() + Color.end
        for k,v in self.parameters.items(): print (Color.ylw + k + ":" + Color.end).ljust(32) + str(v)
        print Color.f_pur + "Processing %i samples" % len(self.samples) + Color.end
        sample_names = [sample.name + '.' + sample.format for sample in self.samples]
        print Color.b_cyn + ' '.join(sample_names) + Color.end

    @property_cached
    def samples(self):
        """Search the raw data directory for samples."""
        sff_files   = glob.glob(self.input_dir + "*.sff")
        fasta_files = glob.glob(self.input_dir + "*.fasta")
        parameters = (self.min_frag_size, self.max_frag_size, self.hard_lag)
        samples = [Sample(f, *parameters) for f in sff_files or fasta_files]
        samples.sort(key=lambda x: natural_sort(x.name))
        return samples

    @property_cached
    def tmp_dir(self):
        """The temporary directory will be removed once the job is done."""
        value = tempfile.gettempdir()
        if not value.endswith('/'): value += '/'
        value += "pyrotrfid/"
        if not os.path.exists(value): os.makedirs(value)
        return value

    @property_cached
    def output_dir(self):
        """Results will fall in here."""
        value = self.input_dir + "pyrotrfid/"
        if not os.path.exists(value): os.makedirs(value)
        return value

    @property_cached
    def input_dir(self):
        """Input files will be searched for in here."""
        value = self.samples_path
        if not value.endswith('/'): value += '/'
        if not os.path.exists(value): raise Exception("Input directory '%s' does not exist." % value)
        return value

    @property_cached
    def wetlab_profile_path(self):
        """Path to original wet lab TRFLP profile"""
        value = self.input_dir + "original_profile.csv"
        if not os.path.exists(value): return False
        return value

    @property_cached
    def qiime_mapping_file(self):
        """Path to 'mapping' file that qiime uses at some moment."""
        value = self.input_dir + "qiime_mapping.txt"
        if self.qiime and not os.path.exists(value): raise Exception("Mapping file '%s' not found." % value)
        return value

    @property_cached
    def reference_path(self):
        """Path to original wet lab TRFLP profile."""
        files = glob.glob(self.input_dir + "*.amb")
        if not files: raise Exception("AMB files not found.")
        return os.path.splitext(files[0])[0]

    @property_cached
    def wetlab_profiles(self):
        """Loads the wetlab peak counts from the appropriate file."""
        # We might be missing it #
        if not self.wetlab_profile_path: return False
        # It is an excel CSV file #
        with open(self.wetlab_profile_path, 'r') as file:
            data = [line.strip().split(';') for line in file]
        # The first line contains the X coordinates #
        x_coordinates = map(int,data.pop(0)[1:])
        # The rest of the file is one sample per line #
        # It starts with the sample name #
        result = {}
        for line in data:
            values = map(float,line[1:])
            peaks = dict(zip(x_coordinates,values))
            peaks = [peaks.get(x, 0.0) for x in range(self.max_frag_size+1)]
            result[line[0]] = peaks
        return result

    @property_cached
    def barcodes(self):
        """A dictionary with a key for every sample and a sequence as every value."""
        with open(self.qiime_mapping_file, 'r') as file:
            data = [line.strip().split('\t') for line in file if not line.startswith("#")]
            return dict(((x[0],x[1]) for x in data))

################################################################################
class Sample(object):
    """Every raw file containing short reads is considered as a separate sample.
    A sample will be converted to FASTAQ format, mapped, sorted and indexed;
    after which it will be virtually digested, the resulting fragments counted
    and grouped by bacteria. Plots will be made comparing these results with the
    results from the wet lab TRFLP.

    :param path: Path to the file containing the short reads. Can be SFF or FASTA.
    :type  path: str
    :param min_frag_size: Restriction fragments shorter than this are ignored.
    :type min_frag_size: int
    :param max_frag_size: Restriction fragments larger than this are ignored.
    :type max_frag_size: int
    :param hard_lag: If not False, the lag will be fixed to this value.
    :type hard_lag: int
    """

    def __init__(self, path, min_frag_size, max_frag_size, hard_lag):
        self.path = path
        self.min_frag_size = min_frag_size
        self.max_frag_size = max_frag_size
        self.hard_lag = hard_lag

    @property
    def format(self):
        """The three letter extension of the short reads file."""
        base, ext = os.path.splitext(self.path)
        return ext[1:]

    @property
    def name(self):
        """The name of the short reads file without the extension."""
        base, ext = os.path.splitext(self.path)
        return os.path.basename(base)

    @property
    def directory(self):
        """The directory in which the short reads file is."""
        return os.path.dirname(self.path)

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)

    def qiime_split_sff_file(self, directory):
        """Convert an SFF sample to FASTA and QUAL formatted files
        in a sub directory in the given directory using QIIME."""
        # Make a sub directory #
        sub_dir = directory + self.name + '_splited/'
        if not os.path.exists(sub_dir): os.mkdir(sub_dir)
        # Put the SFF file alone in a directory #
        shutil.copy(self.path, sub_dir)
        # Run the QIIME script #
        qiime.process_sff(sub_dir)
        # Get the results #
        self.qual_path = sub_dir + self.name + '.qual'
        self.fasta_path = sub_dir + self.name + '.fna'
        self.sff_dump_path = sub_dir + self.name + '.txt'

    def add_barcode(self, directory):
        """Manually adds the short barcode to every sequence
        in the fasta file generated by 'process_sff'."""
        # Make a sub directory #
        sub_dir = directory + self.name + '_barcoded/'
        if not os.path.exists(sub_dir): os.mkdir(sub_dir)
        # The new fasta with barcodes #
        path = sub_dir + self.name + '.fna'
        with open(self.fasta_path, 'r') as input:
            with open(path, 'w') as output:
                for line in input:
                    if line.startswith('>'): output.write(line)
                    else: output.write(self.barcode + line)
        # The new fasta with barcodes #
        self.fasta_path = path

    def qiime_filter(self, directory, qiime_mapping_file, min_read_length, max_read_length):
        """Make a new FASTA file where low quality reads are removed."""
        # Make a sub directory #
        sub_dir = directory + self.name + '_filtered/'
        if not os.path.exists(sub_dir): os.mkdir(sub_dir)
        # Call the filtering #
        qiime.check_id_map(qiime_mapping_file, sub_dir)
        qiime.split_libraries(self.fasta_path, self.qual_path, qiime_mapping_file, len(self.barcode), sub_dir, min_read_length, max_read_length)
        self.filtered_seqs = sub_dir + 'seqs.fna'

    def qiime_denoise_local(self, directory, qiime_mapping_file):
        """Perform a flowgram clustering locally."""
        # Don't make a sub directory #
        sub_dir = directory + self.name + '_clustered/'
        # Keep the future results #
        self.centroids = sub_dir + 'centroids.fasta'
        self.singletons = sub_dir + 'singletons.fasta'
        self.denoiser_mapping = sub_dir + 'denoiser_mapping.txt'
        # Call the denoiser #
        return qiime.denoise_wrapper.parallel(self.sff_dump_path, self.filtered_seqs, qiime_mapping_file, sub_dir)

    def qiime_inflate(self, directory):
        """Make a new FASTA file where clusters containing N
        individuals are wirtten N times to the output."""
        # Make a sub directory #
        sub_dir = directory + self.name + '_denoised/'
        if not os.path.exists(sub_dir): os.mkdir(sub_dir)
        fa_path = sub_dir + self.name + '.fa'
        # Call the inflater #
        qiime.inflate_denoiser_output(self.centroids, self.singletons, self.filtered_seqs, self.denoiser_mapping, self.path, fa_path)
        # The sample becomes a FASTA sample #
        self.path = fa_path

    def convert_to_fastq(self, directory, primer_length):
        """Convert an SFF sample to FASTAQ sample in the given directory.
           The small 8 base pair tag will be removed automatically by Bio.SeqIO.
           In addition to this, we cut off *primer_length* base pairs
           from the start of every read."""
        # Get the paths #
        sff_path = self.path
        fq_path = directory + self.name + '.fq'
        # Open the files #
        sff_handle = open(sff_path, "r")
        fq_handle = open(fq_path, "w")
        # Parse the SFF file with the tag trimming option #
        sff_stream = Bio.SeqIO.parse(sff_handle, "sff-trim")
        # Trim the primer from the reads #
        sff_stream = (record[primer_length:] for record in sff_stream)
        # Write the result #
        self.read_count = Bio.SeqIO.write(sff_stream, fq_handle, "fastq")
        # Close the files #
        sff_handle.close()
        fq_handle.close()
        # The sample becomes a FASTAQ sample #
        self.path = fq_path

    def quality_report(self, directory):
        """Produce the sample sequencing quality report in the given directory."""
        tools.fastqc(self.path, directory)

    def map_reads_local(self, reference, directory):
        """Map the sample to the given reference in the given directory."""
        self.sam_path = directory + self.name + '.sam'
        return tools.bwa_sw.parallel(self.path, reference, self.sam_path)

    def copy_sam_file(self, directory):
        """Copy the sam file to the given directory."""
        new_sam_path = directory + self.name + '.sam'
        shutil.copy(self.sam_path, new_sam_path)
        self.sam_path = new_sam_path

    def create_bam_file(self, directory):
        """Create the sample's bam file from the sam file in the given directory."""
        unsorted_bam = directory + self.name + '_unsorted.bam'
        tools.sam_to_bam(self.sam_path, unsorted_bam)
        tools.sort_bam(unsorted_bam, directory + self.name + '_sorted')
        self.bam_path = directory + self.name + '_sorted.bam'
        tools.index_bam(self.bam_path)

    def report_non_mapping_reads(self, directory):
        """Create a file with all the non mapped reads in the given directory."""
        short_reads = pysam.Samfile(self.bam_path, "rb")
        with open(directory + self.name + "_unmapped_reads.txt", 'w') as file:
            for read in short_reads:
                if read.flag != 0 and read.flag != 16:
                    file.write(read.qname + " " + read.seq + "\n")

    #-------------------------------------------------------------------------#
    def digest(self, enzyme, primer_length, sw_threshold):
        """Virtually digest all the reads in the bam file.
        The documentation for the ``AlignedRead`` object is here:
        http://www.cgat.org/~andreas/documentation/pysam/api.html

        :param enzyme: The enzyme name.
        :type  enzyme: str
        :param primer_length: The amount of base pairs to cut off from every read.
        :type  primer_length: int
        :param sw_threshold: Ignore reads that mapped with a Smith-Waterman score below this.
        :type  sw_threshold: int
        """
        # Get the enzyme from the biopython library #
        try: self.enzyme = getattr(Bio.Restriction, enzyme)
        except AttributeError: raise Exception("Enzyme defined as '%s' is not present in the biopython library." % enzyme)
        # Prepare the peak objects #
        self.peaks = dict([(l,Peak(l)) for l in range(self.max_frag_size+1)])
        # Open the bam file with the pysam tool #
        short_reads = pysam.Samfile(self.bam_path, "rb")
        # Do something for every read #
        for read in short_reads:
            # Skip the ones that didn't map #
            unmapped = True if read.flag != 0 and read.flag != 16 else False
            if unmapped: continue
            # Skip the ones that have a SW score too low #
            sw_score = read.tags[0][1]
            if sw_score < sw_threshold: continue
            # Cut the short read sequence in pieces #
            fragments = self.enzyme.search(Bio.Seq.Seq(read.seq))
            # Skip the ones that don't have a restriction site #
            if not fragments: continue
            # Only the first restriction site is considered #
            fragment_cut_location = fragments[0]
            # Cutting at position 6 will give a fragment of size 5 #
            fragment_length = fragment_cut_location - 1
            # Adjust the size to include the primer that was removed earlier #
            fragment_length += primer_length
            # Ignore those that are too long #
            if fragment_length > self.max_frag_size: continue
            # Retrieve the taxonomy #
            bacteria = Bacteria(short_reads.getrname(read.tid))
            # Create the fragment object #
            fragment = Fragment(fragment_length, bacteria, read.qname, read.rlen, sw_score)
            # Add the fragment to the right peak #
            self.peaks[fragment_length].add_fragment(fragment)

    @property_cached
    def digital_peaks_prop(self):
        """The in silico fragment proportion at each position."""
        digital_peaks_counts = [self.peaks[l].total_fragments for l in sorted(self.peaks)]
        total_counts = sum(digital_peaks_counts)
        return [100 * x / total_counts for x in digital_peaks_counts]

    @property_cached
    def digital_peaks_prop_cut(self):
        """Same thing, but with many low values set to zero."""
        result = [0]*len(self.digital_peaks_prop)
        result[self.min_frag_size:] = self.digital_peaks_prop[self.min_frag_size:]
        return result

    @property_cached
    def digital_peaks_prop_lagged(self):
        """Same thing, but shifted."""
        if not self.wetlab_profile: return None
        return shift_list(self.digital_peaks_prop, self.lag)

    @property_cached
    def digital_peaks_cut_prop_lagged(self):
        """Apply the cutoff and recalculate the proportions."""
        digital_peaks_counts = [l >= self.min_frag_size and self.peaks[l].total_fragments or 0 for l in sorted(self.peaks)]
        total_counts = sum(digital_peaks_counts)
        proporitions = [100 * x / total_counts for x in digital_peaks_counts]
        return shift_list(proporitions, self.lag)

    @property_cached
    def text_digital_peaks(self):
        """A long string containing the relative proportion
        of each fragment peak at each position. Typically the string
        looks something like::

            11_DW-8;0.0;0.0;0.0;0.0;0.245901;0.573770;3.901639;0.0;0.0;0.0;0.0;0.245901;0.225409;
        """
        return delimiter.join([self.name] + map(str, self.digital_peaks_prop))

    @property_cached
    def text_digital_peaks_cut_lagged(self):
        return delimiter.join([self.name] + map(str, self.digital_peaks_cut_prop_lagged))

    @property_cached
    def lag(self):
        if not self.wetlab_profile: return None
        if self.hard_lag: return self.hard_lag
        else: return self.auto_lag

    @property_cached
    def auto_lag(self):
        # Use scipy to correlate #
        correlation = scipy.correlate(self.wetlab_peaks_prop, self.digital_peaks_prop_cut, mode='full')
        # Make a dictionary of possible lags #
        length = len(self.wetlab_peaks_prop)
        lags = xrange(-length+1, length)
        lag_dict = dict(zip(lags,correlation))
        # Apply a hard coded search cutoff #
        for lag in lag_dict:
            if lag < -10 or lag > 10: lag_dict[lag] = 0
        # Best lag #
        return max(lag_dict, key=lag_dict.get)

    @property
    def wetlab_peaks_prop(self):
        """The in vitro fragment proportion at each position."""
        return self.wetlab_profile

    #-------------------------------------------------------------------------#
    def create_annotation_file(self, directory):
        """Create the XLS file containing all the annotations in the given directory."""
        with open(directory + self.name + "_annotations.xls", 'w') as file:
            # Write the first line containing the titles #
            peak_columns = ["Length", "Length shifted", "Count", "Percentage"]
            bacteria_columns = Bacteria.keys
            fragment_columns = ["Normalized SW score", "Normalized SW score", "Short read id"]
            file.write(seperator.join(peak_columns + bacteria_columns + fragment_columns) + '\n')
            # We might have a lag computed #
            lag = self.lag if self.lag else 0
            # Go through every peak and export it #
            generator = (peak.export_all_bacteria(lag) for peak in self.peaks.values() if peak)
            generator = (item for sublist in generator for item in sublist)
            generator = (seperator.join(map(str,line)) + '\n' for line in generator)
            # Each peak can have several bacteria #
            file.writelines(generator)

    def digital_profile_plot(self, directory, file_format):
        """Create the digital profile plot in the given directory in the given format."""
        # Create figure #
        fig = pyplot.figure(figsize=(20,12))
        # Make two axes #
        axes = fig.add_subplot(111)
        fig.subplots_adjust(hspace=0.0, bottom=0.07, top=0.93, left=0.04, right=0.98)
        # Data and source #
        fig.text(0.99, 0.98, time.asctime(), horizontalalignment='right')
        fig.text(0.01, 0.98, __name__.split('.')[0] + ' generated graph', horizontalalignment='left')
        # Vertical ticks every 5 #
        axes.xaxis.set_major_locator(MultipleLocator(5))
        axes.xaxis.set_minor_locator(MultipleLocator(1))
        for tick in axes.xaxis.get_major_ticks(): tick.label.set_rotation('vertical')
        # A dense grid #
        axes.grid(True, which="major", linestyle='-', linewidth=0.25, color='0.5')
        axes.grid(True, which="minor", linestyle=':', linewidth=0.25, color='0.5')
        # Titles #
        fig.suptitle('All peaks are shown. No shift is applied.')
        axes.set_title('pyrotrfid profile plot for sample "%s"' % self.name)
        axes.set_xlabel('Fragment length in base pairs after digital digestion with "%s"' % self.enzyme)
        axes.set_ylabel('In silico fragment proportion [%]')
        # Plotting #
        x_values = sorted(self.peaks)
        axes.bar(x_values, self.digital_peaks_prop, color='k', label='dTRFLP', linewidth=0, width=1.0)
        # Legend #
        axes.legend(prop={'size': 12}, fancybox=True, loc=4)
        # Slightly larger area #
        x_min, x_max, y_min, y_max = axes.axis()
        axes.axis((x_min, x_max, y_min, y_max*1.2))
        # Save it #
        fig.savefig(directory + self.name + "_profile." + file_format)

    def mirror_plot(self, directory, file_format):
        """Create the mirror plot in the given directory in the given format."""
        # Create figure #
        fig = pyplot.figure(figsize=(20,12))
        # Make two axes #
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212, sharex=ax1)
        # Glue them together #
        fig.subplots_adjust(hspace=0.0, bottom=0.07, top=0.93, left=0.04, right=0.98)
        pyplot.setp(ax1.get_xticklabels(), visible=False)
        # Data and source #
        fig.text(0.99, 0.98, time.asctime(), horizontalalignment='right')
        fig.text(0.01, 0.98, __name__.split('.')[0] + ' generated graph', horizontalalignment='left')
        # Vertical ticks every 5 #
        ax2.xaxis.set_major_locator(MultipleLocator(5))
        ax2.xaxis.set_minor_locator(MultipleLocator(1))
        for tick in ax2.xaxis.get_major_ticks(): tick.label.set_rotation('vertical')
        # A dense grid #
        ax1.grid(True, which="major", linestyle='-', linewidth=0.1, color='0.5')
        ax2.grid(True, which="major", linestyle='-', linewidth=0.1, color='0.5')
        ax1.grid(True, which="minor", linestyle=':', linewidth=0.1, color='0.5')
        ax2.grid(True, which="minor", linestyle=':', linewidth=0.1, color='0.5')
        # Titles #
        fig.suptitle('Cuttoff is applied and digital peaks are shifted by %i base pairs' % self.lag)
        ax1.set_title('Mirror plot for sample "%s"' % self.name)
        ax2.set_xlabel('Fragment length in base pairs after digestion with "%s"' % self.enzyme)
        # Plotting #
        x_values = sorted(self.peaks)
        ax1.bar(x_values, self.digital_peaks_cut_prop_lagged, color='c', label='dTRFLP', linewidth=0, width=1.0)
        ax1.set_ylabel('In silico fragment proportion [%]')
        ax2.bar(x_values, self.wetlab_peaks_prop, color='m', label='eTRFLP', linewidth=0, width=1.0)
        ax2.set_ylabel('In vitro fragment proportion [%]')
        # Legend #
        ax1.legend(prop={'size': 12}, fancybox=True, loc=1)
        ax2.legend(prop={'size': 12}, fancybox=True, loc=4)
        # Slightly larger area #
        x_min, x_max, y_min, y_max = ax1.axis()
        ax1.axis((x_min, x_max, y_min, y_max*1.2))
        # Inversion #
        ax2.axis((x_min, x_max, y_min, y_max*1.2))
        ax2.invert_yaxis()
        # Save it #
        fig.savefig(directory + self.name + "_mirror." + file_format)

    def correlation_plot(self, directory, file_format):
        """Create the cross correlation plot in the given directory in the given format."""
        # Create figure #
        fig = pyplot.figure(figsize=(12,8))
        axes = fig.add_subplot(111)
        # Titles #
        fig.suptitle('Given this correlation, the best shift is %i' % self.auto_lag)
        axes.set_title('Cross correlation for sample "%s" with a cutoff at %i' % (self.name, self.min_frag_size))
        axes.set_xlabel('Shift between digital and wet lab profiles in base pairs.')
        axes.set_ylabel('Correlation [no units]')
        # Plotting #
        axes.xcorr(self.wetlab_peaks_prop, self.digital_peaks_prop_cut, normed=True, maxlags=10)
        # Save it #
        fig.savefig(directory + self.name + "_correlation." + file_format)

################################################################################
class Peak(object):
    """All restriction fragments of the same size are regrouped in one peak object.
    Hence a peak object has a unique peak length and contains many fragments.
    Furthermore, inside a peak, fragments coming from the same bacteria are grouped.
    """

    def __init__(self, length):
        self.length = length
        self.fragments = defaultdict(list)

    def __nonzero__(self): return bool(self.fragments)
    def __repr__(self): return '<%s object of size %i>' % (self.__class__.__name__, self.length)

    def add_fragment(self, fragment):
        self.fragments[fragment.bacteria].append(fragment)

    @property_cached
    def total_fragments(self):
        """Return the number of fragments in this peak."""
        return float(sum([1 for bacteria in self.fragments.values() for x in bacteria]))

    @property_cached
    def annotation_string(self):
        """Return the names of the principal bacteria behind this peak."""
        bacterias = sorted(self.fragments.items(), key=lambda x: len(x[1]))[-3:]
        bacterias.reverse()
        names = ['%s (%i)' % (k.name, len(v)) for k,v in bacterias]
        return wrap_string(andify_strings(names), 32)

    def export_all_bacteria(self, lag):
        """Return a generator yielding one text line per bacteria.
        Typically one line looks something like::

            '56 14.8148148148 4 Bacteria Thermi Deinococci 4323 556010_GQ441240.active_nitrogen'

        """
        for bacteria, frags in self.fragments.items():
            peak_info = [self.length, self.length+lag, len(frags), 100*len(frags)/self.total_fragments]
            read_info = [delimiter.join(str(f.sw_score) for f in frags)]
            read_info += [delimiter.join(str(f.sw_norm) for f in frags)]
            read_info += [delimiter.join(f.ref for f in frags)]
            yield peak_info + bacteria.taxon + read_info

################################################################################
class Bacteria(object):
    """A bacteria is described by its taxon. This can be determined by
    parsing an entry in the Green Genes database.
    A bacteria object can be used as a key to a dictionary.
    The annotation_string typically looks something like::

         239015_EU071527.1_and_resistance_microorganisms_European_clean_room_ESTEC_HYDRA_facility_clone_EHFS1_S15a_k__Bacteria;_p__Proteobacteria;_c__Gammaproteobacteria;_o__Xanthomonadales;_f__Xanthomonadaceae;_g__Thermomonas;_Unclassified;_otu_4045

    The abbreviation OTU stands for "Operational Taxonomic Unit".
    """

    keys = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Gender', 'Species', 'OTU', 'Description']
    tags = ['_k__',    '_p__',   '_c__',  '_o__',  '_f__',   '_g__',   '_s__',    '_otu_']

    def __init__(self, annotation_string):
        # A list to contain the result #
        self.taxon = []
        # For every tag, find it in the string #
        for tag in self.tags:
            clade = re.findall('%s([^;]+)' % tag, annotation_string)
            self.taxon.append(clade[0] if clade else '')
        # The 'description' key is special #
        self.taxon.append(annotation_string.split("_k__")[0])

    @property_cached
    def name(self):
        """Return a short and human readable name."""
        for i in [6, 5, 4, 3, 2, 1, 7]:
            if self.taxon[i]: return self.taxon[i]

    def __repr__(self): return '<%s object id %s>' % (self.__class__.__name__, self.taxon[7])
    def __hash__(self): return hash(tuple(self.taxon))
    def __eq__(self, other): return self.taxon == other.taxon

################################################################################
class Fragment(object):
    """Once a short read is digested, a restriction fragment of a certain size is created.
    It has a length, is associated to a specific bacteria, contains the
    reference of the short read from which it came and the mapping
    quality of that short read.

    :param length: Length of the restriction fragment.
    :type  length: int
    :param bacteria: Taxonomic information.
    :type  bacteria: Bacteria
    :param ref: The reference name of the read.
    :type  ref: str
    :param read_length: The length of the read.
    :type  read_length: int
    :param sw_score: The Smith-Watermann score of the alignement of the original read.
    :type  sw_score: int
    """

    def __init__(self, length, bacteria, ref, read_length, sw_score):
        self.length = length
        self.bacteria = bacteria
        self.ref = ref
        self.sw_score = sw_score
        self.read_length = read_length
        self.sw_norm = '%.3g' % (float(sw_score) / read_length)

    def __repr__(self): return '<%s object "%s">' % (self.__class__.__name__, self.name)