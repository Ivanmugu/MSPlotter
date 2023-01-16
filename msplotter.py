# =============================================================================
# BSD 3-Clause License
# 
# Copyright (c) 2022, Ivan Munoz Gutierrez
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# =============================================================================

"""
Make a graphical representation of a blastn alignment.

Multiple Sequence Plotter (MSPlotter) uses GenBank files (.gb) to align the
sequences and plot the genes. To plot the genes, MSPlotter uses the information
from the `CDS` features section. To customize the colors for plotting genes,
you can add a `Color` tag in the `CDS` features with a color in hexadecimal.
For example, to show a gene in green add the tag `/Color="#00ff00"`. To reduce
the manual manipulation of the GenBank file, you can edit the file with
`Geneious` or another software and export the file with the new annotations.

MSPlotter uses `matplotlib`. Therefore, you can modify the parameters in the
`MakeFigure` class to customize your figure.

Usage examples
--------------
To make a figure with default parameters
$ python msplotter -i path/file_1.gb path/file_2.gb path/file_3.gb

To save a figure in pdf format
$ python msplotter -i path/file_1.gb path/file_2.gb path/file_3.gb -f pdf

Notes
-----
1. MSPlotter was designed to plot three sequences with lengths between 8 to 23
   kb. However, the matplotlib parameters can be adjusted for larger, smaller,
   or more sequences.
2. blastn must be installed locally and in the path.

Credits
-------
1. Inspired by easyfig
        Sullivan et al (2011) Bioinformatics 27(7):1009-1010
"""


import argparse
import os
import sys
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord

__version__ = '0.1.1'

class GenBankRecord:
    """Store relevant info from a GenBank file.

    Attributes
    ----------
    name : str
        Sequence name shown next to the `LOCUS` tag.
    accession : str
        Sequence accession number with version.
    description : str
        Description shown next to the `DEFINITION` tag.
    length : int
        Sequence length.
    sequence_start : int
        Start coordinate used for plotting. Default value is zero but changes
        if `alignments_position` of `MakeFigure` class is set to center or
        left.
    sequence_end : int
        End coordinate used for plotting. Default value is zero but changes if
        `alignments_position` of `MakeFigure` class is set to center or left.
    cds : list
        List of `CodingSequence` classes with info of `CDS` tags as product,
        start, end, strand, and color.
    num_cds : int
        Number of CDSs
    """

    def __init__(self, file_name):
        record = SeqIO.read(file_name, 'genbank')
        self.name = record.name
        self.accession = record.id
        self.description = record.description
        self.length = len(record)
        self.sequence_start = 0
        self.sequence_end = self.length
        self.cds = self.parse_gb(record)
        self.num_cds = len(self.cds)

    def parse_gb(self, record):
        """Parse gb file and make a list of `CodingSequence` classes.

        Parameters
        ----------
        record : Bio SeqIO.read object.

        Returns
        -------
        coding_sequences : list
            List of `CodingSequence` classes holding CDSs' information.
        """
        coding_sequences = []
        for feature in record.features:
            if feature.type != 'CDS':
                continue
            product = feature.qualifiers.get('product', None)
            if product != None:
                product = product[0]
            if feature.qualifiers.get('Color', None) is None:
                color = '#ffff00'    # Make yellow default color
            else:
                color = feature.qualifiers['Color'][0]
            strand = feature.strand
            if strand == -1:
                start = feature.location._end
                end = feature.location._start + 1
            else:
                start = feature.location._start + 1
                end = feature.location._end
            # Append cds
            coding_sequences.append(CodingSequence(
                product, start, end, strand, color
            ))
        return coding_sequences

class CodingSequence:
    """Store Coding Sequence (CDS) information from gb file."""
    def __init__(self, product, start, end, strand, color):
        self.product = product
        self.start = int(start)
        self.end = int(end)
        self.strand = int(strand)
        self.color = color

class BlastnAlignment:
    """Store blastn alignment results.

    Attributes
    ----------
    query_name : str
        Name of query sequence.
    hit_name : str
        Name of subject sequence.
    query_len : int
        Length of query sequence.
    hit_len : int
        Length of subject sequence.
    regions : list
        List of `RegionAlignmentResult` classes with info of aligned region as
        query_from, query_to, hit_from, hit_to, and identity.
    """

    def __init__(self, xml_alignment_result):
        with open(xml_alignment_result, 'r') as result_handle:
            blast_record = NCBIXML.read(result_handle)
            self.query_name = blast_record.query
            self.hit_name = blast_record.alignments[0].hit_def
            self.query_len = int(blast_record.query_length)
            self.hit_len = int(blast_record.alignments[0].length)
            self.regions = self.parse_blast_regions(blast_record)

    def parse_blast_regions(self, blast_record):
        """Parse blastn aligned regions to store the information.

        Parameters
        ----------
        blast_record : NCBIXML object
            Harbors blastn alignment results in xml format.

        Returns
        -------
        regions : list
            List of `RegionAlignmentResult` classes with info from alignment.
        """
        regions = []
        for region in blast_record.alignments[0].hsps:
            regions.append(RegionAlignmentResult(
                query_from=int(region.query_start),
                query_to=int(region.query_end),
                hit_from=int(region.sbjct_start),
                hit_to=int(region.sbjct_end),
                identity=int(region.identities),
                positive=int(region.positives),
                align_len=int(region.align_length)
            ))
        return regions

class RegionAlignmentResult:
    """Save blastn results of region that aligned."""
    def __init__(
        self, query_from, query_to, hit_from, hit_to, identity, positive,
        align_len
    ):
        self.query_from = query_from
        self.query_to = query_to
        self.hit_from = hit_from
        self.hit_to = hit_to
        self.identity = identity
        self.positive = positive
        self.align_len = align_len
        self.homology = identity / align_len

def make_fasta_file(gb_files):
    """Make fasta files from GenBank files and save them in local directory.

    Parameters
    ----------
    gb_files : list
        List of GenBank files.
    
    Returns
    -------
    faa_files : list
        List of fasta files names.
    """
    faa_files = []
    for gb_file in gb_files:
        record = SeqIO.read(gb_file, "genbank")
        faa_name = record.name + '.faa'
        new_record = SeqRecord(
            record.seq,
            id=record.id,
            description=record.description
        )
        SeqIO.write(new_record, faa_name, 'fasta')
        faa_files.append(faa_name)
    return faa_files

def run_blastn(faa_files):
    """Run blastn locally and create xml result file(s).

    Parameters
    ----------
    faa_files : list
        List of fasta files.
    
    Returns
    -------
    results : list
        List of xml files' names with blastn results.
    """
    results = []
    for i in range(len(faa_files) - 1):
        output_file = 'result' + str(i) + '.xml'
        blastn_cline = NcbiblastnCommandline(
            query=faa_files[i],
            subject=faa_files[i+1],
            outfmt=5,
            out=output_file)
        stdout, stderr = blastn_cline()
        results.append(output_file)
        print(
            f'BLASTing {faa_files[i]} (query) and {faa_files[i+1]} (subject)'
        )
        print(stdout + '\n' + stderr)
    return results

def delete_files(documents: list) -> None:
    """Delete the files from `documents` list."""
    for document in documents:
        if os.path.exists(document):
            os.remove(document)
        else:
            print(f"File {document} does not exist")

def get_alignment_records(alignment_files: list) -> list:
    """Parse xml alignment files and make list of `BlastnAlignment` classes."""
    alignments = [BlastnAlignment(alignment) for alignment in alignment_files]
    return alignments

def get_gb_records(gb_files: list) -> list:
    """Parse gb files and make list of `GenBankRecord` classes."""
    gb_records = [GenBankRecord(gb_file) for gb_file in gb_files]
    return gb_records

class MakeFigure:
    """Store relevant variables to plot the figure."""
    def __init__(
        self, alignments, gb_records, alignments_position="left",
        add_annotations_genes=False, add_annotations_sequences=False,
        sequence_name="accession", y_separation=10, sequence_color="black",
        sequence_width=3, identity_color="Blues", homology_padding=0.06
    ):
        self.alignments = alignments
        self.gb_records = gb_records
        self.alignments_position = alignments_position
        self.add_annotations_sequences = add_annotations_sequences
        self.add_annotations_genes = add_annotations_genes
        self.sequence_name = sequence_name
        self.y_separation = y_separation
        self.sequence_color = sequence_color
        self.sequence_width = sequence_width
        self.identity_color = identity_color
        self.homology_padding = y_separation * homology_padding
        self.size_longest_sequence = self.get_longest_sequence()

    def get_longest_sequence(self) -> int:
        """Find the longest sequence in gb_records"""
        longest = 0
        for record in self.gb_records:
            if record.length > longest:
                longest = record.length
        return longest

    def adjust_x_right(self, x1, x2):
        """Adjust x coordinates of CDS to alignt figure to the right."""
        difference = self.size_longest_sequence - x2
        return (x1 + difference, x2 + difference)

    def adjust_positions_sequences_right(self):
        """Adjust position of sequences to the right."""
        for record in self.gb_records:
            delta = self.size_longest_sequence - record.length
            record.sequence_start = record.sequence_start + delta
            record.sequence_end = record.sequence_end + delta
            for sequence in record.cds:
                sequence.start = sequence.start + delta
                sequence.end = sequence.end + delta

    def adjust_positions_alignments_right(self):
        """Adjust position of alignments to the right."""
        for alignment in self.alignments:
            delta_query = self.size_longest_sequence - alignment.query_len
            delta_hit = self.size_longest_sequence - alignment.hit_len
            for region in alignment.regions:
                region.query_from = region.query_from + delta_query
                region.query_to = region.query_to + delta_query
                region.hit_from = region.hit_from + delta_hit
                region.hit_to = region.hit_to + delta_hit

    def adjust_positions_sequences_center(self):
        """Adjust position of sequences to the center."""
        pass

    def adjust_position_aligments_center(self):
        """Adjust position if alignmets to the center."""
        pass

    def plot_dna_sequences(self, ax):
        """
        Plot lines that represent DNA sequences.

        Parameters
        ----------
        gb_records : list
            List of GenBankRecord classes
        y_separation : int
            Distance to separate genes in the y axis
        """
        y_distance = len(self.gb_records) * self.y_separation
        # Readjust position sequences to the right if requested
        if self.alignments_position == "right":
            self.adjust_positions_sequences_right()
        for gb_record in self.gb_records:
            x1 = gb_record.sequence_start
            x2 = gb_record.sequence_end
            x_values = np.array([x1, x2])
            y_values = np.array([y_distance, y_distance])
            ax.plot(
                x_values,
                y_values,
                linestyle='solid',
                color='black',
                linewidth=3,
                zorder=1
            )
            y_distance -= self.y_separation

    def truncate_colormap(self, cmap, min_val=0.0, max_val=1.0, n=100):
        name = 'trunc_cmap'
        truncated_cmap = cmap(np.linspace(min_val, max_val, n))
        new_cmap = colors.LinearSegmentedColormap.from_list(
            name,
            truncated_cmap
        )
        return new_cmap

    def plot_homology_regions(self, ax):
        """
        Plot homology regions of aligned sequences.

        Paramenters
        -----------
        alignment : BlastnAlignment class
        """
        y_distance = ((len(self.alignments) + 1) * self.y_separation)
        # make colormap for homology regions
        cmap = plt.colormaps[self.identity_color]
        # Extract subset of colormap for better representation of homologies
        # and update the `self.identity_color` with this new colormap.
        self.identity_color = self.truncate_colormap(
            cmap,
            min_val = 0.0,
            max_val = 0.75,
            n = 100
        )
        # Readjust position of alignment to right if requested
        if self.alignments_position == "right":
            self.adjust_positions_alignments_right()
        # plot regions with homology
        for alignment in self.alignments:
            for region in alignment.regions:
                # get region's coordinates
                x1 = region.query_from
                x2 = region.query_to
                x3 = region.hit_to
                x4 = region.hit_from
                y1 = y_distance - self.homology_padding
                y2 = y_distance - self.homology_padding
                y3 = y_distance - self.y_separation + self.homology_padding
                y4 = y_distance - self.y_separation + self.homology_padding
                xpoints = np.array([x1, x2, x3, x4, x1])
                ypoints = np.array([y1, y2, y3, y4, y1])
                ax.plot(xpoints, ypoints, linewidth=0)
                ax.fill(
                    xpoints,
                    ypoints,
                    facecolor=self.identity_color(region.homology),
                    linewidth=0
                )
            y_distance -= self.y_separation

    def plot_genes(self, ax):
        """
        Plot genes.

        Parameters
        ----------
        gb_records : list
            List of GenBankRecord classes
        y_separation : int
            Distance to separate genes in the y axis
        """
        # Separation of genes of each sequence in the y axis
        y_distance = len(self.gb_records) * self.y_separation
        arrowstyle = mpatches.ArrowStyle(
            "simple", head_width=0.5, head_length=0.2
        )
        # Iterate over GenBankRecords and plot genes
        for gb_record in self.gb_records:
            for gene in gb_record.cds:
                arrow = mpatches.FancyArrowPatch(
                    (gene.start, y_distance),
                    (gene.end, y_distance),
                    arrowstyle=arrowstyle,
                    color=gene.color,
                    mutation_scale=30,
                    zorder=2
                )
                ax.add_patch(arrow)
            y_distance -= self.y_separation

    def annotate_dna_sequences(self, ax):
        """
        Annotate DNA sequences.

        Parameters
        ----------
        gb_records : list
            List of GenBankRecord classes
        y_separation : int
            Distance to separate genes in the y axis
        """
        # Separation of annotations of each sequence in the y axis
        y_distance = len(self.gb_records) * self.y_separation
        for gb_record in self.gb_records:
            if self.sequence_name == 'accession':
                sequence_name = gb_record.accession
            else:
                sequence_name = gb_record.name
            ax.annotate(
                sequence_name,
                xy=(gb_record.sequence_end, y_distance),
                xytext=(10, 0),
                textcoords="offset points"
            )
            y_distance -= self.y_separation

    def annotate_genes(self, ax, gb_record, y_distance):
        """
        Annotate genes of DNA sequence.
        """
        # Annotate genes of first DNA sequence
        for gene in gb_record.cds:
            location_annotation = (gene.start + gene.end) / 2
            ax.annotate(
                gene.product,
                xy=(location_annotation, y_distance),
                xytext=(0, 15),
                textcoords="offset points",
                rotation=90,
                ha="center"
            )

    def make_figure(self, name='figure_1.svg', format='svg'):
        # Change figure size. Default size is 6.4 x 4.8
        fig, ax = plt.subplots(figsize=(8,4.8))
        # Plot DNA sequences
        self.plot_dna_sequences(ax)
        # Annotate DNA sequences
        if self.add_annotations_sequences:
            self.annotate_dna_sequences(ax)
        # Plot homology regions
        self.plot_homology_regions(ax)
        # Plot colorbar
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        fig.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=self.identity_color),
            ax=ax,
            fraction=0.02,
            shrink=0.4,
            pad=0.1,
            # label='Identity (%)',
            orientation='horizontal'
            # boundaries=np.array([0.4, 1])
        )
        # Plot genes
        self.plot_genes(ax)
        # Annotate genes
        if self.add_annotations_genes:
            # Annotate genes first sequence
            self.annotate_genes(
                ax,
                self.gb_records[0],
                self.y_separation * len(self.gb_records)
            )
            # Annotate genes last sequence
            self.annotate_genes(
                ax,
                self.gb_records[len(self.gb_records) - 1],
                self.y_separation
            )
        # # Increase limits in y axis ... I need to fix this part!
        # plt.ylim(
        #     self.y_separation - 5, (len(gb_files) * self.y_separation) + 5
        # )
        # Remove the x-y axis
        plt.axis('off')
        # Adjust the padding between and around subplots.
        plt.tight_layout()
        # Save figure as svg
        plt.savefig(name, format=format)
        # Show plot
        plt.show()

class UserInput:
    """Store information provided by the user via the command line."""
    def __init__(self):
        # Get user input
        info = self.user_input()
        # Store input files
        self.input_files = [Path(document) for document in info.input]
        # Check input files
        self.check_input_files()
        # Store output folder
        if info.output is not None:
            self.output_folder = Path(info.output)
        else:
            self.output_folder = Path('.')
        # Check output folder
        self.check_output_folder()
        # Store figure name
        if info.name is not None:
            self.fig_name = info.name
        else:
            self.fig_name = 'figure.svg'
        # Store figure format
        if info.format is not None:
            self.fig_format = info.format
        else:
            self.fig_format = 'svg'
        # Check that figure name extention matches figure format and make path
        # for output file
        if self.check_figure_extention() == 0:
            self.output_file = self.output_folder / self.fig_name
        else:
            name = self.fig_name + '.' + self.fig_format
            self.output_file = (self.output_folder / name)

    def user_input(self):
        """
        Parse and check command line arguments provided by user.

        Returns
        -------
        info : argparse object
            .input : holds the path to input files
            .output : holds sthe path of output figure
        """
        # Parse arguments and provide help.
        parser = argparse.ArgumentParser(
            prog='msplotter.py',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(__doc__)
        )
        # Make optional arguments.
        parser.add_argument(
            '-i', '--input', required=True, help='path to input files',
            nargs='+'
        )
        parser.add_argument('-o', '--output', help='path to output folder')
        parser.add_argument('-n', '--name', help='name of figure')
        parser.add_argument('-f', '--format', help='format of figure')
        # Parse command line arguments
        info = parser.parse_args()

        return info

    def check_input_files(self) -> bool:
        """Check if path to input files exists."""
        for document in self.input_files:
            if not document.exists():
                sys.exit(f'error: {document} does not exist')
            if not document.is_file():
                sys.exit(f'error: {document} is not a file')
        return True

    def check_output_folder(self) -> bool:
        """Check if output folder exists."""
        if not self.output_folder.exists():
            sys.exit(f'error: {self.output_folder} does not exist')
        elif not self.output_folder.is_dir():
            sys.exit(f'error: {self.output_folder} is not a directory')
        else:
            return True

    def check_figure_extention(self) -> int:
        """Check if figure extention matches figure format."""
        extention = self.fig_name.split('.')
        if len(extention) == 1:
            return 1
        if extention[1] != self.fig_format:
            sys.exit(f'error: file name extention does no match figure format')
        return 0

def main():
    # Get user input
    info = UserInput()
    gb_files = info.input_files
    # gb_files = ['seq1.gb', 'seq2.gb', 'seq3.gb', 'seq4.gb']
    # gb_files = ['CP072204.gb', 'SW4848.gb','CP000645.gb']

    # Create fasta files for BLASTing.
    faa_files = make_fasta_file(gb_files)
    # Run blastn locally.
    xml_results = run_blastn(faa_files)
    # Delete fasta files used for BLASTing.
    delete_files(faa_files)
    # Make a list of `BlastnAlignment` classes from the xml blastn results.
    alignments = get_alignment_records(xml_results)
    # Delete xml documents.
    delete_files(xml_results)
    # Make a list of `GenBankRecord` classes from the gb files.
    gb_records = get_gb_records(gb_files)
    # Make figure.
    figure = MakeFigure(
        alignments,
        gb_records,
        alignments_position="right",
        identity_color="Blues",
        add_annotations_genes=True,
        add_annotations_sequences=False,
        y_separation=10,
        homology_padding=0.08
    )
    figure.make_figure(name=info.output_file, format=info.fig_format)

if __name__ == '__main__':
    main()