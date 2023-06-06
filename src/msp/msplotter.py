"""Make a graphical representation of a blastn alignment.

Multiple Sequence Plotter (MSPlotter) uses GenBank files (.gb) to align the
sequences and plot the genes. To plot the genes, MSPlotter uses the information
from the `CDS` features section. To customize the colors for plotting genes,
you can add a `Color` tag in the `CDS` features with a color in hexadecimal.
For example, to show a gene in green add the tag `/Color="#00ff00"`. To reduce
the manual manipulation of the GenBank file, you can edit the file with
`Geneious` or another software and export the file with the new annotations.

MSPlotter uses the `matplotlib` library. Therefore, you can modify the
parameters in the `MakeFigure` class to customize your figure.

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
        https://mjsull.github.io/Easyfig/

License
-------
This file is part of MSPloter
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
import os
import sys

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord

from msp.arrows import Arrow


class GenBankRecord:
    """Store relevant info from a GenBank file.

    Attributes
    ----------
    file_name : str
        file name.
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
        """
        file_name : Path oject
            Path to file.
        """
        record = SeqIO.read(file_name, 'genbank')
        self.file_name = file_name.stem
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
            if product is not None:
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
            # Append cds.
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
    """Store relevant variables to plot the figure.

    Attributes
    ----------
    alignments : list
        List of `BlastnAlignment` classes created from xml blastn results.
    num_alignments : int
        Number of alignments.
    gb_records : list
        List of `GenBankRecord` classes created from gb files.
    alignments_position : str
        Position of the alignments in the plot (default: `left`).
    add_annotations_genes : Bool
        Annotate genes in plot (default: False).
    add_annotations_sequences : Bool
        Annotate sequences in plot (default: False).
    sequence_name : str
        String to access either `name` or `accession` from GenBankRecord class
        (default: `accession`). Either `name` or `accession` is used to
        annotate the sequence if `add_annotations_sequences` attribute is True.
    y_separation : int
        Distance between sequences in the y-axis (default: 10).
    sequence_color : str
        Color used for lines representing sequences (default: `black`). You can
        use any color name allowed by `Matplotlib`.
    sequence_width : int
        Width of lines representing sequences (default: 3).
    identity_color : str
        Color used to represent regions of homology (default: `Greys`). This
        color represent a `Matplotlib` colormap. Therefore, you should provide
        a valid colormap name.
    homology_padding : float
        Padding between lines representing sequences and regions of homology
        (default: 0.1). The number provided represents a fraction of the
        y_separation.
    save_figure : bool
        Save plotted figure (deault: `True`).
    figure_name : str
        Name to save figure (default: `figure_1.pdf`).
    figure_format : str
        Format to save figure (default: `pdf`).
    """
    def __init__(
        self, alignments, gb_records, alignments_position="left",
        annotate_genes=False, annotate_genes_on_sequence=("top", "bottom"),
        annotate_sequences=False, sequence_name="accession", y_separation=10,
        sequence_color="black", sequence_width=3, identity_color="Greys",
        color_map_range=(0, 0.75), homology_padding=0.1, figure_name=None,
        figure_format=None, use_gui=False
    ):
        self.alignments = alignments
        self.num_alignments = len(alignments)
        self.gb_records = gb_records
        self.alignments_position = alignments_position
        self.annotate_genes = annotate_genes
        self.annotate_genes_on_sequence = annotate_genes_on_sequence
        self.annotate_sequences = annotate_sequences
        self.sequence_name = sequence_name
        self.y_separation = y_separation
        self.sequence_color = sequence_color
        self.sequence_width = sequence_width
        self.identity_color = identity_color
        self.color_map_range = color_map_range
        self.color_map = self.make_colormap(
            identity_color=identity_color,
            min_val=self.color_map_range[0],
            max_val=self.color_map_range[1],
            n=100
        )
        self.homology_padding = y_separation * homology_padding
        self.size_longest_sequence = self.get_longest_sequence()
        self.lowest_homology = int(
            round(self.get_lowest_and_highest_homology()[0] * 100)
        )
        self.highest_homology = int(
            round(self.get_lowest_and_highest_homology()[1] * 100)
        )
        self.figure_name = figure_name
        self.figure_format = figure_format
        self.save_figure = self.check_save_figure()
        self.use_gui = use_gui

    def get_lowest_and_highest_homology(self) -> tuple:
        """Get the lowest and highest homologies in the alignment."""
        lowest = 100
        highest = 0
        for alignment in self.alignments:
            for region in alignment.regions:
                if region.homology < lowest:
                    lowest = region.homology
                if region.homology > highest:
                    highest = region.homology
        return (lowest, highest)

    def get_longest_sequence(self) -> int:
        """Find the longest sequence in gb_records."""
        longest = 0
        for record in self.gb_records:
            if record.length > longest:
                longest = record.length
        return longest

    def adjust_positions_sequences_right(self) -> None:
        """Adjust position of sequences to the right including CDSs."""
        for record in self.gb_records:
            delta = self.size_longest_sequence - record.length
            record.sequence_start = record.sequence_start + delta
            record.sequence_end = record.sequence_end + delta
            for sequence in record.cds:
                sequence.start = sequence.start + delta
                sequence.end = sequence.end + delta

    def adjust_positions_alignments_right(self) -> None:
        """Adjust position of alignments to the right."""
        for alignment in self.alignments:
            delta_query = self.size_longest_sequence - alignment.query_len
            delta_hit = self.size_longest_sequence - alignment.hit_len
            for region in alignment.regions:
                region.query_from = region.query_from + delta_query
                region.query_to = region.query_to + delta_query
                region.hit_from = region.hit_from + delta_hit
                region.hit_to = region.hit_to + delta_hit

    def adjust_positions_sequences_center(self) -> None:
        """Adjust position of sequences to the center including CDSs."""
        for record in self.gb_records:
            shift = (self.size_longest_sequence - record.length) / 2
            record.sequence_start = record.sequence_start + shift
            record.sequence_end = record.sequence_end + shift
            for sequence in record.cds:
                sequence.start = sequence.start + shift
                sequence.end = sequence.end + shift

    def adjust_positions_alignments_center(self) -> None:
        """Adjust position of alignmets to the center."""
        for alignment in self.alignments:
            shift_q = (self.size_longest_sequence - alignment.query_len) / 2
            shift_h = (self.size_longest_sequence - alignment.hit_len) / 2
            for region in alignment.regions:
                region.query_from = region.query_from + shift_q
                region.query_to = region.query_to + shift_q
                region.hit_from = region.hit_from + shift_h
                region.hit_to = region.hit_to + shift_h

    def plot_dna_sequences(self, ax) -> None:
        """Plot lines that represent DNA sequences.

        Parameters
        ----------
        ax : axes, matplotlib object
        """
        y_distance = len(self.gb_records) * self.y_separation
        # Readjust position sequences to the right or center if requested.
        if self.alignments_position == "right":
            self.adjust_positions_sequences_right()
        elif self.alignments_position == 'center':
            self.adjust_positions_sequences_center()
        # Plot lines representing sequences.
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
                linewidth=2,
                zorder=1
            )
            y_distance -= self.y_separation

    def make_colormap(self, identity_color, min_val=0.0, max_val=1.0, n=100):
        """Make color map for homology regions."""
        try:
            cmap = plt.colormaps[identity_color]
        except KeyError:
            sys.exit(
                f"Error: the identity color '{identity_color}' provided is " +
                "not valid.\nUse the help option to find valid colors or " +
                "visit:\n" +
                "https://matplotlib.org/stable/tutorials/colors/colormaps.html\n" +
                "for a complete list of valid options."
            )
        if min_val == 0 and max_val == 1:
            return cmap
        else:
            name = 'trunc_cmap'
            truncated_cmap = cmap(np.linspace(min_val, max_val, n))
            new_cmap = colors.LinearSegmentedColormap.from_list(
                name,
                truncated_cmap
            )
            return new_cmap

    def plot_homology_regions(self, ax):
        """Plot homology regions of aligned sequences.

        Parameters
        ----------
        ax : axes, matplotlib object
        """
        y_distance = ((len(self.alignments) + 1) * self.y_separation)
        # Readjust position of alignment to right or center if requested.
        if self.alignments_position == "right":
            self.adjust_positions_alignments_right()
        elif self.alignments_position == "center":
            self.adjust_positions_alignments_center()
        # Plot regions with homology.
        for alignment in self.alignments:
            for region in alignment.regions:
                # Get region's coordinates.
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
                    facecolor=self.color_map(region.homology),
                    linewidth=0
                )
            y_distance -= self.y_separation

    def plot_colorbar(self, ax):
        """Plot color bar for homology regions.

        Parameters
        ----------
        ax : axes, matplotlib object
        """
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        print(
            "lowest and highest plot colorbar:",
            self.lowest_homology, self.highest_homology
        )
        if self.lowest_homology != self.highest_homology:
            # Place colorbar inside figure.
            # axins = inset_axes(
            #     ax,
            #     width='20%',
            #     height='50%',
            #     loc='center'
            # )
            boundaries = np.linspace(
                self.lowest_homology, self.highest_homology, 100
            )
            colormap = plt.colorbar(
                    plt.cm.ScalarMappable(norm=norm, cmap=self.color_map),
                    ax=ax,                 # Axes to draw colormap
                    shrink=0.25,
                    aspect=10,
                    label='Identity (%)',
                    orientation='horizontal',
                    boundaries=boundaries,
                    ticks=[self.lowest_homology, self.highest_homology]
            )
            colormap.outline.set_visible(False) # Remove colormap frame
            print("Making colormap...")
        else:
            homology_path = mpatches.Patch(
                color=self.color_map(self.highest_homology/100),
                label=f'{self.highest_homology}',
            )
            ax.legend(
                loc="center",
                handles=[homology_path],
                frameon=False,
                title="Identity (%)"
            )

    def plot_genes(self, ax):
        """Plot genes.

        Parameters
        ----------
        ax : axes, matplotlib object
        """
        # Separation of genes of each sequence in the y axis.
        y_distance = len(self.gb_records) * self.y_separation
        arrowstyle = mpatches.ArrowStyle(
            "simple", head_width=0.5, head_length=0.2
        )
        # Iterate over GenBankRecords and plot genes.
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

    def plot_arrows(self, ax):
        """Plot arrows to reprent genes.

        Parameters
        ----------
        ax : axes, matplotlib object
        """
        # Separation of genes of each sequence in the y axis.
        y_distance = len(self.gb_records) * self.y_separation
        # Ratio head_height vs lenght of longest sequence.
        ratio = 0.02
        head_height = self.size_longest_sequence * ratio
        # Iterate over GenBankRecords and plot genes.
        for gb_record in self.gb_records:
            for gene in gb_record.cds:
                arrow = Arrow(
                    x1=gene.start,
                    x2=gene.end,
                    y=y_distance,
                    head_height=head_height
                )
                x_values, y_values = arrow.get_coordinates()
                ax.fill(x_values, y_values, gene.color)
            y_distance -= self.y_separation

    def annotate_dna_sequences(self, ax):
        """Annotate DNA sequences.

        Parameters
        ----------
        ax : axes, matplotlib object
        """
        # Separation of annotations of each sequence in the y axis.
        y_distance = len(self.gb_records) * self.y_separation
        # Annotate sequences.
        for gb_record in self.gb_records:
            # Check if sequence name is valis.
            if self.sequence_name == 'accession':
                sequence_name = gb_record.accession
            elif self.sequence_name == 'name':
                sequence_name = gb_record.name
            elif self.sequence_name == 'fname':
                sequence_name = gb_record.file_name
            else:
                sys.exit(
                    f'Error: invalid sequence name `{sequence_name}` for ' +
                    'annotating sequences.')
            ax.annotate(
                sequence_name,
                xy=(gb_record.sequence_end, y_distance),
                xytext=(10, 0),
                textcoords="offset points"
            )
            y_distance -= self.y_separation

    def annotate_gene_sequences(self, ax):
        """Annotate genes of DNA sequence.
        
        Parameters
        ----------
        ax : axes, matplotlib object
        Note
        ----
        This function annotates genes only of top and bottom sequences.
        """
        # Define dictionaries with parameters to annotate top and bottom
        # sequences. `y_text` indicates how far the annotation is going to be
        # from the position of the gene. `h_alignment` and `v_alignment`
        # indicate the position to rotate the text of the annotations.
        # `gb_record` is a GenBankRecord class that has the coordinates of the
        # genes to annotate. `y_distance` is the position of the sequence in
        # the y-axis.
        top = {
            "y_text": 13,
            "h_alignment": 'left',
            "v_alignment": 'bottom',
            "gb_record": self.gb_records[0],
            "y_distance": self.y_separation * len(self.gb_records)
        }
        bottom = {
            "y_text": -13,
            "h_alignment": 'right',
            "v_alignment": 'top',
            "gb_record": self.gb_records[len(self.gb_records) - 1],
            "y_distance": self.y_separation
        }
        for position in self.annotate_genes_on_sequence:
            if position == 'top':
                parameters = top
            elif position == 'bottom':
                parameters = bottom
            for gene in parameters["gb_record"].cds:
                location_annotation = (gene.start + gene.end) / 2
                ax.annotate(
                    gene.product,
                    xy=(location_annotation, parameters["y_distance"]),
                    xytext=(0, parameters["y_text"]),
                    textcoords="offset points",
                    rotation=90,
                    ha="center",
                    horizontalalignment=parameters["h_alignment"],
                    verticalalignment=parameters["v_alignment"]
                )

    def determine_figure_size(self, num_alignments) -> tuple:
        """Determine figure size."""
        # The Matplotlib default figure size is 6.4 x 4.8.
        width = 6.4
        height = 4.8
        # Increase height after four alignments.
        if num_alignments > 4:
            height = height * (num_alignments / 4)
        return (width, height)

    def make_figure(self):
        """Make figure with matplotlib."""
        # Remove toolbar from plot.
        mpl.rcParams['toolbar'] = 'None'
        # Determine figure size by number of alignments.
        width, height = self.determine_figure_size(self.num_alignments)
        # Change figure size. Matplotlib default size is 6.4 x 4.8.
        if self.lowest_homology == self.highest_homology:
            fig, (ax_1, ax_2) = plt.subplots(
                2, 1,
                figsize=(width, height),
                layout="constrained",
                gridspec_kw={'height_ratios': [15, 1]}
            )
            # Remove axis
            ax_1.set_axis_off()
            ax_2.set_axis_off()
        else:
            fig, ax_1 = plt.subplots(
                figsize=(width, height),
                layout="constrained"
            )
            ax_1.set_axis_off()
        # Plot DNA sequences.
        self.plot_dna_sequences(ax_1)
        # Annotate DNA sequences.
        if self.annotate_sequences:
            self.annotate_dna_sequences(ax_1)
        # Plot homology regions.
        self.plot_homology_regions(ax_1)
        # Plot colorbar.
        if self.lowest_homology == self.highest_homology:
            self.plot_colorbar(ax_2)
        else:
            self.plot_colorbar(ax_1)
        # Plot genes using the Arrow class.
        self.plot_arrows(ax_1)
        # Annotate genes.
        if self.annotate_genes:
            self.annotate_gene_sequences(ax_1)
        return fig

    def check_save_figure(self) -> bool:
        """Check if save figure is True."""
        if (self.figure_format is None) and (self.figure_name is None):
            return False
        else:
            return True

    def save_plot(self) -> None:
        """Save plot."""
        plt.savefig(fname=self.figure_name, format=self.figure_format)

    def display_figure(self) -> None:
        """Display and save figure."""
        # # Adjust the padding between and around subplots.
        # plt.tight_layout()
        # Save figure.
        if self.save_figure and not self.use_gui:
            self.save_plot()
        # Show plot.
        plt.show()

    def close_figure(self) -> None:
        """Close maltplotlib."""
        plt.close()

def app_cli(user_input) -> None:
    """Run msplotter in the command line interface.

    Parameters
    ----------
    user_input : UserInput class object.
    """
    # Get list of input files' paths.
    gb_files = user_input.input_files
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
        alignments_position=user_input.alignments_position,
        identity_color=user_input.identity_color,
        figure_name=user_input.output_path,
        figure_format=user_input.figure_format,
        annotate_genes=user_input.annotate_genes,
        annotate_genes_on_sequence=user_input.annotate_genes_on_sequence,
        annotate_sequences=user_input.annotate_sequences,
        sequence_name=user_input.annotate_sequences_with
    )
    figure.make_figure()
    figure.display_figure()
