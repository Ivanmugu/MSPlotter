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

License
-------
This file is part of MSPloter
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
import os
import sys
from pathlib import Path
from typing import Union

import matplotlib as mpl
from matplotlib.axes import Axes
import matplotlib.colors as colors
from matplotlib.colors import Colormap
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.ticker as ticker
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
            if gene := feature.qualifiers.get('gene', None):
                gene = gene[0]
            if product := feature.qualifiers.get('product', None):
                product = product[0]
            if color := feature.qualifiers.get('Color', None):
                color = feature.qualifiers['Color'][0]
            else:
                color = '#ffff00'    # Make yellow default color
            # Some CDS are composed of more than one parts, like introns, or,
            # in the case of some bacteria, some genes have frameshifts as a
            # regulatory function (some transposase genes have frameshifts as
            # a regulatory function).
            for part in feature.location.parts:
                strand = part._strand
                if strand == -1:
                    start = part._end
                    end = part._start + 1
                else:
                    start = part._start + 1
                    end = part._end
                # Append cds.
                coding_sequences.append(CodingSequence(
                    gene, product, start, end, strand, color
                ))
        return coding_sequences


class CodingSequence:
    """Store Coding Sequence (CDS) information from gb file."""
    def __init__(self, gene, product, start, end, strand, color):
        self.gene = gene
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
            f'BLASTing {faa_files[i]} (query) and {faa_files[i+1]} (subject)\n'
        )
        # print(stdout + '\n' + stderr)
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
    annotate_genes : bool
        Annotate genes in plot (default: False).
    annotate_genes_on_sequence : tuple with the str `top`, and/or `bottom`
        Annotate genes at the top, bottom or both.
    annotate_sequences : bool
        Annotate sequences in plot (defalult: False). If True, the `top` and
        `bottom` sequences are annotated.
    sequence_name : str
        Option: `name`, `fname`, and `accession` (default: `accession`).
        If `name` or `accession` is provided, their values are extracted from
        the GenBankRecord class. If `fname` is provided, the value is obtained
        from the file name. If `add_annotations_sequences` attribute is
        True, either `name`, `fname` or `accession` is used for annotating the
        sequence.
    y_separation : float
        Distance between sequences in the y-axis (default: 10).
    y_limit : float
        Lower limit to show in the plot (default: 5). Use this value to adjust
        the position of the color map and the scale bar.
    sequence_color : str
        Color used for lines representing sequences (default: `black`). You can
        use any color name allowed by `Matplotlib`.
    sequence_width : int
        Width of lines representing sequences (default: 3).
    identity_color : str
        Color used to represent regions of homology (default: `Greys`). This
        color represent a `Matplotlib` colormap. Therefore, you should provide
        a valid colormap name.
    scale_bar : bool
        Show scale bar (default: True).
    color_map_range : tupple
        Tyupple with a min and max value between 0 and 1. The min and max
        values are used in color_map to determine the range of the color_map
        to use (default: (0, 0.75)).
    color_map : matplotlib colormaps object
        Color map to represent homology regions.
    homology_padding : float
        Padding between lines representing sequences and regions of homology
        (default: 0.1). The number provided represents a fraction of the
        y_separation.
    save_figure : bool
        Save plotted figure.
    figure_name : Path
        Output path with figure's name.
    figure_format : str
        Format to save figure.
    figure_width : Union[None, float]
        Width of figure in inches (default: None). If None, MSPlotter
        determines the value.
    figure_height : Union[None, float]
        Height of figure in inches (default: None). If None, MSPlotter
        determines the value.
    dpi : float
        Resolution of figure in dots per inch (default: 300.0).
    user_gui : bool
        Run app in a graphical user interface (default: False).
    """
    def __init__(
        self,
        alignments,
        gb_records,
        alignments_position: str = "left",
        annotate_genes: bool = False,
        annotate_genes_on_sequence: tuple[str] = ("top", "bottom"),
        annotate_genes_from: str = "gene_tag",
        annotate_sequences: bool = False,
        sequence_name: str = "accession",
        y_separation: int = 10,
        y_limit: int = 5,
        sequence_color: str = "black",
        sequence_width: int = 3,
        identity_color: str = "Greys",
        scale_bar: bool = True,
        color_map_range: tuple[float] = (0, 0.75),
        homology_padding: float = 0.1,
        figure_name: Path = (Path.cwd() / 'figure.png'),
        figure_format: str = 'png',
        figure_width: Union[None, float] = None,
        figure_height: Union[None, float] = None,
        dpi: float = 300.0,
        use_gui: bool = False
    ):
        self.alignments = alignments
        self.num_alignments = len(alignments)
        self.gb_records = gb_records
        self.alignments_position = alignments_position
        self.annotate_genes = annotate_genes
        self.annotate_genes_on_sequence = annotate_genes_on_sequence
        self.annotate_genes_from = annotate_genes_from
        self.annotate_sequences = annotate_sequences
        self.sequence_name = sequence_name
        self.y_separation = y_separation
        self.y_limit = y_limit
        self.sequence_color = sequence_color
        self.sequence_width = sequence_width
        self.identity_color = identity_color
        self.scale_bar = scale_bar
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
        self.figure_width = figure_width
        self.figure_height = figure_height
        self.dpi = dpi
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

    def plot_dna_sequences(self, ax: Axes) -> None:
        """Plot lines that represent DNA sequences."""
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

    def draw_scale_bar(self, ax: Axes, bar_position: int = 0) -> None:
        """Draw a horizontal scale bar for DNA length."""
        # Get x_ticks to calculate sequence length.
        x_ticks = ax.get_xticks()
        len_ticks = x_ticks[1] - x_ticks[0]
        # Draw scale_bar.
        ax.plot(
            (0, len_ticks),
            (bar_position, bar_position),
            linestyle='solid',
            color='black',
            linewidth=2,
            zorder=1
        )
        # Annotate scale_bar.
        location_annotation = len_ticks / 2
        ax.annotate(
            f'{self.get_units_size_dna(len_ticks)}',
            fontsize=8,
            xy=(location_annotation, bar_position),
            xytext=(0, -9),
            textcoords='offset points',
            ha='center',
        )

    def get_units_size_dna(self, len_dna: float) -> str:
        """Get the units of the dna length."""
        if len_dna < 1_000_000 and len_dna >= 1_000:
            num = len_dna / 1_000
            return str(round(num, 1)) + ' kbp'
        elif len_dna < 1_000_000_000 and len_dna >= 1_000_000:
            num = len_dna / 1_000_000
            return str(round(num, 1)) + ' Mbp'
        else:
            return str(round(len_dna)) + ' bp'

    def make_colormap(
            self, identity_color: str, min_val: float = 0.0, max_val:
            float = 1.0, n: int = 100
        ) -> Colormap:
        """Make color map for homology regions."""
        try:
            cmap = plt.colormaps[identity_color]
        except KeyError:
            sys.exit(
                f"Error: the identity color '{identity_color}' provided is "
                "not valid.\nUse the help option to find valid colors or "
                "visit:\n"
                "https://matplotlib.org/stable/tutorials/colors/colormaps.html"
                "\nfor a complete list of valid options."
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

    def plot_homology_regions(self, ax: Axes) -> None:
        """Plot homology regions of aligned sequences."""
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

    def draw_colorbar(self, fig, ax: Axes) -> None:
        """Draw color bar for homology regions."""
        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        print(
            "lowest and highest plot colorbar:",
            self.lowest_homology, self.highest_homology
        )
        # if self.lowest_homology != self.highest_homology:
        #     boundaries = np.linspace(
        #         self.lowest_homology, self.highest_homology, 100
        #     )
        #     pos = ax.get_position()
        #     print(pos)
        #     ax.set_position([0, 0, 0.9, 1])
        #     axins = inset_axes(
        #         ax,
        #         width="2%",  # width: 5% of parent_bbox width
        #         height="20%",  # height: 50%
        #         loc="lower left",
        #         bbox_to_anchor=(1.0, 0.01, 1, 1),
        #         bbox_transform=ax.transAxes,
        #         borderpad=0,
        #     )
        #     colormap = fig.colorbar(
        #             plt.cm.ScalarMappable(norm=norm, cmap=self.color_map),
        #             cax=axins,                 # Axes to draw colormap.
        #             # shrink=0.18,
        #             # aspect=10,
        #             orientation='vertical',
        #             boundaries=boundaries,
        #             ticks=[self.lowest_homology, self.highest_homology],
        #     )
        #     colormap.ax.tick_params(labelsize=8)
        #     colormap.set_label(
        #         label="Identity (%)",
        #         size=7,
        #         labelpad=-7,
        #     )
        #     colormap.outline.set_visible(False) # Remove colormap frame.
        #     axins.set_aspect(8)                 # Set aspect of colormap.
        if self.lowest_homology != self.highest_homology:
            boundaries = np.linspace(
                self.lowest_homology, self.highest_homology, 100
            )
            # ax.set_position([0, 0, 0.9, 1])
            axins = inset_axes(
                ax,
                width="15%",  # width: 5% of parent_bbox width
                height="2.5%",  # height: 50%
                loc="lower right",
                bbox_to_anchor=(0., 0.04, 0.95, 0.95),
                bbox_transform=ax.transAxes,
                borderpad=0,
            )
            colormap = fig.colorbar(
                    plt.cm.ScalarMappable(norm=norm, cmap=self.color_map),
                    cax=axins,                 # Axes to draw colormap.
                    orientation='horizontal',
                    boundaries=boundaries,
                    ticks=[self.lowest_homology, self.highest_homology],
            )
            colormap.ax.tick_params(labelsize=8)
            # Format ticks 
            colormap.ax.xaxis.set_major_formatter(
                ticker.FuncFormatter(lambda x, pos: f'{x}%')
            )
            colormap.set_label(
                label="Identity",
                size=7,
                labelpad=-7,
            )
            colormap.outline.set_visible(False) # Remove colormap frame.
            axins.set_aspect(0.1)               # Set aspect of colormap.
        else:
            homology_path = mpatches.Patch(
                color=self.color_map(self.highest_homology/100),
                label=f'{self.highest_homology}%',
            )
            legend = ax.legend(
                loc="center left",
                bbox_to_anchor=(0.8, 0.04),
                handles=[homology_path],
                frameon=False,
                title="Identity",
                fontsize=8
            )
            legend.get_title().set_fontsize(8)

    def plot_genes(self, ax: Axes) -> None:
        """Plot genes."""
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

    def plot_arrows(self, ax: Axes) -> None:
        """Plot arrows to reprent genes."""
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

    def annotate_dna_sequences(self, ax: Axes) -> None:
        """Annotate DNA sequences."""
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
                    f'Error: invalid sequence name `{sequence_name}` for '
                    'annotating sequences.')
            ax.annotate(
                sequence_name,
                xy=(gb_record.sequence_end, y_distance),
                xytext=(10, -4),
                textcoords="offset points"
            )
            y_distance -= self.y_separation

    def annotate_gene_sequences(self, ax: Axes) -> None:
        """Annotate genes of DNA sequence.

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
                if self.annotate_genes_from == "gene_tag":
                    annotation = gene.gene
                else:
                    annotation = gene.product
                ax.annotate(
                    annotation,
                    xy=(location_annotation, parameters["y_distance"]),
                    xytext=(0, parameters["y_text"]),
                    textcoords="offset points",
                    rotation=90,
                    ha="center",
                    horizontalalignment=parameters["h_alignment"],
                    verticalalignment=parameters["v_alignment"]
                )

    def determine_figure_size(
            self, num_alignments: int
        ) -> tuple[float, float]:
        """Determine figure size.

        If user doesn't provide height, adjust height by number of alignments.
        """
        # The Matplotlib default figure size is 6.4 x 4.8 inches.
        if not self.figure_width:
            width = 6.4
        else:
            width = self.figure_width
        if not self.figure_height:
            auto_height = True
            height = 4.8
        else:
            auto_height = False
            height = self.figure_height
        # Adjust height based on four alignments.
        if (auto_height) and (num_alignments > 4):
            height = height * (num_alignments / 4)
        return (width, height)

    def make_figure(self):
        """Make figure with matplotlib."""
        # -- Remove toolbar from plot -----------------------------------------
        mpl.rcParams['toolbar'] = 'None'
        # -- Determine figure size --------------------------------------------
        width, height = self.determine_figure_size(self.num_alignments)
        # -- Change figure size -----------------------------------------------
        fig, ax = plt.subplots(
            # Matplotlib default size is 6.4 x 4.8 inches
            figsize=(width, height),
            layout="constrained"
        )
        # -- Remove figure axis -----------------------------------------------
        ax.set_axis_off()
        # -- Plot DNA sequences -----------------------------------------------
        self.plot_dna_sequences(ax)
        # -- Plot genes using the Arrow class ---------------------------------
        self.plot_arrows(ax)
        # -- Plot homology regions --------------------------------------------
        self.plot_homology_regions(ax)
        # -- Draw colorbar ----------------------------------------------------
        self.draw_colorbar(fig, ax)
        # -- Annotate DNA sequences -------------------------------------------
        if self.annotate_sequences:
            self.annotate_dna_sequences(ax)
        # -- Annotate genes ---------------------------------------------------
        if self.annotate_genes:
            self.annotate_gene_sequences(ax)
        # -- Plot scale bar ---------------------------------------------------
        if self.scale_bar and not self.annotate_genes:
            self.draw_scale_bar(ax, bar_position=self.y_limit)
        # If annotate genes is activated with scale bar, provide more space for
        # annotations.
        elif self.scale_bar and self.annotate_genes:
            self.draw_scale_bar(ax, bar_position=self.y_limit - 5)
        # If neither scale bar nor annotate ges are activated, set the y limit
        # to zero. This will give space for the colorbar.
        else:
            ax.set_ylim(bottom=self.y_limit)

        return fig

    def check_save_figure(self) -> bool:
        """Check if save figure is True."""
        if (self.figure_format is None) and (self.figure_name is None):
            return False
        else:
            return True

    def save_plot(self) -> None:
        """Save plot."""
        plt.savefig(
            fname=self.figure_name, format=self.figure_format, dpi=self.dpi
        )

    def display_figure(self) -> None:
        """Display and save figure."""
        # Adjust the padding between and around subplots.
        plt.tight_layout()
        plt.subplots_adjust(
            left=0, bottom=0, right=1, top=1, wspace=0, hspace=0
        )
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
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
        dpi=user_input.dpi,
        annotate_genes=user_input.annotate_genes,
        annotate_genes_on_sequence=user_input.annotate_genes_on_sequence,
        annotate_sequences=user_input.annotate_sequences,
        sequence_name=user_input.annotate_sequences_with
    )
    figure.make_figure()
    figure.display_figure()
