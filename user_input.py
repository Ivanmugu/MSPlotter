# =============================================================================
# This file is part of MSPloter
#
# BSD 3-Clause License
#
# Copyright (c) 2022, Ivan Munoz Gutierrez
#
# A full description of the license is given in the LICENSE file at:
# https://github.com/ivanmugu/MSPlotter
# =============================================================================

import argparse
import sys
from pathlib import Path


def user_input():
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
        add_help=False,
        prog='msplotter.py',
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Make a graphical representation of a blastn alignment."
        )
    )
    # Make argument groups.
    # required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument(
        '-g', '--gui', help='Run app as graphic user interface.',
        action='store_true'
    )
    optional.add_argument(
        '-i', '--input',
        # required=True,
        help='Path to input files.',
        nargs='+'
    )
    # Optional aguments.
    optional.add_argument(
        '-h', '--help', action='help',
        help='Show this help message and exit.'
    )
    optional.add_argument('-o', '--output', help='Path to output folder.')
    optional.add_argument('-n', '--name', help='Name of figure.')
    optional.add_argument('-f', '--format', help='Format of figure.')
    optional.add_argument(
        '--alignments_position',
        help=(
            'Orientation of the alignments in the plot.\n' +
            'Options: `left`, `center`, and `rigth`.\n' +
            'Default: `left`.'
        )
    )
    optional.add_argument(
        '--identity_color',
        help=(
            'Color map representing homology regions.\n' +
            'For a complete list of valid options visit:\n' +
            'https://matplotlib.org/stable/tutorials/colors/colormaps.html\n' +
            'Some options: `Greys`, `Purples`, `Blues`, and `Oranges`.\n' +
            'Default: `Greys`.'
        )
    )
    optional.add_argument(
        '--annotate_sequences', nargs='?', const='accession',
        help=(
            'Annotate sequences in the plot. If argument is not ' +
            'provided,\n' +
            'sequences will be annotated using `accession` numbers.\n' +
            'Options: `accession`, `name`, and `fname`.\n' +
            '`accession` and `name` are obtained from the `ACCESSION`\n' +
            'and `LOCUS` gb file tags, repectively.\n' +
            '`fname` is the file name.'
        )
    )
    optional.add_argument(
        '--annotate_genes', nargs='?', const='top',
        help=(
            'Annotate genes from top and bottom sequences. If argument is ' +
            'not provided,\n' +
            'the genes at the top of the plot will be annotated.\n' +
            'Options: `top`, `bottom`, and `both`.'
        )
    )
    # Parse command line arguments
    info = parser.parse_args()

    return info


class UserInput:
    """Store information provided user via the command line.

    Attributes
    ----------
    input_files : list
        List of input files' paths as Path objects.
    output_folder : Path object
        Path to output folder (default: current folder).
    figure_name : str
        Name of figure (default: `figure_1.svg`).
    figure_format : str
        Format to make and save figure (default: `svg`).
    output_file : Path object
        Path, including the name and format, of figure.
    alignments_position : str
        Position of the alignments in plot (default: `left`).
    identity_color : str
        Color of shadows representing homology regions (default: `Greys`).
    """
    def __init__(self, user_input):
        """
        Parameters
        ----------
        user_input : argparse object
            Harbors the command line information parsed by argparse.
        """
        self.input_files = self.get_input_files(user_input)
        self.output_folder = self.get_output_folder(user_input)
        self.figure_name = self.get_figure_name(user_input)
        self.figure_format = self.get_figure_format(user_input)
        self.output_file = self.make_output_path(user_input)
        self.alignments_position = self.get_alignments_position(user_input)
        self.identity_color = self.get_identity_color(user_input)
        self.annotate_seq_info = self.get_annotate_sequences_info(user_input)
        self.annotate_sequences = self.annotate_seq_info[0]
        self.sequence_name = self.annotate_seq_info[1]
        self.annotate_genes_info = self.get_annotate_genes_info(user_input)
        self.annotate_genes = self.annotate_genes_info[0]
        self.annotate_genes_on_sequence = self.annotate_genes_info[1]

    def get_input_files(self, user_info) -> list:
        """Get input files and return a list of Path objects."""
        input_files = [Path(document) for document in user_info.input]
        # Check if path to input files exists.
        for document in input_files:
            if not document.exists():
                sys.exit(f'error: {document} does not exist')
            if not document.is_file():
                sys.exit(f'error: {document} is not a file')
        return input_files

    def get_output_folder(self, user_info) -> Path:
        """Get output folder from user input and check if exists."""
        # Get output folder from user
        if user_info.output is not None:
            output_folder = Path(user_info.output)
        else:
            output_folder = Path('.')
        # Check output folder
        if not output_folder.exists():
            sys.exit(f'error: {output_folder} folder does not exist')
        if not output_folder.is_dir():
            sys.exit(f'error: {output_folder} is not a directory')
        return output_folder

    def get_figure_name(self, user_info) -> str:
        """Get figure name from user."""
        if user_info.name is not None:
            figure_name = user_info.name
        else:
            figure_name = 'figure_1.pdf'
        return figure_name

    def get_figure_format(self, user_info) -> str:
        """Get figure format from user."""
        if user_info.format is not None:
            figure_format = user_info.format
        else:
            figure_format = 'pdf'
        return figure_format

    def make_output_path(self, user_info) -> Path:
        """Make output path."""
        if self.check_figure_extention(user_info) == 0:
            output_file = self.output_folder / self.figure_name
        else:
            name = self.figure_name + '.' + self.figure_format
            output_file = self.output_folder / name
        return output_file

    def check_figure_extention(self, user_info) -> int:
        """Check if figure extention matches figure format."""
        # If figure format provided check if name was also provided
        if user_info.name is None and user_info.format is not None:
            sys.exit('error: figure format provided but name not provided')
        # Check if figure name and format match.
        extention = self.figure_name.split('.')
        if len(extention) == 1:
            return 1
        if extention[1] != self.figure_format:
            sys.exit('error: file name extention does no match figure format')
        return 0

    def get_alignments_position(self, user_info) -> str:
        """Check if parameter alignment position is valid."""
        if user_info.alignments_position is not None:
            position = user_info.alignments_position
        else:
            position = 'left'
        # Check that user enter correct parameter.
        if position == "left" or position == "center" or position == "right":
            return position
        else:
            sys.exit(
                f'Error: parameter `alignment_position: {position}` is not ' +
                'valid.\n' +
                'Valid parameters are: `left`, `center`, or `right`.'
            )

    def get_identity_color(self, user_info) -> str:
        """Get color that represent homolgoy regions from user.

        Note
        ----
        The MakeFigure class with the function make_colormap will check color
        input.
        """
        if user_info.identity_color is not None:
            identity_color = user_info.identity_color
        else:
            identity_color = 'Greys'
        return identity_color

    def get_annotate_sequences_info(self, user_info) -> tuple:
        """Get information to annotate sequences in plot."""
        #TODO it looks like fname is off
        annotate = user_info.annotate_sequences
        if annotate is None:
            return (False, '')
        if annotate == 'accession' or annotate == 'name' or annotate == 'fname':
            return (True, annotate)
        else:
            sys.exit(
                f'Error: parameter `annotate_sequence: {annotate}` is not ' +
                'valid.\n' +
                'Valid parameter are: `accession` or `name`.'
            )

    def get_annotate_genes_info(self, user_info) -> tuple:
        """Get information to annotate genes in plot."""
        annotate = user_info.annotate_genes
        if annotate is None:
            return (False, ())
        if annotate == 'top':
            return (True, ('top',))
        elif annotate == 'bottom':
            return (True, ('bottom',))
        elif annotate == 'both':
            return (True, ('top', 'bottom'))
        else:
            sys.exit(
                f'Erro: parameter `annotate_genes: {annotate}` is not ' +
                'valid.\n' +
                'Valid parameters are: `top`, `bottom` or `both`.'
            )
