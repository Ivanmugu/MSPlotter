"""Parse command line arguments provided by user.

License
-------
This file is part of MSPloter
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
import argparse
from argparse import Namespace
import sys
from pathlib import Path
import pkg_resources
from typing import Union

# TODO: test it with app_cli() function.

class UserInput:
    """Store information provided by user via the command line."""
    def __init__(
            self,
            input_files: Union[list[Path], None] = None,
            output_folder: Union[None, Path] = None,
            figure_name: Union[None, str] = None,
            figure_format: Union[None, str] = None,
            dpi: Union[None, float] = None,
            output_path: Union[None, Path] = None,
            alignments_position: Union[None, str] = None,
            identity_color: Union[None, str] = None,
            annotate_sequences: Union[None, bool] = None,
            annotate_sequences_with: Union[None, str] = None,
            annotate_genes: Union[None, bool] = None,
            annotate_genes_on_sequence: Union[None, tuple[str]] = None,
            gui: Union[None, bool] = None,
    ):
        self.input_files = input_files
        self.output_folder = output_folder
        self.figure_name = figure_name
        self.figure_format = figure_format
        self.dpi = dpi
        self.output_path = output_path
        self.alignments_position = alignments_position
        self.identity_color = identity_color
        self.annotate_sequences = annotate_sequences
        self.annotate_sequences_with = annotate_sequences_with
        self.annotate_genes = annotate_genes
        self.annotate_genes_on_sequence = annotate_genes_on_sequence
        self.gui = gui


def parse_command_line_input() -> UserInput:
    """Parse command line arguments provided by user."""
    # Create parser.
    parser = argparse.ArgumentParser(
        add_help=False,
        prog='msplotter',
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "Make a graphical representation of a blastn alignment."
        )
    )
    # Make argument groups.
    helper = parser.add_argument_group('Help')
    required = parser.add_argument_group('Required')
    optional = parser.add_argument_group('Optional')
    gui = parser.add_argument_group('Graphical User Interfase')
    # ================== #
    # Help arguments     #
    # ================== #
    helper.add_argument(
        '-h', '--help', action='help',
        help='Show this help message and exit.'
    )
    prog_version = pkg_resources.get_distribution('msplotter').version
    helper.add_argument(
        '-v', '--version', action='version',
        version=f'%(prog)s {prog_version}',
        help="Show program's version number and exit"
    )
    # ================== #
    # Required arguments #
    # ================== #
    required.add_argument(
        '-i', '--input', nargs='+',
        help=('Path to input files. Provided files must be GenBank files.'), 
    )
    # ================== #
    # Optional arguments #
    # ================== #
    optional.add_argument(
        '-o', '--output',
        help=(
            'Path to output folder.\nDefault: current working directory.'
        )
    )
    optional.add_argument(
        '-n', '--name',
        help=(
            'Name of figure.\nDefault: `figure`.'
        )
    )
    optional.add_argument(
        '-f', '--format',
        help=(
            'Format of figure. Options: pdf, png, and svg.\n'
            'For a complete list of valid options visit:\n'
            'https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html\n'
            'Default: `png`.'
        )
    )
    optional.add_argument(
        '-d', '--dpi',
        help=(
            'Resolution in dots per inch.\n'
            'Default: 300 (high resolution for print).'
        )
    )
    optional.add_argument(
        '--alignments_position',
        help=(
            'Orientation of the alignments in the plot.\n'
            'Options: `left`, `center`, and `rigth`.\n'
            'Default: `left`.'
        )
    )
    optional.add_argument(
        '--identity_color',
        help=(
            'Color map representing homology regions.\n'
            'For a complete list of valid options visit:\n'
            'https://matplotlib.org/stable/tutorials/colors/colormaps.html\n'
            'Some options: `Greys`, `Purples`, `Blues`, and `Oranges`.\n'
            'Default: `Greys`.'
        )
    )
    optional.add_argument(
        '--annotate_sequences', nargs='?', const='accession',
        help=(
            'Annotate sequences in the plot.\n'
            'Options: `accession`, `name`, and `fname`.\n'
            '`accession` and `name` are obtained from the `ACCESSION`\n'
            'and `LOCUS` gb file tags, repectively. `fname` is the file\n'
            'name.\n'
            'If the flag is provided without argument, the sequences will\n'
            'be annotated using `accession` numbers.'
        )
    )
    optional.add_argument(
        '--annotate_genes', nargs='?', const='top',
        help=(
            'Annotate genes from top and bottom sequences.\n'
            'Options: `top`, `bottom`, and `both`.\n'
            'If the flag is provided without argument, only the genes at\n'
            'the top of the plot will be annotated.'
        )
    )
    # ================== #
    # GUI arguments      #
    # ================== #
    gui.add_argument(
        '-g', '--gui', help='Run app in a graphical user interface.',
        action='store_true'
    )

    # Parse command line arguments
    command_line_info = parser.parse_args()
    # Get command line arguments
    user_input = get_command_line_arguments(command_line_info)

    return user_input


def get_command_line_arguments(command_line_info: Namespace) -> UserInput:
    """Store the command line input into a UserInput class."""
    # Initiate UserInput class.
    user_input = UserInput()
    # =========================================================================
    # If GUI, check if user didn't provide extra arguments.
    # =========================================================================
    if gui := command_line_info.gui:
        user_input.gui = gui
        check_if_user_provided_only_gui_argument(command_line_info)
        return user_input
    # =========================================================================
    # If not gui, store and check all provided arguments.
    # =========================================================================
    # input_files is the only mandatory argument, unless user requests gui.
    if input_files := command_line_info.input:
        user_input.input_files = check_input_files(input_files)
    else:
        sys.exit('Error: the `--input` argument is mandatory.')
    # If user doesn't provide output folder, use current working directory.
    if output_folder := command_line_info.output:
        user_input.output_folder = check_output_folder(output_folder)
    else:
        user_input.output_folder = Path.cwd()
    # If user doesn't provide figure_name, use `figure` as name.
    if figure_name := command_line_info.name:
        user_input.figure_name = figure_name
    else:
        user_input.figure_name = 'figure'
    # If user doesn't provide figure format, use `png`.
    if figure_format := command_line_info.format:
        user_input.figure_format = figure_format
    else:
        user_input.figure_format = 'png'
    # If user provided figure_name and/or figure_format check correctness.
    check_figure_name_and_format(user_input)
    # If user doesn't provide dpi, make 300 as default.
    if dpi := command_line_info.dpi:
        user_input.dpi = float(dpi)
    else:
        user_input.dpi = 300.0
    # Create output_path using output_folder, figure_name, and figure_format
    user_input.output_path = make_output_path(user_input)
    # If user didn't provide alignments position, use 'left'
    if alignments_position := command_line_info.alignments_position:
        user_input.alignments_position = check_alignments_position(
                alignments_position
            )
    else:
        user_input.alignments_position = 'left'
    # The MakeFigure class with the function make_colormap of the msplotter
    # module will check the correctness of identity_color.
    if identity_color := command_line_info.identity_color:
        user_input.identity_color = identity_color
    else:
        user_input.identity_color = 'Greys'
    # If annotate_sequences has any value, then user_input.annotate_sequence is
    # set to True, and user_input.annotate_sequences_with gets the value.
    # Otherwise, set user_input.annotate_sequences to False.
    if annotate_sequences := command_line_info.annotate_sequences:
        user_input.annotate_sequences = True
        user_input.annotate_sequences_with = check_annotate_sequences(
                annotate_sequences
            )
    else:
        user_input.annotate_sequences = False
    # If annotate_genes has any value, then user_input.annotate_genes is set to
    # True and user_input.annotate_genes_on_sequence gets the value.
    # Otherwise, set user_input.annotate_genes to False
    # The get_annotate_genes_info function checks the arguments.
    if annotate_genes := command_line_info.annotate_genes:
        user_input.annotate_genes = True
        user_input.annotate_genes_on_sequence = get_annotate_genes_info(
                annotate_genes
            )
    else:
        user_input.annotate_genes = False

    return user_input


def check_if_user_provided_only_gui_argument(
        command_line_info: Namespace
    ) -> None:
    """Make sure that user do not provide extra flags when activating GUI."""
    arguments = vars(command_line_info)
    counter = 0
    for _, value in arguments.items():
        if value:
            counter += 1
    # if there are more than one arguments, then exit with error.
    if counter > 1:
        sys.exit(
            'Error: too many arguments provided.\n' +
            'If you want to activate the GUI, only provide the `--gui` flag.'
        )


def check_mandatory_arguments(user_input: Namespace) -> None:
    """Check mandatory arguments and conflicts with gui."""
    if not user_input.input_files and not user_input.gui:
        sys.exit('Error: the `--input` argument is mandatory.')
    if (
        (user_input.input_files and user_input.gui)
    ):
        sys.exit(
            'Error: too many arguments.\nIf you want to activate the GUI '
            'you only need to provide the `--gui` flag.\n'
            'Otherwise, no need to provide the `--gui` flag.'
        )


def check_input_files(input_files: list[str]) -> list[Path]:
    """Check input files."""
    # Convert list of input files into list of paths.
    input_files = [Path(infile) for infile in input_files]
    # Check if paths exist and if they are paths to files.
    for document in input_files:
        if not document.exists():
            sys.exit(f'Error: `{document}` does not exist')
        if not document.is_file():
            sys.exit(f'Error: `{document}` is not a file')
    return input_files


def check_output_folder(output_folder: str) -> Path:
    """Check output folder."""
    # Convert output output folder into Path
    output_folder = Path(output_folder)
    # Check output folder.
    if not output_folder.is_dir():
        sys.exit(f'Error: `{output_folder}` is not a directory')
    return output_folder


def check_figure_name_and_format(user_input: Namespace) -> None:
    """Check if figure extention matches figure format."""
    # Check if figure name has extension.
    extension = user_input.figure_name.split('.')
    # If figure name does not have extension add it with format.
    if len(extension) == 1:
        user_input.figure_name = (
            user_input.figure_name + '.' + user_input.figure_format
        )
    # If figure name has extension check if it matches figure format.
    elif extension[-1] != user_input.figure_format:
        sys.exit(
            'Error: file name extension does no match figure format.\n'
            'The default format is `png`. If you want a different format '
            'provide it using the `--format` flag.'
        )


def make_output_path(user_input: Namespace) -> Path:
    """Make output path."""
    if len(user_input.figure_name.split('.')) == 0:
        figure_name = user_input.figure_name + '.' + user_input.figure_format
        output_path = user_input.output_folder / figure_name
    else:
        output_path = user_input.output_folder / user_input.figure_name
    return output_path


def check_alignments_position(alignments_position: str) -> Union[None, str]:
    """Check position to align the alignments."""
    # Check that user enter correct parameter.
    if (
        alignments_position == "left"
        or alignments_position == "center"
        or alignments_position == "right"
    ):
        return alignments_position
    else:
        sys.exit(
            f'Error: parameter `alignment_position: {alignments_position}` '
            'is not valid.\n'
            'Valid parameters are: `left`, `center`, or `right`.'
        )


def check_annotate_sequences(annotate_sequences: str) -> Union[None, str]:
    """Check correct input for annotate sequences."""
    if (
        annotate_sequences == 'accession'
        or annotate_sequences == 'name'
        or annotate_sequences == 'fname'
    ):
        return annotate_sequences
    else:
        sys.exit(
            f'Error: parameter `annotate_sequence: {annotate_sequences}` is '
            'not valid.\n'
            'Valid parameter are: `accession`, `name`, or `fname`.'
        )


def get_annotate_genes_info(annotate_genes: str) -> Union[None, tuple[str]]:
    """Get information to annotate genes in plot."""
    if annotate_genes == 'top':
        return ('top')
    elif annotate_genes == 'bottom':
        return ('bottom')
    elif annotate_genes == 'both':
        return ('top', 'bottom')
    else:
        sys.exit(
            f'Error: parameter `annotate_genes: {annotate_genes}` is not '
            'valid.\n'
            'Valid parameters are: `top`, `bottom` or `both`.'
        )
