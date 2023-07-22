<p align="center">
   <img src="./src/msp/images/logo.png" alt="MSPlotter" width="350">
</p>

# Make a graphical representation of a blastn alignment

Multiple Sequence Plotter (MSPlotter) uses GenBank files (.gb) to align the
sequences and plot the genes. MSPlotter uses the information from the `CDS`
features section to plot the genes. To customize the colors for plotting genes,
you can add a `Color` tag in the `CDS` features with color in hexadecimal.
For example, add the tag `/Color="#00ff00"` to show a gene in green. To avoid
direct manual manipulation of the GenBank file, you can edit the file with
`Geneious` or another software and export the file with the new annotations.

MSPlotter is an easy-to-use option for people with little coding knowledge or
problems running the classic app `easyfig`. The program offers a graphical user
interface (GUI) and a command line interface (CLI).

If the user knows Python, MSPlotter uses `matplotlib`. Therefore, to customize
your figure, the user can modify the parameters in the `MakeFigure` class of
the `msplotter` module or make any necessary change in any part of the code.

I am developing MSPlotter in my free time, therefore, if you find a bug, it may
take me some time to fix it, but I will try to do my best to fix the problems
as soon as possible. Also, if you have any suggestions, let me know, and I will
try to implement them.

## Requirements

- [Python](https://www.python.org/) 3.11 or later
- [biopython](https://biopython.org/) 1.81 or later
- [customtkinter](https://customtkinter.tomschimansky.com/) 5.1 or later
- [matplotlib](https://matplotlib.org/) 3.7 or later
- [blastn](https://www.ncbi.nlm.nih.gov/books/NBK569861/) must be installed
  locally and in the path

MSPlotter has been tested in macOS and Windows.

## Installation

First, create a virtual environment with `conda` or `venv`. Then, install
msplotter using pip as follows:

```bash
pip install msplotter
```

## Usage and options

To run the GUI type:

```bash
msplotter --gui
```

Output GUI:

<p align="center">
   <img src="./src/msp/images/MSPlotter_gui.png" alt="MSPlotter" width="350">
</p>

To view all the options run:

```bash
msplotter --help
```

Partial output CLI:

```console
usage: msplotter [-h] [-v] [-i INPUT [INPUT ...]] [-o OUTPUT] [-n NAME] [-f FORMAT] [-d DPI]
                 [--alignments_position ALIGNMENTS_POSITION] [--identity_color IDENTITY_COLOR]
                 [--annotate_sequences [ANNOTATE_SEQUENCES]] [--annotate_genes [ANNOTATE_GENES]] [-g]

Make a graphical representation of a blastn alignment.

Help:
  -h, --help            Show this help message and exit.
  -v, --version         Show program's version number and exit

Required:
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Path to input files. Provided files must be GenBank files.

Optional:
  -o OUTPUT, --output OUTPUT
                        Path to output folder.
                        Default: current working directory.
  -n NAME, --name NAME  Name of figure.
                        Default: `figure`.
  -f FORMAT, --format FORMAT
                        Format of figure. Options: pdf, png, and svg.
                        For a complete list of valid options visit:
                        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html
                        Default: `png`.
  -d DPI, --dpi DPI     Resolution in dots per inch.
                        Default: 300 (high resolution for print).
  --alignments_position ALIGNMENTS_POSITION

```

## Usage examples CLI

To make a figure with default parameters:

```bash
msplotter -i path/file_1.gb path/file_2.gb path/file_3.gb
```

To save a figure in pdf format:

```bash
msplotter -i path/file_1.gb path/file_2.gb path/file_3.gb -f pdf
```

## Notes

I started this project to make a figure paper with three sequences with lengths
between 8 to 23 kb. However, the matplotlib parameters can be adjusted for
larger, smaller, or more sequences.

## Credits

Inspired by easyfig: Sullivan et al (2011) Bioinformatics 27(7):1009-1010

## License

BSD 3-Clause License
