<p align="center">
   <img src="./src/msp/images/logo.png" alt="MSPlotter" width="350">
</p>

# Make a graphical representation of a blantn alignment

Multiple Sequence Plotter (MSPlotter) uses GenBank files (.gb) to align the
sequences and plot the genes. To plot the genes, MSPlotter uses the information
from the `CDS` features section. To customize the colors for plotting genes,
you can add a `Color` tag in the `CDS` features with a color in hexadecimal.
For example, to show a gene in green add the tag `/Color="#00ff00"`. To avoid
direct the manual manipulation of the GenBank file, you can edit the file with
`Geneious` or another software and export the file with the new annotations.

MSPlotter uses `matplotlib`. Therefore, to customize your figure, you can
modify the parameters in the `MakeFigure` class of the `msplotter` module.

## Requirements

- [Python](https://www.python.org/) 3.11 or later
- [biopython](https://biopython.org/) 1.81 or later
- [customtkinter](https://customtkinter.tomschimansky.com/) 5.1 or later
- [matplotlib](https://matplotlib.org/) 3.7 or later
- [blastn](https://www.ncbi.nlm.nih.gov/books/NBK569861/) must be installed locally and in the path

## Installation

Create a virtual environment and install msplotter using pip as follows:

```bash
pip intall msplotter
```

## Usage and options

To view all the options run:

```bash
msplotter --help
```

Output:

```console
usage: msplotter [-h] [-v] [-i INPUT [INPUT ...]] [-o OUTPUT] [-n NAME] [-f FORMAT]
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
                        Format of figure.
                        Default: `pdf`.
  --alignments_position ALIGNMENTS_POSITION
                        Orientation of the alignments in the plot.
                        Options: `left`, `center`, and `rigth`.
                        Default: `left`.
  --identity_color IDENTITY_COLOR
                        Color map representing homology regions.
                        For a complete list of valid options visit:
                        https://matplotlib.org/stable/tutorials/colors/colormaps.html
                        Some options: `Greys`, `Purples`, `Blues`, and `Oranges`.
                        Default: `Greys`.
  --annotate_sequences [ANNOTATE_SEQUENCES]
                        Annotate sequences in the plot.
                        Options: `accession`, `name`, and `fname`.
                        `accession` and `name` are obtained from the `ACCESSION`
                        and `LOCUS` gb file tags, repectively. `fname` is the file
                        name.
                        If the flag is provided without argument, the sequences will
                        be annotated using `accession` numbers.
  --annotate_genes [ANNOTATE_GENES]
                        Annotate genes from top and bottom sequences.
                        Options: `top`, `bottom`, and `both`.
                        If the flag is provided without argument, only the genes at
                        the top of the plot will be annotated.

Graphic User Interfase:
  -g, --gui             Run app in a graphic user interface.
```

## Usage examples

To make a figure with default parameters:

```bash
msplotter -i path/file_1.gb path/file_2.gb path/file_3.gb
```

To save a figure in pdf format:

```bash
msplotter -i path/file_1.gb path/file_2.gb path/file_3.gb -f pdf
```

If you don't like the terminal and prefer a graphical user interface:

```bash
msplotter --gui
```

## Notes

I started this project to make a figure paper with three sequences with lengths
between 8 to 23 kb. However, the matplotlib parameters can be adjusted for
larger, smaller, or more sequences.

## Credits

Inspired by easyfig: Sullivan et al (2011) Bioinformatics 27(7):1009-1010

## License

BSD 3-Clause License
