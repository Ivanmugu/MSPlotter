# msplotter

Make a graphical representation of a blantn alignment

Multiple Sequence Plotter (MSPlotter) uses GenBank files (.gb) to align the
sequences and plot the genes. To plot the genes, MSPlotter uses the information
from the `CDS` features section. To customize the colors for plotting genes,
you can add a `Color` tag in the `CDS` features with a color in hexadecimal.
For example, to show a gene in green add the tag `/Color="#00ff00"`. To reduce
the manual manipulation of the GenBank file, you can edit the file with
`Geneious` or another software and export the file with the new annotations.

MSPlotter uses `matplotlib`. Therefore, you can modify the parameters in the
`MakeFigure` class to customize your figure.

## Usage examples

To make a figure with default parameters

```bash
python msplotter -i path/file_1.gb path/file_2.gb path/file_3.gb
```

To save a figure in pdf format

```bash
$ python msplotter -i path/file_1.gb path/file_2.gb path/file_3.gb -f pdf
```

## Notes

1. MSPlotter was designed to plot three sequences with lengths between 8 to 23
   kb. However, the matplotlib parameters can be adjusted for larger, smaller,
   or more sequences.
2. blastn must be installed locally and in the path.

## Credits

1. Inspired by easyfig
   Sullivan et al (2011) Bioinformatics 27(7):1009-1010

## License

BSD 3-Clause License
