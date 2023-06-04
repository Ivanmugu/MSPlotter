"""Class to plot alignments.

License
-------
This file is part of MSPloter
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
from tkinter import filedialog
import customtkinter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt


class Plot(customtkinter.CTkToplevel):
    """Plot alignment."""
    def __init__(self, matplotlib_figure, msplotter_figure):
        """
        matplotlib_figure : matplotlib Figure object class
        msplotter_figure : msplotter Figure object class
        """
        super().__init__()
        # Set plot canvas and variables
        self.fig = matplotlib_figure      # matplotlib object.
        self.figure = msplotter_figure    # msplotter object.
        # Set canvas for plot
        self.title("Graphic represenation of alignments")
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.plot = self.canvas.get_tk_widget()
        self.plot.pack(side='top', fill='both', expand=True, padx=10, pady=10)

        # Save button
        self.save_button = customtkinter.CTkButton(
            self, text='Save', command=self.save_figure
        )
        self.save_button.pack(
            side='bottom', pady=10
        )

    def save_figure(self):
        """Save plot."""
        f = filedialog.asksaveasfilename(
            initialdir='.',
            title='Save file as',
            filetypes=(
                ('Encapsulated Postcript', '.eps'),
                ('Joint Photographic Experts Group', '.jpg'),
                ('Joint Photographic Experts Group', '.jpeg'),
                ('Portable Document Format', '.pdf'),
                ('PGF code for LaTeX', '.pgf'),
                ('Portable Network Graphics', '.png'),
                ('Postscript', '.ps'),
                ('Raw RGBA bitmap', '.raw'),
                ('Raw RGBA bitmap', '.rgba'),
                ('Scalable Vector Graphics', '.svg'),
                ('Scalable Vector Graphics', '.svgz'),
                ('Tagged Image File Format', '.tif'),
                ('Tagged Image File Format', '.tiff'),
                ('WevP Image Format', '.webp')
            )
        )
        self.figure.figure_name = f
        self.figure.figure_format = f.split('.')[-1]
        self.figure.save_plot()
