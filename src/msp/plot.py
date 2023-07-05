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
                ('Portable Document Format', '.pdf'),
                ('Portable Network Graphics', '.png'),
                ('Scalable Vector Graphics', '.svg'),
            )
        )
        self.figure.figure_name = f
        self.figure.figure_format = f.split('.')[-1]
        self.figure.save_plot()
