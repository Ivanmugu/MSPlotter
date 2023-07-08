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
from pathlib import Path


class Plot(customtkinter.CTkToplevel):
    """Plot alignment."""
    def __init__(self, matplotlib_figure, msplotter_figure):
        """
        matplotlib_figure : matplotlib Figure object class
        msplotter_figure : msplotter Figure object class
        """
        super().__init__()
        # Set plot canvas and variables
        self.figure_plt = matplotlib_figure      # matplotlib object.
        self.figure_msp = msplotter_figure    # msplotter object.
        # Set layout parameters
        self.title("Graphic represenation of alignments")
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        # Create plot frame
        self.plot_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.plot_frame.grid(row=0, sticky='nsew')
        # Set canvas for plot
        self.canvas = FigureCanvasTkAgg(self.figure_plt, self.plot_frame)
        self.canvas.draw()
        self.plot = self.canvas.get_tk_widget()
        self.plot.pack(expand=True)

        # Create save frame
        self.save_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.save_frame.grid(row=1, sticky='nswe')
        self.save_frame.grid_columnconfigure(4, weight=1)
        # File type label
        self.file_type_label = customtkinter.CTkLabel(
            self.save_frame, text='File type:',
        )
        self.file_type_label.grid(row=0, column=0, padx=(10,5), pady=20)
        # Variable to store file_type_menu.
        self.file_type_var = customtkinter.StringVar(self, 'png')
        # File type menu
        self.file_type_menu = customtkinter.CTkOptionMenu(
            self.save_frame,
            values=['pdf', 'png', 'svg'],
            variable=self.file_type_var
        )
        self.file_type_menu.grid(
            row=0, column=1, padx=5, pady=20
        )
        # Resolution label
        self.resolution_label = customtkinter.CTkLabel(
            self.save_frame, text='Resolution:'
        )
        self.resolution_label.grid(
            row=0, column=2, padx=5, pady=20
        )
        # Variable to store resolution_menu
        self.resolution_var = customtkinter.StringVar(self, '300 DPI')
        # Resolution menu
        self.resolution_menu = customtkinter.CTkOptionMenu(
            self.save_frame,
            values=['72 DPI', '150 DPI', '300 DPI', '600 DPI'],
            variable=self.resolution_var
        )
        self.resolution_menu.grid(
            row=0, column=3, padx=5, pady=20
        )
        # Save button
        self.save_button = customtkinter.CTkButton(
            self.save_frame, text='Save', command=self.save_figure
        )
        self.save_button.grid(
            row=0, column=4, padx=5, pady=20
        )

    def save_figure(self):
        """Save plot."""
        f = filedialog.asksaveasfilename(
            initialdir='.',
            title='Save file as',
        )
        self.figure_msp.figure_name = Path(f + '.' + self.file_type_var.get())
        # self.figure.figure_format = f.split('.')[-1]
        self.figure_msp.figure_format = self.file_type_var.get()
        self.figure_msp.dpi = float(self.resolution_var.get().split()[0])
        self.figure_msp.save_plot()
