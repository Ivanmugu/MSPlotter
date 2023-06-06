"""Class to make the GUI.

License
-------
This file is part of MSPloter
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
from pathlib import Path

from tkinter import filedialog
import customtkinter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import msp.msplotter as msp
from msp.colormap_picker import ColormapPicker
from msp.plot import Plot


class App(customtkinter.CTk):
    """msplotter GUI."""
    def __init__(self):
        super().__init__()
        # Variables for BLASTing and plotting.
        self.fig = None                        # matplotlib object.
        self.gb_files: list = None
        self.figure = None                     # msplotter object.
        self.identity_color: str = "Greys"
        self.colormap_range: tuple = (0, 0.75)
        self.annotate_sequences: bool = False
        self.annotate_genes: bool = False
        # Set layout parameters
        self.geometry('600x520')
        self.title('MSPlotter')
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(1, weight=1)
        # Protocol to close app, including plot if exist.
        self.protocol('WM_DELETE_WINDOW', self.on_closing)
        # Variable to store the ColormapPicker class used in the
        # launch_colormap_picker function.
        self.colormap_app = None

        # ################ #
        # Navigation frame #
        # ################ #
        # Create navigation frame.
        self.navigation_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.navigation_frame.grid(row=0, column=0, columnspan=3, sticky='nsew')
        self.navigation_frame.grid_columnconfigure(0, weight=1)
        # Logo label
        self.logo_label = customtkinter.CTkLabel(
            self.navigation_frame, text='MSPloter',
            font=customtkinter.CTkFont(size=20, weight='bold')
        )
        self.logo_label.grid(row=0, column=0, padx=(10,5), pady=20)
        # Select button
        self.select_button = customtkinter.CTkButton(
            self.navigation_frame, text='Select',
            command=self.get_files_path
        )
        self.select_button.grid(row=0, column=1, padx=5, pady=20)
        # Clear button
        self.clear_button = customtkinter.CTkButton(
            self.navigation_frame, text='Clear',
            command=self.clear_input,
            state='disabled'
        )
        self.clear_button.grid(row=0, column=2, padx=5, pady=20)
        # Plot button
        self.plot_button = customtkinter.CTkButton(
            self.navigation_frame, text='Plot',
            command=self.plot_figure,
            state='disabled'
        )
        self.plot_button.grid(row=0, column=3, padx=(5,10), pady=20)

        # ################ #
        # Appearance frame #
        # ################ #
        self.appearance_frame = customtkinter.CTkFrame(
            self, corner_radius=5,
        )
        self.appearance_frame.grid(
            row=1, column=0, padx=(10, 5), pady=10, sticky='nsew'
        )
        self.appearance_frame.grid_rowconfigure(9, weight=1)
        # Appearance label
        self.system_appearance_mode = customtkinter.get_appearance_mode()
        if self.system_appearance_mode == 'Dark':
            self.appearance_fg_color = '#333333'
        else:
            self.appearance_fg_color = 'Gray80'
        self.appearance_label = customtkinter.CTkLabel(
            self.appearance_frame,
            text='Appearance',
            font=customtkinter.CTkFont(size=16),
            width=120,
            height=50,
            corner_radius=5,
            fg_color=self.appearance_fg_color,
        )
        self.appearance_label.grid(
            row=0, column=0, pady=(0,10), sticky='nswe'
        )
        # Homology color label
        self.homology_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Homology color:',
        )
        self.homology_label.grid(row=1, column=0, pady=(10,0))
        # Colormap picker
        self.format_button = customtkinter.CTkButton(
            self.appearance_frame,
            text='Choose colormap',
            command=self.launch_colormap_picker
        )
        self.format_button.grid(row=2, column=0, padx=20, pady=(0, 10))
        # Align plot label
        self.align_plot_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Align plot:'
        )
        self.align_plot_label.grid(row=3, column=0, pady=(10, 0))
        # Variable to store align_plot menu selection
        self.align_plot_var = customtkinter.StringVar(self, 'Left')
        # Align plot menu
        self.align_plot = customtkinter.CTkOptionMenu(
            self.appearance_frame,
            values=['Left', 'Center', 'Right'],
            variable=self.align_plot_var,
        )
        self.align_plot.grid(row=4, column=0, padx=20, pady=(0, 10))
        # Annotate sequences label
        self.annotate_seq_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Annotate sequences:'
        )
        self.annotate_seq_label.grid(row=5, column=0, pady=(10, 0))
        # Variable to store annotate_seq menu selection
        self.annotate_seq_var = customtkinter.StringVar(self, 'No')
        # Annotate sequences menu
        self.annotate_seq_menu = customtkinter.CTkOptionMenu(
            self.appearance_frame,
            values=['No', 'Yes'],
            variable=self.annotate_seq_var,
            command=lambda _:self.update_annotate_seq_var()
        )
        self.annotate_seq_menu.grid(row=6, column=0, pady=(0,10))
        # Annotate genes label
        self.annotate_genes_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Annotate genes:'
        )
        self.annotate_genes_label.grid(row=7, column=0, pady=(10, 0))
        # Variable to store annotate_genes menu selection
        self.annotate_genes_var = customtkinter.StringVar(self, 'No')
        # Annotate genes menu
        self.annotate_genes_menu = customtkinter.CTkOptionMenu(
            self.appearance_frame,
            values=['No', 'Yes'],
            variable=self.annotate_genes_var,
            command=lambda _:self.update_annotate_genes_var()
        )
        self.annotate_genes_menu.grid(row=8, column=0, pady=(0,10))
        # Reset button
        self.reset_button = customtkinter.CTkButton(
            self.appearance_frame, text='Reset', command=self.reset_appearance
        )
        self.reset_button.grid(row=9, column=0, pady=(20, 10))

        # ############# #
        # Display frame #
        # ############# #
        # Create display frame
        self.display_frame = customtkinter.CTkFrame(
            self, corner_radius=0,
            fg_color='transparent'
        )
        self.display_frame.grid(row=1, column=1, columnspan=2, sticky='nsew')
        self.display_frame.grid_rowconfigure(0, weight=1)
        self.display_frame.grid_columnconfigure(0, weight=1)
        self.display_window = customtkinter.CTkTextbox(
            self.display_frame, wrap='word'
        )
        self.display_window.grid(row=0, column=0, padx=(5, 10), pady=10,
            sticky='nsew'
        )
        self.display_window.insert(
            'end',
            '\n\n\n\nWelcome to MSPlotter!\n\n' +
            'Select your GenBank sequences, plot, and save.\n\n'
            "If you don't like the default parameters, change the appearance."
        )
        self.display_window.configure(state='disabled')

    # =========================================================================
    # Functionality
    # =========================================================================
    def append_paths_gb_files(self, input_files: tuple) -> None:
        """Append paths of gb files into self.gb_files list."""
        print(input_files)
        if self.gb_files is None:
            self.gb_files = [Path(element) for element in input_files]
        else:
            for element in input_files:
                self.gb_files.append(Path(element))

    def get_files_path(self) -> None:
        """Get paths of gb files and append them to self.gb_files list.

        The names of the gb files that are going to be analyzed are printed in
        the display frame.
        """
        self.append_paths_gb_files(filedialog.askopenfilenames(
            initialdir=".",
            title="Select a GenBank file",
            filetypes=(('GenBank files', '*.gb'), ('All files', '*.*'))
            ))
        self.display_window.configure(state='normal')
        self.display_window.delete('1.0', 'end')
        # Print files to be analyzed
        self.display_window.insert(
            'end', 'Files are going to be BLASTed in the next order:\n')
        for i, file_path in enumerate(self.gb_files):
            self.display_window.insert(
                'end', f'{i+1} --> {file_path.name}\n')
        self.display_window.configure(state='disabled')
        self.clear_button.configure(state='normal')
        self.plot_button.configure(state='normal')

    def clear_input(self):
        """Clean input data store in self.gb_files."""
        self.gb_files = None
        self.display_window.configure(state='normal')
        self.display_window.delete('1.0', 'end')
        self.display_window.configure(state='disabled')
        self.plot_button.configure(state='disabled')
        self.save_button.configure(state='disabled')
        self.clear_button.configure(state='disabled')

    def launch_colormap_picker(self):
        """Pick colormap and range for homology regions."""
        if self.colormap_app is None or not self.colormap_app.winfo_exists():
            self.colormap_app = ColormapPicker(self.get_colormap_data)
        else:
            self.colormap_app.focus()

    def get_colormap_data(self, cmap_name, cmap_range):
        """Get colormap data from colormap_picker."""
        self.identity_color = cmap_name
        self.colormap_range = (
            round(cmap_range[0])/100, round(cmap_range[1])/100
        )
        self.display_window.configure(state='normal')
        self.display_window.delete('1.0', 'end')
        self.display_window.insert(
            'end',
            f'The homology regions are going to be shown in `{cmap_name}` ' +
            f'with a range of colors between `{round(cmap_range[0])}-'
            f'{round(cmap_range[1])}`.'
        )
        self.display_window.configure(state='disabled')
        # print(f'Selected cmap: {self.identity_color}')
        # print(f'Selected cmap range: {self.colormap_range}')

    def update_annotate_seq_var(self):
        if self.annotate_seq_var.get() == 'No':
            self.annotate_sequences = False
        else:
            self.annotate_sequences = True

    def update_annotate_genes_var(self):
        if self.annotate_genes_var.get() == 'No':
            self.annotate_genes = False
        else:
            self.annotate_genes = True

    def reset_appearance(self):
        """Reset appearance parameters."""
        # Reset color map
        self.identity_color = "Greys"
        self.colormap_range = (0, 0.75)
        # Reset align plot
        self.align_plot_var.set('Left')
        self.align_plot.configure(variable=self.align_plot_var)
        # Reset annotate sequences
        self.annotate_seq_var.set('No')
        self.annotate_seq_menu.configure(variable=self.annotate_seq_var)
        # Reset Annotate genes
        self.annotate_genes_var.set('No')
        self.annotate_genes_menu.configure(variable=self.annotate_genes_var)

    def plot_figure(self):
        """Plot alignments using msplotter."""
        # Create fasta files for BLASTing.
        faa_files = msp.make_fasta_file(self.gb_files)
        # Run blastn locally.
        xml_results = msp.run_blastn(faa_files)
        # Delete fasta files used for BLASTing.
        msp.delete_files(faa_files)
        # Make a list of `BlastnAlignment` classes from the xml blastn results.
        alignments = msp.get_alignment_records(xml_results)
        # Delete xml documents.
        msp.delete_files(xml_results)
        # Make a list of `GenBankRecord` classes from the gb files.
        gb_records = msp.get_gb_records(self.gb_files)
        # Make figure.
        self.figure = msp.MakeFigure(
            alignments,
            gb_records,
            alignments_position=self.align_plot_var.get().lower(),
            identity_color=self.identity_color,
            color_map_range=self.colormap_range,
            # figure_name=self.gui_input.output_path,
            # figure_format=self.gui_input.figure_format,
            annotate_sequences=self.annotate_sequences,
            annotate_genes=self.annotate_genes,
            # sequence_name=info.sequence_name
            use_gui=True
        )
        # # Plot figure using msplotter
        # self.figure.make_figure()
        # self.figure.display_figure()
        # Plot figure using the Plot class
        self.fig = self.figure.make_figure()
        # Plot figure
        Plot(self.fig, self.figure)

    def on_closing(self):
        if self.fig is not None:
            self.figure.close_figure()
        self.quit()
        self.destroy()


def app_gui():
    app = App()
    app.mainloop()


if __name__ == '__main__':
    app_gui()
