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

import msp.msplotter as msp
from msp.colormap_picker import ColormapPicker
from msp.plot import Plot
from msp.spinbox import FloatSpinbox

class App(customtkinter.CTk):
    """msplotter GUI."""
    def __init__(self):
        super().__init__()
        # -- Variables for BLASTing and plotting ------------------------------
        self.figure_plt = None                 # matplotlib object.
        self.gb_files: list = None
        self.figure_msp = None                 # msplotter object.
        self.identity_color: str = "Greys"
        self.colormap_range: tuple = (0, 0.75)
        self.y_limit: float = 0                # to adjust position of colormap
        self.scale_bar: bool = False
        self.annotate_sequences: bool = False
        # if self.annotate_sequences is True, self.sequence_name is used for
        # annotating the sequences.
        self.sequence_name: str = "accession"
        self.annotate_genes: bool = False
        self.annotate_genes_from: str = "gene_tag"

        # -- Set layout parameters --------------------------------------------
        self.geometry('750x600')
        self.title('MSPlotter')
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure((2,3), weight=1)
        self.grid_columnconfigure((0,1), weight=0)
        # Protocol to close app, including plot if exist.
        self.protocol('WM_DELETE_WINDOW', self.on_closing)
        # Variable to store the ColormapPicker class used in the
        # launch_colormap_picker function.
        self.colormap_app = None

        # -- Navigation frame -------------------------------------------------
        # Create navigation frame.
        self.navigation_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.navigation_frame.grid(row=0, column=0, columnspan=4, sticky='nsew')
        self.navigation_frame.grid_columnconfigure(0, weight=1)
        # Logo label
        self.logo_label = customtkinter.CTkLabel(
            self.navigation_frame, text='MSPloter',
            font=customtkinter.CTkFont(size=24, weight='bold')
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
            state='disabled',
        )
        self.clear_button.grid(row=0, column=2, padx=5, pady=20)
        # Plot button
        self.plot_button = customtkinter.CTkButton(
            self.navigation_frame, text='Plot',
            command=self.plot_figure,
            state='disabled',
            fg_color='#b300b3',
            hover_color='#800080',
        )
        self.plot_button.grid(row=0, column=3, padx=(5,10), pady=20)

        # -- Appearance frame -------------------------------------------------
        self.appearance_frame = customtkinter.CTkFrame(
            self, corner_radius=5,
        )
        self.appearance_frame.grid(
            row=1, column=0, columnspan=2, padx=(10, 5), pady=10, sticky='nsew'
        )
        self.appearance_frame.grid_rowconfigure(9, weight=1)
        self.appearance_frame.grid_columnconfigure((0,1), weight=0)
        # Appearance label
        self.appearance_label = customtkinter.CTkLabel(
            self.appearance_frame,
            text='Appearance',
            font=customtkinter.CTkFont(size=18, weight='bold'),
            width=120,
            height=50,
            corner_radius=5,
            # fg_color=self.appearance_fg_color,
        )
        self.appearance_label.grid(
            row=0, column=0, columnspan=2, pady=(10,10), sticky='nswe'
        )
        # Align plot label
        self.align_plot_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Align plot:'
        )
        self.align_plot_label.grid(row=1, column=0, pady=(10, 0))
        # Variable to store align_plot menu selection
        self.align_plot_var = customtkinter.StringVar(self, 'Left')
        # Align plot menu
        self.align_plot = customtkinter.CTkOptionMenu(
            self.appearance_frame,
            values=['Left', 'Center', 'Right'],
            variable=self.align_plot_var,
        )
        self.align_plot.grid(row=2, column=0, padx=20, pady=(0, 10))
        # Annotate sequences label
        self.annotate_seq_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Annotate sequences:'
        )
        self.annotate_seq_label.grid(row=3, column=0, pady=(10, 0))
        # Variable to store annotate_seq menu selection
        self.annotate_seq_var = customtkinter.StringVar(self, 'No')
        # Annotate sequences menu
        self.annotate_seq_menu = customtkinter.CTkOptionMenu(
            self.appearance_frame,
            values=['No', 'From acc number', 'From file name'],
            variable=self.annotate_seq_var,
            command=lambda _:self.update_annotate_seq()
        )
        self.annotate_seq_menu.grid(row=4, column=0, pady=(0,10))
        # Annotate genes label
        self.annotate_genes_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Annotate genes:'
        )
        self.annotate_genes_label.grid(row=5, column=0, pady=(10,0))
        # Variable to store annotate_genes menu selection
        self.annotate_genes_var = customtkinter.StringVar(self, 'No')
        # Annotate genes menu
        self.annotate_genes_menu = customtkinter.CTkOptionMenu(
            self.appearance_frame,
            values=['No', 'From gene tag', 'From product tag'],
            variable=self.annotate_genes_var,
            command=lambda _:self.update_annotate_genes()
        )
        self.annotate_genes_menu.grid(row=6, column=0, pady=(0,10))
        # Homology color label
        self.homology_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Homology color:',
        )
        self.homology_label.grid(row=1, column=1, pady=(10,0))
        # Colormap picker
        self.format_button = customtkinter.CTkButton(
            self.appearance_frame,
            text='Choose colormap',
            command=self.launch_colormap_picker
        )
        self.format_button.grid(row=2, column=1, padx=20, pady=(0, 10))
        # Scale bar label
        self.scale_bar_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Scale bar:'
        )
        self.scale_bar_label.grid(row=3, column=1, pady=(10,0))
        # Variable to store scale_bar menu selection
        self.scale_bar_var = customtkinter.StringVar(self, 'No')
        # Scale bar menu
        self.scale_bar_menu = customtkinter.CTkOptionMenu(
            self.appearance_frame,
            values=['No', 'Yes'],
            variable=self.scale_bar_var,
            command=lambda _:self.update_scale_bar()
        )
        self.scale_bar_menu.grid(row=4, column=1, pady=(0,10))
        # Color map position label
        self.cmap_position_label = customtkinter.CTkLabel(
            self.appearance_frame, text='Position colormap:'
        )
        self.cmap_position_label.grid(row=5, column=1, pady=(10, 0))
        # Color map position spinbox
        self.cmap_position_spinbox = FloatSpinbox(
            self.appearance_frame,
            width=140,
            height=28,
            command=self.update_y_limit
        )
        self.cmap_position_spinbox.grid(row=6, column=1, pady=(0, 10))
        # Create a validation function to check for float input in entry boxes
        self.validation = self.register(self.validate_input)
        # Change width check box variable
        self.width_check_var = customtkinter.StringVar(value='off')
        # Change width check box
        self.width_check = customtkinter.CTkCheckBox(
            self.appearance_frame,
            text='Change figure width:',
            variable=self.width_check_var,
            onvalue='on',
            offvalue='off',
            command=self.change_figure_width_event
        )
        self.width_check.grid(
            row=7, column=0, sticky='w', padx=(10,0), pady=(20, 10))
        # Change width entry box
        self.width_entry = customtkinter.CTkEntry(
            self.appearance_frame,
            validate='key',
            validatecommand=(self.validation, '%P'),
        )
        self.width_entry.grid(
            row=7, column=1, sticky='we', padx=(0,10), pady=(20, 10))
        self.width_entry.configure(state='disabled')
        # Change height check box variable
        self.height_check_var = customtkinter.StringVar(value='off')
        # Change height check box
        self.height_check = customtkinter.CTkCheckBox(
            self.appearance_frame,
            text='Change figure height:',
            variable=self.height_check_var,
            onvalue='on',
            offvalue='off',
            command=self.change_figure_height_event
        )
        self.height_check.grid(
            row=8, column=0, sticky='w', padx=(10,0), pady=(10))
        # Change height entry box
        self.height_entry = customtkinter.CTkEntry(
            self.appearance_frame,
            validate='key',
            validatecommand=(self.validation, '%P'),
        )
        self.height_entry.grid(
            row=8, column=1, sticky='we', padx=(0,10), pady=(10))
        self.height_entry.configure(state='disabled')
        # Reset button
        self.reset_button = customtkinter.CTkButton(
            self.appearance_frame, text='Reset', fg_color='transparent',
            text_color=('gray10', '#DCE4EE'), border_width=2,
            command=self.reset_appearance
        )
        self.reset_button.grid(row=9, column=0, columnspan=2, pady=(10, 10))

        # -- Display frame ---------------------------------------------------
        # Create display frame
        self.display_frame = customtkinter.CTkFrame(
            self, corner_radius=0,
            fg_color='transparent'
        )
        self.display_frame.grid(row=1, column=2, columnspan=2, sticky='nsew')
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

    def update_annotate_seq(self):
        if (annotate := self.annotate_seq_var.get()) == 'No':
            self.annotate_sequences = False
            return
        else:
            self.annotate_sequences = True
        if annotate == "From acc number":
            self.sequence_name = "accession"
        else:
            self.sequence_name = "fname"

    def update_annotate_genes(self):
        if (annotate := self.annotate_genes_var.get()) == 'No':
            self.annotate_genes = False
            return
        else:
            self.annotate_genes = True
        if annotate == "From gene tag":
            self.annotate_genes_from = "gene_tag"
        else:
            self.annotate_genes_from = "product_tag"

    def update_scale_bar(self):
        if self.scale_bar_var.get() == 'No':
            self.scale_bar = False
        else:
            self.scale_bar = True

    def update_y_limit(self):
        self.y_limit = self.cmap_position_spinbox.get()

    def change_figure_width_event(self):
        if self.width_check_var.get() == 'on':
            self.width_entry.configure(state='normal')
        else:
            self.width_entry.configure(state='disabled')

    def change_figure_height_event(self):
        if self.height_check_var.get() == 'on':
            self.height_entry.configure(state='normal')
        else:
            self.height_entry.configure(state='disabled')

    def validate_input(self, value):
        if value == "":
            return True
        try:
            float(value)
            return True
        except ValueError:
            return False

    def get_figure_size(self) -> tuple[float, float]:
        width = self.width_entry.get()
        height = self.height_entry.get()
        if not width:
            width = None
        else:
            width = float(width)
        if not height:
            height = None
        else:
            height = float(height)
        return (width, height)

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
        # Reset annotate genes
        self.annotate_genes_var.set('No')
        self.annotate_genes_menu.configure(variable=self.annotate_genes_var)
        # Reset scale bar
        self.scale_bar_var.set('No')
        self.scale_bar_menu.configure(variable=self.scale_bar_var)
        # Reset y_limit
        self.y_limit = 0
        self.cmap_position_spinbox.set(0.0)
        # Reset figure width and height functionality
        self.width_check.deselect()
        self.height_check.deselect()
        self.width_entry.configure(state='normal')
        self.width_entry.delete(0, 'end')
        self.width_entry.configure(state='disabled')
        self.height_entry.configure(state='normal')
        self.height_entry.delete(0, 'end')
        self.height_entry.configure(state='disabled')

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
        # Get figure size
        width, height = self.get_figure_size()
        # Make figure.
        self.figure_msp = msp.MakeFigure(
            alignments,
            gb_records,
            figure_width=width,
            figure_height=height,
            alignments_position=self.align_plot_var.get().lower(),
            identity_color=self.identity_color,
            color_map_range=self.colormap_range,
            annotate_sequences=self.annotate_sequences,
            sequence_name=self.sequence_name,
            annotate_genes=self.annotate_genes,
            annotate_genes_from=self.annotate_genes_from,
            scale_bar=self.scale_bar,
            y_limit=self.y_limit,
            use_gui=True
        )
        # Plot figure using the Plot class
        self.figure_plt = self.figure_msp.make_figure()
        # Plot figure
        Plot(self.figure_plt, self.figure_msp)

    def on_closing(self):
        if self.figure_plt is not None:
            self.figure_msp.close_figure()
        self.quit()
        self.destroy()


def app_gui():
    app = App()
    app.mainloop()


if __name__ == '__main__':
    app_gui()
