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

import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from pathlib import Path
import msplotter as msp


class App:
    def __init__(self, root):
        self.main = self.main_window(root)
        self.display_window = self.make_display()
        self.gb_files = None
        self.figure = None
        self.select_button = self.make_select_button()
        self.clear_button = self.make_clear_button()
        self.plot_button = self.make_plot_button()
        self.save_button = self.make_save_button()

    def main_window(self, root):
        root.title('MSPlotter')
        root.columnconfigure(0, weight=1)
        root.geometry('600x400')
        main = ttk.Frame(root, padding=(30, 30))
        main.grid()
        return main

    def make_display(self):
        display_window = tk.Text(self.main, width=40, height=20, state='disabled')
        display_window.grid(row=0, column=0, columnspan=3, sticky='EW')
        return display_window

    def make_select_button(self):
        # Select button
        select_button = tk.Button(
            self.main, height=1, width=10, text='Select',
            command=self.get_files_path
        )
        select_button.grid(row=1, column=0)
        return select_button

    def make_clear_button(self):
        clear_button = tk.Button(
            self.main, height=1, width=10, text='Clear',
            command=self.clear_input,
            state='disabled'
        )
        clear_button.grid(row=1, column=1)
        return clear_button

    def make_plot_button(self):
        plot_button = tk.Button(
            self.main, height=1, width=10, text='Plot',
            command=self.plot_figure,
            state='disabled'
        )
        plot_button.grid(row=1, column=2)
        return plot_button

    def make_save_button(self):
        save_button = tk.Button(
            self.main, height=1, width=10, text='Save',
            command=self.save_figure,
            state='disabled'
        )
        save_button.grid(row=2, column=1)
        return save_button

    def add_files(self, input_files):
        if self.gb_files is None:
            self.gb_files = [Path(element) for element in input_files]
        else:
            for element in input_files:
                self.gb_files.append(Path(element))

    def get_files_path(self):
        self.add_files(filedialog.askopenfilenames(
            initialdir=".",
            title="Select a GenBank file",
            filetypes=(('GenBank files', '*.gb'), ('All files', '*.*'))
            ))
        self.display_window.config(state='normal')
        self.display_window.delete('1.0', 'end')
        # Print files to be analyzed
        self.display_window.insert(
            'end', 'Files are going to be BLASTed in the next order:\n')
        for i, file_path in enumerate(self.gb_files):
            self.display_window.insert(
                'end', f'{i+1} --> {file_path.name}\n')
        self.display_window.config(state='disabled')
        self.clear_button.config(state='normal')
        self.plot_button.config(state='normal')

    def clear_input(self):
        self.gb_files = None
        self.display_window.config(state='normal')
        self.display_window.delete('1.0', 'end')
        self.display_window.config(state='disabled')
        self.plot_button.config(state='disabled')
        self.save_button.config(state='disabled')

    def plot_figure(self):
        self.save_button.config(state='normal')
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
            # alignments_position=info.alignments_position,
            # identity_color=info.identity_color,
            figure_name='test_gui.pdf',
            figure_format='pdf',
            # annotate_sequences=info.annotate_sequences,
            # sequence_name=info.sequence_name
            use_gui=True
        )
        self.figure.make_figure()
        self.figure.display_figure()

    def save_figure(self):
        self.figure.save_plot()

def run_gui():
    root = tk.Tk()
    app = App(root)
    app.main.mainloop()

if __name__ == '__main__':
    run_gui()