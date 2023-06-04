"""Class to make a colormap picker.

License
-------
This file is part of MSPloter
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""

import os
from pathlib import Path
from PIL import Image

import customtkinter

from msp.slider_widget import Slider


class ColormapPicker(customtkinter.CTkToplevel):

    cmaps = [
        'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrBr',
        'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu',
        'PuBuGn', 'BuGn', 'YlGn'
    ]

    def __init__(self, get_colormap_data=None):
        """
        Parameters
        ----------
        get_colormap_data : function
            Function that links ColormapPicker and App to transfer the selected
            colormap and colormap's range.
        """
        super().__init__()
        # Connect colormap_picker through the get_colormap_data function.
        self.get_colormap_data = get_colormap_data
        self.cmap = None
        self.cmap_range = None

        # set grid layout 1x2
        self.geometry('590x150')
        self.resizable(False, False)
        self.title(
            'Pick a colormap and a range of colors for the homology regions'
        )
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

        # =====================================================================
        # Create navigation frame
        # =====================================================================
        self.nav_frame = customtkinter.CTkFrame(self)
        self.nav_frame.grid(row=0, column=0, sticky='nsew')
        self.nav_frame.grid_rowconfigure(0, weight=1)
        self.nav_frame.grid_rowconfigure(1, weight=1)
        # ############################################## #
        # Create combobox and variable to save selection #
        # ############################################## #
        self.optionmenu_var = customtkinter.StringVar(self, value='Greys')
        self.combobox = customtkinter.CTkOptionMenu(
            master=self.nav_frame,
            values=self.cmaps,
            command=lambda _:self.update_colormap(),
            variable=self.optionmenu_var
        )
        self.combobox.grid(row=0, column=0, padx=20, pady=(20, 0))
        self.combobox.grid_rowconfigure(0, weight=1)

        # ###################### #
        # Create select button #
        # ###################### #
        self.select = customtkinter.CTkButton(
            self.nav_frame, text='Select', command=self.transfer_data
        )
        self.select.grid(row=1, column=0, padx=20, pady=(0, 20))
        self.select.grid_rowconfigure(0, weight=1)

        # =====================================================================
        # Colormap_frame
        # =====================================================================
        self.colormap_frame = customtkinter.CTkFrame(self, corner_radius=0,
            fg_color='transparent'
        )
        self.colormap_frame.grid(row=0, column=1, sticky='nsew')
        self.colormap_frame.grid_rowconfigure(0, weight=1)
        self.colormap_frame.grid_rowconfigure(1, weight=1)
        self.colormap_frame.grid_columnconfigure(0, weight=1)
        # ############# #
        # Make colormap #
        # ############# #
        # Initialized selected colormap and range_colormap.
        self.selected_colormap = None
        self.range_colormap = None
        # Get path to current directory
        self.current_dir = Path(os.path.dirname(os.path.realpath(__file__)))
        # Show selected colormap.
        self.update_colormap()

        # ################# #
        # show range slider #
        # ################# #
        self.range_slider = Slider(
            self.colormap_frame,
            width=390,
            height=60,
            min_val=0,
            max_val=100,
            init_lis=[0, 75],
            show_value=True,
        )
        self.range_slider.grid(
            row=1, column=0, padx=(0, 0), pady=(0, 30)
        )

    def update_colormap(self):
        """Show and update selected colormap."""
        print('updating selected colormap...')
        self.selected_colormap = self.optionmenu_var.get()
        colormap_with_extention = self.selected_colormap + '.png'
        colormap_path = (
            self.current_dir / 'images' / colormap_with_extention
        )
        colormap_image = customtkinter.CTkImage(
            Image.open(colormap_path), size=(364.032, 14.57)
        )
        self.show_colormap = customtkinter.CTkLabel(
            self.colormap_frame, text='', image=colormap_image
        )
        self.show_colormap.grid(row=0, column=0, padx=(0, 0), pady=(40, 0))
        print('optionmenu dropdown clicked:', self.selected_colormap)

    def transfer_data(self):
        """Transfer colormap data values to main window."""
        self.range_colormap = self.range_slider.getValues()
        self.get_colormap_data(self.selected_colormap, self.range_colormap)
        self.destroy()

if __name__ == '__main__':
    app = ColormapPicker()
    app.mainloop()
