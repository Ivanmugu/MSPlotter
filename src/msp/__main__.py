"""MSPlotter main function.

License
-------
This file is part of MSPloter
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
from msp.user_input import parse_command_line_input
from msp.msplotter import app_cli
from msp.gui import app_gui


def main():
    user_input = parse_command_line_input()
    if user_input.gui:
        app_gui()
    else:
        app_cli(user_input)


if __name__ == "__main__":
    main()