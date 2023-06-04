"""MSPlotter main function.

License
-------
This file is part of MSPloter
BSD 3-Clause License
Copyright (c) 2023, Ivan Munoz Gutierrez
"""
from msp.user_input import user_input, UserInput
from msp.msplotter import app_cli


def main():
    info = user_input()
    app_cli(UserInput(info))


if __name__ == "__main__":
    main()