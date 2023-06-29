# HopPyBar GRANTA Functionality

The HopPyBar program uses the files in this folder to add metadata to the analyzed file. This is used to streamline GRANTA import, specifically, though it could probably be adapted to fit import schema for any arbitrary database.

# Usage
Since GRANTA relies on a number of discrete values for some attributes, we've separated them here into individual files which match the options (exactly) in the GRANTA database that will eventually receive the test data. Unless otherwise noted, options are added one-per-line in the corresponding files. These values are then read into the granta_input_ui.py file to populate dropdown boxes for user selection.

# Future Work
Planned to add the option of direct import to GRANTA via the MI Python STK toolkit, which works for Windows and Linux, but not for macOS. 