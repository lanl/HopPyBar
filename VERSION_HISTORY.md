# HopPyBar Version History
| Version | Notes |
| ------ | ------ |
| HopPyBar.py | analyzes data, little clunky |
| HopPyBar_align | first attempt at allowing to align pulses after initial selection |
| HopPyBar_align2 | switch to a different method of selection & alignment. Instead of selecting windows, we select the start of each pulse. Still needs to put everything together and resample the window.|
| HopPyBar_align_20180919 | I think I hammered out the alignment stuff
| HopPyBar_align_20180920 | Tried to clean up code, including the equations for strain, and then commensurate changes to the strainrate and stress equations.|
| HopPyBar_align_20181030 | <ul><li>Tried cleaning up by substituting new equations from Addessio, Blumenthal, & Gray;</li><li>switching completely to metric units;</li><li>split bar wavespeeds and moduli to account for potentially different bars.</li><li>Added Excel export of all data, and changed true strain rate formula to do instantaneous differences of true strain and time.</li></ul>|
| HopPyBar_20181114 | Further tweaks, and set window to maximum size (can be scaled down to arbitrary number manually, as usual).|
| HopPyBar_LANL_HEbar_20181128 | Branch to read in data from the HE bar at LANL. Hopefully we can pull some of the constants from the headers, too, to save some hassle.|
| HopPyBar_20181128 | Added option of having different incident & transmitted bars. |
| HopPyBar_20190508 | <ul><li>Added output of timeshifted pulses to ensure it was doing what it was supposed to;</li><li> BUGFIX: Fixed bug in in2mm def;</li><li>Added output for variance, so see the total error from comparing incident with trans and refl pulses.</li></ul>|
| HopPyBar_20190529 | <ul><li>First attempt at adding aspects for dispersion correction.</li><li>Added export of input_parameters for easier reference.</li></ul>|
| HopPyBar_20190625 | Smallish tweaks. If no background is selected, set background to zero instead of crashing. |
| HopPyBar_20190913 | Attempting to add a dead-reckoning time of arrival estimate for the pulses. For now, we will just leave it in the background as extra data rather than break anything that's working now. We could calculate sound speed based on final timing, too. Maybe add later. |
| HopPyBar_20200123 | <ul><li>Adding better hooks for Sierra Peaks Kolsky program (including exporting an input parameters in a format that looks similar to the first two columns of the Sierra Peaks output file).</li><li>Added dead-reckoning for reflected pulse. </li></ul> |
| HopPyBar_20200124 | Same as previous version with many unused snippets removed for clarity. May still need implementation of features from previous versions. |
| HopPyBar_20200309 | Started moving things out to modules (e.g. unit conversions)|
| HopPyBar_20200827	| <ul><li>Incorporated edits by Ginny Euser to move to pandas dataframes in many more places. Now uses DataFrames pretty extensively, which looks much nicer and eliminates some of the weird things we were doing to get the data lengths to match up.</li><li>Some changes to be compatible with interactive python (e.g. block=True for plotting);</li><li>Changed way that strain plateau of incident pulse is defined to get striker velocity. Still to-do: implement dispersion correction.</li></ul>|
| HopPyBar_20200901	|	<ul><li>Added dispersion correction via Bragov method</li><li>Expanded import/export options</li><li> Moved from tkinter to pyQt5 for GUI elements</li></ul> |
| HopPyBar_20201006 | Moved to gitlab @ LANL. Cleaned up some things. Plan to continue to develop dispersion correction, increase export options (e.g. Granta format) |
| HopPyBar_20210506 | <ul><li>Implementing dispersion correction into final calculations</li><li> Updating output to include GRANTA format</li><li> Save dispersion corrections to the HopPyBar run directory so we don't have to redo these so often</li><li> Added variance output on alignment of pulses. Added variance display to pulse alignment widget</li><li> Patched the input file for the mini-Kolsky so it works again (broke due to unit mismatch with dispersion correction).</li></ul> |
| HopPyBar_20210608	| Transmitted shift now uses either user-selected window or the tt_dr_shift calculated value if none is selected. |
| HopPyBar_20210710 | <ul><li>Added some snippets for an attempt at a non-folder-based test name (so far unsuccessful)</li><li> Added output for bar forces for friction test analysis</li><li> Added Granta clean file output. Changed Granta output to root folder so we don't have to dig it out of the Analysis folder. |
| HopPyBar_20210721	|	Added ability to save preferences for input options and to suppress output, plus small bugfixes. |
| HopPyBar_20210728 | <ul><li> Added .raw reconstruction for formats that don't start with .raw as input</li><li> Read in timestep for each of these instead of using fixed, to account for variable nature of scopes. </li><li>Bonus: more robust, even if standard gets changed for some reason. </li><li>Possible BUG: FLAG might have variable timestep, which would screw with things if .raw file were used.</li></ul> |
| HopPyBar_20210922 | Added multiplier for the transmitted wave alignment to make it easier to align things |
| HopPyBar_20211117	| <ul><li>Effort to remove globals from GUIs and return explicit values. Seems to work. Might not be totally robust.</li><li>BUGFIX: changed dispersion correction from bar_i_id to calibration_id_i (and same for trans) because the bars are probably more generically named than the calibrations. </li><li>Confirmed the code works in Python3.10.</li></ul> |
| HopPyBar_20211220 | <ul><li>Fixed some more bugs with redux of GUIs (to remove globals). Should be okay now. </li><li> Possible BUG: FLAG might have variable timestep, which would screw with things if .raw file were used. </li><li> Started moving to functions for major components (like plotting, 1-wave and 2-wave analysis) to clean the code further.</li></ul>|
| HopPyBar_20221005 | Minor bug fixes with pandas (drop(0) & deprecation warning for append) |
| HopPyBar_20221213 | Minor fixes to clear some deprecation warnings in index handling and ExcelWriter | 
| HpoPyBar_20230628 | <ul><li>First open-source release!!! Added banners and updated licensing info.</li><li> Added dispersion correction method select (including "None"). </li><li> Organized GRANTA options in separate folder, no more hard-coding dropbox values.</li></ul>