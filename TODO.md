# TODO for HopPyBar Python Program
### _version = '20260526'_

## IN PROGRESS
- Moving things towards functionalizing major operations and having multiple files for large chunks of code to make things easier to read/follow. So far, mostly COMPLETED output plots and output for Excel, but plan to do wherever possible to clean up code.
- General cleaning of variables, once things are more settled.
- Add more options for dispersion correction method. Currently have Bragov method and None. Might try to add some other popular options in the future, or at least make it easier to hook them in. (e.g. Shin et al., 2022)
- Add option to output both True and Engineering values simultaneously.
- Add option to direct import to GRANTA on Windows/Linux.
- Add direct analysis of PDV traces
- Add REL as an import option. Currently there is a separate program that reconstructs a LANL-style .raw file to use, which is fine but inconvenient.

## PLANNED
+ Add "confirm" and "redo timing" buttons to the 1/2/3-wave plot. Use while flag = False to stay in the loop until breakout.
+ Add ability to calibrate bars. Read from oscilloscope, perform calculations, write calibration files for future use. There's an initial stub for reading in data directly from oscilloscope, so this could be added as an import style in the future.

## POSSIBLE BUGS
+ FLAG might have variable timestep, which would screw with things if .raw file were used. Check whether this is the case.

## COMPLETE
* Add reconstruction of .raw file for formats that don't start with it (minibar, APS, FLAG). This should streamline GRANTA upload, file storage, and subsequent reanalysis
* Add option for suppressing outputs/plots
* Add option for suppressing GRANTA output
* Add option for engineering stress-strain output instead of true
* Fix instructions on transmitted selection screen to indicate no-selection will default to calculated values.
* Add "gain" adjustment to alignment curves to multiply signals for easier alignment; add better ticks to graph to aid alignment.
