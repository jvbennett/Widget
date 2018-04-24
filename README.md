# the Widget

Author: Jake Bennett

April 2018

The Widget is a tool for performing and testing the dE/dx reconstruction and calibration procedure
at BESIII and Belle II. Its modular design is intended to be useful for various tests as well as for
obtaining the actual calibration constants for a data sample.

To compile, make sure the $ROOTSYS environment variable points to your local ROOT installation. Then 
move into the Widget directory and run "make". Next move to the WidgetExe directory and run "make" again.

To perform a calibration, copy the run/template directory and change the parameters as needed. The
hadron.template.txt configuration file is to be used with the HadronCalibration executable. To see
the available parameters for the executable, use the -h flag (HadronCalibration -h).

A more complete documentation is available in the doc directory.
