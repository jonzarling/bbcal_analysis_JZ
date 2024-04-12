# jz_bbcal_analysis
Analysis code for data from the "Baby BCAL" (or BBCAL) used for benchmarking the Barrel Imaging Calorimeter (BIC) at the upcoming Electron Ion Collider.

[comment]: <> (This is a comment, it will not be included in MD file)
[comment]: <> (https://pandao.github.io/editor.md/en.html used for markdown viewer)

### Contents
- **hd_root_plugin/** contains plugin that converts from raw datastream (.evio extension files) to ROOT tree formatted output, using Hall D software (see setup below).
- Example ROOT file output from plugin included.
- Input evio files are ~ 10 GB in size and not included. Location on Jefferson Lab computing network available upon request.
- **offline_analysis/** contains code for further analysis of BBCAL output
    - **gain_calib/** contains code for determining calibration constants (i.e. "gains") to convert from ADC channel readout to convert to energy deposited in a cell.
    - **data_quality_monitoring/** contains code to perform quality checks on data within a run

*More code for resolution extraction, NPE (number of photoelectrons), and cosmic analysis to be added soon.*

## Processing Raw Data (hd_root plugin)
This plugin takes inputs from raw evio input files and stores relevant information for the BBCAL in a ROOT file. Some branch variables and code are Hall D specific, such as Pair Spectrometer variables (which begin with `PS_` in their string name)

###### Setup:
1. [Instructions on setting up GlueX (Hall D) software container](https://halldweb.jlab.org/wiki/index.php/HOWTO_use_the_GlueX_Singularity_Container)
2. `cd offline_analysis`
3. `scons install` to compile plugin
4. See `./run_plugin_testfile` for example of how to run plugin on command line.
5. See output file `bbcal_121216.root` for example of one output file.

## Energy Calibration
###### From raw ADC channel readouts (i.e. fADC pulse integral or amplitude) one can determine conversion factors to energy given a few caveats.

Reference: [see starting on slide 20](https://indico.bnl.gov/event/20800/contributions/83132/attachments/50849/86924/jz_babyBCAL_11.7.23.pdf)

###### Assumptions:
- EM shower energy is known independently of the BBCAL system on an event-by-event basis.
- Sufficient statistics and occupancy on all channels
- Clean event sample with little/no background

###### To run:
1. cd `offline_analysis/gain_calib`
2. `root -b -q 'MakeGainsAndTrees.cxx+("[filename.root]")'`
3. See example input file in `offline_analysis/gain_calib/input_plugin_files` and output file in  `offline_analysis/gain_calib/output_Ecalib_files` for example of output file.

## Data Quality Checking
Data quality checking code is written using PyROOT. If you run into errors while importing modules, run the following first (for CSH terminal, use export if using `bash` instead):
`setenv PYTHONPATH $ROOTSYS/lib:$PYTHONPATH`


###### To run:
1. cd `offline_analysis/data_quality_monitoring`
2. `BBCALQualityCheckSingleRun.py [filename.root] [output_tagname]`
3. Output files
4. See example file in `offline_analysis/gain_calib/output_root_files/bbcal_trees_Ecalib_121216.root`.


##Other Resources

#### [Baby BCAL Wiki Page](https://halldweb.jlab.org/wiki/index.php/Baby_BCAL_at_JLab)

#### [Run Log Spreadsheet](https://docs.google.com/spreadsheets/d/1Rz7lqmvchuP_fqpP4ckbNFOBRCjJ-Ccz-_Y8l-p5JJU/edit?usp=sharing) for BBCAL in Hall D
(spreadsheet tracks PS runs in March '23 and cosmic runs in fall)
