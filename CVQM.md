# CVQM Overview
CVQM performs automated processing on a pair of video files. One contains an original video sequence (e.g., straight from the camera) and the other contains a processed video sequence (e.g., after coding and transmission and decoding). 

All video sequences must be in uncompressed AVI files, in either the UYVY or RGB color space. Original and processed video sequences must contain the same image size and frames per second. CVQM compares the original video sequence to the processed video sequence (i.e., the sequence that has been processed by the video system under test). 

Every original / processed video sequence pair is run through two main steps. First, the requested calibration is run, and results are saved in a file named after the processed video file, with '_calibration.txt' appended. Second, the requested model is run, and results are saved in a file named after the processed video file, with '_model.txt' appended.

Any errors that occur are recorded in a file named after the processed video file, with '_errors.txt' appended. If that file is absent, no errors occurred. Each line on this file will start with a number, indicating error type, and then a string describing the error. 
Type 1 errors are fatal data input/output issues (e.g., operation cannot continue due to invalid input argument, file read error, or file write error). Type 2 errors are calibration warnings (e.g., still sequence, temporal registration failure, spatial shift beyond search limits). Type 3 errors are non-fatal warnings.

## VQM Model Options
The following VQM models options are available:

* 'none'- No model will be calculated - run calibration only (see below).
* 'general' - NTIA General Model, as given in the ANSI T1.801.03-2003, ITU-T J.144 (03/04), ITU-R BT.1683 (06/04), and [TR-02-392](https://its.ntia.gov/publications/details?pub=2423). 
* 'developers' - Developer's model, a fast running variant of the General Model, described in [TR-02-392](https://its.ntia.gov/publications/details?pub=2423)
* 'lowbw' - Low Bandwidth Model, described [here](https://its.ntia.gov/publications/details?pub=2575)
* 'fastlowbw' - Fast Low Bandwidth Model, a fast running variant of the low bandwidth model 

The CVQM software provides a clean implementation of these VQM models. 
It is intended to be ported into other programming languages or to be called from other software. 

## Calibration Options
The following calibration options are available. 
Technical details are available in [TR-02-392](https://its.ntia.gov/publications/details?pub=2423) and [TR-06-433b](https://its.ntia.gov/publications/details?pub=2465).

### 'none'
No calibration will be performed. Assume that the first frame in original and processed video files align temporally. Run model with default calibration values which assumes that the processed video file is perfectly calibrated.

### 'manual' 
Read the calibration file created on a previous run. The values on the beginning of each line may be manually modified.

### 'rrcal' 
Perform reduced reference calibration as given in ntia_tr_06_433a.pdf (except assume no spatial scaling). These algorithms use random processes, which may yield slightly different results from one run to another. 
It is highly recommended that rrcal results be median filtered across 7 or more different video sequences that have been sent through the same video system (see Calibration Note below).

### 'rrcalscale' 
Perform reduced reference calibration as given in ntia_tr_06_433a.pdf, including estimating spatial scaling (e.g., stretch). These algorithms use random processes, which may yield slightly different results from one run to another. 
It is highly recommended that rrcalscale results be median filtered across 7 or more different video sequences that have been sent through the same video system (see Calibration Note below).

### 'rrcal2' 
Improved version of 'rrcal' - version 2, as specified in ITU-T Recommendation J.244 and ntia_tr_08_433b.pdf (see also 'rrcal2scale'). Includes estimate of Cb and Cr gain and offset; and slightly improved luminance gain & offset algorithm.

### 'rrcal2scale' 
Improved version of 'rrcalscale' - version 2, as specified in ITU-T Recommendation J.244 and ntia_tr_08_433b.pdf (see also 'rrcal2'). Includes estimate of Cb and Cr gain and offset; and slightly improved luminance gain & offset algorithm.

### 'frcal' 
Perform full reference bandwidth calibration as given in ANSI T1.801.03-2003, ITU-T J.144 (03/04), and ITU-R BT.1683 (06/04). 
Preferably, results should be median filtered across several different video sequences that have been sent through the same video system (see Calibration Note below).

### 'frtime' 
Performs full reference temporal registration and valid region estimation as given in ANSI T1.801.03-2003, ITU-T J.144 (03/04), and ITU-R BT.1683 (06/04). No other calibration will be performed. 
Run model with default calibration values. Median filtering is not necessary. Suitable for video systems that are known to never shift, scale, or change the luminance levels.

### 'rrtime' 
Performs reduced reference temporal registration and valid region estimation only. No other calibration will be performed. Run model with default calibration values. Median filtering is not necessary. 
Suitable for video systems that are known to never shift, scale, or change the luminance levels.

### 'frtimemanual' 
Read the calibration file created on a previous run. The values on the beginning of each line may be manually modified. 
Then, ignore delay specification and perform full reference temporal registration as given in ANSI T1.801.03-2003, ITU-T J.144 (03/04), and ITU-R BT.1683 (06/04).

### 'rrtimemanual' 
Read the calibration file created on a previous run. The values on the beginning of each line may be manually modified. Then, ignore delay specification and perform reduced reference temporal registration from 'rrcal' and 'rrcal2'.

## Calibration Search Ranges
The following calibration search ranges are presumed for 'rrcal' and 'rrcalscale:

* Maximum temporal registration uncertainty (1 second)
* Maximum spatial shift search (+/- 4 pixels for QCIF, +/- 8 pixels for CIF, and +/- 20 pixels for NTSC & PAL & HDTV).
* Maximum spatial scaling (6% for QCIF & CIF, 10% for NTSC & PAL & HDTV)

## Calibration Note:
Increased calibration accuracy may be obtained by median filtering results over multiple video sequences. In this case, several video sequences must be run through the exact same system, with all system parameters held constant. 
Run with 'none' specified for the model, and examine all calibration results produced by these runs. The following values can be median filtered for increased accuracy: Horizontal Shift, Vertical Shift, Luminance Gain, Luminance Offset, Horizontal Scale, and Vertical Scale.

Median Filtering Example: For each of the above calibration quantities, form a list of the values obtained from each video sequence (e.g., Horizontal Shift = [10 10 11 10 -13 10 9 10 17] - notice the erroneous horizontal shift values of "-13" and "17").
Calculate the median value by sorting the values and returning the 50% percentile value (e.g., 10 in the previous example). Then, write that value to each processed video sequence's calibration file (e.g., "10 Horizontal Shift") on the correct line. 
If CVQM is then run with 'manual' specified for calibration, the median filtered calibration values will be used.
