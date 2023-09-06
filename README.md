# Video Quality Metric (VQM) Archival Repository

This repository provides an archival record of video quality metric (VQM) research performed by ITS between 1989 and 2011. 
This research focused on full reference (FR) and reduced reference (RR) methods, which compare a processed video sequence with an original, high quality video sequence. 
Thus, quality estimation for a processed video sequence requires access to the original, high quality video sequence. 

This software has not been maintained since 2013 and is provided for archival and research purposes only. 
ITS does not have the funding to perform troubleshooting. 
All souce code was developed for MATLAB® R2013b (8.2). Future versions of MATLAB may cause the code not to run.

These tools are available free of charge for any commercial or non-commercial use, in accordance with the terms described in the [LICENSE](https://github.com/NTIA/vqm/blob/master/LICENSE.md).

More information on these tools, including publications that fully document the VQM metrics and INLSA algorithm, can be found on the [ITS website](https://its.ntia.gov/research-topics/video-quality-research/software/).

| Tool | Description |
| ---- | ----- |
| `bvqm` | Batch Video Quality Metric |
| `cvqm` | Command Line Video Quality Metric |
| `fdf` | Estimate HRC Fraction Dropped Frames using RR and NR VQMs |
| `psnr` | Peak Signal to Noise Ratio Search and PSNR Variable Frame Delay |
| `vfd` | Variable Frame Delays |
| `inlsa` | Iterated Nested Least Squares Algorithm, for combining multiple datasets |
| `screening` | Methods to screen subjects in a subjective test |
| `si` | Spatial Information Filter |


## Video Quality Metrics: BVQM, CVQM, FDF, PSNR, and VFD

The VQM software tools provide both standardized and non-standardized methods for measuring the video quality of digital video systems. 
Each VQM tool estimates how people perceive video quality. 
These VQM algorithms compare the processed video (output) with the original video (input). 
These models are suitable when the original video is good quality or better. 
Quality problems from the camera and original production are not considered. 

Although ITS released several metrics, the most impactful metrics for video quality assessment are as follows:

- NTIA General Model, referred to as 'VQM' in literature, as given in the ANSI T1.801.03-2003, ITU-T J.144 (03/04), ITU-R BT.1683 (06/04), [this report](https://its.ntia.gov/publications/details.aspx?pub=2423), and [this journal article](https://its.ntia.gov/publications/details.aspx?pub=2576). Finalized in 2001.
- NTIA Developers Model, a fast running variant of the General Model, described in [this report](https://its.ntia.gov/publications/details.aspx?pub=2423). Finalized in 2001.
- NTIA Video Quality Model for Variable Frame Delay (VQM_VFD), described in [this report](https://its.ntia.gov/publications/details.aspx?pub=2556). Finalized in 2011.

The Batch Video Quality Metric (BVQM) Software performs out-of-service, lab bench testing. 
BVQM can processing and analyses of multiple video scenes and multiple video systems at once. 
BVQM reads video sequences from files, and reports results to the screen. 
BVQM includes a variety of calibration options, quality models, and graphical presentation of results. 
For an overview of the BVQM tool, click [here](https://its.ntia.gov/publications/details.aspx?pub=2558).

The Command Line Video Quality Metric (CVQM) Software is a command line program for performing out-of-service, lab bench testing. 
The video quality calibration and metric algorithms within CVQM are basically identically to those offered by BVQM. 
The two differences are that (1) CVQM is called from a command line (e.g., Windows Accessory "Command Prompt"), and (2) CVQM cannot compare the results from multiple video clips to improve calibration accuracy.

CVQM runs on one pair of video files at a time, and writes results to files. 
CVQM demonstrates how to split the metric calculation for an RR workflow.
For an overview of the CVQM tool, click [here](cvqm.md).

FDF and VFD implement algorithms to calculate variable frame delay; these functions are also available within metrics supplied by BVQM and CVQM. 

File `psnr.m` calculates peak signal to noise ratio (PSNR) according to ITU-T Rec. J.340. 
It also computes PSNR with variable frame delays removed (PSNR-VFD). These metrics are also available in the BVQM and CVQM software.

For each software package, there are three download options: the MATLAB® source, a 32-bit compiled version, and a 64-bit compiled version. In each case, the relevant documentation is provided. 
Each metric is stored in this repository such that for the given version of a given `<name>` from the above table:
``` 
+ <name>
|- <name>_source_v123   # directory contains the source for that metric version.
|- <name>_pc32_v123     # directory contains the 32-bit executable for that metric version. 
|- <nmae>_pc64_v123     # directory contains the 64-bit executable for that metric version.
```

You may also download zipped versions of the `pc32` and `pc64` directories through [Releases](https://github.com/NTIA/vqm/releases).

Readmes on how to use each metric are available under the name `<name>_pc_readme.pdf`.

Note that this distribution does not include the MATLAB Runtime, which must first be installed before running any of the compiled versions of these metrics. The runtime can be installed using `MCRInstaller.exe`, located [here](https://www.mathworks.com/products/compiler/matlab-runtime.html). Be aware that the compiled versions of the metrics were compiled in MATLAB R2013b (8.2). Future versions of MATLAB may cause the code not to run. This issue applies to both the compiled code and the source code. 

## Iterated Nested Least Squares Algorithm (INLSA)

The inlsa directory contains an iterated nested least-squares algorithm (INLSA) for fitting multiple data sets

INLSA is an algorithm that allows multiple subjective datasets to be fitted to a single subjective scale. 
INLSA computes the fit from a common set of objective metrics. 
Files inlsa.m and pars_inlsa.m implement INLSA.
File inlsa_demo.m creates made-up data for three (3) experiments and plots that data before and after running INLSA. The user can actually see what INLSA does to the data. Also the user gets a concrete example of how to call INLSA.  The user can just replace the made-up data with real data and the use INLSA.

- Click [here](https://its.ntia.gov/publications/2428.aspx) for a reference. This document fully describes the algorithm.
- Click [here](https://its.ntia.gov/publications/2578.aspx) for a comparison of INLSA and subjective mapping. This paper demonstrates why data from different subjective tests must be fitted.
- Click [here](https://its.ntia.gov/publications/2494.aspx) for an example of how multiple subjective datasets can be combined using overlapping subjective datasets.

## Subject Screening
The screening directory contains four methods for screening subjects during a subjective test.

## Spatial Information (SI) Filter

The si directory contains MATLAB code implementing the spatial information (SI) filter. 
This filter detects long edges and estimates edge angle. 
The filter size can be adjusted to any odd number size (e.g., 13 by 13 pixels or 21 by 21 pixels).

## Contact

For questions, contact Margaret Pinson, (720) 601-7314, <a href="mailto:mpinson@ntia.gov">mpinson@ntia.gov</a>
