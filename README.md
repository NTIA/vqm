# Video Quality Metric (VQM)

This repository provides an archival record of video quality metric (VQM) research performed by ITS between 1989 and 2011. 
This research focused on full reference (FR) and reduced reference (RR) methods, which compare a processed video sequence with an original, high quality video sequence. 
Thus, quality estimation for a processed video sequence requires access to the original, high quality video sequence. 

This software has not been maintained since 2013 and is provided for archival and research purposes only. 
ITS does not have the funding to perform troubleshooting. 
All souce code was developed for MATLAB® R2013b (8.2). Future versions of MATLAB may cause the code not to run.

These tools are available free of charge for any commercial or non-commercial use, in accordance with the terms described in the [LICENSE](https://github.com/NTIA/vqm/blob/master/LICENSE.md).

More information on these tools, including publications that fully document the VQM metrics and INLSA algorithm, can be found on the [ITS website](https://its.ntia.gov/research-topics/video-quality-research/software/).

| Name | Title |
| ---- | ----- |
| `bvqm` | Batch Video Quality Metric |
| `cvqm` | Command Line Video Quality Metric |
| `fdf` | Estimate HRC Fraction Dropped Frames using RR and NR VQMs |
| `psnr` | Peak Signal to Noise Ratio Search and PSNR Variable Frame Delay |
| `vfd` | Variable Frame Delays |
| `inlsa` | Iterated Nested Least Squares Algorithm, for combining multiple datasets |
| `screening` | Methods to screen subjects in a subjective test |
| `si` | Spatial Information Filter |


## Video Quality Metrics

The following VQM tools use FR and RR techniques. 
BVQM provides a graphical user interface, batch processing, and all VQMs released by ITS between 2000 and 2013. 
CVQM provides a command line interface for the same VQMs for a single video sequence (i.e., processed video sequence and associated original sequence). 
CVQM demonstrates how to split the metric calculation for an RR workflow.
FDF and VFD implement algorithms to calculate variable frame delay; these functions are also available within metrics supplied by BVQM and CVQM. 
PSNR implements the renown metric. 

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

The `inlsa` directory contains MATLAB code implementing INLSA. 
This algorithm that allows multiple subjective datasets to be fitted to a single subjective scale. 
INLSA computes the fit from a common set of objective metrics.

## Subject Screening
The `subject_screening` directory contains four methods for screeing subjects.

## Spatial Information (SI) Filter

The `si` directory contains MATLAB code implementing the spatial information (SI) filter. 
This filter detects long edges and estimates edge angle. 
The filter size can be adjusted to any odd number size (e.g., 13 by 13 pixels, 21 by 21 pixels).

## Contact

For questions, contact Margaret Pinson, (720) 601-7314, <a href="mailto:mpinson@ntia.gov">mpinson@ntia.gov</a>
