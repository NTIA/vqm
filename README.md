# VQM
## Video Quality Metrics
This repository contains the Video Quality Metric (VQM) software developed by ITS between 2002 and 2010. 
The VQM algorithms estimate how people perceive video quality. These VQM algorithms compare the processed video (output) with the original video (input). These models are suitable when the original video is good quality or better. Quality problems from the camera and original production are not considered. 
The VQM tools are available free of charge for any commercial or non-commercial use, in accordance with the terms described in the [LICENSE](https://github.com/NTIA/vqm/blob/master/LICENSE.md).

For each software package, there are three download options: the MATLAB® source, a 32-bit compiled version, and a 64-bit compiled version. In each case, the relevant documentation is provided. 

The available software packages are:

| Name | Title |
| ---- | ----- |
| `bvqm` | Batch Video Quality Metric |
| `cvqm` | Command Line Video Quality Metric |
| `fdf` | Estimate HRC Fraction Dropped Frames using RR and NR VQMs |
| `psnr` | Peak Signal to Noise Ratio Search and PSNR Variable Frame Delay |
| `vfd` | Variable Frame Delays |

These VQM tools use full reference (FR) and reduced reference (RR) techniques. Thus, quality estimation for a processed video sequence requires access to the original, high quality video sequence. BVQM provides a graphical user interface, batch processing, and all VQMs developed by ITS from 2000 through 2013. CVQM provides a command line interface for the same VQMs for a single video sequence (i.e., processed video sequence and associated original sequence). FDF and VFD implement algorithms to calculate variable frame delay; these functions are also available within metrics supplied by BVQM and CVQM. PSNR implements the renown metric.  

## Structure

Each metric is stored in this repository such that for the given version of a given `<name>` from the above table:
``` 
+ <name>
|- <name>_source_v123   # directory contains the source for that metric version.
|- <name>_pc32_v123     # directory contains the 32-bit executable for that metric version. 
|- <name>_pc64_v123     # directory contains the 64-bit executable for that metric version.
```

Readmes on how to use each metric are available under the name `<name>_pc_readme.pdf`.

Note that this distribution does not include the MATLAB Runtime, which must first be installed before running any of the compiled versions of these metrics. The runtime can be installed using `MCRInstaller.exe`, located [here](https://www.mathworks.com/products/compiler/matlab-runtime.html). Be aware that the compiled versions of the metrics were compiled in MATLAB R2013b (8.2). Future versions of MATLAB may cause the code not to run. This issue applies to both the compiled code and the source code. 

## Usage

### Batch Video Quality Metric (BVQM) Software 
BVQM is a Windows program for performing out-of-service / lab bench testing. 
BVQM can perform processing and analyses of multiple video scenes and multiple video systems at once. 
BVQM reads video sequences from files, and reports results to the screen. 
BVQM includes a variety of calibration options, quality models, and graphical presentation of results. 

* BVQM_pc32_v20 has BVQM version 2.0, compiled for 32-bit Windows operating systems.
* BVQM_pc64_v20 has BVQM version 2.0, compiled for 64-bit Windows operating systems.
* BVQM_src_v20 has MATLAB source code for BVQM 2.0.
* [Documentation](https://its.ntia.gov/publications/details.aspx?pub=2558)

BVQM contains unnecessary functionality and residual code remaining from the metric development process. 

### Command Line Video Quality Metric (CVQM) Software
CVQM is a Windows command line program for performing out-of-service / lab bench testing. The video quality calibration and metric algorithms within CVQM are basically identical to those offered by BVQM. 
The two differences are that (1) CVQM is called from a command line (e.g., Windows Accessory "Command Prompt"), and (2) CVQM cannot compare the results from multiple video clips to improve calibration accuracy.  
CVQM runs on one pair of video files at a time, and writes results to files. 

* CVQM_pc32_v30 has CVQM version 3.0 compiled for the 32-bit Windows operating system.
* CVQM_pc64_v30 has CVQM version 3.0 compiled for the 64-bit Windows operating system.
* CVQM_src_v30 has MATLAB source code for CVQM version 3.0.\
* [Documentation](CVQM.md)

CVQM is a clean implementation of the VQM metrics that are available in BVQM. 
The software functionality is split to support an in-service RR system (i.e., live calculations on the upstream and downstream video feeds, communicating via a secondary communication channel). 

### PSNR
This command-line MATALAB software calculates peak signal to noise ratio (PSNR) according to ITU-T Rec. J.340. Also computes PSNR with variable frame delays removed (PSNR-VFD). These metrics are also available in the BVQM and CVQM software.
The algorithm is described [here](https://its.ntia.gov/publications/details.aspx?pub=2500).

### VFD
This command-line MATLAB software calculates variable frame delay (VFD) between an original video sequence and a processed video sequence. This algorithm and code are used by BVQM, CVQM, and PSNR software packages.
The algorithm is described [here](https://its.ntia.gov/publications/2500.aspx).

### FDF
This command-line MATLAB software calculates fraction dropped frames (FDF) from VFD statistics. This algorithm and code are used by BVQM, CVQM, and PSNR software packages.
The algorithm is described [here](https://its.ntia.gov/publications/2493.aspx).

### IVQM - Obsolete - Not Available
The In-Service Video Quality Metric (IVQM) was a "proof of concept" software implementation.
This software demonstrates that the reduced reference (RR) algorithms provided within BVQM and CVQM can actually run in-service. 

The IVQM software is no longer available. 
The video capture code is not compatible with any current video capture device.
The remainder of the code is provided by the CVQM software package. 

### VQM - Obsolete - Not Available
VQM_pc (compiled for Windows) and VQM_lx (compiled for Linux) have been replaced by BVQM. 
This software was developed under a cooperative research and development agreement (CRADA) by ITS and Intel. 
This was the first version of the ITS video quality models that was suitable for use by people outside of our laboratory. 
The source code cannot be redistributed.
The code cannot be compiled or run on modern operating systems.

## Contact

For questions, contact Margaret Pinson, <a href="mailto:mpinson@ntia.gov">mpinson@ntia.gov</a>
