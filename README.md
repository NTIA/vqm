# VQM
## Video Quality Metrics
This repository contains the Video Quality Metric (VQM) software developed by ITS. A full description of each of the files provided here can be found <a href="https://www.its.bldrdoc.gov/resources/video-quality-research/guides-and-tutorials/description-of-vqm-tools.aspx" target="_blank"> here on the ITS Web Site</a>. The VQM tools are available free of charge for any commercial or non-commercial use, in accordance with the terms described in the [LICENSE](https://github.com/NTIA/vqm/blob/master/LICENSE.md).

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

## Usage

Each metric is stored in this repository such that for the given version of a given `<name>` from the above table:
``` 
+ <name>
|- <name>_source_v123   # directory contains the source for that metric version.
|- <name>_pc32_v123     # directory contains the 32-bit executable for that metric version. 
|- <nmae>_pc64_v123     # directory contains the 64-bit executable for that metric version.
```

You may also download zipped versions of the `pc32` and `pc64` directories through [Releases](https://github.com/NTIA/vqm/releases).

## Contact

For questions, contact Margaret Pinson, (303) 497-3579, <a href="mailto:mpinson@ntia.gov">mpinson@ntia.gov</a>

Readmes on how to use each metric are available under the name `<name>_pc_readme.pdf`.

Note that this distribution does not include the MATLAB Runtime, which must first be installed before running any of the compiled versions of these metrics. The runtime can be installed using `MCRInstaller.exe`, located [here](https://www.mathworks.com/products/compiler/matlab-runtime.html). Be aware that the compiled versions of the metrics were compiled in MATLAB R2013b (8.2). Future versions of MATLAB may cause the code not to run. This issue applies to both the compiled code and the source code. 
