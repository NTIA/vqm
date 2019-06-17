# VQM
## Video Quality Metrics
This repository contains the Video Quality Metrics developed by ITS.

For each metric, the source, a 32-bit compiled version, a 64-bit compiled version, and relevant documentation are provided. 

The available metrics are:

| Name | Title |
| ---- | ----- |
| `bvqm` | Batch Video Quality Metric |
| `cvqm` | Command Line Video Quality Metric |
| `fdf` | Estimate HRC Fraction Dropped Frames using RR and NR VQMs |
| `psnr` | Peak Signal to Noise Ratio Search and PSNR Variable Frame Delay |
| `vfd` | Variable Frame Delays |

## Usage

Each metric is stored in this repository such that for the given version of a given `<metric>`:
``` 
+ <metric>
|- <metric>_source_v123   # directory contains the source for that metric version.
|- <metric>_pc32_v123     # directory contains the 32-bit executable for that metric version. 
|- <metric>_pc64_v123     # directory contains the 64-bit executable for that metric version.
```

Readmes on how to use each metric are available under the name `<metric>_pc_readme.pdf`.