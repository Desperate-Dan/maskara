# maskara
A tool to create a coverage mask from BAM files to apply during consensus generation, in particular creates a tsv than can be used with bcftools consensus.
## Installation
```
git clone https://github.com/Desperate-Dan/maskara.git && cd maskara
pip install .
```

You can now install maskara directly from pip!

```
pip install maskara
```
## Usage
```
Maskara v1.1.7

usage: maskara.py [-h] [-d DEPTH] [-r REF_NAME] [-o OUTPUT_NAME] [-m FASTA_TO_MASK] [-i] [-q QUALITY] [--mmm] [--reads READS] [-v] input_file

Creates a coverage mask to apply to your lovely consensus fasta.

positional arguments:
  input_file            Path to the BAM file you want to create a mask for

options:
  -h, --help            show this help message and exit

Optional:
  -d DEPTH, --depth DEPTH
                        If coverage is below this it will be masked. Default = 20
  -r REF_NAME, --ref-name REF_NAME
                        Name of ref the bam files were aligned to. Default = "MN908947.3"
  -o OUTPUT_NAME, --output-name OUTPUT_NAME
                        Prefix for the output. Default = "depth_mask"
  -m FASTA_TO_MASK, --mask FASTA_TO_MASK
                        Mask a consensus sequence with your newly produced mask
  -i, --inverse         Return bed file of positions EQUAL OR ABOVE the chosen depth
  -q QUALITY, --quality QUALITY
                        Choose the minimum base quality for consideration in coverage counting. Default = 20
  --mmm                 Multi-Map Mode: Get a depth mask for all references in your bam file with at least X reads
  --reads READS         Set the read limit for Multi-Map Mode. If a reference has at least this many reads it will have a depth mask made (note this is not the same as depth)
  -v, --version         Return Maskara version
```
