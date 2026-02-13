# HRPC_Becker_group
Data analysis workflows for HRPC

## bedRMod conversion

We convert the mapping produced by ThermoFisher BioPharma Finder to the bedRMod format to account for the RNA modifications that were found. Please see the following table for a description of how the bedRMod fields were used.

| bedRMod Field | Type | Description |
|--|--|--|
| `chrom` | string | Chromosome name |
| `chromStart` | integer | Start position of the modification |
| `chromEnd` | integer | End position of the modification |
| `name` | string | MODOMICS modification name |
| `score`| float | BioPharma Finder Confidence Score |
| `strand`| string | Feature strand (Default: ".") |
| `thickStart` | integer | Start position of the modification |
| `thickEnd` | integer | End position of the modification |
| `itemRgb` | integer | Display color (Default: 0) |
| `coverage` | integer | Number of fragments with valid MS/MS data (error <10ppm, ASR 1.0, Conf. Score 100), where modification has occurred at this position |
| `frequency`| float | `coverage` divided by number of total valid fragments detected |



> Written with [StackEdit](https://stackedit.io/).
