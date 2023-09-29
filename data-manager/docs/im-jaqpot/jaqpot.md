# The Data Manager jaqpot Job documentation

This describes how to run the `jaqpot` job in the `im-jaqpot` collection.

## What the job does

This job predicts molecular properties for the input structures with models from Jaqpot package. Available endpoints are:
- Aqueous solubility model
- hERG model
- AMES model
- CYP2C9 inhibition model
- CYP3A4 inhibition model
- CYP2C9 substrate model
- CYP2D6 inhibition model
- Lipophilicity model
- PPBR model
- HIA model
- CYP2D6 substrate model
- Bioavailability model
- Clearance Microsome model
- LD50 model
- CYP3A4 Substrate CarbonMangels  model
- CaCO2 Wang model
- DILI model
- VDss Lombardo model
- Clearance Hepatocyte model
- Half Life Obach model
- Blood Brain Barrier model
- PGP model


## How to run the job

### Inputs

* **Molecules to predict**: input molecules in sdf or SMILES format.
* **Models**: models to use for prediction, several can be given simultaneously.

### Options
* **Output file**: name of the returned file. Default: `result.sdf`.
* **Delimiter**: delimiter to use in output file. Default: tab. Ignored in sdf output
* **ReadHeader**: Read header from the input file. Default: False. Ignored in sdf output
* **WriteHeader**: Write header line to output file. Default: True. Ignored in sdf output

## Related topics

* [Virtual screening](https://github.com/InformaticsMatters/virtual-screening)
