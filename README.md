# Mesa-Stella-SESNe
A set of Python scripts for creating and running MESA+Stella model grids for stripped-envelope supernovae.

## How it works
The main component is within ```MesaStellaCore.py```.  This script reads the configuration file ```SetupConfig.cfg``` and an input simlist (```InputFiles/simlist.csv``` by default), then creates a set of MESA and Stella simulation grids with the parameters specified within the simlist.  It then runs the MESA simulations with a user-defined thread count.  Once the MESA sims finish, it runs the Stella component of the simulations in parallel, each on its own thread.  Stella has limited parallelization, so this really speeds up the process.  Output data is held within ```ModelGrids```, though CSVs of some output data are exported to ```DataExports```.

## Getting Started

### Installation

Make sure you have Docker installed on your system.  Read the installation guidelines for your system here: https://docs.docker.com/engine/install/

Once you have Docker installed, clone this repo to a directory of your choice with ```git clone https://github.com/amartas/Mesa-Stella-SESNe.git```.  Then, create a docker container with ```docker build -t <name>:latest```.  Docker will run through the Dockerfile and build the container.  You can enter the container with ```docker run -it <name>:latest```, and get to your simulations!  The primary Mesa directory is held within ```~/MESA``` alongside my scripts for grid-building.

### How to create model grids

#### Photometry
As an input, you need to give your photometry points in either g, r, or i.  (I may update this to include more at a later date).  This input file is a .csv that *must* meet the following specifications:

You must have:
- A column ```jd``` that contains the Julian date of the photometry point, in float form
- A column ```mag``` that contains the *apparent* magnitude of the photometry point, in float form
- A column ```magerr``` that contains the *apparent* magnitude error of the photometry point, in float form
- A column ```filter``` that contains the filter used for that photometry point as a string (no quotes).  This can be any filter, but the core script will only read ones labeled "g", "r", or "i", for Gaia *gri* bandpasses.

There can be no NaN or empty values in the CSV.  I suggest using Pandas to create this file.

#### Simlist
As another input, you need the actual parameters of the models you wish to create.  Currently, the core script supports:

- Progenitor ZAMS mass ($M_\odot$)
- Explosion energy (in $10^{50}$ ergs)
- CSM velocity (km/s)
- CSM mass loss rate ($M_\odot$/year)
- Duration of the mass loss (years)
- The Dutch wind scaling factor $\eta$ (This is zero by default - if you want to strip envelopes, increase this value)
- Whether the simulation has CSM (Boolean: 1 True, 0 False)

The core script reads another CSV with this information.  It must meet the following specifications:

You must have columns:
- ```mass```: float, ZAMS mass ($M_\odot$)
- ```energy```: float, explosion energy ($10^{50}$ ergs)
- ```windscalar```: float, Dutch wind scaling factor $\eta$
- ```metallicity```: float, metal mass fraction
- ```hefrac```: float, helium mass fraction
- ```ni56```: float, Ni56 mass ($M_\odot$)
- ```csmvelo```: float, CSM velocity (km/s)
- ```csmrate```: float, CSM mass loss rate ($M_\odot$/yr)
- ```csmtime```: float, CSM mass loss duration (yr)
- ```progoptimize```: logical 1 or 0, enables progenitor optimization
- ```csmoptimize```: logical 1 or 0, enables CSM optimization
- ```gridtag```: string, identifier for different sets of models; exported data will be saved under this name in ```DataExports```

There can be no NaN or empty values in the CSV.

##### ProgOptimize

For cases where the progenitor parameters (mass, $\eta$, metallicity, helium mass fraction) are the same, one can skip the pre-core-collapse modeling if a model has already been built for those parameters.  Each pre-core-collapse model you build is saved to ```ProgOptimize``` by default.  If this is enabled (set to 1), ```MesaStellaCore.py``` will assume that there's already a model built with these parameters, and will use the saved model rather than generating a new, identical one.  This hugely saves on time for grids involving variable $\eta$.

##### CSMOptimize

If you fix the progenitor and explosion properties & only vary CSM parameters, you can completely eliminate modeling prior to MESA's construction of the CSM in ```shock_part_5```.  If you enable this, you need to place the prior step's model (```shock_part_4.mod```) from a given post-core-collapse model (```PostCC```) into ```InputFiles```.  You can find this ```.mod``` file from a given model by looking at ```mesa-24.08.1/ModelGrids/<MODEL_DIR>/PostCC/shock_part_4.mod```.  **BEWARE: This will make all simulations with CSMOptimize enabled use this model - do not run multiple progenitors in the same grid with CSMOptimize enabled.**

#### Config

Open ```SetupConfig.cfg``` with the text editor of your choice and fill in your user parameters and object parameters.  Don't skate over anything - make sure you've read through all the config options and made sure they're what you want - you don't want to have spent all that computation time building sims that don't even apply to your research.

### Running the models

You should have everything set up now!  All you need to do now is run ```MesaStellaCore.py``` in the Python environment from earlier, and it'll start chugging along!  Some data is exported in the form of a CSV to ```DataExports```, but all the output data is stored in subdirectories within ```mesa-24.08.1/ModelGrids```.  Go read the MESA documentation to learn to read it!  Make sure to move these sims somewhere else *outside* the parent directory, as ```MesaStellaCore.py``` will *not* overwrite these sims if you are rerunning with identical input parameters, throwing an error.

