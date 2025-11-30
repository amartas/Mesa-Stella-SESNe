# Mesa-Stella-SESNe
A set of Python scripts for creating and running MESA+Stella model grids for stripped-envelope supernovae.

## How it works
The main logic behind the gridding is ```MesaStellaCore.py```.  This script reads the configuration file ```SetupConfig.cfg``` and an input simlist (```InputFiles/simlist.csv``` by default), then creates a set of MESA and Stella simulation grids with the parameters specified within the simlist.  ```Scheduler.py``` then runs the MESA simulations.  After completing the MESA sims, it runs the Stella component of the simulations in parallel, each on its own thread.  Stella has limited parallelization, so this really speeds up the process.  The full output data is held within ```mesa-24.08.1/ModelGrids```, though photometry and various explosion profiles are saved to ```DataExports```.

## Getting Started

### Installation

Make sure you have Docker installed on your system.  Read the installation guidelines for your system here: https://docs.docker.com/engine/install/

Once you have Docker installed, clone this repo to a directory of your choice with ```git clone https://github.com/amartas/Mesa-Stella-SESNe.git```.  Then, create a docker container with ```docker build -t <name>:latest .```.  Docker will run through the Dockerfile and build the container.  You can enter the container with ```docker run -it <name>:latest```, and get to your simulations!  The primary Mesa directory is held within ```~/MESA``` alongside my scripts for grid-building.

### How to create model grids

The primary text editor is nano.  Use ```nano <file name>``` to open the editor.

#### Simlist
To define your grid, you need the actual parameters of the models you wish to create, formatted as a CSV.  The CSV must have the following columns:

- ```mass```: float, ZAMS mass ($M_\odot$)
- ```energy```: float, explosion energy ($10^{50}$ ergs)
- ```windscalar```: float, Dutch wind scaling factor $\eta$
- ```metallicity```: float, metal mass fraction
- ```hefrac```: float, helium mass fraction
- ```ni56```: float, Ni56 mass ($M_\odot$)
- ```alpha_MLT```: float, $\alpha_\mathrm{MLT}$, which is a parameter used in mixing length calculations
- ```alpha_semiconv```: float, $\alpha_\mathrm{sc}$, which is a parameter used for semiconvective mixing
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

You should have everything set up now!  All you need to do now is run ```MesaStellaCore.py``` in the Python environment from earlier, and it'll start chugging along!  Explosion profiles and the light curves are exported to ```DataExports```, but all the output data is stored in subdirectories within ```mesa-24.08.1/ModelGrids```.  Go read the MESA documentation to learn to read it!  Make sure to move these sims somewhere else *outside* the parent directory, as ```MesaStellaCore.py``` will *not* overwrite these sims if you are rerunning with identical input parameters.

