import numpy as np
import os
import pandas as pd
import shutil
import configparser
import subprocess
import signal
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import logging
import time



### Set up logging

Logfiles = ["Logs/Latest.log", "Logs/MESA.log", "Logs/Stella.log"]

# Check if any of the log files exist
if any(os.path.exists(log) for log in Logfiles):
    # Create a timestamp in year-month-day_hms format
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    archive_dir = f"Logs/Archive/{timestamp}"
    os.makedirs(archive_dir, exist_ok=True)

    # Move existing logs into the archive directory
    for log in Logfiles:
        if os.path.exists(log):
            shutil.move(log, archive_dir)

# Primary logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO) # Minimum level passed to console

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

file_handler = logging.FileHandler("Logs/Latest.log", mode="w")
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# MESA and Stella logs

MesaLogger = logging.getLogger("MESA")
MesaLogger.propagate = False
MesaLogger.setLevel(logging.INFO)
MesaFileHandler = logging.FileHandler("Logs/MESA.log", mode="w")
MesaFileHandler.setLevel(logging.INFO)
MesaFileHandler.setFormatter(formatter)
MesaLogger.addHandler(MesaFileHandler)

StellaLogger = logging.getLogger("Stella")
StellaLogger.propagate = False
StellaLogger.setLevel(logging.INFO)
StellaFileHandler = logging.FileHandler("Logs/Stella.log", mode="w")
StellaFileHandler.setLevel(logging.INFO)
StellaFileHandler.setFormatter(formatter)
StellaLogger.addHandler(StellaFileHandler)

###

# Define main directories
MainDir = os.path.dirname(os.path.realpath(__file__))
GridDir = os.path.join(MainDir, "mesa-24.08.1/ModelGrids/")
ProgOptimizeDir = os.path.join(MainDir, "ProgOptimize")
SourceDir_12M = os.path.join(GridDir, "000_Source_12M")
SourceDir_20M = os.path.join(GridDir, "000_Source_20M")
DataDir = os.path.join(MainDir, "DataExports")
InputDir = os.path.join(MainDir, "InputFiles")

# Load user-set information from config file
config = configparser.ConfigParser(inline_comment_prefixes="#")
config.read("SetupConfig.cfg")

# MAIN
User = config["MAIN"]["User"]
MesaSDKDir = config["MAIN"]["MesaSDK_Dir"]
NumThreads = eval(config["MAIN"]["NumThreads"])
SimlistName = config["MAIN"]["SimlistName"]
TimeoutTime = eval(config["MAIN"]["TimeoutTime"])

# Define observed data globals
obsinfo = {}
obsinds = {}
obsmags = {}
obstimes = {}
dmags = {}
ll0 = {}

def BandConv(data, bandout):
    out = {}
    
    for band in ["U", "B", "V", "R", "I"]:
        out[band] = data[f"M{band}"]
    
    # These conversions are based upon Lupton (2005).  These equations are from solving the system with Mathematica for the ugriz components.
    out["u"] = 7.30606e-11 * (-5.66804e9 + 4.31686e10*data["MB"] + 77341.9 * np.sqrt(3.58517e10 - 3.81558e10*data["MB"] + 1.453e11*data["MB"]**2 - 3.768e11*data["MV"]))
    out["g"] = 3.44116e-8 * (-1.90779e6 + 1.453e7*data["MB"] - 38.1182 * np.sqrt(3.58517e10 - 3.81558e10*data["MB"] + 1.453e11*data["MB"]**2 - 3.768e11*data["MV"]))
    out["r"] = 4.75955e-11 * (1.14344e9 - 7.65731e9*data["MB"] + 20088.3 * np.sqrt(3.58517e10 - 3.81558e10*data["MB"] + 1.453e11*data["MB"]**2 - 3.768e11*data["MV"]) + 3.6325e10*data["MV"])
    out["i"] = 1.29688e-13 * (2.76958e12 + 6.7614e12*data["MB"] + 2.6263e13*data["MR"] - 1.7738e7 * np.sqrt(3.58517e10 - 3.81558e10*data["MB"] + 1.453e11*data["MB"]**2 - 3.768e11*data["MV"]) - 3.2075e13*data["MV"])
    out["z"] = 2.05854e-15 * (6.85196e14 + 4.25968e14*data["MB"] + 1.65457e15*data["MR"] - 1.11749e9 *np.sqrt(3.58517e10 - 3.81558e10*data["MB"] + 1.453e11*data["MB"]**2 - 3.768e11*data["MV"]) - 2.02072e15*data["MV"])

    return out[bandout]

class InvalidSimType(Exception):
    pass

class Sim:
    def __init__(self, mass, energy, ni56, windscalar, metallicity, HeFrac, csmtime, csmrate, csmvelo, CSMOptimize, ProgOptimize, gridtag):
        # Non-CSM parameters
        self.mass = mass
        self.energy = energy
        self.ni56 = ni56
        self.windscalar = windscalar
        self.metallicity = metallicity
        self.HeFrac = HeFrac
        self.ProgOptimize = ProgOptimize
        
        # CSM parameters
        self.csmtime = csmtime
        self.csmrate = csmrate
        self.csmvelo = csmvelo
        self.CSMOptimize = CSMOptimize
        
        # Others
        self.GridTag = gridtag
        
        # Select pre-CC model based upon mass
        if self.mass <= 18:
            self.TheSourceDir = SourceDir_12M
        else:
            self.TheSourceDir = SourceDir_20M
        
        dirname =(
        f"M{self.mass}_" # Mass
        f"E{self.energy}_" # Energy
        f"Ni{self.ni56}_" # Nickel
        f"Z{self.metallicity}_" # Metallicity
        f"He{self.HeFrac}_" # Helium mass fraction
        f"Eta{self.windscalar}_" # Wind scaling factor
        f"WT{self.csmtime}_" # Mass loss duration
        f"WR{self.csmrate}_" # Mass loss rate
        f"WV{self.csmvelo}" # Velocity of lost mass
        )
        self.dirname = dirname
        self.simdir = os.path.join(GridDir, dirname)
        
        self.premodname =(
        f"M{self.mass}_" # Mass
        f"Z{self.metallicity}_" # Metallicity
        f"He{self.HeFrac}_" # Helium mass fraction
        f"Eta{self.windscalar}" # Wind scaling factor    
        ".mod"
        )
        
        logger.info(f"Created Sim instance with name {dirname}")
        logger.info(f"CSM optimization is: {self.CSMOptimize}")
        
    def MakeSource(self):
        """Creates source shell scripts for running the make and run files later on"""
        
        def UpdateBlock(fp, block, start_marker = "# BEGIN BLOCK", end_marker = "# END BLOCK"):
            with open(fp, "r", encoding="utf-8") as f:
                lines = f.readlines()
        
            # Find indices of the start and end markers
            start_index = next(i for i, line in enumerate(lines) if line.strip() == "# BEGIN BLOCK")
            end_index   = next(i for i, line in enumerate(lines) if line.strip() == "# END BLOCK")
        
            # Remove the old block (including markers)
            del lines[start_index:end_index + 1]
        
            # Create the new block with markers
            blocklines = block.split("\n")
            new_block = [
                start_marker + "\n",
                *[line + "\n" for line in blocklines],
                end_marker + "\n"
            ]
        
            # Insert the new block in place of the old block
            lines[start_index:start_index] = new_block
        
            # Write updated lines back to the file
            with open(fp, "w", encoding="utf-8") as f:
                f.writelines(lines)
        
        mesadir = os.path.join(MainDir, "mesa-24.08.1")
    
        if User == "root":  
            env = (
            f'export MESA_DIR="{mesadir}"\n'
            f'export OMP_NUM_THREADS={NumThreads}\n'
            f'export MESASDK_ROOT="/root/mesasdk"\n'
            'source "$MESASDK_ROOT/bin/mesasdk_init.sh"\n'
            'export PATH="$PATH:$MESA_DIR/scripts/shmesa"'
            )
        else:
            env = (
            f'export MESA_DIR="{mesadir}"\n'
            f'export OMP_NUM_THREADS={NumThreads}\n'
            f'export MESASDK_ROOT="{MesaSDKDir}"\n'
            'source "$MESASDK_ROOT/bin/mesasdk_init.sh"\n'
            'export PATH="$PATH:$MESA_DIR/scripts/shmesa"'
            )
        
        for simtype in ["PreCC", "PostCC"]:
            fp = os.path.join(self.TheSourceDir, f"{simtype}/run_mesa.sh")
            UpdateBlock(fp, env)
            
        # Update the optimize shell script in case CSM optimization is enabled
        fp = os.path.join(self.TheSourceDir, "PostCC/run_mesa_optimized.sh")
        UpdateBlock(fp, env)
        
        if User == "root":
            env2 = (
            "\n"
            f'export MESA_DIR="{mesadir}"\n'
            f'export OMP_NUM_THREADS=1\n'
            f'export MESASDK_ROOT="{MesaSDKDir}"\n'
            'source "$MESASDK_ROOT/bin/mesasdk_init.sh"\n'
            'export PATH="$PATH:$MESA_DIR/scripts/shmesa"'
            )
        else:
            env2 = (
            "\n"
            f'export MESA_DIR="{mesadir}"\n'
            f'export OMP_NUM_THREADS=1\n'
            f'export MESASDK_ROOT="{MesaSDKDir}"\n'
            'source "$MESASDK_ROOT/bin/mesasdk_init.sh"\n'
            'export PATH="$PATH:$MESA_DIR/scripts/shmesa"'
            )
        
        # Update the shell script for Stella
        fp = os.path.join(self.TheSourceDir, "PostCC/stella/run_stella.sh")
        UpdateBlock(fp, env2)
        
        logger.info("Updated source shell scripts")

    def CreateSim(self):
        """Creates simulation directory with MESA and Stella"""
        
        def ConfigInlist(filepath, line, value):
            """Replaces placeholders in inlists with requested values"""
            with open(filepath, 'r') as file:
                lines = file.readlines()
            
            lines[line] = lines[line].replace("PLACEHOLDER", str(value))
            
            # Write the lines back to the inlist
            with open(filepath, 'w') as file:
                file.writelines(lines)
            logger.info(f"Set line {line} in file '{filepath}' to {value}")
        
        # Copy source to simdir
        shutil.copytree(self.TheSourceDir,
                        self.simdir
                        )        
        # Check if it's a CSM sim or not to configure MESA for optimized speed during directory setup
        if self.CSMOptimize == True:
            
            shutil.copy(os.path.join(InputDir, "PreCSM.mod"),
                        os.path.join(self.TheSourceDir, "PostCC/shock_part4.mod")
                        )
            
            logger.info("Copied CSM acclerator model")
        # If progenitor optimization is true, check if a progenitor model has already been built.  If so, use it
        if self.ProgOptimize == True:
            premod_path = os.path.join(ProgOptimizeDir, self.premodname)
            if os.path.exists(premod_path) == True:
                shutil.copy(premod_path,
                            os.path.join(self.simdir, "PostCC/pre_ccsn.mod")
                            )
                logger.info(f"Copied progenitor acclerator model '{self.premodname}'")


        ### Configure inlist(s) for the pre-CC MESA model
        
        # Metallicity
        ConfigInlist(
            os.path.join(self.simdir, "PreCC/inlist_mass_Z_wind_rotation"),
            3,
            self.metallicity
            )
        ConfigInlist(
            os.path.join(self.simdir, "PreCC/inlist_mass_Z_wind_rotation"),
            9,
            self.metallicity
            )
        ConfigInlist(
            os.path.join(self.simdir, "PreCC/inlist_mass_Z_wind_rotation"),
            14,
            self.metallicity
            )
        
        # Helium fraction
        ConfigInlist(
            os.path.join(self.simdir, "PreCC/inlist_mass_Z_wind_rotation"),
            15,
            self.HeFrac
            )
        
        # Mass
        ConfigInlist(
            os.path.join(self.simdir, "PreCC/inlist_mass_Z_wind_rotation"),
            13,
            self.mass
            )
        
        # Dutch wind scaling factor
        ConfigInlist(
            os.path.join(self.simdir, "PreCC/inlist_mass_Z_wind_rotation"),
            16,
            self.windscalar
            )
        
        ### Configure inlists for the post-CC MESA model
        
        # Explosion energy
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_edep"),
            49,
            str(self.energy) + "d+50"
            )
        
        # Mass
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_edep"),
            37,
            self.mass
            )
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_mass_Z"),
            12,
            self.mass
            )
        
        # Metallicity
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_mass_Z"),
            3,
            self.metallicity
            )
        
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_mass_Z"),
            7,
            self.metallicity
            )
        
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_mass_Z"),
            11,
            self.metallicity
            )
        
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_edep"),
            38,
            self.metallicity
            )
        
        # Ni56
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_shock_part3"),
            27,
            self.ni56
            )
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_shock_part5"),
            56,
            self.ni56
            )
        
        # CSM parameters
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_stella"),
            9,
            self.csmtime
            )
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_stella"),
            10,
            self.csmrate
            )
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_stella"),
            11,
            self.csmvelo
            )
            
        # Adds CSM cells if CSM is present
        if self.csmrate != 0:
            ConfigInlist(
                os.path.join(self.simdir, "PostCC/inlist_stella"),
                4,
                40
                )
        else:
            ConfigInlist(
                os.path.join(self.simdir, "PostCC/inlist_stella"),
                4,
                0
                )
        
        logger.info("Updated inlists")
        
    def RunSim(self, simtype):
        """Runs a simulation of a given type (PreCC, PostCC, Stella)"""
        # so hip to be square
        def RunShellWithMESA(filename):
            process = subprocess.Popen(f"./{filename}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    
            # This prints the output as it comes from the script
            MesaLogger.info("------------- Beginning MESA simulation -------------")
            for line in process.stdout:
                MesaLogger.info(line)
            
            for line in process.stderr:
                MesaLogger.error(line)
                
            MesaLogger.info("------------- Finished MESA simulation -------------")
            
        def RunShellWithStella(filename):
            StellaLogger.info("------------- Beginning Stella simulation -------------")
            process = subprocess.Popen(f"./{filename}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, shell=True)
            # This prints the output as it comes from the script
            for line in process.stdout:
                StellaLogger.info(line)
                
            
            for line in process.stderr:
                StellaLogger.error(line)
            
            StellaLogger.info("------------- Finished Stella simulation -------------")
            
        if simtype == "PreCC":
            # Run the shell script
            logger.info("Beginning pre-core-collapse simulation")
            os.chdir(os.path.join(self.simdir, "PreCC")) # This is necessary for MESA's rn script to run properly
            RunShellWithMESA("run_mesa.sh")
            logger.info("Finished pre-core-collapse simulation")
            # This copies the output from the pre-CC model to the post-CC model
            
            shutil.copyfile(
                os.path.join(self.simdir, "PreCC/final.mod"),
                os.path.join(self.simdir, "PostCC/pre_ccsn.mod")
                )
            
            logger.info("Copied pre-core-collapse model to the post-core-collapse simulation")
            
            shutil.copyfile(
                os.path.join(self.simdir, "PreCC/final.mod"),
                os.path.join(ProgOptimizeDir, self.premodname)
                )
            
            logger.info("Saved pre-core-collapse model to ProgOptimize in case we can use it later")

        elif simtype == "PostCC":
            # Choose to run optimized method or not
            if self.CSMOptimize == True:
                logger.info("Beginning post-core-collapse simulation with optimization")
                os.chdir(os.path.join(self.simdir, "PostCC"))
                RunShellWithMESA("run_mesa_optimized.sh")
            else:
                logger.info("Beginning post-core-collapse simulation without optimization")
                os.chdir(os.path.join(self.simdir, "PostCC"))
                RunShellWithMESA("run_mesa.sh")
            
            logger.info("Finished post-core-collapse simulation")
                
        elif simtype == "Stella":
            
            logger.info("Beginning Stella simulation")
            os.chdir(os.path.join(self.simdir, "PostCC"))
            
            # Copy MESA's output to Stella
            shutil.copyfile("mesa.abn", os.path.join(self.simdir, "PostCC/stella/modmake/mesa.abn"))
            shutil.copyfile("mesa.hyd", os.path.join(self.simdir, "PostCC/stella/modmake/mesa.hyd"))
            
            # Run the shell script
            os.chdir(os.path.join("stella"))
            
            RunShellWithStella("run_stella.sh")
            
            logger.info("Finished Stella simulation")
    
        else:
            logger.error(f"Invalid simulation type '{simtype}' was passed to RunSim")
            raise InvalidSimType("Invalid simulation type - is it 'PreCC', 'PostCC', or 'Stella'?")
    
    def ExportData(self):
        
        datapath = os.path.join(self.simdir, "PostCC/stella/res/mesa.tt")
        
        # This searches mesa.tt for the CSV header, finds it, and then only interprets it and the subsequent lines as a data table.
        csvstr = "time           Tbb         vFe        Teff      Rlast_sc   R(tau2/3)    Mbol     MU      MB      MV      MI      MR   Mbolavg  gdepos"
        with open(datapath, "r") as file:
            for linenum, line in enumerate(file):
                if csvstr in line:
                    found = True
                    break
        if found:
            data = pd.read_csv(datapath,
                               skiprows=linenum,
                               skipinitialspace=True,
                               header=0,
                               delimiter=" ",
                               )
            
            for band in ["u", "g", "r", "i", "z"]:
                data[band] = BandConv(data, band)
            
            GridDir = os.path.join(DataDir, self.GridTag)
            
            if os.path.exists(GridDir) == False:
                os.mkdir(GridDir)
            
            fpfinal = os.path.join(GridDir, f"Data_{self.dirname}.csv")
            data.to_csv(fpfinal)
            logger.info(f"Exported simulation data from '{self.simdir}'")
        else:
            logger.error("An error occured within ExportData; the data header may not have been found")

class TimeoutException(Exception):
    pass

def signal_handler(signum, frame):
    raise TimeoutException("")


# Import params from simlist
Simlist = pd.read_csv(os.path.join(InputDir, SimlistName))

Simarr = np.array([])

logger.info("Imported simlist")

# Iterate over every simulation parameter set in the simlist
for index, row in Simlist.iterrows():
    
    # Non-CSM parameters
    mass = row["mass"]
    energy = row["energy"]
    Ni56 = row["ni56"]
    metallicity = row["metallicity"]
    HeFrac = row["hefrac"]
    windscalar = row["windscalar"]
    ProgOptimize = True if row["progoptimize"] == 1 else False
    
    # CSM parameters
    csmvelo = row["csmvelo"]
    csmrate = row["csmrate"]
    csmtime = row["csmtime"]
    CSMOptimize = True if row["csmoptimize"] == 1 else False
    
    # Grid tag
    GridTag = row["gridtag"]
    
    signal.signal(signal.SIGALRM, signal_handler)
    
    # Timeout after the timeout time
    signal.alarm(TimeoutTime)
    
    try:
        sim1 = Sim(mass, energy, Ni56, windscalar, metallicity, HeFrac, csmtime, csmrate, csmvelo, CSMOptimize, ProgOptimize, GridTag)
        
        Simarr = np.append(Simarr, sim1)
        
        sim1.MakeSource()
            
        sim1.CreateSim()
        logger.info(f"Created simulation with index {index}")
        
        # If CSM optimization is off:
        if CSMOptimize != True:
            # And if progenitor optimization is off, run PreCC.  Otherwise, skip it since we're optimizing with CSM or the progenitor
            if ProgOptimize != True:
                logger.info("------------- Running pre-core-collapse model -------------")
                sim1.RunSim("PreCC")
                logger.info("------------- Finished pre-core-collapse model -------------")
        
        if ProgOptimize == True:
            logger.info("Progenitor optimization is true.  Skipping pre-CC modeling.")
        logger.info("------------- Running post-core-collapse model -------------")
        sim1.RunSim("PostCC")
        logger.info("------------- Finished pre-core-collapse model -------------")
    except Exception as err:
        logger.error(f"An exception occured while running simulation with index {index}; Exception: {err}")
    finally:
        # Disable the alarm
        signal.alarm(0)
        del sim1

logger.info("------------- Finished MESA simulations -------------")

# Run the Stella sims in parallel
results = []
Executor = ProcessPoolExecutor

logger.info("------------- Beginning Stella sl------------")

with Executor(max_workers=NumThreads) as executor:
    # Submit tasks for each instance
    futures = {executor.submit(partial(sim.RunSim, "Stella")) for sim in Simarr}

    for future in as_completed(futures):
        try:
            result = future.result()
            logger.info("Finished Stella simulation")
            results.append(result)
        except Exception as E:
            logger.info(f"Stella generated an exception: {E}")
            logger.info("Continuing with next Stella simulation...")

for sim in Simarr:
    try:
        sim.ExportData()
    except Exception as err:
        logger.info(f"Data exporting threw an exception: {err}")


logger.info("------------- Finished Stella simulations.  Done! -------------")














