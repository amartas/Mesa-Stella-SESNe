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
import glob
from io import StringIO
import re
import h5py
from scipy.optimize import curve_fit as curve_fit

###

# Define main directories
MainDir = os.path.dirname(os.path.realpath(__file__))
GridDir = os.path.join(MainDir, "ModelGrids/")
ProgOptimizeDir = os.path.join(MainDir, "ProgOptimize")
SourceDir = os.path.join(GridDir, "000_Source")
#SourceDir_20M = os.path.join(GridDir, "000_Source_20M")
DataDir = os.path.join(MainDir, "DataExports")
InputDir = os.path.join(MainDir, "InputFiles")

# Load user-set information from config file
config = configparser.ConfigParser(inline_comment_prefixes="#")
config.read("SetupConfig.cfg")

# MAIN
User = config["MAIN"]["User"]
MesaSDKDir = config["MAIN"]["MesaSDK_Dir"]
NumThreads = eval(config["MAIN"]["NumThreads"])
SimThreads = eval(config["MAIN"]["SimThreads"])
SimlistName = config["MAIN"]["SimlistName"]
SimMemory = eval(config["MAIN"]["SimMemory"])
TotMemory = eval(config["MAIN"]["TotMemory"])
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
    def __init__(self,
        mass,
        windscalar,
        metallicity,
        hefrac,
        energy,
        ni56,
        alpha_noH,
        alpha_H,
        alpha_semiconv,
        CO_mthd,
        CO_f,
        CO_f0,
        csmtime,
        csmrate,
        csmvelo,
        CSMOptimize,
        ProgOptimize,
        DoPre,
        gridtag
        ):
        
        # Base stellar parameters
        self.mass = mass # Mass
        self.energy = energy # Energy
        self.ni56 = ni56 # Nickel produced in explosion
        self.windscalar = windscalar # Dutch wind scaling factor
        self.metallicity = metallicity # Metallicity
        self.hefrac = hefrac # Helium mass fraction
        
        # MLT parameters + semiconvection
        self.alpha_noH = alpha_noH # MLT mixing alpha for H-rich regions
        self.alpha_H = alpha_H # MLT mixing alpha for H-poor regions
        self.alpha_semiconv = alpha_semiconv # Semiconvection mixing parameter
        
        # Convective overshoot
        self.CO_mthd = CO_mthd # CO method (none, step, exponential)
        self.CO_f = CO_f # CO parameter f
        self.CO_f0 = CO_f0 # CO parameter f0
        
        if CO_mthd == "none":
            self.mthdint = 0
        elif CO_mthd == "step":
            self.mthdint = 1
        elif CO_mthd == "exponential":
            self.mthdint = 2
        else:
            raise ValueError("Convective overshoot method invalid (none, step, exponential)")
        
        # CSM parameters
        self.csmtime = csmtime # Mass loss duration
        self.csmrate = csmrate # Mass loss rate
        self.csmvelo = csmvelo # Velocity of lost mass
        
        # Optimization
        self.CSMOptimize = CSMOptimize
        self.ProgOptimize = ProgOptimize
        self.DoPre = DoPre
        
        # Others
        self.GridTag = gridtag # Name for exporting data
        self.TheSourceDir = SourceDir # Which source directory
        
        dirname =(
        f"M{self.mass}_"
        f"E{self.energy}_"
        f"Ni{self.ni56}_"
        f"Z{self.metallicity}_"
        f"He{self.hefrac}_"
        f"Eta{self.windscalar}_"
        f"AlH{self.alpha_H}_"
        f"AlC{self.alpha_noH}"
        f"AlSC{self.alpha_semiconv}"
        f"COf{self.CO_f}_"
        f"COff{self.CO_f0}_"
        f"COMTH{self.mthdint}_"
        f"WT{self.csmtime}_"
        f"WR{self.csmrate}_"
        f"WV{self.csmvelo}"
        )
        self.dirname = dirname
        self.simdir = os.path.join(GridDir, dirname)
        
        self.premodname =(
        f"M{self.mass}_"
        f"Z{self.metallicity}_"
        f"He{self.hefrac}_"
        f"Eta{self.windscalar}_"
        f"AlH{self.alpha_H}_"
        f"AlC{self.alpha_noH}_"
        f"AlSC{self.alpha_semiconv}_"
        f"COf{self.CO_f}_"
        f"COff{self.CO_f0}_"
        f"COMTH{self.mthdint}"
        ".mod"
        )
        
        ### Set up logging

        Logfiles = [f"Logs/Sim_{self.dirname}.log", f"Logs/MESA_{self.dirname}.log", f"Logs/Stella_{self.dirname}.log"]

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
        self.SimLogger = logging.getLogger(__name__)
        self.SimLogger.setLevel(logging.INFO)

        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO) # Minimum level passed to console

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self.SimLogger.addHandler(console_handler)

        file_handler = logging.FileHandler(f"Logs/Sim_{self.dirname}.log", mode="w")
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        self.SimLogger.addHandler(file_handler)

        # MESA and Stella logs

        self.MesaLogger = logging.getLogger("MESA")
        self.MesaLogger.propagate = False
        self.MesaLogger.setLevel(logging.INFO)
        MesaFileHandler = logging.FileHandler(f"Logs/MESA_{self.dirname}.log", mode="w")
        MesaFileHandler.setLevel(logging.INFO)
        MesaFileHandler.setFormatter(formatter)
        self.MesaLogger.addHandler(MesaFileHandler)

        self.StellaLogger = logging.getLogger("Stella")
        self.StellaLogger.propagate = False
        self.StellaLogger.setLevel(logging.INFO)
        StellaFileHandler = logging.FileHandler(f"Logs/Stella_{self.dirname}.log", mode="w")
        StellaFileHandler.setLevel(logging.INFO)
        StellaFileHandler.setFormatter(formatter)
        self.StellaLogger.addHandler(StellaFileHandler)
        
        self.SimLogger.info(f"Created Sim instance with name {dirname}")
        self.SimLogger.info(f"CSM optimization is: {self.CSMOptimize}")
        
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
        
        mesadir = os.path.join(MainDir, "mesa-25.12.1")
    
        if User == "root":  
            env = (
            f'export MESA_DIR="{mesadir}"\n'
            f'export OMP_NUM_THREADS={SimThreads}\n'
            f'export MESASDK_ROOT="/root/mesasdk"\n'
            'source "$MESASDK_ROOT/bin/mesasdk_init.sh"\n'
            'export PATH="$PATH:$MESA_DIR/scripts/shmesa"'
            )
        else:
            env = (
            f'export MESA_DIR="{mesadir}"\n'
            f'export OMP_NUM_THREADS={SimThreads}\n'
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
        
        self.SimLogger.info("Updated source shell scripts")

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
            self.SimLogger.info(f"Set line {line} in file '{filepath}' to {value}")
        
        # Copy source to simdir
        shutil.copytree(self.TheSourceDir,
                        self.simdir
                        )        
        # Check if it's a CSM sim or not to configure MESA for optimized speed during directory setup
        if self.CSMOptimize == True:
            
            shutil.copy(os.path.join(InputDir, "PreCSM.mod"),
                        os.path.join(self.simdir, "PostCC/shock_part4.mod")
                        )
            
            self.SimLogger.info("Copied CSM accelerator model")
        # If progenitor optimization is true, check if a progenitor model has already been built.  If so, use it
        if self.ProgOptimize == True:
            premod_path = os.path.join(ProgOptimizeDir, self.premodname)
            if os.path.exists(premod_path) == True:
                shutil.copy(premod_path,
                            os.path.join(self.simdir, "PostCC/pre_ccsn.mod")
                            )
                self.SimLogger.info(f"Copied progenitor accelerator model '{self.premodname}'")
            else:
                self.SimLogger.error(f"Progenitor accelerator model '{self.premodname}' does not exist.  Failing this model and continuing...")
                raise ValueError("Progenitor accelerator model DNE")


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
            self.hefrac
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
        
        # Mixing lengths for MLT + semiconvection
        ConfigInlist(
            os.path.join(self.simdir, "PreCC/inlist_common"),
            105,
            self.alpha_H
            )
        ConfigInlist(
            os.path.join(self.simdir, "PreCC/inlist_common"),
            106,
            self.alpha_noH
            )
        ConfigInlist(
            os.path.join(self.simdir, "PreCC/inlist_common"),
            115,
            self.alpha_semiconv
            )
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_common"),
            97,
            self.alpha_H
            )
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_common"),
            98,
            self.alpha_noH
            )
        ConfigInlist(
            os.path.join(self.simdir, "PostCC/inlist_common"),
            92,
            self.alpha_semiconv
            )
        
        # Convective overshoot
        f_mthd_lines = [144, 152, 160]
        f_lines = [148, 156, 164]
        f0_lines = [149, 157, 165]
        
        for f in f_lines:
            ConfigInlist(
                os.path.join(self.simdir, "PreCC/inlist_common"),
                f,
                self.CO_f
                )
            
        for f0 in f0_lines:
            ConfigInlist(
                os.path.join(self.simdir, "PreCC/inlist_common"),
                f0,
                self.CO_f0
                )
            
        for mthd in f_mthd_lines:
            ConfigInlist(
                os.path.join(self.simdir, "PreCC/inlist_common"),
                mthd,
                self.CO_mthd
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
        
        self.SimLogger.info("Updated inlists")
        
    def RunSim(self, simtype):
        """Runs a simulation of a given type (PreCC, PostCC, Stella)"""
        # so hip to be square
        def RunShellWithMESA(filename):
            process = subprocess.Popen(f"./{filename}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
    
            # This prints the output as it comes from the script
            self.MesaLogger.info(f"{'~'*10} Beginning MESA simulation {'~'*10}")
            for line in process.stdout:
                self.MesaLogger.info(line)
            
            for line in process.stderr:
                self.MesaLogger.error(line)
                
            self.MesaLogger.info(f"{'~'*10} Finished MESA simulation {'~'*10}")
            
        def RunShellWithStella(filename):
            self.StellaLogger.info(f"{'~'*10} Beginning Stella simulation {'~'*10}")
            process = subprocess.Popen(f"./{filename}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, shell=True)
            # This prints the output as it comes from the script
            for line in process.stdout:
                self.StellaLogger.info(line)
                
            
            for line in process.stderr:
                self.StellaLogger.error(line)
            
            self.StellaLogger.info(f"{'~'*10} Finished Stella simulation {'~'*10}")
            
        if simtype == "PreCC":
            # Run the shell script
            self.SimLogger.info("Beginning pre-core-collapse simulation")
            os.chdir(os.path.join(self.simdir, "PreCC")) # This is necessary for MESA's rn script to run properly
            RunShellWithMESA("run_mesa.sh")
            self.SimLogger.info("Finished pre-core-collapse simulation")
            # This copies the output from the pre-CC model to the post-CC model
            
            shutil.copyfile(
                os.path.join(self.simdir, "PreCC/final.mod"),
                os.path.join(self.simdir, "PostCC/pre_ccsn.mod")
                )
            
            self.SimLogger.info("Copied pre-core-collapse model to the post-core-collapse simulation")
            
            shutil.copyfile(
                os.path.join(self.simdir, "PreCC/final.mod"),
                os.path.join(ProgOptimizeDir, self.premodname)
                )
            
            self.SimLogger.info("Saved pre-core-collapse model to ProgOptimize in case we can use it later")

        elif simtype == "PostCC":
            # Choose to run optimized method or not
            if self.CSMOptimize == True:
                self.SimLogger.info("Beginning post-core-collapse simulation with optimization")
                os.chdir(os.path.join(self.simdir, "PostCC"))
                RunShellWithMESA("run_mesa_optimized.sh")
            else:
                self.SimLogger.info("Beginning post-core-collapse simulation without optimization")
                os.chdir(os.path.join(self.simdir, "PostCC"))
                RunShellWithMESA("run_mesa.sh")
            
            self.SimLogger.info("Finished post-core-collapse simulation")
                
        elif simtype == "Stella":
            
            self.SimLogger.info("Beginning Stella simulation")
            os.chdir(os.path.join(self.simdir, "PostCC"))
            
            # Copy MESA's output to Stella
            shutil.copyfile("mesa.abn", os.path.join(self.simdir, "PostCC/stella/modmake/mesa.abn"))
            shutil.copyfile("mesa.hyd", os.path.join(self.simdir, "PostCC/stella/modmake/mesa.hyd"))
            
            # Run the shell script
            os.chdir(os.path.join("stella"))
            
            RunShellWithStella("run_stella.sh")
            
            self.SimLogger.info("Finished Stella simulation")
    
        else:
            self.SimLogger.error(f"Invalid simulation type '{simtype}' was passed to RunSim")
            raise InvalidSimType("Invalid simulation type - is it 'PreCC', 'PostCC', or 'Stella'?")
    
    def ExportProfile(self, homologous=False):
        # Columns in the Stella model snapshot
        cols = [
                "index",
                "cell_mass",
                "cell_center_mass",
                "radius",
                "velocity",
                "density",
                "rad_pressure",
                "temperature",
                "rad_temperature",
                "opacity",
                "tau",
                "outer_mass",
                "outer_radius",
                "h1",
                "he3",
                "he4",
                "c12",
                "n14",
                "o16",
                "ne20",
                "na23",
                "mg24",
                "si28",
                "s32",
                "ar36",
                "ca40",
                "ti44",
                "cr48",
                "cr60",
                "fe52",
                "fe54",
                "fe56",
                "co56",
                "ni56",
                "luminosity",
                "n_bar",
                "n_e"
                ]
        def ReadProfile(fn):
            with open(fn, "r") as f:
                lines = f.readlines()
            
            # find header
            for i, line in enumerate(lines):
                if "mass of cell" in line:
                    header_line_index = i
                    break

            data_lines = lines[header_line_index + 1:]
            df = pd.read_csv(
                StringIO("".join(data_lines)),
                sep="\s+",
                names=cols,
                engine="python"
            )
            df.set_index("index", inplace=True)
            return df
        
        resdir = os.path.join(self.simdir, "PostCC/stella/res")        
        os.chdir(resdir)
    
        # Get profiles and corresponding days
        files = np.array(glob.glob("mesa.day*_post_Lbol_max.data"))
        days = np.array([])
        for f in files:
            m = re.search(r"day(\d+(?:\.\d+)?)_post", f)
            if m:
                val = float(m.group(1))
                days = np.append(days, val)
        
        # Sort for convenience
        sind = np.argsort(days)
        days = days[sind]
        files = files[sind]
        
        nzones = len(ReadProfile(files[0]))
        mname = f"Model_H_{self.dirname}.h5" if homologous else f"Model_{self.dirname}.h5"

        ExportDir = os.path.join(DataDir, self.GridTag)

        if os.path.exists(ExportDir) == False:
            os.mkdir(ExportDir)

        fn = os.path.join(ExportDir, mname)
        
        if os.path.exists(fn):
            self.SimLogger.warning(f"The h5 file {fn} already exists, I will skip and not overwrite")
            return
        
        f = h5py.File(fn, mode="a")
        
        # actual time and the homologous time
        times = f.create_dataset("time", (len(days),), dtype="f")
        times[:] = days
        times_homologous = f.create_dataset("time_homologous", (len(days),), dtype="f")
        homologous_dev = f.create_dataset("homologous_deviation", (len(days),), dtype="f")
        
        # create datasets
        for col in cols[1:]: # Ignore index column
            f.create_dataset(col, (len(days), nzones), dtype="f")
        
        # add the data from the Stella model snapshots
        for i in range(len(files)):
            file = files[i]
            day = days[i]
            df = ReadProfile(file)
            rad = df["radius"].to_numpy()
            velo = df["velocity"].to_numpy()
            
            def HExpansion(r, t):
                return r / t
            popt, _ = curve_fit(HExpansion, rad, velo, p0=(day*3600*24))
            time_homol = popt[0]
            times_homologous[i] = time_homol / 3600 / 24 # to days
            
            # Enforce homologous expansion if requested
            if homologous:
                v_homo = rad / time_homol
                
                # Find the accumulated deviation from homologous expansion
                accerr = np.trapezoid(abs(velo - v_homo), x=rad)
                acc = np.trapezoid(v_homo, x=rad)
                dev = accerr / acc
                
                homologous_dev[i] = dev
                df["velocity"] = v_homo
            
            for col in cols[1:]:
                f[col][i, :] = df[col].to_numpy()
        
        self.SimLogger.info(f"Exported {self.simdir} to {fn}")
        
    def ExportPhotometry(self):
        datapath = os.path.join(self.simdir, "PostCC/stella/res/mesa.tt")
        
        # This searches mesa.tt for the CSV header, finds it, and then only interprets it and the subsequent lines as a data table
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
            
            ExportDir = os.path.join(DataDir, self.GridTag)
            
            if os.path.exists(ExportDir) == False:
                os.mkdir(ExportDir)
            
            fpfinal = os.path.join(ExportDir, f"Phot_{self.dirname}.csv")
            data.to_csv(fpfinal)
            self.SimLogger.info(f"Exported photometry from '{self.simdir}'")
        else:
            self.SimLogger.error("An error occured within ExportPhotometry; the data header may not have been found")

class TimeoutException(Exception):
    pass

def signal_handler(signum, frame):
    raise TimeoutException("Stage timed out")
