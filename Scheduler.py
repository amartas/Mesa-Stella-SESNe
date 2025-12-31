import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import numpy as np
from MesaStellaCore import Sim, logger, InputDir, SimlistName, NumThreads, SimThreads, SimMemory, TotMemory
from concurrent.futures import ProcessPoolExecutor
import time
import os

class ResourceScheduler:
    def __init__(self, mem_limit, thread_limit, MaxWorkers=None, min_start_interval=1.0):
        self.mem_limit = float(mem_limit)
        self.thread_limit = int(thread_limit)
        
        if MaxWorkers is None:
            MaxWorkers = self.thread_limit
        
        self.executor = ThreadPoolExecutor(max_workers=MaxWorkers)
        self.lock = threading.Lock()
        self.cond = threading.Condition(self.lock)
        self.used_memory = 0.0
        self.used_threads = 0
        self.futures = set()
        self.next_start_time = 0.0
        self.min_start_interval = float(min_start_interval)
    
    def AcquireResources(self, mem_cost, thread_cost):
        with self.cond:
            # wait until resources are available
            while ((self.used_memory + mem_cost > self.mem_limit) or
                   (self.used_threads + thread_cost > self.thread_limit)):
                self.cond.wait()
            # alloc
            self.used_memory += mem_cost
            self.used_threads += thread_cost
    
    def ReleaseResources(self, mem_cost, thread_cost):
        with self.cond:
            self.used_memory -= mem_cost
            self.used_threads -= thread_cost
            if self.used_memory < 0:
                self.used_memory = 0
            if self.used_threads < 0:
                self.used_threads = 0
            self.cond.notify_all()
            
    def submit(self, func, *, mem_cost, thread_cost, **kwargs):
        self.AcquireResources(mem_cost, thread_cost)

        def wrapper():
            try:
                # stagger starts so make doesn't freak out
                if self.min_start_interval > 0:
                    with self.cond:
                        now = time.monotonic()
                        if now < self.next_start_time:
                            delay = self.next_start_time - now
                        else:
                            delay = 0.0
                        # reserve the next start time
                        self.next_start_time = max(self.next_start_time, now) + self.min_start_interval

                    if delay > 0:
                        time.sleep(delay)

                return func(**kwargs)
            finally:
                self.ReleaseResources(mem_cost, thread_cost)

        future = self.executor.submit(wrapper)
        with self.lock:
            self.futures.add(future)
        future.add_done_callback(lambda f: self.futures.discard(f))

        return future
        
    def shutdown(self):
        self.executor.shutdown(wait=True)


def BuildSims():
    simlist_path = os.path.join(InputDir, SimlistName)
    df = pd.read_csv(simlist_path)

    sims = []

    for _, row in df.iterrows():
        mass = row["mass"]
        metallicity = row["metallicity"]
        hefrac = row["hefrac"]
        windscalar = row["windscalar"]
        
        energy = row["energy"]
        Ni56 = row["ni56"]
        
        alpha_H = row["alpha_H"]
        alpha_noH = row["alpha_noH"]
        alpha_semiconv = row["alpha_semiconv"]
        
        CO_mthd = row["CO_mthd"]
        CO_f = row["CO_f"]
        CO_f0 = row["CO_f0"]
        
        csmvelo = row["csmvelo"]
        csmrate = row["csmrate"]
        csmtime = row["csmtime"]
        
        ProgOptimize = bool(row["progoptimize"])
        CSMOptimize = bool(row["csmoptimize"])

        GridTag = row["gridtag"]

        sim = Sim(mass,
                  windscalar,
                  metallicity,
                  hefrac,
                  energy,
                  Ni56,
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
                  GridTag
                  )

        sim.MakeSource()
        sim.CreateSim()
        sims.append(sim)

    logger.info(f"Created {len(sims)} simulation directories")
    return sims

def RunMesaScheduled(
    total_memory_gb,
    per_sim_memory_gb,
    total_threads=None,
    per_sim_threads=1,
):
    if total_threads is None:
        total_threads = NumThreads

    sims = BuildSims()

    scheduler = ResourceScheduler(
        mem_limit=total_memory_gb,
        thread_limit=total_threads,
        MaxWorkers=total_threads,
        min_start_interval = 2.0,
    )

    all_futures = []

    def run_stage(sim, stage):
        logger.info(f"Starting {stage} for {sim.dirname}")
        sim.RunSim(stage)
        logger.info(f"Finished {stage} for {sim.dirname}")

    for sim in sims:
        # per sim resource consumption
        sim_mem = getattr(sim, "mesa_mem_gb", per_sim_memory_gb)
        sim_threads = getattr(sim, "mesa_threads", per_sim_threads)

        if not sim.CSMOptimize and not sim.ProgOptimize:
            # Need PreCC, then PostCC
            pre_future = scheduler.submit(
                run_stage, mem_cost=sim_mem, thread_cost=sim_threads,
                sim=sim, stage="PreCC"
            )
            all_futures.append(pre_future)

            def schedule_postcc(pre_fut, sim=sim, sim_mem=sim_mem, sim_threads=sim_threads):
                if pre_fut.exception() is not None:
                    logger.error(f"PreCC failed for {sim.dirname}, skipping PostCC.")
                    return
                post_future = scheduler.submit(
                    run_stage, mem_cost=sim_mem, thread_cost=sim_threads,
                    sim=sim, stage="PostCC"
                )
                all_futures.append(post_future)

            pre_future.add_done_callback(schedule_postcc)

        else:
            post_future = scheduler.submit(
                run_stage, mem_cost=sim_mem, thread_cost=sim_threads,
                sim=sim, stage="PostCC"
            )
            all_futures.append(post_future)

    # wait for all MESA futures to finish before running Stella
    for fut in as_completed(all_futures):
        try:
            _ = fut.result()
        except Exception as exc:
            logger.error(f"MESA stage raised exception: {exc}")

    scheduler.shutdown()

    logger.info("All MESA stages completed, ready to run Stella")
    return sims

def RunStella(sim):
    sim.RunSim("Stella")
    sim.ExportPhotometry()
    sim.ExportProfile()

def RunExplosion(sims):
    logger.info(f"{'~'*10} Beginning Stella simulations {'~'*10}")

    results = []
    with ProcessPoolExecutor(max_workers=NumThreads) as executor:
        futures = {executor.submit(RunStella, sim): sim for sim in sims}
        for fut in as_completed(futures):
            sim = futures[fut]
            try:
                fut.result()
                logger.info(f"Finished Stella + exports for {sim.dirname}")
            except Exception as exc:
                logger.error(f"Stella/extract failed for {sim.dirname}: {exc}")
                results.append(exc)

    logger.info(f"{'~'*10} Finished Stella simulations. Done! {'~'*10}")



sims = RunMesaScheduled(
        total_memory_gb=TotMemory,
        per_sim_memory_gb=SimMemory,
        total_threads=NumThreads,
        per_sim_threads=SimThreads,
    )
RunExplosion(sims)


