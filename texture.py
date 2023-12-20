from matflow import load_workflow
from pathlib import Path
from damask import VTK
import damask ; from damask import Grid
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm
import defdap
from defdap.quat import Quat
import h5py
import re
import numpy as np
import os
import shutil # "
import matplotlib


def write_quat_txts(workflow_path, phase):
    """Extract quaternions from a workflow result dir to a numpy array, write array to .txt file for MTEX"""
    
    workflow = load_workflow(workflow_path)
    for tasknum in range(0, len(workflow.tasks)):
        if "simulate_volume_element_loading" in str(workflow.tasks[tasknum].unique_name):
            print(f"Extracting quaternions from task_{str(tasknum+1)}_{workflow.tasks[tasknum].name}...")
            Path(workflow_path+"/task_"+str(tasknum+1)+"_"+workflow.tasks[tasknum].unique_name+"/"+phase+"_oris/").mkdir(parents=True, exist_ok=True)
        
            oris = workflow.tasks[tasknum].elements[0].outputs.volume_element_response['phase_data'][phase+'_orientations']['data']['quaternions']
    
            n_incs = oris.shape[0] # oris=(n_incs, n_quats, n_quatcomps)
            for inc in tqdm(range(0, n_incs)):
                quat_txt = open(workflow_path+"/task_"+str(tasknum+1)+"_"+workflow.tasks[tasknum].unique_name+"/"+phase+"_oris/"+phase+"_inc"+str(inc)+"_oris.txt",'w')
                quat_txt.write("xyzw\n")
        
                n_quats = oris.shape[1]
                for quat_num in range(0, n_quats):
                    quat = oris[inc, quat_num]
            
                    quat_txt.write("   %.15f   %.15f   %.15f   %.15f\n"\
                    %(quat[0], quat[1], quat[2], quat[3]))
            print(f"Saved to {workflow_path}task_{tasknum+1}_{workflow.tasks[tasknum].name}/{phase}_oris/")