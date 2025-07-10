#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 09:52:31 2022

@author: yutasuzuki
"""

import numpy as np
import json
import glob
import os
import datetime
import matplotlib
import seaborn as sns        
from multiprocessing import Pool
from tqdm import tqdm
import scipy.io

from PupilAnalysisToolbox.lib.pre_processing_cls import re_sampling
from PupilAnalysisToolbox.lib.asc2array_cls import asc2array_cls
import shutil

sns.set()
sns.set_style("whitegrid")
sns.set_palette("Set2")

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


dt_now = datetime.datetime.now()
date = dt_now.strftime("%Y%m%d")

#%% path and condition name

rootFolder = "./results/N-back"

condition_list0 = {"audio":[],
                  "light":[],
                  "Nback":[]
                  }

mmName1 = ["0_ipRGC-","1_ipRGC+","2_light flux"]
mmName2 = ["N=1","N=2"]

condition_name = list(condition_list0.keys())


#%% experimental parameters

cfg={   
     "TIME_START":-1,
     "TIME_START_RUN":-10,
     "TIME_END":2,
     "WID_ANALYSIS":2,
     "WID_BP_ANALYSIS":2,
     "usedEye":1,
     "WID_FILTER":[],
     "RESAMPLING_RATE":250,
     # "mmFlag":False,
     "mmFlag":True,
     "normFlag":False,
     "s_trg":[],
     "visualization":True,
     "MS":False,
     "DOT_PITCH":0.278,   
     "VISUAL_DISTANCE":100,
     "acceptMSRange":4.5,
     "MIN_DURATION":0.1,
     "SCREEN_RES":[960*2, 540*2],
     "rejectFlag":[]
     }

center = np.array(cfg['SCREEN_RES'])/2
t2a = asc2array_cls(cfg)
  

#%% analysis batch
def run(iSub):

    condition_list = condition_list0.copy()

    cfg["numOfZero_original"]=[]
    cfg["numOfZero_th"]=[]
    
    datHash={"PDR":[],
             "gazeX":[],
             "gazeY":[],
             "ipRGC":[],
             "Run":[],
             "order":[],
             "response":[],
             "target":[],
             "hourCls":[],
             "RT":[],
             "sub":[]
              }
    
    datHash_run={"PDR":[],
                 "Run":[],
                 "sub":[],
                 "order":[],
                 "hourCls":[],
                 "light2":[],
                 "Sleepiness":[],
                 "Fatigue":[],
                 }
    
    for c in list(condition_list.keys()):
        datHash[c] = []
        datHash_run[c] = []
    
    cfg["rejectFlag"] = []

    #% ------------------ load json file (eye data) -------------------
    print("Processing --> " + iSub + "...")
    
    ascFileName = glob.glob(rootFolder + "/" + iSub + "/*.asc")
    ascFileName.sort()
    
    #%%
    for iRun,ascFile in tqdm(enumerate(ascFileName)):
        
        cfg["fName"] = "Run"+str(iRun)
        # print(1)
        # fDate = matFileName[iRun].split('_')[-2][:8]
        fTime = ascFileName[iRun].split('_')[-2][-6:]
        hour = int(fTime[:2])
        if 10 <= hour < 13:
            hourCls = 0
        elif 13 <= hour < 16:
            hourCls = 1
        elif 16 <= hour < 19:
            hourCls = 2
       
        matFileName = glob.glob(rootFolder + "/" + iSub + "/*"+fTime+"*.mat")
        matFileName.sort()
        f = matFileName[0]
        
        print(scipy.io.loadmat(f)["cfg"][0,0]["res_summary"][0,0]["Run1"][0,0]["N1back"][0,0]["Light"][0,0])
        
        print(scipy.io.loadmat(f)["cfg"][0,0]["res_summary"][0,0]["Run1"][0,0]["N1back"][0,0]["Light"][0,0])
        
        #% ------------------ load json file (eye data) -------------------
        dat = t2a.dataExtraction(ascFile)
    
        #% ------------------ load json file (eye data) -------------------
        eyeData,events,initialTimeVal,fs = t2a.dataParse(dat)
        ave,sigma = t2a.getAve(eyeData["pupilData"])
        
        #% ------------------ load json file (eye data) -------------------
        eyeData = t2a.blinkInterp(eyeData)
        
        if cfg["normFlag"]:
            pupilData = t2a.pupilNorm(eyeData["pupilData"], ave, sigma).reshape(-1)
        else:
            pupilData = np.mean(eyeData["pupilData"],axis=0)
            pupilData = pupilData.reshape(-1)
        
        if eyeData["usedEye"] == "L":
            iEyes=0
        else:
            iEyes=1

        zeroArray = np.repeat(eyeData["zeroArray"][iEyes].copy(),4)
        zeroArray = zeroArray[:len(pupilData)]
                
        st = np.argwhere( np.diff(zeroArray) == 1).reshape(-1)
        ed = np.argwhere( np.diff(zeroArray) == -1).reshape(-1)
        if len(ed) < len(st):
            ed = np.r_[ed,len(zeroArray)]
        
        # gazeX = np.mean(eyeData["gazeX"],axis=0)
        # gazeY = np.mean(eyeData["gazeY"],axis=0)
        
        cfg["rejectFlag"].append(eyeData["rejectFlag"])
       
        cfg["SAMPLING_RATE"] = int(fs)
      
        #% ------------------ response -------------------
        seq_onset  = [[int((int(e[0])- initialTimeVal) * (fs/1000)),e[1]] for e in events["MSG"] if e[1] == "onset_stim"]       

        # onset_fix  = [[int((int(e[0])- initialTimeVal) * (fs/1000)),e[1]] for e in events["MSG"] if e[1] == "onset_fixation"]       
        
        for c in list(condition_list.keys()):
            condition_list[c] = [int(e[2]) for e in events["MSG"] if e[1] == "Condition_" + c]           
       
        key_press = [int(int(e[0])- initialTimeVal) for e in events["MSG"] if "key_response" in e[1]]
        
        response = []
        RT = []
        for i in np.arange(len(seq_onset)):
            
            tmp_RT = []
            if i == (len(seq_onset)-1): # edge process
                for k in key_press:
                    if seq_onset[-1][0] < k:
                        tmp_RT.append(k-seq_onset[i][0])
            else: 
                if condition_list["Nback"][i-1] != condition_list["Nback"][i]:
                    response.append(0)
                    RT.append(0)
                    continue
                
                for k in key_press:
                    if seq_onset[i][0] <= k and seq_onset[i+1][0] > k:
                       tmp_RT.append(k-seq_onset[i][0])
                        
            if len(tmp_RT) > 0:
                if tmp_RT[-1] > 300:
                    response.append(1)
                    RT.append(tmp_RT[-1])
                else:
                    response.append(0)
                    RT.append(0)
                
            else:
                response.append(0)
                RT.append(0)
       
        datHash["response"] = np.r_[datHash["response"],np.array(response)]
        datHash["RT"] = np.r_[datHash["RT"],np.array(RT)]
        
        ind = np.argwhere(np.diff(np.r_[2, np.array(condition_list["Nback"])]) != 0).reshape(-1)
        ind = np.r_[ind, len(seq_onset)]
        
        if "_7.asc" in ascFile:
            tmp_runCount=4
        else:
            tmp_runCount=1
        
        for i in np.arange(len(ind)-1):
        
            datHash_run["PDR"].append(pupilData[np.arange(int(seq_onset[ind[i]][0]+fs*cfg["TIME_START_RUN"]),int(seq_onset[ind[i+1]-1][0]))].tolist())
            datHash_run["Run"] = np.r_[datHash_run["Run"],int(tmp_runCount)]
            datHash_run["Nback"] = np.r_[datHash_run["Nback"],condition_list["Nback"][ind[i]]]
            datHash_run["light"] = np.r_[datHash_run["light"],condition_list["light"][ind[i]]]
            datHash_run["light2"] = np.r_[datHash_run["light2"],
                                          scipy.io.loadmat(f)["cfg"][0,0]["res_summary"][0,0]["Run"+str(tmp_runCount)][0,0]["N1back"][0,0]["Light"][0,0]]
            datHash_run["sub"].append(int(iSub[1:3]))     
            datHash_run["order"].append(fTime)
            datHash_run["hourCls"].append(int(hourCls))
            
            datHash_run["Fatigue"].append(int(scipy.io.loadmat(f)["cfg"][0,0]["res_summary"][0,0]["Run"+str(tmp_runCount)][0,0]["tiredness"][0,0]))
            datHash_run["Sleepiness"].append(int(scipy.io.loadmat(f)["cfg"][0,0]["res_summary"][0,0]["Run"+str(tmp_runCount)][0,0]["drowsiness"][0,0]))
                
            if i%2 == 1:
                tmp_runCount+=1
                
                
        #% ------------------ target -------------------
        target = [0]
        for i in np.arange(1,len(condition_list["audio"])):
            if condition_list["Nback"][i-1] != condition_list["Nback"][i]:
                target.append(0)
                continue
            
            if condition_list["Nback"][i] == 1:
                if condition_list["audio"][i-1] == condition_list["audio"][i]:
                    target.append(1)
                else:
                    target.append(0)
           
            else:
                if condition_list["audio"][i-2] == condition_list["audio"][i]:
                    target.append(1)
                else:
                    target.append(0)
                    
        datHash["target"] = np.r_[datHash["target"],np.array(target)]

        #% ------------------ Run ------------------
        
        tmp  = [int(e[2]) for e in events["MSG"] if e[1] == "Condition_Nback"]       
        
        if "_7.asc" in ascFile:
            tmp_runCount=4
        else:
            tmp_runCount=1
        
        for r in np.r_[np.diff(tmp),0]:
            datHash["Run"].append(tmp_runCount)            
            if r == -1:
                tmp_runCount+=1
        
        #% ------------------ PDR ------------------
        for s in seq_onset:
            tmp = np.arange(int(s[0]+fs*cfg["TIME_START"]),int(s[0]+fs*cfg["TIME_END"]))
            if len(pupilData) < tmp[-1]:
                datHash["PDR"].append(np.r_[pupilData[tmp[0]:],np.zeros(tmp[-1]-len(pupilData)+1)].tolist())
            else:
                datHash["PDR"].append(pupilData[tmp].tolist())

            datHash["order"].append(int(fTime))
            datHash["hourCls"].append(int(hourCls))

            # datHash["Run"] = np.r_[datHash["Run"],np.ones(len(seq_onset))*(iRun+1)]
        
        #% ------------------ sub & condition ----------------------- 
        datHash["sub"] = np.r_[datHash["sub"],np.ones(len(seq_onset))*int(iSub[1:3])]      
        
        for c in list(condition_list.keys()):
            datHash[c] = np.r_[datHash[c],np.array(condition_list[c])]

    #%%    
    cfg["numOfZero_th"] = np.mean(cfg["numOfZero_th"])
    
    cfg["SAMPLING_RATE"] = fs
    cfg["conditionName"] = list(condition_list0.keys())
    cfg["conditionName1"] = mmName1
    cfg["conditionName2"] = mmName2
    
    for mm in list(datHash.keys()):
        if not isinstance(datHash[mm],list):
            datHash[mm] = datHash[mm].tolist()

    for mm in list(datHash_run.keys()):
        if not isinstance(datHash_run[mm],list):
            datHash_run[mm] = datHash_run[mm].tolist()

    for mmName in ["PDR"]: #,"gazeX","gazeY"
        datHash[mmName] = re_sampling(datHash[mmName],
                                      int(len(datHash[mmName][0])*(cfg["RESAMPLING_RATE"]/cfg["SAMPLING_RATE"])))        
        datHash_run[mmName] = re_sampling(datHash_run[mmName],
                                          int(len(datHash_run[mmName][0])*(cfg["RESAMPLING_RATE"]/cfg["SAMPLING_RATE"])))
        
    datHash["cfg"] = cfg
    datHash_run["cfg"] = cfg

    with open("./data/"+date+"/"+iSub+"_trial.json","w") as f:
        json.dump(datHash,f)

    with open("./data/"+date+"/"+iSub+"_run.json","w") as f:
        json.dump(datHash_run,f)
    
    return datHash

#%%    
if __name__ == '__main__':
 
    folderName = glob.glob("./data/*")
    folderName.sort()
    
    if len(folderName) == 0:
        os.mkdir("./data/" + date)
    else:
        folderName = folderName[-1]
    
    folderList=[]
    for filename in os.listdir(rootFolder):
        if os.path.isdir(os.path.join(rootFolder, filename)): 
            folderList.append(filename)
    folderList.sort()
    # folderList=folderList[-2:]
    
    log = glob.glob(folderName+"/s*")
    log.sort()

    # if os.name == 'nt':
    #     log = [f.split('\\')[-1][:3] for f in log]
    # else:
    #     log = [f.split('/')[-1][:3] for f in log]
    
    # folderList = [f for f in folderList if not f in log]
        
    if not os.path.exists("./data/"+date):
        os.mkdir("./data/"+date)
    
    # for f in glob.glob(folderName+"/*"):
    #     if not os.path.exists("./data/"+date+"/"+f.split('/')[-1]):
    #         shutil.move(f, "./data/"+date+"/")
    
    with Pool(6) as p:
        tmp_datHash = p.map(run, folderList)
    