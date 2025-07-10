#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 17:33:15 2023

@author: yutasuzuki
"""

import numpy as np
import matplotlib.pyplot as plt
import json
import glob
import os
import pandas as pd
import seaborn as sns        
import datetime
# from statsmodels.stats.anova import AnovaRM
# import pyper
import scipy
from tqdm import tqdm

# from LightSource import LightSource,showSpectra
from PupilAnalysisToolbox.lib.pre_processing_cls import re_sampling
# from rejectBlink_PCA import rejectBlink_PCA

from PupilAnalysisToolbox.lib.pre_processing_cls import ov,fft_ave,hanning,getfft,zero_padding

sns.set(font_scale=1.2)
sns.set_style("whitegrid")
sns.set_palette("Set2")

dt_now = datetime.datetime.now()
date = dt_now.strftime("%Y%m%d")

if not os.path.exists("./figure/"+date):
    os.makedirs("./figure/"+date,exist_ok=True)
   
ext='.png'

#%% time course


datHashRun = {}

cfg={
    "windowL":20,
    "WID_BASELINE":[[-1,0]],
    "WID_FILTER":np.array([]),
    "METHOD":1, #subtraction
    # "METHOD":2, #subtraction
    "FLAG_LOWPASS":False,
    "visualization":False,
    "THRES_DIFF":0.02
    }


folderName = glob.glob("./data/*")
folderName.sort()
folderName = folderName[-1]

log = glob.glob(folderName+"/s*_run.json")
log.sort()

# reject due to the less response rate 
log = [l for l in log if not "s24_run.json" in l]

for subName in log:
    f = open(subName)
    tmp_datHash = json.load(f)
    f.close()

    for mmName in list(tmp_datHash.keys()):
        if mmName != "cfg":
            if not mmName in list(datHashRun.keys()):
                datHashRun[mmName] = []
            datHashRun[mmName] = datHashRun[mmName] + tmp_datHash[mmName]

cfg.update(tmp_datHash["cfg"])

numOfTime = [len(f) for f in datHashRun["PDR"]]
y1_targeted = [f[:min(numOfTime)] for f in datHashRun["PDR"]]
y1_targeted = np.array(y1_targeted)

cfg["TIME_START"] = cfg["TIME_START_RUN"]
cfg["TIME_END"] = int(y1_targeted.shape[1] / cfg["RESAMPLING_RATE"])

x = np.linspace(cfg["TIME_START"],cfg["TIME_END"],y1_targeted.shape[1])

tmp=[]
for i in np.arange(y1_targeted.shape[0]):
    tmp_y = y1_targeted[i].copy()
    tmp_y[abs(np.gradient(tmp_y)) > cfg["THRES_DIFF"]] = np.nan
    tmp.append(tmp_y)

y1_targeted = np.array(tmp).copy()

mmName1 = ["low ipRGC","high ipRGC"]
mmName2 = ["1-back","2-back"]
mmName3 = ["1st","2nd"]

#%% ################# FFT ##########################
fs = 250
overlap = 80
frame = y1_targeted.shape[1]

zeropad = 2 ** 16
zeropad = 0

# datHash['FFT'] = []
amp = []
for iTrial in tqdm(np.arange(y1_targeted.shape[0])):
    # aa
    tmp_y = y1_targeted[iTrial,:]
    
    time_array, N_ave = ov(tmp_y, fs, overlap, frame)
    time_array, acf = hanning(time_array, N_ave, frame)    
    
    time_array = zero_padding(np.array(time_array.copy()), zeropad)
    
    spectrum, fft_array, phase, freq = getfft(time_array, fs)
   
    amp.append(fft_array.mean(axis=0))
  

y1_targeted = np.array(re_sampling(y1_targeted,int(y1_targeted.shape[1]/10)))
x = np.linspace(cfg["TIME_START"],cfg["TIME_END"],y1_targeted.shape[1])

amp = np.array(amp)[:,:y1_targeted.shape[1]]
freq = freq[:y1_targeted.shape[1]]

#%% make dataframe

df_timeCourse = pd.DataFrame()
df_timeCourse["Pupil[mm]"] = y1_targeted.reshape(-1)
df_timeCourse["Amplitude"] = amp.reshape(-1)
df_timeCourse["freq"] = np.tile(freq, y1_targeted.shape[0])
df_timeCourse["Time[s]"] = np.tile(x, y1_targeted.shape[0])
df_timeCourse["vel"] = np.c_[np.zeros((y1_targeted.shape[0],1)),np.diff(y1_targeted)].reshape(-1)
df_timeCourse["ipRGC"] = np.tile(np.array(datHashRun["light"]).reshape(np.array(datHashRun["Run"]).shape[0],1), y1_targeted.shape[1]).reshape(-1)
df_timeCourse["sub"] = np.tile( np.array(datHashRun["sub"]).reshape(np.array(datHashRun["sub"]).shape[0],1), y1_targeted.shape[1]).reshape(-1)
df_timeCourse["hourCls"] = np.tile( np.array(datHashRun["hourCls"]).reshape(np.array(datHashRun["hourCls"]).shape[0],1), y1_targeted.shape[1]).reshape(-1)
df_timeCourse["Run"] = np.tile( np.array(datHashRun["Run"]).reshape(np.array(datHashRun["Run"]).shape[0],1), y1_targeted.shape[1]).reshape(-1)
df_timeCourse["Fatigue"] = np.tile( np.array(datHashRun["Fatigue"]).reshape(np.array(datHashRun["Fatigue"]).shape[0],1), y1_targeted.shape[1]).reshape(-1)
df_timeCourse["Fatigue"]  = df_timeCourse["Fatigue"] - 1

df_timeCourse["Sleepiness"] = np.tile( np.array(datHashRun["Sleepiness"]).reshape(np.array(datHashRun["Sleepiness"]).shape[0],1), y1_targeted.shape[1]).reshape(-1)
df_timeCourse["Sleepiness"] = df_timeCourse["Sleepiness"] - 1
  
for c in cfg["conditionName"][1:]:
    df_timeCourse[c] = np.tile( np.array(datHashRun[c]).reshape(np.array(datHashRun[c]).shape[0],1), y1_targeted.shape[1]).reshape(-1)

df_timeCourse["light"] = [mmName1[int(i-1)] for i in df_timeCourse["light"].values]
df_timeCourse["Nback"] = [mmName2[int(i-1)] for i in df_timeCourse["Nback"].values]

# df_timeCourse_plot = df_timeCourse.groupby(["sub","Time[s]","light","Nback"],as_index=False).agg('mean', numeric_only=True)

df_timeCourse_plot_ave = df_timeCourse.groupby(["sub","Time[s]","Run","Nback","light"],as_index=False).agg('mean', numeric_only=True)

df_timeCourse_plot_ave = df_timeCourse_plot_ave.dropna(how="any")
df_timeCourse_plot_ave.to_json(folderName + "/data_timeCourse.json")


df_timeCourse_plot_ave = df_timeCourse[df_timeCourse["Time[s]"]>0]

df_timeCourse_plot_ave = df_timeCourse_plot_ave.dropna(how="any")
df_ave = df_timeCourse_plot_ave.groupby(["sub","Run","Nback","light"],as_index=False).agg('mean', numeric_only=True)
df_ave.to_json(folderName + "/data_timeCourse_RunAve.json")


df_timeCourse_plot_ave = df_timeCourse_plot_ave.dropna(how="any")
df_ave = df_timeCourse_plot_ave.groupby(["sub","Nback","light"],as_index=False).agg('mean', numeric_only=True)
df_ave.to_json(folderName + "/data_timeCourse_ave.json")


#%% show time-course

tmp_plot = df_timeCourse.groupby(["sub","freq","Nback","light"],as_index=False).agg('mean', numeric_only=True)

plt.figure()
g = sns.FacetGrid(tmp_plot,col="Nback", hue="light")
g.map(sns.lineplot, "freq", "Amplitude", errorbar='se',alpha=0.4)
g.set(xlim=[0,1],ylim=[0,0.05])
plt.legend()

tmp_plot = tmp_plot[(tmp_plot["freq"]>0.4)&(tmp_plot["freq"]<0.6)]

tmp_plot = tmp_plot.groupby(["sub","Nback","light"],as_index=False).agg('mean', numeric_only=True)

plt.figure()
g = sns.FacetGrid(tmp_plot,col="Nback")
g.map(sns.pointplot, "light", "Amplitude", errorbar='se').add_legend()
# plt.legend()
# sns.move_legend(g,'upper right',title = '',bbox_to_anchor=(1, 0.9),borderaxespad=0,frameon=False)

#%% show time-course

df_timeCourse_plot_ave = df_timeCourse.groupby(["sub","freq","Nback","light"],as_index=False).agg('mean', numeric_only=True)

tmp_plot = df_timeCourse_plot_ave[df_timeCourse_plot_ave["Time[s]"]>0]
# tmp_plot = tmp_plot[tmp_plot["Time[s]"]<30]

plt.figure()
grid = sns.FacetGrid(tmp_plot, row="sub",col="light", hue="Nback")
grid.map(sns.lineplot, "Time[s]", "Pupil[mm]", errorbar='se')
plt.legend(bbox_to_anchor=(1.5, 1), loc='upper right', borderaxespad=0)
plt.savefig("./figure/"+date+"/pupil_timecourse_run"+ext, bbox_inches="tight",pad_inches=0.2, dpi=400)

#%% Sleepiness

df = pd.DataFrame()
for mmName in ["light","Nback","order","Run","sub","Sleepiness","Fatigue","hourCls"]:
    df[mmName] = datHashRun[mmName]

df["PDR"] = [np.mean(p) for p in datHashRun["PDR"]]

df["Sleepiness"]  = df["Sleepiness"] - 1
df["Fatigue"]  = df["Fatigue"] - 1

df["light"] = [mmName1[int(i-1)] for i in df["light"].values]
df["Nback"] = [mmName2[int(i-1)] for i in df["Nback"].values]

df["orderNum"]=0

for iSub in df["sub"].unique():
    order = df[df["sub"]==iSub]["order"].unique()
    if int(order[0]) < int(order[1]):
        df.loc[(df["sub"]==iSub)&(df["order"]==order[0]),"orderNum"]=0
        df.loc[(df["sub"]==iSub)&(df["order"]==order[1]),"orderNum"]=1
    else:
        df.loc[(df["sub"]==iSub)&(df["order"]==order[1]),"orderNum"]=0
        df.loc[(df["sub"]==iSub)&(df["order"]==order[0]),"orderNum"]=1
        
df = df.sort_values(["sub","Nback","light","order","Run"])

# %%

# plt.figure()
# sns.lmplot(data=df,
#             x="PDR",
#             y="Sleepiness",
#             y="Fatigue",
#             col="sub",
#             hue="light",
#             col_wrap=6).add_legend()

# plt.figure()
# sns.lmplot(data=df,
#             x="PDR",
#             # y="Sleepiness",
#             y="Fatigue",
#             # col="sub",
#             hue="sub",
#             # hue="light"
#             ).add_legend()


#%% save questionaire data

df_eval = df.copy()
df_eval = df_eval[df_eval["Nback"]=="1-back"]


for iSub in df_eval["sub"].unique():
    for i,immName2 in enumerate(['Sleepiness','Fatigue']):
        df_eval.loc[df_eval["sub"]==iSub,immName2+"_norm"] = scipy.stats.zscore(df_eval[df_eval["sub"]==iSub][immName2])

tmp_Sleepiness = df_eval.copy()
tmp_Sleepiness["ans"] = tmp_Sleepiness["Sleepiness"]
tmp_Sleepiness["ans_norm"] = tmp_Sleepiness["Sleepiness_norm"]
tmp_Sleepiness['ans_norm'].fillna(0, inplace = True)

tmp_Sleepiness["Task"] = "Sleepiness"
del tmp_Sleepiness["Sleepiness_norm"],tmp_Sleepiness["Fatigue_norm"]

tmp_tired = df_eval.copy()
tmp_tired["ans"] = tmp_tired["Fatigue"]
tmp_tired["ans_norm"] = tmp_tired["Fatigue_norm"]

tmp_tired['ans_norm'].fillna(0, inplace = True)

tmp_tired["Task"] = "Fatigue"
del tmp_tired["Sleepiness_norm"],tmp_tired["Fatigue_norm"]

df_eval_all = pd.concat([tmp_Sleepiness,tmp_tired])

# df_eval_all = df_eval_all.dropna(how="any")

tmp_plot = df_eval_all.groupby(["sub","Run","light","Task"],as_index=False).agg('mean', numeric_only=True)

tmp_plot.to_json(folderName+"/questionaire.json")

# %%

tmp_plot = df.groupby(["sub","hourCls","orderNum","light"],as_index=False).agg('mean', numeric_only=True)


