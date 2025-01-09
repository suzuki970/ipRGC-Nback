#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 10:08:48 2024

@author: yutasuzuki
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import json
import glob
import os
import pandas as pd
import seaborn as sns        
import datetime
import scipy.stats
# from statsmodels.stats.anova import AnovaRM
  
from PupilAnalysisToolbox.lib.pre_processing_cls import rejectedByOutlier,getNearestValue
# from rejectBlink_PCA import rejectBlink_PCA
from PupilAnalysisToolbox.lib.LightSource import LightSource,showSpectra

# from pixel_size import pixel_size

sns.set(font_scale=1.5)
sns.set_style("whitegrid")
sns.set_palette("Set2")


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


dt_now = datetime.datetime.now()
date = dt_now.strftime("%Y%m%d")

if not os.path.exists("./figure/"+date):
    os.makedirs("./figure/"+date,exist_ok=True)

ext = ".pdf"

folderName = glob.glob("./data/*")
folderName.sort()
folderName = folderName[-1]

rootFolder = "./results/monitorTest/"

#%% data loading

cfg = {
    'numOfLEDs': [1,2,3,4,5,6],
    'winlambda' :[380,780],
    # 'Y'    :3,
    'Y'    :3,
    'maxoutPw' : 255, # 100%
    'rod' : np.linspace(0,0.01,10, endpoint=True),
    # 'ipRGC' : np.linspace(0,1000,1000, endpoint=True),
    'ipRGC' : np.linspace(1,1000,2000, endpoint=True),
    "visualization":False,
    "VISUAL_ANGLE":2,
    "projector":True,
    # "plus":True
    "plus":False
    }

cfg["outPw"] = 120
cfg["x"] = 0.40053
cfg["y"] = 0.25057

LEDCube = LightSource(cfg)

ipRGC = LEDCube.getipRGCfunc()

fileName_ipRGC = glob.glob(rootFolder+"light.csv")
fileName = glob.glob(rootFolder+"control.csv")

# csv_input = pd.read_csv(filepath_or_buffer = fileName[0],  sep=",", header=None)
csv_input = pd.read_csv(filepath_or_buffer = fileName[0],  sep=",", header=0)

# csv_input = csv_input.drop(12,axis=0)
csv_input = csv_input.reset_index(drop=True)

csv_input_original = csv_input.copy()

csv_input = csv_input.drop(csv_input.columns[np.arange(18).tolist()], axis=1)
csv_input = csv_input.drop(csv_input.columns[np.arange(-25,-1).tolist()], axis=1)
csv_input = csv_input.drop(csv_input.columns[[401]], axis=1)

dr1 = np.arange(380,781)
# np.array([int(l[:3]) for l in csv_input.iloc[0,:].values])

# csv_input = np.array(csv_input.iloc[1:,:].values, dtype='float')
# csv_input[abs(csv_input)<(10**(-4)*2)]=0
# csv_input[csv_input<0]=0

df = []
for iLight in np.arange(csv_input.shape[0]):
    tmp_df = pd.DataFrame()
    tmp_df["wavelength[nm]"] = dr1 
    tmp_df["power"] = csv_input.iloc[iLight,:].values
    tmp_df["ipRGC"] = csv_input.iloc[iLight,:].values*ipRGC["ipRGC"]
    
    # if iLight == csv_input.shape[0]-1:
    #     tmp_df["Light"] = "high ipRGC"
    # else:
    tmp_df["Light"] = "low ipRGC{:02}".format(iLight+1)
    
    for iLoc in ['x','y','X','Y','Z','L*','a*','b*']:
        tmp_df[iLoc] = csv_input_original[iLoc][iLight]

    df.append(tmp_df)
    
csv_input = pd.read_csv(filepath_or_buffer = fileName_ipRGC[0],  sep=",", header=0)
csv_input = csv_input.reset_index(drop=True)
csv_input_original = csv_input.copy()

csv_input = csv_input.drop(csv_input.columns[np.arange(18).tolist()], axis=1)
csv_input = csv_input.drop(csv_input.columns[np.arange(-25,-1).tolist()], axis=1)
csv_input = csv_input.drop(csv_input.columns[[401]], axis=1)


iLight=13
tmp_df = pd.DataFrame()
tmp_df["wavelength[nm]"] = dr1 
tmp_df["power"] = csv_input.iloc[iLight,:].values
tmp_df["ipRGC"] = csv_input.iloc[iLight,:].values*ipRGC["ipRGC"]

tmp_df["Light"] = "high ipRGC"

for iLoc in ['x','y','X','Y','Z','L*','a*','b*']:
    tmp_df[iLoc] = csv_input_original[iLoc][iLight]

df.append(tmp_df)
    
df = pd.concat(df)

df_ipRGC = df.groupby(["Light"],as_index=False).agg('sum', numeric_only=True)

#%% plot CIE xy 

plt.figure()
# g = sns.FacetGrid(df,hue="Light")
# g.map(sns.lineplot, "wavelength[nm]", "power", errorbar='se',alpha=0.8)

g=sns.lineplot(x="wavelength[nm]", y="power",hue= "Light",data=df, errorbar='sd')
sns.move_legend(g,'upper right',title = '',bbox_to_anchor=(1.5, 1),borderaxespad=0,frameon=False)

df=df.reset_index()
df.to_json(folderName+"/spectra_all.json")

df_ave = df.groupby(["Light"],as_index=False).agg('mean', numeric_only=True)
df_ave["ipRGC"] = df_ipRGC["ipRGC"]

plt.figure()
# g=sns.scatterplot(x='x',y='y',hue="Light",palette="Set1",data=df_ave)
g=sns.scatterplot(x='x',y='y',style="Light",hue="Light",palette="Set1",data=df_ave,s=100)
# plt.legend(g,loc='upper right')
# sns.move_legend(g, "center right", bbox_to_anchor=(1, 1))
sns.move_legend(g,'upper right',title = '',bbox_to_anchor=(1.3, 1),borderaxespad=0,frameon=False)

#%% data loading

folderList = glob.glob("./results/Brightness/*")
folderList.sort()

folderList = [l for l in folderList if not "s24" in l]

df_brightness = []
for iSub,subName in enumerate(folderList):
   
    #% ------------------ load json file (eye data) -------------------
    print("Processing --> " + subName[-3:] + "...")

    matFileName = glob.glob(subName + "/*_param.json")

    for iRun in np.arange(len(matFileName)):

        f = open(matFileName[iRun])
        param = json.load(f)
        f.close()

    tmp_df = pd.DataFrame()
    tmp_df["ans"] =  [c["selected"] for c in param["coeff"]]
    tmp_df["sub"] = int(subName[-2:])
    # tmp_df["numOfLight"] = int(folderList[iSub][1:3])
   
    df_brightness.append(tmp_df)
    
df_brightness = pd.concat(df_brightness)

df_brightness = df_brightness.reset_index(drop=True)

rejectNum = rejectedByOutlier(df_brightness,[],"ans")

df_brightness = df_brightness.drop(df_brightness.index[rejectNum])

df_brightness = df_brightness.groupby(["sub"],as_index=False).agg(np.nanmean)

#%% data loading 

df_behavior = []
for ff in ["./results/stimTuning/"]:
   
    folderList=[]
    for filename in os.listdir(ff):
        if os.path.isdir(os.path.join(ff, filename)): 
            folderList.append(filename)
    folderList.sort()
    
    folderList = [l for l in folderList if not "s24" in l]


    for iSub,subName in enumerate(folderList):
    
        ################ load json file (eye data) ################
        print("Processing --> " + subName[-3:] + "...")
       
        f_cfg = glob.glob(ff + "/" + subName + "/*[!param].json")
        f_cfg.sort()
        f_param = glob.glob(ff + "/" + subName + "/*_param.json")
        f_param.sort()
       
        for iRun in np.arange(len(f_cfg)):
            f = open(f_param[iRun])
            param = json.load(f)
            f.close()
        
            f = open(f_cfg[iRun])
            cfg = json.load(f)
            f.close()
            
            tmp_df = pd.DataFrame()
            tmp_df["Color"] = [c["Color"] for c in cfg["condition_frame"]]
            tmp_df["Locs"] = [c["Locs"] for c in cfg["condition_frame"]]
        
            tmp_df["ans"] = param["ans"]
            
            # tmp_df["sub"] = 's'+str(int(folderList[iSub][1:3]))
            tmp_df["sub"] = int(folderList[iSub][1:3])
            
            
            # if len(np.argwhere(tmp_df["Color"].values == 0)) > 0:
            #     # x = np.r_[spectra_plus["Yxy"][0][1], np.array(spectra_minus["Yxy"])[:,1]]
            #     # y = np.r_[spectra_plus["Yxy"][0][2], np.array(spectra_minus["Yxy"])[:,2]]
            #     x = np.r_[spectra["iprgc+"]["Yxy"][0][1], np.array(spectra["control1_0"]["Yxy"])[:,1]]
            #     y = np.r_[spectra["iprgc+"]["Yxy"][0][2], np.array(spectra["control1_0"]["Yxy"])[:,2]]
        
            #     tmp_df["x"] = np.round(x[tmp_df["Color"]],4)
            #     tmp_df["y"] = np.round(y[tmp_df["Color"]],4 )
    
            # else:
            #     # x = np.array(spectra_minus["Yxy"])[:,1]
            #     # y = np.array(spectra_minus["Yxy"])[:,2]
            #     x = np.array(spectra["control1_0"]["Yxy"])[:,1]
            #     y = np.array(spectra["control1_0"]["Yxy"])[:,2]
       
            #     tmp_df["x"] = np.round(x[tmp_df["Color"]-1],4)
            #     tmp_df["y"] = np.round(y[tmp_df["Color"]-1],4 )
        
            df_behavior.append(tmp_df)

df_behavior = pd.concat(df_behavior)

df_behavior["adjusted"]=0
df_behavior["selected"]=0
df_behavior["Y"]=0

for ff in ["./results/Brightness/"]:
   
    folderList=[]   
    for filename in os.listdir(ff):
        if os.path.isdir(os.path.join(ff, filename)): 
            folderList.append(filename)
    folderList.sort()
    
    folderList = [l for l in folderList if not "s24" in l]

    for iSub,subName in enumerate(folderList):
       
        ################ load json file (eye data) ################       
        f_cfg = glob.glob(ff + "/" + subName + "/*[!param].json")
        f_cfg.sort()
        f_param = glob.glob(ff + "/" + subName + "/*_param.json")
        f_param.sort()
       
        for iRun in np.arange(len(f_cfg)):
            f = open(f_param[iRun])
            param = json.load(f)
            f.close()
        
            f = open(f_cfg[iRun])
            cfg = json.load(f)
            f.close()
            
        df_behavior.loc[df_behavior["sub"]== int(folderList[iSub][1:3]),"adjusted"] = param["coeff_adjusted"]
        df_behavior.loc[df_behavior["sub"]== int(folderList[iSub][1:3]),"selected"] = param["selected"]
        # df_behavior.loc[df_behavior["sub"]== int(folderList[iSub][1:3]),"Y"] = df_ave.iloc[param["selected"]]["Y"]


df_behavior_same = df_behavior[df_behavior["Color"] == 0]
df_behavior = df_behavior[df_behavior["Color"] != 0]

df_behavior_ave = df_behavior.groupby(["sub","Color"],as_index=False).agg(np.nanmean)
df_behavior_ave["norm"]=0

df_max = []
for iSub in np.unique(df_behavior_ave["sub"]):
    
    tmp = df_behavior_ave[df_behavior_ave["sub"]==iSub]
    tmp.loc[tmp["sub"]==iSub,"norm"] = np.array(scipy.stats.zscore(tmp["ans"]))
    df_behavior_ave.loc[df_behavior_ave["sub"]==iSub,"norm"] = np.array(scipy.stats.zscore(tmp["ans"]))
    
    tmp_sub = tmp[(tmp["ans"]==tmp.max()["ans"])]
    
    df_max.append(tmp_sub.iloc[0:1,:])
    
df_max = pd.concat(df_max)

df_max["Y"] = df_ave["Y"][df_max["Color"].values].values
df_max.to_json(folderName+"/df_max.json")

#%% probablity of yes at the same condition

tmp_plot = df_behavior_same.groupby(["sub"],as_index=False).agg(np.mean)

plt.figure()
g=sns.boxplot(data=tmp_plot, x="Color", y="ans")
# sns.move_legend(g,'upper right',title = '',
#                 bbox_to_anchor=(1.2, 1),
#                 borderaxespad=0,frameon=False)

tmp_plot.to_json(folderName+"/stimSelect_same.json")


#%%
# plt.figure()
# sns.pointplot(data=df_behavior_ave, x="Color", y="norm", ci='sd')

plt.figure()
# grid = sns.FacetGrid(df_ave, col="sub",col_wrap=5)
# grid.map(sns.pointplot, "Color", "ans", color=np.array([0.5,0.5,0.5]),fill=True,alpha=0.8)
# g=sns.pointplot(data=df_behavior_ave, x="Color", y="norm", hue="sub", ci='sd')
g=sns.pointplot(data=df_behavior_ave, x="Color", y="ans", hue="sub", ci='sd')
# grid.set(xlabel='Color', ylabel='Prob. of the "same"')
sns.move_legend(g,'upper right',title = '',bbox_to_anchor=(1.2, 1),borderaxespad=0,frameon=False)

# df_max = df_ave.groupby(["sub"],as_index=False).agg(np.max)

for iSub in df_max["sub"]:
    g.scatter(df_max[df_max["sub"]==iSub]["Color"].values-1, df_max[df_max["sub"]==iSub]["ans"].values, marker="o", color="black", s=100)

plt.savefig("./figure/"+date+"/stimulus_Tuning"+ext,bbox_inches="tight",pad_inches=0.2)

df_behavior_ave.to_json(folderName+"/df_behavior_ave.json")


    # g.scatter(df_max[df_max["sub"]==iSub]["Color"].values-1, df_max[df_max["sub"]==iSub]["norm"].values, marker="o", color="black", s=120)

    # if df_behavior_ave[(df_behavior_ave["sub"]==iSub)&
    #           (df_behavior_ave["Color"]==col)]["norm"].values == df_max[df_max["sub"]==iSub]["norm"].values:
        # df_max.loc[df_max["sub"]==iSub,"Color"] = col
    
        # g.scatter(col-1, df_max[df_max["sub"]==iSub]["norm"].values, marker="o", color="black", s=120)

    # if df_behavior_ave[(df_behavior_ave["sub"]==iSub)&
    #           (df_behavior_ave["Color"]==col)]["ans"].values == df_max[df_max["sub"]==iSub]["ans"].values:
        
        # df_max.loc[df_max["sub"]==iSub,"Color"] = col
    
    

#%%

df_ave["size"]=0
for iLight in np.arange(len(df_ave)-1)+1:
    if len(df_max[df_max["Color"]==iLight])!=0:
        df_ave.loc[iLight,"size"] = len(df_max[df_max["Color"]==iLight])

plt.figure(figsize=(6,6))
g=sns.scatterplot(x='x',y='y',
                  # style="Light",
                  # hue="Light",
                  hue="size",
                  palette="OrRd",
                  size="size",
                  sizes=(50, 400),
                  data=df_ave)
g.set_aspect('equal')
# g.scatterplot(x=df_ave["x"],y=df_ave["y"],color="black")
g.scatter(df_ave[df_ave["size"]==0]["x"].values, 
          df_ave[df_ave["size"]==0]["y"].values, marker="o", color="gray", s=30)
plt.plot(df_ave[df_ave["Light"]=="high ipRGC"]["x"],
         df_ave[df_ave["Light"]=="high ipRGC"]["y"],'ro')

sns.move_legend(g,'upper right',title = '',
                bbox_to_anchor=(1.2, 1),borderaxespad=0,frameon=False)

plt.savefig("./figure/"+date+"/xy_bubble"+ext,bbox_inches="tight",pad_inches=0.2)

df_ave.to_json(folderName+"/spectra.json")


#%% data loading

l = [0,52,104]
rgb = ["R","G","B"]
df_spectra=[]

for i,iProjector in enumerate(["blueFilt_","yellowFilt_"]):
    # aa
    fileName = glob.glob("./PupilAnalysisToolbox/lib/LEDdata/projector/Spectrum/"+iProjector+"R.csv")

    csv_input = pd.read_csv(filepath_or_buffer = fileName[0],  sep=",", header=0)
    
    x = csv_input["x"]
    y = csv_input["y"]
    
    csv_input = csv_input.drop(csv_input.columns[np.arange(18).tolist()], axis=1)
    csv_input = csv_input.drop(csv_input.columns[np.arange(-25,-1).tolist()], axis=1)
    csv_input = csv_input.drop(csv_input.columns[[401]], axis=1)

    # plt.plot(csv_input.T)
    
    for irgb in np.arange(3):
        tmp_df = pd.DataFrame()
        tmp_df["wavelength[nm]"] = dr1 
        tmp_df["power"] = csv_input.iloc[l[irgb],:].values
        tmp_df["x"] = round(x[l[irgb]],3)
        tmp_df["y"] = round(y[l[irgb]],3)
   
        tmp_df["Light"] = iProjector+rgb[irgb]
        
        df_spectra.append(tmp_df)

df_spectra = pd.concat(df_spectra)

plt.figure()
g=sns.lineplot(x='wavelength[nm]',y='power',hue="Light",palette="Set1",data=df_spectra)
sns.move_legend(g,'upper right',title = '',bbox_to_anchor=(1.3, 1),borderaxespad=0,frameon=False)

plt.savefig("./figure/"+date+"/spectra"+ext,bbox_inches="tight",pad_inches=0.2)

df_spectra = df_spectra.reset_index(drop=True)

df_spectra.to_json(folderName+"/spectra_projector.json")

df_spectra = df_spectra.groupby(["Light"],as_index=False).agg('mean', numeric_only=True)

#%%


ipRGC = LEDCube.getipRGCfunc()

fileName = glob.glob(rootFolder+"lightFlux.csv")

csv_input = pd.read_csv(filepath_or_buffer = fileName[0],  sep=",", header=0)

# csv_input = csv_input.drop(12,axis=0)
# csv_input = csv_input.reset_index(drop=True)

csv_input_original = csv_input.copy()

csv_input = csv_input.drop(csv_input.columns[np.arange(18).tolist()], axis=1)
csv_input = csv_input.drop(csv_input.columns[np.arange(-25,-1).tolist()], axis=1)
csv_input = csv_input.drop(csv_input.columns[[401]], axis=1)

dr1 = np.arange(380,781)
# np.array([int(l[:3]) for l in csv_input.iloc[0,:].values])

# csv_input = np.array(csv_input.iloc[1:,:].values, dtype='float')
# csv_input[abs(csv_input)<(10**(-4)*2)]=0
# csv_input[csv_input<0]=0

df = []
for iLight in np.arange(csv_input.shape[0]):
    tmp_df = pd.DataFrame()
    tmp_df["wavelength[nm]"] = dr1 
    tmp_df["power"] = csv_input.iloc[iLight,:].values
    tmp_df["ipRGC"] = csv_input.iloc[iLight,:].values*ipRGC["ipRGC"]
    tmp_df["sub"] = iLight+1
    
    if iLight == csv_input.shape[0]-1:
        tmp_df["Light"] = "high ipRGC"
    else:
        tmp_df["Light"] = "low ipRGC{:02}".format(iLight+1)
    
    for iLoc in ['x','y','X','Y','Z','L*','a*','b*']:
        tmp_df[iLoc] = csv_input_original[iLoc][iLight]

    df.append(tmp_df)
    
df = pd.concat(df)

df = df.reset_index(drop=True)

df.to_json(folderName+"/spectra_flux.json")

df_ave = df.groupby(["sub"],as_index=False).agg('mean', numeric_only=True)

# A=69.7
# a=0.88
# B=28
# b=0.922
# C=-14.9
# c=1.104
# D=0.674

# x = np.arange(390,781)
# x_max=485
# x = x_max/x

# sx = 1/(np.exp(A*(a-x))+np.exp(B*(b-x))+np.exp(C*(c-x))+D)


#%%

csv_input = pd.read_csv(filepath_or_buffer = rootFolder+"rgb50.csv",  sep=",", header=0)

x = csv_input["x"]
y = csv_input["y"]

csv_input = csv_input.drop(csv_input.columns[np.arange(18).tolist()], axis=1)
csv_input = csv_input.drop(csv_input.columns[np.arange(-25,-1).tolist()], axis=1)
csv_input = csv_input.drop(csv_input.columns[[401]], axis=1)

np.round(np.array([0.353668,0.391385]),3)

# f = open("/Users/yutasuzuki/GoogleDrive/PupilAnalysisToolbox/python/preprocessing/lib/LEDdata/photoreceptors.json")
# s = json.load(f)
# f.close()

# plt.plot(np.arange(380,781),s["energyFundamentals"])
# plt.plot(np.arange(380,781),ipRGC["ipRGC"])

