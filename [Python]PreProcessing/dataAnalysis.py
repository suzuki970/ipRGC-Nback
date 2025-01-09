#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 08:43:51 2022

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
from statsmodels.stats.anova import AnovaRM

import scipy as sp
import scipy

from PupilAnalysisToolbox.lib.pre_processing_cls import pre_processing,rejectDat,getNearestValue,SDT
from PupilAnalysisToolbox.lib.pre_processing_cls import rejectedByOutlier
from PupilAnalysisToolbox.lib.rejectBlink_PCA import rejectBlink_PCA

sns.set(font_scale=2)
sns.set_style("whitegrid")
sns.set_palette("Set2")

def annotate(data, **kwargs):
    r, p = sp.stats.pearsonr(data[kwargs["var1"]], data[kwargs["var2"]])
    ax = plt.gca()
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
            transform=ax.transAxes)

dt_now = datetime.datetime.now()
date = dt_now.strftime("%Y%m%d")

if not os.path.exists("./figure/"+date):
    os.mkdir("./figure/"+date)
    
ext='.png'

mmName1 = ["low ipRGC","high ipRGC"]
mmName2 = ["1-back","2-back"]

#%% ------------------ data loading  -----------------------------------------   

cfg={
    "windowL":[],
    "WID_ANALYSIS":1.5,
    "WID_BASELINE":[[-0.2,0]],
    "WID_FILTER":np.array([]),
    "METHOD":1, #subtraction
    # "METHOD":2, #subtraction
    "FLAG_LOWPASS":False,
    "visualization":False,
    # "THRES_DIFF":0.03
    "THRES_DIFF":0.05
    }

folderName = glob.glob("./data/*")
folderName.sort()
folderName = folderName[-1]

log = glob.glob(folderName+"/s*_trial.json")
log.sort()

datHash = {}
zeroArrays = []
for subName in log:
    f = open(subName)
    tmp_datHash = json.load(f)
    f.close()

    for mmName in list(tmp_datHash.keys()):
        if mmName != "cfg":
            if not mmName in list(datHash.keys()):
                datHash[mmName] = []
            datHash[mmName] = datHash[mmName] + tmp_datHash[mmName]
        else:
            zeroArrays.append(tmp_datHash["cfg"]["numOfZero_th"])

cfg.update(tmp_datHash["cfg"])

# print("ave="+str(round(np.mean(zeroArrays),3))+",SD="+str(round(np.std(zeroArrays),4)))


#%% plot time series

reject = {}

pp = pre_processing(cfg)
  
y1_targeted,reject["PDR"] = pp.pre_processing(np.array(datHash["PDR"]).copy())

x = np.linspace(cfg["TIME_START"],cfg["TIME_END"],y1_targeted.shape[1])

# pca,reject["PCA"]= rejectBlink_PCA(y1_targeted[:,getNearestValue(x,cfg["WID_BASELINE"][0][0]):getNearestValue(x,cfg["WID_ANALYSIS"])])

reject = reject["PDR"]
reject = list(set(reject))
reject.sort()

datHash, y1_targeted = rejectDat(datHash, reject, y1_targeted)

y_bp = np.array(datHash["PDR"])[:,getNearestValue(x,cfg["WID_BASELINE"][0][0]):getNearestValue(x,0)]
y_bp = np.nanmean(y_bp,axis=1)

y_bp_task = np.array(datHash["PDR"])[:,getNearestValue(x,0):]
y_bp_task = np.nanmean(y_bp_task,axis=1)


#%% make df time course

df_timeCourse = pd.DataFrame()
df_timeCourse["Pupil[mm]"] = y1_targeted.reshape(-1)
df_timeCourse["Time[s]"] = np.tile(x, y1_targeted.shape[0])
df_timeCourse["ipRGC"] = np.tile(np.array(datHash["light"]).reshape(np.array(datHash["Run"]).shape[0],1), y1_targeted.shape[1]).reshape(-1)
df_timeCourse["sub"] = np.tile( np.array(datHash["sub"]).reshape(np.array(datHash["sub"]).shape[0],1), y1_targeted.shape[1]).reshape(-1)
df_timeCourse["Run"] = np.tile( np.array(datHash["Run"]).reshape(np.array(datHash["Run"]).shape[0],1), y1_targeted.shape[1]).reshape(-1)

for c in cfg["conditionName"][1:]:
    df_timeCourse[c] = np.tile( np.array(datHash[c]).reshape(np.array(datHash[c]).shape[0],1), y1_targeted.shape[1]).reshape(-1)


df_timeCourse["light"] = [mmName1[int(i-1)] for i in df_timeCourse["light"].values]
df_timeCourse["Nback"] = [mmName2[int(i-1)] for i in df_timeCourse["Nback"].values]

df_timeCourse_plot = df_timeCourse.groupby(["sub","Time[s]","light","Nback"],as_index=False).agg('mean', numeric_only=True)
df_timeCourse_plot = df_timeCourse_plot[(df_timeCourse_plot["Time[s]"]>=cfg["WID_BASELINE"][0][0])&
                                        (df_timeCourse_plot["Time[s]"]<=1.5)]

# grid = sns.FacetGrid(df_timeCourse_plot, col="light", hue="Nback")
grid = sns.FacetGrid(df_timeCourse_plot, col="Nback", hue="light")
grid.map(sns.lineplot, "Time[s]", "Pupil[mm]", errorbar='sd')

grid.fig.set_figwidth(10)
grid.fig.set_figheight(6)

plt.legend(bbox_to_anchor=(1.6, 1), loc='upper right', borderaxespad=0)

df_ave = df_timeCourse_plot.groupby(["sub","light","Nback"],as_index=False).agg('mean', numeric_only=True)

plt.figure()
sns.barplot(x="Nback", y="Pupil[mm]",hue="light", data=df_ave)
plt.legend(bbox_to_anchor=(1.5, 1), loc='upper right', borderaxespad=0)

#%%

tmp_plot = df_timeCourse_plot.groupby(["sub","Time[s]","light","Nback"],as_index=False).agg('mean', numeric_only=True)

tmp_plot = tmp_plot.dropna(how="any")
tmp_plot.to_json(folderName + "/data_timeCourse_each.json")

plt.figure()
g = sns.FacetGrid(tmp_plot,col="Nback", hue="light")
g.map(sns.lineplot, "Time[s]", "Pupil[mm]", errorbar='se',alpha=0.4).add_legend()


#%% SDT

df = pd.DataFrame()
for mmName in ["sub","target","response","RT","order","Run"]:
    df[mmName] =  np.array(datHash[mmName])

df["Pupil[a.u.]"] = np.nanmean(y1_targeted,axis=1)
df["BP"] =  y_bp
df["BP_task"] =  y_bp_task

for c in cfg["conditionName"]:  
    df[c] = np.array(datHash[c])

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
  
    
for mmName in ["hit","FA","CR","miss"]:
    df[mmName] = 0

df.loc[(df["target"]==1) & (df["response"]==1),"hit"]  = 1
df.loc[(df["target"]==0) & (df["response"]==1),"FA"]   = 1
df.loc[(df["target"]==0) & (df["response"]==0),"CR"]   = 1
df.loc[(df["target"]==1) & (df["response"]==0),"miss"] = 1


tmp = df.groupby(["sub",'light','Nback'],as_index=False,sort=False).agg('mean', numeric_only=True)

df_target = df[df["target"]==1].copy()

tmp = df_target.groupby(["sub",'light','Nback'],as_index=False,sort=False).agg('mean', numeric_only=True)

tmp = df_target[df_target["RT"] > 0]

#########################################################
q75 = tmp.groupby(["sub"],as_index=False).quantile(.75, numeric_only=True)
q25 = tmp.groupby(["sub"],as_index=False).quantile(.25, numeric_only=True)

IQR = q75 - q25

lower = q25 - IQR*3
upper = q75 + IQR*3

df["lower"] = 0
df["upper"] = 0
for iSub in np.unique(df["sub"]):
    df.loc[df["sub"]==iSub,"lower"] = float(lower["RT"][lower["sub"]==iSub].values)
    df.loc[df["sub"]==iSub,"upper"] = float(upper["RT"][lower["sub"]==iSub].values)

df = df[(df["RT"] == 0) |
        (df["RT"] > df["lower"]) &
        (df["RT"] < df["upper"])]

#########################################################

df_target = df[df["target"]==1].copy()

df_target["numOftrials[trial]"] = 0
for iSub in np.unique(df_target["sub"]):
    for iNback in np.unique(df_target["Nback"]):
        for iLight in np.unique(df_target["light"]):
            tmp = len(df_target[(df_target["sub"] == iSub) &
                     (df_target["Nback"] == iNback) &
                     (df_target["light"] == iLight)]["numOftrials[trial]"])
            df_target.loc[(df_target["sub"] == iSub) &
                           (df_target["Nback"] == iNback) &
                           (df_target["light"] == iLight),"numOftrials[trial]"] = np.arange(1,tmp+1)

# tmp_df_target = df_target.dropna(subset=['Pupil[a.u.]'])
tmp_df_target = df_target.dropna(subset=['BP'])
tmp_df_target = tmp_df_target[df_target["RT"]!=0]

# sns.lmplot(data=tmp_df_target,
#             x="BP",y="RT",hue="light",col="sub",col_wrap=5)

############## STD ##############

df_sum = df.groupby(["sub","light","Nback"],as_index=False).agg(np.sum)

sdt_data = pd.DataFrame()

t = SDT(df_sum.at[0,"hit"], df_sum.at[0,"miss"], df_sum.at[0,"FA"], df_sum.at[0,"CR"])

for iSTD in list(t.keys()):
    df_sum[iSTD] = 0

for i, item in df_sum.iterrows():
    n1 = item["hit"] + item["miss"]
    n2 = item["FA"] + item["CR"]
    
    # for mmName in ["hit","miss"]:
    #     df_sum.at[i,mmName] = df_sum.at[i,mmName] / n1
        
    # for mmName in ["CR","FA"]:
    #     df_sum.at[i,mmName] = df_sum.at[i,mmName] / n2
    t = SDT(df_sum.at[i,"hit"], df_sum.at[i,"miss"], df_sum.at[i,"FA"], df_sum.at[i,"CR"])
    
    for mmName in list(t.keys()):
        df_sum.at[i,mmName] = t[mmName].values


############## STD each Run ##############

df_sum_run = df.groupby(["sub","Run","light","Nback"],as_index=False).agg(np.sum)
df_sum_run_target = df[df["target"]==1].groupby(["sub","Run","light","Nback"],as_index=False).agg(np.mean)

sdt_data = pd.DataFrame()

t = SDT(df_sum_run.at[0,"hit"], df_sum_run.at[0,"miss"], df_sum_run.at[0,"FA"], df_sum_run.at[0,"CR"])

for iSTD in list(t.keys()):
    df_sum_run_target[iSTD] = 0

for i, item in df_sum_run.iterrows():
    n1 = item["hit"] + item["miss"]
    n2 = item["FA"] + item["CR"]    
    # for mmName in ["hit","miss"]:
    #     df_sum.at[i,mmName] = df_sum.at[i,mmName] / n1
        
    # for mmName in ["CR","FA"]:
    #     df_sum.at[i,mmName] = df_sum.at[i,mmName] / n2
    t = SDT(df_sum_run.at[i,"hit"], df_sum_run.at[i,"miss"], df_sum_run.at[i,"FA"], df_sum_run.at[i,"CR"])
    
    for mmName in list(t.keys()):
        df_sum_run_target.at[i,mmName] = t[mmName].values

#%% reject sub

tmp = df_sum.groupby(["sub"],as_index=False,sort=False).agg('sum', numeric_only=True)

tmp["resRate"] = tmp["response"] / tmp["target"]

sigma = np.std(tmp["resRate"] ) 

lower = np.mean(tmp["resRate"] ) - sigma*3
upper = np.mean(tmp["resRate"] ) + sigma*3

rejectSub = tmp[(tmp["resRate"]<lower)|(tmp["resRate"]>upper)]["sub"]

for iSub in rejectSub:    
    df_target = df_target[df_target["sub"]!=iSub]
    df_sum = df_sum[df_sum["sub"]!=iSub]
    df_sum_run = df_sum_run[df_sum_run["sub"]!=iSub]
    df_sum_run_target = df_sum_run_target[df_sum_run_target["sub"]!=iSub]
    df = df[df["sub"]!=iSub]
    

print(str(rejectSub) + " was rejected because of less number of responses.")

############## Number of ipRGC+ first ##############

tmp = df_sum.groupby(["sub",'light','Nback'],as_index=False,sort=False).agg('mean', numeric_only=True)
tmp = tmp[(tmp["Nback"]=="1-back")&(tmp["light"]=="high ipRGC")]

print("Number of ipRGC+ first = " + str(len(tmp[tmp["orderNum"]==0])) + " out of " + str(len(tmp)))

df_sum_run.to_json(folderName + "/df_run.json")
df_target.to_json(folderName + "/df_targetAll.json")
df_sum_run_target.to_json(folderName + "/df_SDT.json")



#%% RT

# tmp_df = df_target.copy()
tmp_df = df.copy()

# tmp_df["numOfTrialAll"] = df_target["numOftrials[trial]"]

tmp_df = tmp_df[tmp_df["RT"]!=0]

tmp_df2 = tmp_df.copy()
tmp_df2 = tmp_df2.reset_index(drop=True)
tmp_df = tmp_df.reset_index(drop=True)

# t_std = tmp_df2.groupby(["sub",'Nback'],as_index=False,sort=False).agg('std', numeric_only=True)
# t_ave = tmp_df2.groupby(["sub",'Nback'],as_index=False,sort=False).agg('mean', numeric_only=True)

# tmp_df2["lower"]=0
# tmp_df2["upper"]=0

# for iSub in np.unique(tmp_df2["sub"]):
#     for iNback in np.unique(tmp_df2["Nback"]):
        
#         tmp_df2.loc[(tmp_df2["sub"]==iSub)&(tmp_df2["Nback"]==iNback),"lower"] = \
#             t_ave[(t_ave["sub"]==iSub)&(t_ave["Nback"]==iNback)]["RT"].values[0] - t_std[(t_std["sub"]==iSub)&(t_std["Nback"]==iNback)]["RT"].values[0]*3
        
#         tmp_df2.loc[(tmp_df2["sub"]==iSub)&(tmp_df2["Nback"]==iNback),"upper"] = \
#             t_ave[(t_ave["sub"]==iSub)&(t_ave["Nback"]==iNback)]["RT"].values[0] + t_std[(t_std["sub"]==iSub)&(t_std["Nback"]==iNback)]["RT"].values[0]*3

# sd = tmp_df2.groupby(["sub","Nback"],as_index=False).agg("std", numeric_only=True)
# ave = tmp_df2.groupby(["sub","Nback"],as_index=False).agg("mean", numeric_only=True)

# tmp_df2["lower"] = 0
# tmp_df2["upper"] = 0
# for iSub in np.unique(tmp_df2["sub"]):
#     for iNback in np.unique(tmp_df2["Nback"]):
#         tmp_df2.loc[(tmp_df2["sub"]==iSub)&(tmp_df2["Nback"]==iNback),"lower"] = \
#         ave[(ave["sub"]==iSub)&(ave["Nback"]==iNback)]["RT"].values[0] - \
#             sd[(sd["sub"]==iSub)&(sd["Nback"]==iNback)]["RT"].values[0]*3
        
#         tmp_df2.loc[(tmp_df2["sub"]==iSub)&(tmp_df2["Nback"]==iNback),"upper"] = \
#         ave[(ave["sub"]==iSub)&(ave["Nback"]==iNback)]["RT"].values[0] + \
#             sd[(sd["sub"]==iSub)&(sd["Nback"]==iNback)]["RT"].values[0]*3
        
# reject = tmp_df2[(tmp_df2["RT"] < tmp_df2["lower"]) | (tmp_df2["RT"] > tmp_df2["upper"])].index.values

# reject = rejectedByOutlier(tmp_df2,[],"RT",False)
# tmp_df2 = tmp_df2.drop(reject)

tmp_df.to_json(folderName + "/df_RT.json")
tmp_df2.to_json(folderName + "/df_RT2.json")

tmp = tmp_df2.groupby(["sub",'light','Nback'],as_index=False,sort=False).agg('median', numeric_only=True)

for iNback in ["1-back","2-back"]:
    d1 = tmp[(tmp["Nback"]==iNback)&
             (tmp["light"]=="low ipRGC")]["RT"].values
    d2 = tmp[(tmp["Nback"]==iNback)&
             (tmp["light"]=="high ipRGC")]["RT"].values
     
    print(scipy.stats.ttest_rel(d1, d2).pvalue)

#%% hit rate

df["nonTarget"]=1
df.loc[df["target"]==1,"nonTarget"]=0

df_notarget = df[df["target"]==0].copy()
df_notarget_mean = df_notarget.groupby(["sub",'light','Nback'],as_index=False,sort=False).agg('sum', numeric_only=True)

df_notarget_mean["FA"] = df_notarget_mean["FA"] / df_notarget_mean["nonTarget"]
df_notarget_mean["CR"] = df_notarget_mean["CR"] / df_notarget_mean["nonTarget"]
df_notarget_mean.to_json(folderName + "/df_nonTarget.json")


df_notarget_mean = df_notarget.groupby(["sub",'light','Nback'],as_index=False).agg('mean', numeric_only=True)

for iSTD in ["FA","CR"]:
    plt.figure()
    g = sns.FacetGrid(df_notarget_mean, 
                      # x="Nback",
                     # y=iSTD,
                     col="Nback",
                     # col_wrap=5,
                     # hue="light"
                     )
    g.map(sns.pointplot,
             "light",iSTD,
             errorbar="se").add_legend()
    
    sns.move_legend(g,'upper right',title = '',
                bbox_to_anchor=(0.6, 0.9),borderaxespad=0,frameon=False)

    g.fig.set_figheight(7)
    g.fig.set_figwidth(7)
    g.set(xlabel="",ylabel="Accuracy [% of hit]")
    g.set_titles("{col_name}")
    
#%%

df_notarget_mean = df_notarget.groupby(["sub",'light','Nback'],as_index=False).agg('mean', numeric_only=True)

df_target_mean = df_target.groupby(["sub","light","Nback"],as_index=False).agg('mean', numeric_only=True)

df_target_mean.to_json(folderName + "/df_target_mean.json")

for iSTD in ["hit","miss"]:
    plt.figure()

    g = sns.FacetGrid(df_target_mean, 
                      # x="Nback",
                     # y=iSTD,
                     col="Nback",
                     # col_wrap=5,
                     # hue="light"
                     )
    g.map(sns.pointplot,
             "light",iSTD,
             errorbar="sd").add_legend()
    
    sns.move_legend(g,'upper right',title = '',
                bbox_to_anchor=(0.6, 0.9),borderaxespad=0,frameon=False)

    g.fig.set_figheight(7)
    g.fig.set_figwidth(7)
    g.set(xlabel="",ylabel="Accuracy [% of hit]")
    g.set_titles("{col_name}")
        
    # grid.fig.subplots_adjust(top=0.85)
    g.fig.suptitle(iSTD)
    # plt.legend(bbox_to_anchor=(1.8, 1), loc='upper right', borderaxespad=0)
    
    # plt.xlabel("")
    plt.savefig("./figure/"+date+"/hit"+ext,bbox_inches="tight",pad_inches=0.2)
    
    print("###################")
    print("####### "+iSTD+" #######")
    print("###################")
    
    print(AnovaRM(data=df_target_mean, depvar=iSTD,
                  subject='sub', within=['light','Nback']).fit())
    
    for iNback in ["1-back","2-back"]:
        d1 = df_target_mean[(df_target_mean["Nback"]==iNback)&
                            (df_target_mean["light"]=="low ipRGC")][iSTD].values
        d2 = df_target_mean[(df_target_mean["Nback"]==iNback)&
                            (df_target_mean["light"]=="high ipRGC")][iSTD].values
        
        print(scipy.stats.ttest_rel(d1, d2).pvalue)
        
df_target_mean["hit+CR"] = df_target_mean["hit"]/3 + df_notarget_mean["CR"]/7
df_target_mean["miss+FA"] = df_target_mean["miss"]/3 + df_notarget_mean["FA"]/7
    

for iSTD in ["hit+CR","miss+FA"]:
    plt.figure()
    sns.pointplot(data=df_target_mean,
          x="Nback",y=iSTD,
          hue="light",
          errorbar="sd",
          dodge=True)

    print("###################")
    print("####### "+iSTD+" #######")
    print("###################")
    print(AnovaRM(data=df_target_mean, depvar="hit+CR",
                  subject='sub', within=['light','Nback']).fit())
     
    for iNback in ["1-back","2-back"]:
        d1 = df_target_mean[(df_target_mean["Nback"]==iNback)&
                            (df_target_mean["light"]=="low ipRGC")][iSTD].values
        d2 = df_target_mean[(df_target_mean["Nback"]==iNback)&
                            (df_target_mean["light"]=="high ipRGC")][iSTD].values
         
        print(scipy.stats.ttest_rel(d1, d2).pvalue)


#%% time-course of RT

# df_target_plot = df_target.copy()
# df_target_plot[df_target_plot["RT"]==0] = np.nan

# plt.figure()
# grid = sns.FacetGrid(df_target_plot,
#                      col = "Nback",
#                      hue = "light")

# grid.map(sns.lineplot, "numOftrials[trial]", "RT")
# # grid.map(sns.pointplot, "numOftrials[trial]", "RT")

# # , errorbar='sd', dodge=True)

# # sns.pointplot(x="numOftrials[trial]", y="RT", hue = "light", data=df_target_plot, errorbar='sd', dodge=True)
# grid.fig.subplots_adjust(top=0.85)
# # grid.fig.suptitle('RT')
# grid.fig.set_figwidth(10)
# grid.fig.set_figheight(6)

# plt.legend(bbox_to_anchor=(1.8, 1), loc='upper right', borderaxespad=0)

# plt.savefig("./figure/"+date+"/RT_timecourse"+ext,bbox_inches="tight",pad_inches=0.2)

            
#%% throughput

# tmp_rt = df_target.copy()
# tmp_rt = tmp_rt[(tmp_rt["hit"]==1)]

# tmp_rt = tmp_rt.groupby(["sub","light","Nback"],as_index=False).agg('mean', numeric_only=True)

# tmp_hit = df_target.groupby(["sub","light","Nback"],as_index=False).agg('mean', numeric_only=True)

# tmp_hit["RT"] = tmp_rt["RT"]
# tmp_hit["throughput"] = tmp_hit["hit"] * (1/tmp_hit["RT"])*1000

# plt.figure(figsize=(3,8))
# sns.pointplot(data=tmp_hit,
#               x="Nback",y="throughput",
#               hue="light",
#               errorbar="sd", dodge=True)
# grid.fig.subplots_adjust(top=0.85)
# grid.fig.suptitle('throughput')
# plt.legend()

# plt.savefig("./figure/"+date+"/throughput"+ext,bbox_inches="tight",pad_inches=0.2)

# # grid = sns.FacetGrid(df_target_mean, 
# #                       row="sub",
# #                       col = "Nback",
# #                       hue = "Nback")
# # grid.map(sns.barplot, "light", "hit", errorbar='sd')


# print("###################")
# print("#### Throughput ###")
# print("###################")

# for iNback in ["1-back","2-back"]:
#     d1 = tmp_hit[(tmp_hit["Nback"]==iNback)&
#                         (tmp_hit["light"]=="low ipRGC")]["throughput"].values
#     d2 = tmp_hit[(tmp_hit["Nback"]==iNback)&
#                         (tmp_hit["light"]=="high ipRGC")]["throughput"].values
     
#     print(scipy.stats.ttest_rel(d1, d2).pvalue)
    
# print(AnovaRM(data=tmp_hit, depvar='throughput',
#               subject='sub', within=['light','Nback']).fit())

#%% RT


df_target_mean = df_target.copy()
df_target_mean = df_target_mean[(df_target_mean["hit"]==1)]

# df_target_mean = df_target_mean.groupby(["sub","light","Nback"],as_index=False).agg('mean', numeric_only=True)

df_target_mean = df_target_mean.reset_index(drop=True)

reject = rejectedByOutlier(df_target_mean,[],"RT")

df_target_mean = df_target_mean.drop(reject)

plt.figure(figsize=(3,8))
sns.pointplot(data=df_target_mean,
              x="Nback",
              y="RT",
              hue="light",errorbar="sd", dodge=True)


# sns.pointplot(data=df_target_mean,
#               x="Nback",y="RT",
#               hue="light")

grid.fig.subplots_adjust(top=0.85)
grid.fig.suptitle('RT')
plt.legend()
plt.savefig("./figure/"+date+"/RT"+ext,bbox_inches="tight",pad_inches=0.2)


print("###################")
print("####### RT ########")
print("###################")

# print(AnovaRM(data=df_target_mean, depvar='RT',
#               subject='sub', within=['light','Nback']).fit())

# for mmName in ["1-back","2-back"]:
#     d1 = df_target_mean[(df_target_mean["Nback"]==mmName)&
#                         (df_target_mean["light"]=="low ipRGC")]["RT"].values
#     d2 = df_target_mean[(df_target_mean["Nback"]==mmName)&
#                         (df_target_mean["light"]=="high ipRGC")]["RT"].values
    
#     print(scipy.stats.ttest_rel(d1, d2).pvalue)

