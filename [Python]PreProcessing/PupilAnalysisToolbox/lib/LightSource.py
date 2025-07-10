#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:19:11 2021

@author: yutasuzuki
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import os
import glob
import json
import seaborn as sns        
import warnings
import datetime
from scipy import interpolate
import pickle
from tqdm import tqdm
import matplotlib.patches as patches
import math
from scipy.optimize import minimize

from pre_processing import getNearestValue

# from scipy.interpolate import griddata
# from scipy.optimize import fmin_bfgs
# from sklearn.datasets import make_classification
# from sklearn.model_selection import train_test_split
# import itertools
# from pre_processing_cls import getNearestValue

warnings.simplefilter('ignore')

dt_now = datetime.datetime.now()
date = dt_now.strftime("%Y%m%d")

#%%
def f_optLEDs(x,iLED,lamb,base,LEDs,norm=False):    
   
    out = {'LED1':[],
           'LED2':[],
           'LED3':[],
           'LED4':[]}
  
    for i,iiLED in enumerate(iLED):
        for wavL,_ in enumerate(lamb):
            out["LED"+str(i+1)].append(np.float(LEDs[iiLED-1][wavL](x[i])))
    
    for i,iiLED in enumerate(iLED):
        out["LED"+str(i+1)] = np.array(out["LED"+str(i+1)]).max()

    if norm:
        tmp = []
        for i,b in enumerate(base):
            tmp.append(out["LED"+str(i+1)])
        
        m = np.array(tmp).max()
        tmp = np.array(tmp) / m
        
        res = 0
        for i,(b,t) in enumerate(zip(base,tmp)):
            res = res + abs(b - t)
    
    else:
        res = 0
        for i,b in enumerate(base):
            res = res + abs(b - out["LED"+str(i+1)])
        # y0 = out["LED1"].max()+out["LED2"]+out["LED3"]+out["LED4"]
        # y0 = y0 / max(y0)
    
        # return np.sum(abs(base-y0))
    # plt.figure()
    # plt.plot(np.arange(4),base,color="black")
    # plt.plot(np.arange(0),res,color="black")
    
    # plt.plot(lamb,out["LED1"]+out["LED2"]+out["LED3"]+out["LED4"],color="red")
    # print(res)
    return res
                  

    # return np.sum(abs((out["LED1"]+out["LED2"]+out["LED3"]+out["LED4"])-base))
   
    
#%%
xyzTolms = [[0.3897,0.6890,0.0787],
            [0.2298,1.1834,0.0464],
            [0.0000,0.0000,1.0000]]

xyzToRGB = [[3.2410, -1.5374,-0.4986],
            [-0.9692, 1.8760, 0.0416],
             [0.0556, -0.2040, 1.0507]] 

XYZ_d50 = {'X':0.9642,
           'Y':1,
           'Z':0.8251}

winWaveLen = [380, 780]

class LightSource:    
    def __init__(self, cfg):
        
        """
        
        Input:    
        cfg                - dict of parameters for analysis
        - "numOfLEDs"      - list of used number of LEDs (LED cube)
        - "winlambda"      - list of wave length winfow
        - "maxoutPw"       - int of maximum output values (1000 for LED cube)
        - "rod"            - int or list of rod candidate (unused)
        - "ipRGC"          - int or list of ipRGC candidate
        - "visualization"  - boolean for figs
        - "VISUAL_ANGLE"   - int of visual angle for model (2 or 10[degrees])
        
        - "x"              - float of x values in Yxy color space for target light 
        - "y"              - float of y values in Yxy color space for target light 
        - "Y"              - float of Y values in Yxy color space for target light 
        
        - ourPw            - 
       
        Exapmple:
        cfg = {
                "numOfLEDs": [1,2,3,4,5,6,7,9,10],
                "winlambda" :[380,780],
                "maxoutPw" : 1000, # 100%
                "rod" : np.linspace(0,0.01,10, endpoint=True),
                "ipRGC" : np.arange(0,1000,10),
                "visualization":True,
                "VISUAL_ANGLE":10
               }

    
        cfg["x"]  = 0.30619
        cfg["y"] = 0.25572
        cfg["Y"] = 3
        
        cfg["outPw"] = 500

        LEDCube = LightSource(cfg)
        res = LEDCube.seekCombinations()
        
        
        """

        self.cfg = cfg
        
        # self.cfg['ipRGC_outPw'] = np.array(self.cfg['ipRGC'])*self.cfg['outPw']
        self.cfg['ipRGC_outPw'] = np.array(self.cfg['ipRGC'])
        
        # self.cfg['rod_outPw'] = np.array(self.cfg['rod'])*self.cfg['outPw']
        self.cfg['rod_outPw'] = np.array(self.cfg['rod'])
        
        self.XYZ_d65 = {'X':0.95047,
                        'Y':1,
                        'Z':1.08883}
        
        x = cfg['x']
        y = cfg['y']
        
        Y = cfg['Y'] # depends on the used devices
        
        S = Y/y
        X = S*x
        Z = S-X-Y
     
        self.XYZ_d65['X'] = X
        self.XYZ_d65['Y'] = Y
        self.XYZ_d65['Z'] = Z
          
        self.filePath = glob.glob("PupilAnalysisToolbox/lib/LEDdata/")[0]
        
        self.peakLength = [420,450,475,505,525,540,555,595,610,635,660,705] # LED Cube
 
        # print(str(X) + "," + str(Y) + "," + str(Z))

    def getD65Spectra(self,step=1):
        
        ############ data loading (D65 and color matching function) ############
        Y_d65 = pd.read_csv(filepath_or_buffer= self.filePath + "spectrum_D65_1.csv", sep=",")
        
        tmp_Y_d65 = interpolate.PchipInterpolator(Y_d65["lambda"].values, Y_d65["Y"])
        
        dr = np.arange(Y_d65["lambda"].values[0],Y_d65["lambda"].values[-1]+1,step)
        
        new_Y_d65 = pd.DataFrame()
        new_Y_d65["lambda"] = dr
        new_Y_d65["Y"] = tmp_Y_d65(dr)

        new_Y_d65 = new_Y_d65[(new_Y_d65['lambda'] >= self.cfg['winlambda'][0]) & (new_Y_d65['lambda'] <= self.cfg['winlambda'][1])]
        
        self.Y_d65 = new_Y_d65
        return new_Y_d65
    
    def getRodfunc(self,step=1):
        ############ data loading (D65 and color matching function) ############
        # rod = pd.read_csv(filepath_or_buffer= self.filePath + "Luminous_efficiency_functions_scotopic.csv", sep=",")
        if self.cfg["VISUAL_ANGLE"] == 2:
            # rod = pd.read_csv(filepath_or_buffer = self.filePath + "Luminous_efficiency_functions_photopic.csv", sep=",")
            # rod = pd.read_csv(filepath_or_buffer = "/Users/yutasuzuki/GoogleDrive/PupilAnalysisToolbox/python/preprocessing/lib/LEDdata/linCIE2008v2e_1.csv", sep=",")
            rod = pd.read_csv(filepath_or_buffer = self.filePath + "linCIE2008v2e_1.csv", sep=",")
        else:
            rod = pd.read_csv(filepath_or_buffer = self.filePath + "linCIE2008v10e_1.csv", sep=",")
            
        tmp_rod = interpolate.PchipInterpolator(rod["lambda"].values, rod["sensitivity"].values)
        
        dr = np.arange(rod["lambda"].values[0],rod["lambda"].values[-1]+1,step)
        
        new_rod = pd.DataFrame()
        new_rod["lambda"] = dr
        new_rod["rod"] = tmp_rod(dr)

        new_rod = new_rod[(new_rod['lambda'] >= self.cfg['winlambda'][0]) & (new_rod['lambda'] <= self.cfg['winlambda'][1])]
                
        self.rod = new_rod
        return new_rod
    
    def getXYZfunc(self,step=1):
        
        # cl_func = pd.read_csv(filepath_or_buffer="./data/color-matching_function_wl-xyz.csv", sep=",")
        if self.cfg["VISUAL_ANGLE"] == 2:
            cl_func = pd.read_csv(filepath_or_buffer = self.filePath + "lin2012xyz2e_1_7sf.csv", sep=",")
        else:
            cl_func = pd.read_csv(filepath_or_buffer = self.filePath + "lin2012xyz10e_1_7sf.csv", sep=",")
            
        tmp_cl_func = []
        for mmName in ["X","Y","Z"]:
            tmp_cl_func.append(interpolate.PchipInterpolator(cl_func["lambda"].values, cl_func[mmName]))

        dr = np.arange(cl_func["lambda"].values[0],cl_func["lambda"].values[-1]+1,step)
        # dr = np.arange(winWaveLen[0],winWaveLen[1]+1,step)
        
        
        new_cl_func = pd.DataFrame()
        new_cl_func["lambda"] = dr
        for i,mmName in enumerate(["X","Y","Z"]):
            new_cl_func[mmName] = tmp_cl_func[i](dr)
        
        # cl_func = pd.read_csv(filepath_or_buffer="./data/whitelight/sbrgb10w.csv", sep=",")
        new_cl_func = new_cl_func[(new_cl_func['lambda'] >= self.cfg['winlambda'][0]) & (new_cl_func['lambda'] <= self.cfg['winlambda'][1])]
        self.cl_func = new_cl_func
      
        return  new_cl_func
   
    def getLMSfunc(self,step=1):
        
        if self.cfg["VISUAL_ANGLE"] == 2:
            lms = pd.read_csv(filepath_or_buffer = self.filePath + "linss2_10e_5.csv", encoding="ms932", sep=",")
        else:
            lms = pd.read_csv(filepath_or_buffer = self.filePath + "linss10e_5.csv", encoding="ms932", sep=",")
        
        lms['S'][np.isnan(lms['S'])] = 0
        
        tmp_lms = []
        for mmName in ["L","M","S"]:
            tmp_lms.append(interpolate.PchipInterpolator(lms["lambda"].values, lms[mmName]))

        dr = np.arange(lms["lambda"].values[0],lms["lambda"].values[-1]+1,step)
        
        new_lms = pd.DataFrame()
        new_lms["lambda"] = dr
        for i,mmName in enumerate(["L","M","S"]):
            new_lms[mmName] = tmp_lms[i](dr)

        new_lms = new_lms[(new_lms['lambda'] >= self.cfg['winlambda'][0]) & (new_lms['lambda'] <= self.cfg['winlambda'][1])]
      
        # lms = lms[(lms['lambda'] >= self.cfg['winlambda'][0]) & (lms['lambda'] <= self.cfg['winlambda'][1])]
        self.lms = new_lms
        
        if self.cfg["visualization"]:
            
            tmp = new_lms.melt(id_vars = "lambda",var_name="LEDs",value_name="Radiance")

            plt.figure()
            sns.lineplot(x="lambda", y="Radiance", hue = "LEDs", data=tmp)
            plt.xlabel("Wavelength[nm]")
            plt.ylabel("Illuminant")
            plt.savefig("./LMS.pdf")

            
        return new_lms
    
    def getipRGCfunc(self,step=1):
        
        f = open(self.filePath + "photoreceptors.json")
        s = json.load(f)
        f.close()

        # plt.plot(np.arange(380,781),s["energyFundamentals"])
        # plt.plot(np.arange(380,781),ipRGC["ipRGC"])


        # ipRGC = pd.read_csv(filepath_or_buffer = self.filePath + "ipRGCSpectrum.csv", sep=",")

        # tmp_ipRGC = interpolate.PchipInterpolator(ipRGC["lambda"].values, ipRGC["ipRGC"])
        
        # dr = np.arange(ipRGC["lambda"].values[0],ipRGC["lambda"].values[-1]+1,step)
        dr = np.arange(s["nomogram"]["S"][0],s["nomogram"]["S"][0]+s["nomogram"]["S"][2])
        
        new_ipRGC = pd.DataFrame()
        new_ipRGC["lambda"] = dr
        # new_ipRGC["ipRGC"] = tmp_ipRGC(dr)
        new_ipRGC["ipRGC"] = s["energyFundamentals"]

        # new_ipRGC = new_ipRGC[(new_ipRGC['lambda'] >= self.cfg['winlambda'][0]) & (new_ipRGC['lambda'] <= self.cfg['winlambda'][1])]
        # new_ipRGC['ipRGC'] = new_ipRGC['ipRGC'] /max(new_ipRGC['ipRGC'])
        
        self.ipRGC = new_ipRGC
        
        if self.cfg["visualization"]:
            
            tmp = new_ipRGC.melt(id_vars = "lambda",var_name="LEDs",value_name="Radiance")

            plt.figure()
            sns.lineplot(x="lambda", y="Radiance", hue = "LEDs", data=tmp)
            plt.xlabel("Wavelength[nm]")
            plt.ylabel("Illuminant")
            plt.savefig("./ipRGC.pdf")
            
        return new_ipRGC


    # def plotSpectrum(self):
    #     for ixyz,ilms,irgb in zip(['X','Y','Z'],['L','M','S'],['r','g','b']):
    #         plt.subplot(3,1,2); 
    #     plt.plot(cl_func['lambda'],cl_func[ixyz],irgb)
    
    #     plt.subplot(3,1,3); 
    #     plt.plot(lms['lambda'],lms[ilms],irgb)
    
    #     plt.plot(lms['lambda'],ipRGC['ipRGC'],'k')
    #     # plt.savefig("./color_matching_func.pdf")
    #     plt.show()
# plt.figure(figsize=(6,10))
# plt.subplot(3,1,1);plt.plot(Y_d65['lambda'],Y_d65['Y'])

    def getLEDspectra(self):
        
        ############ data loading (12LEDs embed in LED cube) Output power = 100 (max) ############        
        dat_LED_spectrum = pd.DataFrame()
        
        if self.cfg["visualization"]:
            plt.figure()
            
        for iLED in np.arange(1,max(self.cfg['numOfLEDs'])+1):
            csv_input = pd.read_csv(filepath_or_buffer = self.filePath + "LED/LED" + str(iLED) + ".csv", encoding="ms932", sep=",")
            csv_input = csv_input.drop(['Unnamed: 1','Unnamed: 2','0'], axis=1)
        
            lam = csv_input.drop(np.arange(24), axis=0)
            
            dr = lam['メモ'].tolist()
            dr = [int(i[0:3]) for i in dr]
            dr = np.array(dr)
            
            
            csv_input = csv_input.set_index('メモ')
        
            dat_LED_spectrum['lambda'] = dr
            dat_LED_spectrum['LED' + str(iLED)] = np.array(np.float64(lam[str(self.cfg['outPw'])]))
            
            t = np.array([float(csv_input[str(self.cfg['maxoutPw'])]["X"]),
                          float(csv_input[str(self.cfg['maxoutPw'])]["Y"]),
                          float(csv_input[str(self.cfg['maxoutPw'])]["Z"])])
            t = t/max(t)

            rgb = np.round(np.dot(np.array(xyzToRGB),t),4)*0.5
            rgb[rgb<0] = 0
            rgb[rgb>1] = 1
        
            if self.cfg["visualization"]:
                plt.plot(dat_LED_spectrum['lambda'], np.array(np.float64(lam[str(self.cfg['outPw'])])),color=(rgb))
                # print(rgb)
                
                ind = np.argmax(np.array(np.float64(lam[str(self.cfg['outPw'])])))
                plt.text(dat_LED_spectrum['lambda'][ind], np.array(np.float64(lam[str(self.cfg['outPw'])]))[ind], "LED"+  str(iLED))
            
                # plt.savefig("./LEDs_power.pdf")
                # plt.show()
       
        self.dat_LED_spectrum = dat_LED_spectrum
        
        return dat_LED_spectrum

    def getLEDNTT(self, dat_LED_spectrum):
              
        csv_input = pd.read_csv(filepath_or_buffer = self.filePath + "LED/LED1.csv", encoding="ms932", sep=",")
        csv_input = csv_input.drop(['Unnamed: 1','Unnamed: 2','0'], axis=1)
    
        lam = csv_input.drop(np.arange(24), axis=0)
        dr1 = lam['メモ'].tolist()
        dr1 = [int(i[0:3]) for i in dr1]
        dr1 = np.array(dr1)
            
        if self.cfg["visualization"]:
            plt.figure()
          
        for iLED in np.arange(1,max(self.cfg["numOfLEDs"])+1):
           
            tmp_dat_LED_spectrum = pd.DataFrame()
                        
            folderList=[]
            for filename in os.listdir(self.filePath + "LED/LED"+str(iLED)):
                folderList.append(filename)
            folderList.sort()

            dat = []
            # for pw in np.arange(0,110,10):
            for pw in folderList:
                # print(self.filePath + "LED/LED"+str(iLED)+"/LED"+str(iLED)+'_'+str(pw)+".csv")
            
                with open(self.filePath + "LED/LED"+str(iLED)+"/"+pw, 'r') as temp_f:
                    # get No of columns in each line
                    col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
                    
                column_names = [i for i in range(0, max(col_count))]
                column_names[0] = "lambda"
                column_names[1] = "Y"
                tmp = pd.read_csv(filepath_or_buffer = self.filePath + "LED/LED"+str(iLED)+"/"+pw, sep=",", header=None, names=column_names)
                tmp = tmp[:110]
                dr5 = np.array(tmp["lambda"].astype(int))
                  
                tmp = tmp[(dr5 >= winWaveLen[0]) & (dr5 <= winWaveLen[1])]
                  
                dat.append(tmp['Y'].values)
            
            # tmp_dat_LED_spectrum['lambda'] = dr1
            # tmp_dat_LED_spectrum['LED1'] = np.array(np.float64(lam[str(self.cfg['maxoutPw'])]))
                
            dat = np.array(dat) 
            # dat = dat / np.max(dat) * max(tmp_dat_LED_spectrum['LED1'].values)
            
            dr5 = dr5[(dr5 >= winWaveLen[0]) & (dr5 <= winWaveLen[1])]
            
            y = []
            for i,pw in enumerate(np.arange(len(folderList))):
                  yy = interpolate.PchipInterpolator(dr5, dat[i,:])
                  y.append(yy(np.arange(winWaveLen[0],winWaveLen[1]+1)))
            
            y = np.array(y)
            
            y2 = pd.DataFrame()
            tmp_y = []
            for i,dr in enumerate(dr1):
                yy = interpolate.PchipInterpolator(np.arange(0,1010,100), y[:,i])
                tmp_y.append(yy(np.arange(0,1010,10)))
         
            y2 = pd.DataFrame(tmp_y,columns = np.arange(0,1010,10))
            y2["lambda"] = dr1
            y2 = y2.set_index('lambda')
            
            dat_LED_spectrum["LED"+str(iLED)] = y2[self.cfg["outPw"]].values
            # self.dat_LED_spectrum = dat_LED_spectrum
               
            if self.cfg["visualization"]:
                plt.plot(dat_LED_spectrum['lambda'], y2[self.cfg["outPw"]].values)
                # print(rgb)
                 
                    # ind = np.argmax(y2[self.cfg["outPw"]].values)
                    # plt.text(dat_LED_spectrum['lambda'][ind],y2[self.cfg["outPw"]].values[ind], "LED"+  str(iLED))
    
                ############################################################################
                # LUT = []
                # for wavL in dr1:
                #     tmp = y2.loc[wavL,].values
                #     LUT.append(interpolate.PchipInterpolator(np.arange(0,1010,10),tmp))
                  
                # pickle.dump(LUT, open("LED"+str(iLED)+".pkl", 'wb'))
                ############################################################################
            
        return dat_LED_spectrum
    
    def getProjectorNTT(self):
              
        path = self.filePath + "projector/Spectrum/dat_LED_spectrum_"+str(self.cfg["outPw"])+".json"
        if os.path.isfile(path):
            dat_LED_spectrum = pd.read_json(path)
            
        else:  
            dat_LED_spectrum = pd.DataFrame()
            LEDcount = 1
            for iProjector in ["blueFilt","yellowFilt"]:
                    # aa
                    fileName = glob.glob(self.filePath + "projector/Spectrum/"+iProjector+"_R.csv")
                    fileName = glob.glob("/Users/yutasuzuki/GoogleDrive/PupilAnalysisToolbox/python/preprocessing/lib/LEDdata/projector/Spectrum/"+iProjector+"_R.csv")
                    # measuredWaveLen = np.unique([int(c[-9:-6]) for c in fileName])
                    
                    dat = []
                
                    csv_input = pd.read_csv(filepath_or_buffer = fileName[0],  sep=",", header=None)
                    
                    Y = np.array(csv_input.iloc[1:,2].values, dtype='float')
                    
                    csv_input = csv_input.drop(csv_input.columns[np.arange(18).tolist()], axis=1)
                    csv_input = csv_input.drop(csv_input.columns[np.arange(-25,-1).tolist()], axis=1)
                    csv_input = csv_input.drop(csv_input.columns[[401]], axis=1)
                    
                    dr1 = np.array([int(l[:3]) for l in csv_input.iloc[0,:].values])
                    
                    csv_input = np.array(csv_input.iloc[1:,:].values, dtype='float')
                    # csv_input[abs(csv_input)<(10**(-4)*2)]=0
                    csv_input[csv_input<0]=0
                
                    for irgb in np.arange(3):
                        
                        ind =  np.array(np.arange(irgb*(255/5+1),irgb*(255/5+1)+(255/5+1)), dtype='int')
                        tmp_Y = Y[ind]
                        tmp_dat = csv_input[ind,:]
                        ind_sorted=np.argsort(np.c_[np.arange(255,128,-5),np.arange(0,128,5)].reshape(-1))
                        
                        tmp_dat = tmp_dat[ind_sorted]
                        tmp_Y = tmp_Y[ind_sorted]
                        
                        
                    # import matplotlib.cm as cm
                    # plt.figure()
                    # for i in range(tmp_dat.shape[0]):
                    #     plt.plot(dr1, tmp_dat[i,:], linestyle='solid' , color=cm.Blues(i/tmp_dat.shape[0]))
                        
                    # plt.figure()
                    # for i in range(tmp_Y.shape[0]):
                    #     plt.plot(i, tmp_Y[i], 'o', color=cm.Blues(i/tmp_Y.shape[0]))

                   
                        y2 = pd.DataFrame(tmp_dat.T,columns = np.arange(0,256,5))
                        y2["lambda"] = dr1
                        y2 = y2.set_index('lambda')
                        
                        dat_LED_spectrum["LED"+str(LEDcount)] = y2[self.cfg["outPw"]].values
                        # dat_LED_spectrum["LED"+str(LEDcount)] = y2.loc[:,cfg["outPw"]].values
                               
                        dat_LED_spectrum["lambda"] = dr1
                        
                        ###########################################################################
                        
                        path_pickle = self.filePath + "projector/Spectrum/LED"+str(LEDcount)+".pkl"
                        
                        if not os.path.isfile(path_pickle):
                            print(path_pickle)
                            LUT = []
                            for i,dr in enumerate(dr1):
                                yy = interpolate.PchipInterpolator(np.arange(0,256,5), tmp_dat[:,i])
                                # LUT.append(yy(np.arange(1,256)))
                                LUT.append(yy)
                              
                            pickle.dump(LUT, open(path_pickle, 'wb'))
                            
                        ############################################################################
                        
                        LEDcount+=1
                
                    
            # # for iProjector in np.arange(2)+1:
            # for iProjector in ["_y","_notch"]:
            # # for iProjector in [2,3]:
            #     # for irgb in ["red","green","blue"]:
            #     for irgb in ["r","g","b"]:
                     
            #         fileName = glob.glob(self.filePath + "projector/"+irgb+str(iProjector)+"/*.csv")
            #         measuredWaveLen = np.unique([int(c[-9:-6]) for c in fileName])
                    
            #         dat = []
            #         for iirgb in measuredWaveLen:
                        
            #             if iirgb < 10:
            #                 fileName = glob.glob(self.filePath + "projector/"+irgb+str(iProjector)+"/*_00"+str(iirgb)+"_*.csv")
            #             elif (iirgb >= 10) and (iirgb < 100):
            #                 fileName = glob.glob(self.filePath + "projector/"+irgb+str(iProjector)+"/*_0"+str(iirgb)+"_*.csv")
            #             else:
            #                 # fileName = glob.glob("/Users/yutasuzuki/GoogleDrive/PupilAnalysisToolbox/python/preprocessing/lib/LEDdata/projector/blue1/*.csv")
            #                 fileName = glob.glob(self.filePath + "projector/"+irgb+str(iProjector)+"/*_"+str(iirgb)+"_*.csv")
            #             # print(fileName) 
            #             fileName.sort()
                    
            #             tmp_dat = []
            #             for f in fileName[1:]:
                        
            #                 csv_input = pd.read_csv(filepath_or_buffer = f, encoding="ms932", sep=",", header=None)
            #                 dr1 = csv_input[0].values
                            
            #                 tmp_dat.append(csv_input[1].values[(dr1 >= winWaveLen[0]) & (dr1 <= winWaveLen[1])])
    
            #             dat.append(np.array(tmp_dat).mean(axis=0))
                    
            #         dr1 = dr1[(dr1 >= winWaveLen[0]) & (dr1 <= winWaveLen[1])] 
                    
            #         # if irgb == "green":
            #         #     return dat,dr1
            #         y = []
            #         for i,pw in enumerate(dat):
            #             # print(i)
            #             yy = interpolate.PchipInterpolator(dr1, pw)
            #             y.append(yy(np.arange(winWaveLen[0],winWaveLen[1]+1)))
                    
            #         y = np.array(y)
                    
            #         # y2 = pd.DataFrame()
            #         tmp_y = []
            #         for i,dr in enumerate(dr1):
            #             yy = interpolate.PchipInterpolator(measuredWaveLen, y[:,i])
            #             tmp_y.append(yy(np.arange(1,256)))
                 
            #         y2 = pd.DataFrame(tmp_y,columns = np.arange(1,256))
            #         y2["lambda"] = dr1
            #         y2 = y2.set_index('lambda')
                    
            #         dat_LED_spectrum["LED"+str(LEDcount)] = y2[self.cfg["outPw"]].values
                    
            #         dat_LED_spectrum["lambda"] = dr1
                
                
            #         ###########################################################################
                    
            #         path_pickle = self.filePath + "projector/LED"+str(LEDcount)+".pkl"
                  
            #         if not os.path.isfile(path_pickle):
            #             print(path_pickle)
            #             LUT = []
            #             for wavL in dr1:
            #                 tmp = y2.loc[wavL,].values
            #                 LUT.append(interpolate.PchipInterpolator(np.arange(1,256),tmp))
                          
            #             pickle.dump(LUT, open(path_pickle, 'wb'))
                        
            #         ############################################################################
            #         LEDcount+=1
                
            dat_LED_spectrum.to_json(path)
    
            # print(rgb)
         
            # ind = np.argmax(y2[self.cfg["outPw"]].values)
            # plt.text(dat_LED_spectrum['lambda'][ind],y2[self.cfg["outPw"]].values[ind], "LED"+  str(iLED))
    
        return dat_LED_spectrum    

    def getXYZvalues(self):

        if self.cfg["projector"]:
            dat_LED_spectrum = self.getProjectorNTT()
        else:
            dat_LED_spectrum = self.getLEDspectra()
            dat_LED_spectrum = self.getLEDNTT(dat_LED_spectrum)
        
        cl_func = self.getXYZfunc()
        ipRGC   = self.getipRGCfunc()
        rod   = self.getRodfunc()
        Y_d65 = self.getD65Spectra()
        
        ############ XYZ values of 12LEDs ############
        k = 100/sum(Y_d65['Y'].values * cl_func['Y'].values)
        
        dat_LED_XYZ = {'X':[],'Y':[],'Z':[],'rod':[],'ipRGC':[]}
        
        for iLED in np.arange(1,max(self.cfg['numOfLEDs'])+1):
            dat_LED_XYZ['X'].append(k * sum(dat_LED_spectrum['LED' + str(iLED)].values * cl_func['X'].values))
            dat_LED_XYZ['Y'].append(k * sum(dat_LED_spectrum['LED' + str(iLED)].values * cl_func['Y'].values))
            dat_LED_XYZ['Z'].append(k * sum(dat_LED_spectrum['LED' + str(iLED)].values * cl_func['Z'].values))
            dat_LED_XYZ['rod'].append(sum(dat_LED_spectrum['LED' + str(iLED)].values * rod['rod'].values))
            dat_LED_XYZ['ipRGC'].append(sum(dat_LED_spectrum['LED' + str(iLED)].values * ipRGC['ipRGC'].values))
            
            # dat_LED_XYZ['X'].append(k * sum(dat_LED_spectrum['LED' + str(iLED)].values * lms['L'].values))
            # dat_LED_XYZ['Y'].append(k * sum(dat_LED_spectrum['LED' + str(iLED)].values * lms['M'].values))
            # dat_LED_XYZ['Z'].append(k * sum(dat_LED_spectrum['LED' + str(iLED)].values * lms['S'].values))
            
        return dat_LED_XYZ
    
        ############ XYZ values of 12LEDs ############
        
        # XYZ_d65 = {'X':sum(Y_d65['Y'].values * lms['L'].values),
        #             'Y':sum(Y_d65['Y'].values * lms['M'].values),
        #             'Z':sum(Y_d65['Y'].values * lms['S'].values)}
         
        # XYZ_d65['X'] = XYZ_d65['X'] / XYZ_d65['Y']
        # XYZ_d65['Z'] = XYZ_d65['Z'] / XYZ_d65['Y']
        # XYZ_d65['Y'] = XYZ_d65['Y'] / XYZ_d65['Y']
        
        # plt.plot((XYZ_d65['X']/sum(XYZ_d65.values())),(XYZ_d65['Y']/sum(XYZ_d65.values())),'o')
            
    def seekCombinations(self):

        dat_LED_XYZ = self.getXYZvalues()
        if self.cfg["projector"]:
            dat_LED_spectrum = self.getProjectorNTT()
        else:
            dat_LED_spectrum = self.getLEDspectra()
            dat_LED_spectrum = self.getLEDNTT(dat_LED_spectrum)
            
        cl_func = self.getXYZfunc()
        ipRGC   = self.getipRGCfunc()
        lms_func = self.getLMSfunc()
        
        res = {'lambda':[],
               'LEDs':[],
               'spectrum':[],
               'coeff':[],
               'XYZ':[],'Yxy':[],'LMS':[],'ipRGC':[],'ipRGC_pw':[],
               'RGB':[],
               
               'corrected_coeff':[],
               'corrected_LMS':[],'corrected_XYZ':[],'corrected_ipRGC':[],'corrected_Yxy':[],
               'corrected_RGB':[],
               'corrected_spectrum':[]
               }
        
        ############ seek combination ############
        for iipRGC in self.cfg['ipRGC_outPw']:
            for comb in list(itertools.combinations(np.array(self.cfg['numOfLEDs'])-1,4)):
                # if 8 in comb:
                #     continue
                
                dat = np.matrix([[dat_LED_XYZ['X'][comb[0]], dat_LED_XYZ['X'][comb[1]], dat_LED_XYZ['X'][comb[2]], dat_LED_XYZ['X'][comb[3]]],
                                 [dat_LED_XYZ['Y'][comb[0]], dat_LED_XYZ['Y'][comb[1]], dat_LED_XYZ['Y'][comb[2]], dat_LED_XYZ['Y'][comb[3]]],
                                 [dat_LED_XYZ['Z'][comb[0]] ,dat_LED_XYZ['Z'][comb[1]] ,dat_LED_XYZ['Z'][comb[2]], dat_LED_XYZ['Z'][comb[3]]],
                                 [dat_LED_XYZ['ipRGC'][comb[0]] ,dat_LED_XYZ['ipRGC'][comb[1]] ,dat_LED_XYZ['ipRGC'][comb[2]], dat_LED_XYZ['ipRGC'][comb[3]]]])
                
                dat_inv = dat**(-1)
            
                coeff = np.dot(np.array(dat_inv), np.array([self.XYZ_d65['X'],self.XYZ_d65['Y'],self.XYZ_d65['Z'],iipRGC])).reshape(-1)
                
                # print(coeff)
                # print(any((x < 0 for x in coeff.tolist())))
                # print(any((x > self.cfg['maxoutPw'] for x in coeff.tolist())))
                
                if not (any((x < 0 for x in coeff.tolist()))) | (any((x > self.cfg['maxoutPw'] for x in coeff.tolist()))):
                    # spectrum_simulated_d65 = [dat_LED_spectrum['LED' + str(comb[0]+1)].values * (coeff[0]/self.cfg['maxoutPw']),
                    #                           dat_LED_spectrum['LED' + str(comb[1]+1)].values * (coeff[1]/self.cfg['maxoutPw']),
                    #                           dat_LED_spectrum['LED' + str(comb[2]+1)].values * (coeff[2]/self.cfg['maxoutPw']),
                    #                           dat_LED_spectrum['LED' + str(comb[3]+1)].values * (coeff[3]/self.cfg['maxoutPw'])]
                    
                    # divided by the outPw because LED peak was set as "outPw" 
                    spectrum_simulated_d65 = [dat_LED_spectrum['LED' + str(comb[0]+1)].values * (coeff[0]/self.cfg['outPw']),
                                              dat_LED_spectrum['LED' + str(comb[1]+1)].values * (coeff[1]/self.cfg['outPw']),
                                              dat_LED_spectrum['LED' + str(comb[2]+1)].values * (coeff[2]/self.cfg['outPw']),
                                              dat_LED_spectrum['LED' + str(comb[3]+1)].values * (coeff[3]/self.cfg['outPw'])]
                    # print(np.dot(np.array(dat), np.array(coeff)))
                          
                    # spectrum_simulated_d65 = [dat_LED_spectrum['LED' + str(comb[0]+1)].values * (coeff[0]),
                    #                           dat_LED_spectrum['LED' + str(comb[1]+1)].values * (coeff[1]),
                    #                           dat_LED_spectrum['LED' + str(comb[2]+1)].values * (coeff[2]),
                    #                           dat_LED_spectrum['LED' + str(comb[3]+1)].values * (coeff[3])]
                    
                    tmp = spectrum_simulated_d65.copy()
                    spectrum_simulated_d65 = np.sum(spectrum_simulated_d65,axis=0)
                          
                    # print(any((x < 0 for x in spectrum_simulated_d65.tolist())))
                    if not any((x < 0 for x in spectrum_simulated_d65.tolist())):
                        
                        t = [sum(spectrum_simulated_d65 * cl_func['X'].values),
                             sum(spectrum_simulated_d65 * cl_func['Y'].values),
                             sum(spectrum_simulated_d65 * cl_func['Z'].values)]
                        
                        t = np.array(t)/t[1]
                        
                        lms = [sum(spectrum_simulated_d65 * lms_func['L']),
                               sum(spectrum_simulated_d65 * lms_func['M']),
                               sum(spectrum_simulated_d65 * lms_func['S'])]
                    
                        res['LEDs'].append((np.array(comb)+1).tolist())
                        # res['spectrum'].append(np.array(tmp).tolist())
                        res['spectrum'].append(np.array(tmp).tolist())
                        res['lambda'].append(np.array(dat_LED_spectrum["lambda"]).tolist())
                        res['coeff'].append(np.round(coeff).tolist())
                        res['XYZ'].append(np.round(t,3).tolist())
                        
                        res['Yxy'].append(np.round([t[1],(t[0]/sum(t)),(t[1]/sum(t))],5).tolist())
                        
                        # res['LMS'].append(np.dot(np.array(xyzTolms),np.array(t)).tolist())
                        
                        res['LMS'].append(np.round(lms,5).tolist())
                        
                        res['ipRGC'].append(np.round(sum(spectrum_simulated_d65 * ipRGC['ipRGC'].values),3))
                        res['ipRGC_pw'].append(float(iipRGC))
                        
                        t = t*0.2
                        rgb = np.round(np.dot(np.array(xyzToRGB),t),4)
                        rgb[rgb<0] = 0
                        rgb[rgb>1] = 1
                        rgb = rgb * 255
                        res['RGB'].append(rgb.tolist())
                  
        print("We found " +str(len(res["LEDs"])))
        return res
    
    def seekCombinations5LEDs(self,dat_LED_XYZ,dat_LED_spectrum,cl_func,ipRGC,rod):

        res = {'lambda':[],
               'LEDs':[],
               'spectrum':[],
               'coeff':[],
               'XYZ':[],'Yxy':[],'LMS':[],'ipRGC':[],'rod':[],
               'RGB':[],
               
               'corrected_coeff':[],
               'corrected_LMS':[],'corrected_XYZ':[],'corrected_ipRGC':[],'corrected_Yxy':[],
               'corrected_RGB':[],
               'corrected_spectrum':[]
               }
        
        ############ seek combination ############
        for irod in self.cfg['rod_outPw']:
            for iipRGC in self.cfg['ipRGC_outPw']:
                for comb in list(itertools.combinations(np.array(self.cfg['numOfLEDs'])-1,5)):
                    M = np.matrix([[dat_LED_XYZ['X'][comb[0]], dat_LED_XYZ['X'][comb[1]], dat_LED_XYZ['X'][comb[2]], dat_LED_XYZ['X'][comb[3]], dat_LED_XYZ['X'][comb[4]]],
                                   [dat_LED_XYZ['Y'][comb[0]], dat_LED_XYZ['Y'][comb[1]], dat_LED_XYZ['Y'][comb[2]], dat_LED_XYZ['Y'][comb[3]], dat_LED_XYZ['Y'][comb[4]]],
                                   [dat_LED_XYZ['Z'][comb[0]], dat_LED_XYZ['Z'][comb[1]] ,dat_LED_XYZ['Z'][comb[2]], dat_LED_XYZ['Z'][comb[3]], dat_LED_XYZ['Z'][comb[4]]],
                                   [dat_LED_XYZ['rod'][comb[0]], dat_LED_XYZ['rod'][comb[1]], dat_LED_XYZ['rod'][comb[2]], dat_LED_XYZ['rod'][comb[3]], dat_LED_XYZ['rod'][comb[4]]],
                                   [dat_LED_XYZ['ipRGC'][comb[0]], dat_LED_XYZ['ipRGC'][comb[1]], dat_LED_XYZ['ipRGC'][comb[2]], dat_LED_XYZ['ipRGC'][comb[3]], dat_LED_XYZ['ipRGC'][comb[4]]]
                                   ])
                    
                    # dat = np.matrix([[dat_LED_XYZ['X'][comb[0]], dat_LED_XYZ['X'][comb[1]], dat_LED_XYZ['X'][comb[2]], dat_LED_XYZ['X'][comb[3]]],
                    #                  [dat_LED_XYZ['Y'][comb[0]], dat_LED_XYZ['Y'][comb[1]], dat_LED_XYZ['Y'][comb[2]], dat_LED_XYZ['Y'][comb[3]]],
                    #                  [dat_LED_XYZ['Z'][comb[0]], dat_LED_XYZ['Z'][comb[1]] ,dat_LED_XYZ['Z'][comb[2]], dat_LED_XYZ['Z'][comb[3]]],
                    #                  [dat_LED_XYZ['rod'][comb[0]],   dat_LED_XYZ['rod'][comb[1]],   dat_LED_XYZ['rod'][comb[2]],   dat_LED_XYZ['rod'][comb[3]]]
                    #                  ])
                    
                    M_inv = M**(-1)
                    # M_inv2 = np.linalg.inv(M)
                
                    coeff = np.dot(np.array(M_inv), np.array([self.XYZ_d65['X'],self.XYZ_d65['Y'],self.XYZ_d65['Z'],irod,iipRGC])).reshape(-1)
                    
                    # coeff = np.dot(np.array(M_inv),np.array([XYZ_d65['X'],XYZ_d65['Y'],XYZ_d65['Z'],irod])).reshape(-1)
                    
                    # coeff = np.dot(np.array([XYZ_d65['X'],XYZ_d65['Y'],XYZ_d65['Z'],irod,iipRGC]),np.array(dat_inv)).reshape(-1)
                    # coeff = np.dot(np.array(M_inv),np.array([XYZ_d65['X'],XYZ_d65['Y'],XYZ_d65['Z'],irod,iipRGC])).reshape(-1)
                    # coeff = np.dot(np.array(M_inv),np.array([XYZ_d65['X'],XYZ_d65['Y'],XYZ_d65['Z'],0.1,iipRGC])).reshape(-1)
                    
                    # test = np.dot(np.array(dat),coeff)
                    
                    # return coeff
                
                    if not (any((x < 0 for x in coeff.tolist()))) | (any((x > self.cfg['maxoutPw'] for x in coeff.tolist()))):
                        spectrum_simulated_d65 = [dat_LED_spectrum['LED' + str(comb[0]+1)].values * (coeff[0]/self.cfg['outPw']),
                                                  dat_LED_spectrum['LED' + str(comb[1]+1)].values * (coeff[1]/self.cfg['outPw']),
                                                  dat_LED_spectrum['LED' + str(comb[2]+1)].values * (coeff[2]/self.cfg['outPw']),
                                                  dat_LED_spectrum['LED' + str(comb[3]+1)].values * (coeff[3]/self.cfg['outPw']),
                                                  dat_LED_spectrum['LED' + str(comb[4]+1)].values * (coeff[4]/self.cfg['outPw'])]
                        
                        tmp = spectrum_simulated_d65.copy()
                        spectrum_simulated_d65 = np.sum(spectrum_simulated_d65,axis=0)
                              
                        if not any((x < 0 for x in spectrum_simulated_d65.tolist())):
                            
                            t = [sum(spectrum_simulated_d65 * cl_func['X'].values),
                                 sum(spectrum_simulated_d65 * cl_func['Y'].values),
                                 sum(spectrum_simulated_d65 * cl_func['Z'].values)]   
                        
                            res['LEDs'].append((np.array(comb)+1).tolist())
                            res['spectrum'].append(np.array(tmp).tolist())
                            res['lambda'].append(np.array(dat_LED_spectrum["lambda"]).tolist())
                            res['coeff'].append(np.round(coeff).tolist())
                            res['XYZ'].append(t)
                            res['Yxy'].append([t[1],(t[0]/sum(t)),(t[1]/sum(t))])
                            
                            res['LMS'].append(np.dot(np.array(xyzTolms),np.array(t)).tolist())
                            
                            # res['LMS'].append([sum(spectrum_simulated_d65 * lms['L'].values),
                            #                    sum(spectrum_simulated_d65 * lms['M'].values),
                            #                    sum(spectrum_simulated_d65 * lms['S'].values)])
                            
                            res['ipRGC'].append(sum(spectrum_simulated_d65 * ipRGC['ipRGC'].values))
                            
                            res['rod'].append(sum(spectrum_simulated_d65 * rod['rod'].values))
                            
                            t = t/max(t)
                            rgb = np.round(np.dot(np.array(xyzToRGB),t),4)
                            rgb[rgb<0] = 0
                            rgb[rgb>1] = 1
                            res['RGB'].append(rgb.tolist())
                        
        return res
    
    
    def optimizeLight(self, spectra, y_measured,norm=False):
        
        if self.cfg["projector"]:
            
            LEDs = []
            for iLED in np.arange(1,max(self.cfg['numOfLEDs'])+1):
                LEDs.append(pickle.load(open(self.filePath + "projector/LED"+str(iLED)+".pkl", 'rb')))
            
            waveL = np.arange(spectra["lambda"][0][0],spectra["lambda"][0][-1]+1)
            
            tmp = []
            for iCoeff,iLED,iSpectrum in tqdm(zip(spectra['coeff'],spectra['LEDs'],spectra["spectrum"])):
             
                result = minimize(f_optLEDs, np.array(iCoeff), 
                                  method="Nelder-Mead", 
                                  bounds=((0, 255),(0, 255),(0, 255),(0, 255)),
                                  args=(iLED,waveL,np.array(iSpectrum).max(axis=1),LEDs))
                
                t = np.round(result.x)
                tmp.append(t.tolist())
            
            # res["optimized_coeff"] = tmp
            
        # path = self.filePath + spectrometer + "/LED_summary.json"
    
        # if os.path.isfile(path):
        #     f = open(os.path.join(path))
        #     lut_spec = json.load(f)
        #     f.close()
            
        else:
            LEDs = []
            for iLED in np.arange(1,max(self.cfg['numOfLEDs'])+1):
                LEDs.append(pickle.load(open(self.filePath + "LED/LED"+str(iLED)+".pkl", 'rb')))
            
            waveL = np.arange(spectra["lambda"][0],spectra["lambda"][-1]+1)
            
            dat = []
            for coeff in spectra["LEDs"]:
                dr = int(np.argwhere(waveL == self.peakLength[coeff-1]).reshape(-1))
            
                dat.append(y_measured[(dr-5):(dr+5)].max())
                
            # print(dat)
            
            result = minimize(f_optLEDs, np.array(spectra["corrected_coeff"]), 
                              method="Nelder-Mead",
                              bounds=((0, 1000),(0, 1000),(0, 1000),(0, 1000)),
                              args=(spectra["LEDs"], waveL, np.array(dat), LEDs, norm))
        
        return np.round(result.x).tolist()
    
    def optimizeLEDs(self,res,spectrometer = "cl-500a"):
        
        dat_LED_XYZ = self.getXYZvalues()
        
        if self.cfg["projector"]:
            dat_LED_spectrum = self.getProjectorNTT()
        else:
            dat_LED_spectrum = self.getLEDspectra()
            dat_LED_spectrum = self.getLEDNTT(dat_LED_spectrum)
        
        cl_func = self.getXYZfunc()
        ipRGC   = self.getipRGCfunc()
        lms_func = self.getLMSfunc()
        
        
        if self.cfg["projector"]:
            
            LEDs = []
            for iLED in np.arange(1,max(self.cfg['numOfLEDs'])+1):
                LEDs.append(pickle.load(open(self.filePath + "projector/LED"+str(iLED)+".pkl", 'rb')))
            
            waveL = np.arange(res["lambda"][0][0],res["lambda"][0][-1]+1)
            
            tmp = []
            for iCoeff,iLED,iSpectrum in tqdm(zip(res['coeff'],res['LEDs'],res["spectrum"])):
             
                result = minimize(f_optLEDs, np.array(iCoeff), 
                                  method="Nelder-Mead", 
                                  bounds=((0, 255),(0, 255),(0, 255),(0, 255)),
                                  args=(iLED,waveL,np.array(iSpectrum).max(axis=1),LEDs))
                
                t = np.round(result.x)
                tmp.append(t.tolist())
            
            res["optimized_coeff"] = tmp
            
            
        else:
            LEDs = []
            for iLED in np.arange(1,max(self.cfg['numOfLEDs'])+1):
                LEDs.append(pickle.load(open(self.filePath + "LED/LED"+str(iLED)+".pkl", 'rb')))
            
            waveL = np.arange(res["lambda"][0][0],res["lambda"][0][-1]+1)
            
            tmp = []
            for iCoeff,iLED,iSpectrum in tqdm(zip(res['coeff'],res['LEDs'],res["spectrum"])):
                
                # yy = interpolate.PchipInterpolator(res["lambda"][0], np.array(iSpectrum).sum(axis=0))
             
                result = minimize(f_optLEDs, np.array(iCoeff), 
                                  method="Nelder-Mead", 
                                  bounds=((0, 1000),(0, 1000),(0, 1000),(0, 1000)),
                                  args=(iLED,waveL,np.array(iSpectrum).max(axis=1),LEDs))
                
                t = np.round(result.x)
                tmp.append(t.tolist())
            
            res["optimized_coeff"] = tmp
            
            
            
        for mmName in ["optimized_spectrum","optimized_XYZ","optimized_Yxy","optimized_LMS","optimized_ipRGC","optimized_RGB"]:
            res[mmName] = []


        ################# validation using measured spectrum #################
            
        # if spectrometer == "projector":
        
            
        if spectrometer == "cl-500a":
            
            ################# validation using measured spectrum #################
            path = self.filePath + "LED//.json"
            if os.path.isfile(path):
                f = open(os.path.join(path))
                lut_spec = json.load(f)
                f.close()
    
            else:
                lut_spec = {}
                st = np.arange(1,1001)
                
                LEDs = []
                for iLED in np.arange(1,13):
                    LEDs.append(pickle.load(open(self.filePath + "LED/LED"+str(iLED)+".pkl", 'rb')))
              
                for iLED in  np.unique( res['LEDs']):
                    tmp=[]
                    for i,dr in enumerate(res["lambda"][0]):
                        tmp2=[]
                        for pw in st:
                            tmp2.append(LEDs[iLED-1][i](pw))
                        tmp.append(np.array(tmp2))
                    
                    lut_spec['LED' + str(iLED)] = np.array(tmp).T
                    lut_spec['LED' + str(iLED) + '_peak'] = np.array(tmp).max(axis=0)
                        
                    for mm in list(lut_spec.keys()):
                        if not isinstance(lut_spec[mm],list):
                            lut_spec[mm] = lut_spec[mm].tolist()
        
                with open(os.path.join(self.filePath + "LED/LED_summary.json"),"w") as f:
                    json.dump(lut_spec,f)
                    
        else:
            lut_spec = {}
            st = np.arange(1,1001)
            
            for iLED in  np.unique( res['LEDs']):
                
                csv_input = pd.read_csv(filepath_or_buffer = self.filePath + "LED/LED" + str(iLED) + ".csv", encoding="ms932", sep=",")
                csv_input = csv_input.drop(['Unnamed: 1','Unnamed: 2','0'], axis=1)
            
                csv_input = csv_input.drop(np.arange(24), axis=0)
                dr = csv_input['メモ'].tolist()
                x = np.arange(10,1010,10)
                
                lut_spec['LED' + str(iLED)] = []
                for ilambda in np.arange(len(dr)):
                    
                    y = np.float64(csv_input.values[ilambda,1:])
                    a1, a2, b = np.polyfit(x, y, 2)
                    
                    lut_spec['LED' + str(iLED)].append(np.array(a1 * st**2 + a2 * st + b))
                    
                lut_spec['LED' + str(iLED) + '_peak'] = np.max(np.array(lut_spec['LED' + str(iLED)]),axis=0)
                lut_spec['LED' + str(iLED)] = np.array(lut_spec['LED' + str(iLED)])
        
        
        for i,(iLight,iCoeff) in enumerate(zip(res['LEDs'],res['optimized_coeff'])):
            
            spectrum_simulated_d65 = []
            for iLED,ic in zip(iLight,iCoeff):
                spectrum_simulated_d65.append(lut_spec['LED' + str(iLED)][int(ic)])
            
            spectrum_simulated_d65 = np.sum(spectrum_simulated_d65,axis=0)
            # spectrum_simulated_d65 = spectrum_simulated_d65[(res["lambda"] >= self.cfg['winlambda'][0]) & (res["lambda"] <= self.cfg['winlambda'][1])]
            
            t = [sum(spectrum_simulated_d65 * cl_func['X'].values),
                 sum(spectrum_simulated_d65 * cl_func['Y'].values),
                 sum(spectrum_simulated_d65 * cl_func['Z'].values)]
            
            lms = [sum(spectrum_simulated_d65 * lms_func['L']),
                   sum(spectrum_simulated_d65 * lms_func['M']),
                   sum(spectrum_simulated_d65 * lms_func['S'])]
            
            res['optimized_spectrum'].append(spectrum_simulated_d65.tolist())
            res['optimized_XYZ'].append(t)
            res['optimized_Yxy'].append([t[1],(t[0]/sum(t)),(t[1]/sum(t))])
            
            # res['optimized_LMS'].append(np.dot(np.array(xyzTolms),np.array(t)).tolist())
            
            res['optimized_LMS'].append(lms)
            res['optimized_ipRGC'].append(sum(spectrum_simulated_d65 * ipRGC['ipRGC'].values))
    
            t = t/max(t)
            rgb = np.round(np.dot(np.array(xyzToRGB),t),4)
            rgb[rgb<0] = 0
            rgb[rgb>1] = 1
            res['optimized_RGB'].append(rgb.tolist())
            
        return res

    def validation(self,res,spectrometer="cl-500a"):
     
        lms_func = self.getLMSfunc()
        if spectrometer == "projector":
            
            path = self.filePath + "projector/Spectrum/LED_summary.json"
            
            if os.path.isfile(path):
                f = open(os.path.join(path))
                lut_spec = json.load(f)
                f.close()
    
            else:
                lut_spec = {}
                # st = np.arange(1,256)
                st = np.arange(255)
                
                flieNameLEDs = glob.glob(self.filePath + "/projector/Spectrum/*.pkl")
                # flieNameLEDs = glob.glob("/Users/yutasuzuki/GoogleDrive/PupilAnalysisToolbox/python/preprocessing/lib/LEDdata/projector/Spectrum/*.pkl")
                flieNameLEDs.sort()
                # print(flieNameLEDs)
                # flieNameLEDs =  glob.glob("./LEDdata/projector/*.pkl")
                LEDs = []
                for f in flieNameLEDs:
                    LEDs.append(pickle.load(open(f, 'rb')))
                    
                for iLED in self.cfg['numOfLEDs']:
                    tmp=[]
                    for i,dr in enumerate(np.arange(self.cfg["winlambda"][0],self.cfg["winlambda"][1]+1)):
                        tmp2=[]
                        for pw in st:
                            tmp2.append(LEDs[iLED-1][i](pw))

                        tmp.append(np.array(tmp2))
                    
                    lut_spec['LED' + str(iLED)] = np.array(tmp).T
                    lut_spec['LED' + str(iLED) + '_peak'] = np.array(tmp).max(axis=0)
                        
                    for mm in list(lut_spec.keys()):
                        if not isinstance(lut_spec[mm],list):
                            lut_spec[mm] = lut_spec[mm].tolist()
        
                # lut_spec.to_json(path)
                 
                with open(os.path.join(self.filePath + "projector/Spectrum/LED_summary.json"),"w") as f:
                    json.dump(lut_spec,f)
            
        elif spectrometer == "cl-500a":
            ################# validation using measured spectrum #################
            path = self.filePath + "LED/LED_summary.json"
            if os.path.isfile(path):
                f = open(os.path.join(path))
                lut_spec = json.load(f)
                f.close()
    
            else:
                lut_spec = {}
                st = np.arange(0,1001)
                
                LEDs = []
                for iLED in np.arange(1,13):
                    LEDs.append(pickle.load(open(self.filePath + "LED/LED"+str(iLED)+".pkl", 'rb')))
              
                for iLED in np.unique(res['LEDs']):
                    tmp=[]
                    for i,dr in enumerate(res["lambda"][0]):
                        tmp2=[]
                        for pw in st:
                            tmp2.append(LEDs[iLED-1][i](pw))
                        tmp.append(np.array(tmp2))
                    
                    lut_spec['LED' + str(iLED)] = np.array(tmp).T
                    lut_spec['LED' + str(iLED) + '_peak'] = np.array(tmp).max(axis=0)
                        
                    for mm in list(lut_spec.keys()):
                        if not isinstance(lut_spec[mm],list):
                            lut_spec[mm] = lut_spec[mm].tolist()
        
                with open(os.path.join(self.filePath + "LED/LED_summary.json"),"w") as f:
                    json.dump(lut_spec,f)
        
        else:
            for iLED in  np.unique( res['LEDs']):
                
                csv_input = pd.read_csv(filepath_or_buffer = self.filePath + "LED/LED" + str(iLED) + ".csv", encoding="ms932", sep=",")
                csv_input = csv_input.drop(['Unnamed: 1','Unnamed: 2','0'], axis=1)
            
                csv_input = csv_input.drop(np.arange(24), axis=0)
                dr = csv_input['メモ'].tolist()
                x = np.arange(10,1010,10)
                
                lut_spec['LED' + str(iLED)] = []
                for ilambda in np.arange(len(dr)):
                    
                    y = np.float64(csv_input.values[ilambda,1:])
                    a1, a2, b = np.polyfit(x, y, 2)
                    
                    lut_spec['LED' + str(iLED)].append(np.array(a1 * st**2 + a2 * st + b))
                    
                lut_spec['LED' + str(iLED) + '_peak'] = np.max(np.array(lut_spec['LED' + str(iLED)]),axis=0)
                lut_spec['LED' + str(iLED)] = np.array(lut_spec['LED' + str(iLED)])
        
        ################# validation using measured spectrum #################
        
        for iLight,ipickedLED in zip(res['spectrum'],res['LEDs']):
            corrected_coeff = []
            for iLED,numLED in zip(iLight,ipickedLED):
                corrected_coeff.append(int(getNearestValue(lut_spec['LED' + str(numLED) + '_peak'],max(iLED))))
                a = int(getNearestValue(lut_spec['LED' + str(numLED) + '_peak'],max(iLED)))
                
                # plt.plot(lut_spec['LED' + str(numLED)][a])
                # plt.plot(iLED)
                
            res['corrected_coeff'].append(corrected_coeff)
        
        
        for i,(iLight,iCoeff) in enumerate(zip(res['LEDs'],res['corrected_coeff'])):
            spectrum_simulated_d65 = []
            for iLED,ic in zip(iLight,iCoeff):
                spectrum_simulated_d65.append(lut_spec['LED' + str(iLED)][int(ic)])
            
            spectrum_simulated_d65 = np.array(spectrum_simulated_d65).sum(axis=0)
            
            # spectrum_simulated_d65 = spectrum_simulated_d65[(res["lambda"][0] >= self.cfg['winlambda'][0]) & (res["lambda"][0]<= self.cfg['winlambda'][1])]
        
            # spectrum_simulated_d65 = spectrum_simulated_d65[np.arange(self.cfg['winlambda'][0],self.cfg['winlambda'][1]+1) % 5 == 0]
            
            t = [sum(spectrum_simulated_d65 * self.cl_func['X'].values),
                  sum(spectrum_simulated_d65 * self.cl_func['Y'].values),
                  sum(spectrum_simulated_d65 * self.cl_func['Z'].values)]
            
            lms = [sum(spectrum_simulated_d65 * lms_func['L']),
                    sum(spectrum_simulated_d65 * lms_func['M']),
                    sum(spectrum_simulated_d65 * lms_func['S'])]
            
            # t = np.array(t) / t[1]
            res['corrected_spectrum'].append(spectrum_simulated_d65.tolist())
            res['corrected_XYZ'].append(np.round(t,3).tolist())
            res['corrected_Yxy'].append(np.round([t[1],(t[0]/sum(t)),(t[1]/sum(t))],5).tolist())
            
            # res['corrected_LMS'].append(np.dot(np.array(xyzTolms),np.array(t)).tolist())
            
            res['corrected_LMS'].append(np.round(lms,5).tolist())
            
            res['corrected_ipRGC'].append(np.round(sum(spectrum_simulated_d65 * self.ipRGC['ipRGC'].values),3))
    
            t = t/max(t)
            rgb = np.round(np.dot(np.array(xyzToRGB),t),4)
            rgb[rgb<0] = 0
            rgb[rgb>1] = 1
            res['corrected_RGB'].append(rgb.tolist())

        return res
        # return lut_spec
    
    def rejectOutlier(self,res,varName = "corrected_Yxy"):

        
        ################# reject outlier #################
        reject = {"coeff":[],
                  "Y":[],
                  "xy":[]
                  }
        
        #%% coeff
        for i, coeff in enumerate(res["corrected_coeff"]):
            if self.cfg["maxoutPw"]-1 in coeff:
                reject["coeff"].append(i)
            if 0 in coeff:
                reject["coeff"].append(i)
        
        
        #%% xy
        
        if "optimized" in varName:
            x = np.array(res[varName])[:,1]
            y = np.array(res[varName])[:,2]
           
            x_norm = x - x.mean()
            y_norm = y - y.mean()
            
            stdX = x_norm.std()*0.5
            stdY = y_norm.std()*0.2
            
            tmp_x = x_norm**2
            tmp_y = y_norm**2
            
            P = (tmp_x/stdX**2)+(tmp_y/stdY**2)-1
            
            reject["xy"] = np.argwhere(P > 0).reshape(-1).tolist()
                
            fig = plt.figure()
            ax = plt.axes()
            e = patches.Ellipse(xy=(0,0), width=stdX*2, height=stdY*2, fill=False, ec='r')
            ax.add_patch(e)
            plt.plot(x_norm,y_norm,'o')
            plt.axis('scaled')
            plt.xlim([-stdX*1.5,stdX*1.5])
            plt.ylim([-stdY*1.5,stdY*1.5])
    
    
        #%% Y
        y = np.array(res[varName])[:,0]
    
        reject["Y"] = np.argwhere((y > (y.mean()+y.std()*0.5)) | (y < (y.mean()-y.std()*0.5))).reshape(-1).tolist()
        
        plt.hist(y)
        plt.vlines(y.mean()+y.std()*0.5, 0, 200,"red")
        plt.vlines(y.mean()-y.std()*0.5, 0, 200,"red")
        
        rejectNum = np.unique(reject["coeff"] + reject["Y"] + reject["xy"])
        
        for mm in res.keys():
            res[mm] = [d for i,d in enumerate(res[mm]) if not i in rejectNum]
        
        
        return res
   
    
    def getMinMax(self,res, varName="optimized_ipRGC"):
    
        ################# choose min/max and middle level ipRGC #################
        
        # filt = (res['ipRGC'] == max(res['ipRGC'])) | (res['ipRGC'] == min(res['ipRGC']))
        
        mid = getNearestValue(res[varName],np.mean(res[varName]))
        
        
        filt = [np.argmin(res[varName]),np.argmax(res[varName]),mid]
        # filt = [np.argmin(res[varName]),np.argmax(res[varName])]
        # filt = [np.argmax(res[varName])]
        # filt = [mid]
        
        # for mm in list(res.keys()):
        #     res[mm] = [d for i,d in enumerate(res[mm]) if filt[i]]
       
        for mm in list(res.keys()):
            res[mm] = [d for i,d in enumerate(res[mm]) if i in filt]
        
        print("Target / Control ratio = " + str(res[varName][-1]/res[varName][0]))
        
        return res
    
    def plot(self,res,dat_LED_spectrum):
        
        ################# data plot #################
        
        plt.figure(figsize=(8,8))
        plt.style.use('ggplot')
        
        text = " L = {}, M = {}, /n S = {}, ipRGC = {}"
        
        for i,(sp,d_lms,d_Yxy,rgb) in enumerate(zip(res['spectrum'],res['LMS'],res['Yxy'],res['RGB'])):
            
            plt.subplot(2, 3, 1)
            plt.plot(dat_LED_spectrum['lambda'],np.sum(np.array(sp),axis=0),color=(np.array(rgb)))
            # plt.ylim([0,0.02])
            plt.xlabel("wavelength")
            plt.ylabel("Power")
            
            plt.text(500, 0.01+(i*0.003), 
                      text.format(np.round(d_lms[0],4),
                                  np.round(d_lms[1],4),
                                  np.round(d_lms[2],4),
                                  np.round(res['ipRGC'][i],4)),
                                  fontsize=10)
            plt.subplot(2, 3, 2)
            plt.plot(d_Yxy[1],d_Yxy[2],'o',color=(np.array(rgb)))
            # plt.xlim([0.28,0.34])
            # plt.ylim([0.3,0.4])
            plt.xlabel("x")
            plt.ylabel("y")
        
        plt.subplot(2, 3, 3)
        plt.hist(d_Yxy[0])
        plt.hist(np.array(res['Yxy'])[:,0])
        # plt.xlim([0.28,0.34])
        # plt.ylim([0.3,0.4])
        plt.xlabel("x")
        plt.ylabel("y")
        
        for i,(sp,d_lms,d_Yxy,rgb) in enumerate(zip(res['corrected_spectrum'],res['corrected_LMS'],res['corrected_Yxy'],res['corrected_RGB'])):
            
            plt.subplot(2, 3, 4)
            plt.plot(dat_LED_spectrum['lambda'],sp,color=(np.array(rgb)))
            # plt.ylim([0,0.02])
            plt.xlabel("wavelength")
            plt.ylabel("Power")
            
            plt.text(500, 0.01+(i*0.003), 
                      text.format(np.round(d_lms[0],4),
                                  np.round(d_lms[1],4),
                                  np.round(d_lms[2],4),
                                  np.round(res['corrected_ipRGC'][i],4)),
                                  fontsize=10)  
            plt.subplot(2, 3, 5)
            plt.plot(d_Yxy[1],d_Yxy[2],'o',color=(np.array(rgb)))
            # print(rgb)
            # plt.xlim([0.28,0.34])
            # plt.ylim([0.3,0.4])
            plt.xlabel("x")
            plt.ylabel("y")
            
        plt.subplot(2, 3, 6)
        plt.hist(np.array(res['corrected_Yxy'])[:,0])
        # plt.xlim([0.28,0.34])
        # plt.ylim([0.3,0.4])
        plt.xlabel("x")
        plt.ylabel("y")
        # plt.savefig("./calcuratedOut.pdf")
        
        # plt.savefig("res.eps")
        plt.show()
        

    def saveData(self,res,fName="test"):
    
                
        if self.cfg["plus"]:
            fName = "x"+str(int(self.cfg['x']*100)) + "y"+str(int(self.cfg['y']*100))+ "ipRGC"+str(int(res["ipRGC"][-1]*100))+"v"+str(int(self.cfg["VISUAL_ANGLE"]))
            
        else:
            fName = "x"+str(int(self.cfg['x']*100)) + "y"+str(int(self.cfg['y']*100))+ "ipRGC"+str(int(res["ipRGC"][0]*100))+"v"+str(int(self.cfg["VISUAL_ANGLE"]))
        
        print("The file is saved at ./data_LEDCube_"+date+fName+".json")
        with open(os.path.join("./data_LEDCube_"+date+fName+".json"),"w") as f:
            json.dump(res,f)


    def showProjectorSpectra(self, numLEDs, pw=0):
    
        led = pickle.load(open(self.filePath + "projector/Spectrum/LED"+str(numLEDs)+".pkl", 'rb'))
        
        # tmp = []
        # for dr in np.arange(len(led)):
        #     tmp.append(np.float(led[dr](pw)))

    
        # csv_input = pd.read_csv(filepath_or_buffer = self.filePath + "LED/LED1.csv", encoding="ms932", sep=",")
        # csv_input = csv_input.drop(['Unnamed: 1','Unnamed: 2','0'], axis=1)
    
        # lam = csv_input.drop(np.arange(24), axis=0)
        # dr1 = lam['メモ'].tolist()
        # dr1 = [int(i[0:3]) for i in dr1]
        # dr1 = np.array(dr1)
           
        # plt.plot(dr1,tmp)

        return led
    
      
        # dat_LED_spectrum['lambda'] = dr
        # dat_LED_spectrum['LED' + str(numLEDs)] = np.array(np.float64(lam[str(pw)]))
        
        # plt.plot(dat_LED_spectrum['lambda'], np.array(np.float64(lam[str(pw)])))
        
        # ind = np.argmax(np.array(np.float64(lam[str(pw)])))
        # plt.text(dat_LED_spectrum['lambda'][ind], np.array(np.float64(lam[str(pw)]))[ind], "LED" + str(numLEDs))
        
    def showLEDSpectra(self, numLEDs, pw=0):
    
        led = pickle.load(open(self.filePath + "LED/LED"+str(numLEDs)+".pkl", 'rb'))
        
        tmp = []
        for dr in np.arange(len(led)):
            tmp.append(np.float(led[dr](pw)))

    
        csv_input = pd.read_csv(filepath_or_buffer = self.filePath + "LED/LED1.csv", encoding="ms932", sep=",")
        csv_input = csv_input.drop(['Unnamed: 1','Unnamed: 2','0'], axis=1)
    
        lam = csv_input.drop(np.arange(24), axis=0)
        dr1 = lam['メモ'].tolist()
        dr1 = [int(i[0:3]) for i in dr1]
        dr1 = np.array(dr1)
           
        plt.plot(dr1,tmp)

        return tmp
    
      
        # dat_LED_spectrum['lambda'] = dr
        # dat_LED_spectrum['LED' + str(numLEDs)] = np.array(np.float64(lam[str(pw)]))
        
        # plt.plot(dat_LED_spectrum['lambda'], np.array(np.float64(lam[str(pw)])))
        
        # ind = np.argmax(np.array(np.float64(lam[str(pw)])))
        # plt.text(dat_LED_spectrum['lambda'][ind], np.array(np.float64(lam[str(pw)]))[ind], "LED" + str(numLEDs))
        
    def showMeasLightSummary(self,spectrometer="cl-500a"):
        
        rod   = self.getRodfunc()
        ipRGC   = self.getipRGCfunc()
        
        if spectrometer == "projector":
            path = self.filePath + spectrometer + "/dat_LED_spectrum_"+str(self.cfg["outPw"])+".json"

        if os.path.isfile(path):
            dat_LED_spectrum = pd.read_json(path)
            
        tmp = dat_LED_spectrum.melt(id_vars = "lambda",var_name="LEDs",value_name="Radiance")

        plt.figure()
        sns.lineplot(x="lambda", y="Radiance", hue = "LEDs", data=tmp)
        plt.savefig("./output.pdf")

        figCount = 1
        

        # figCount+=2
        
        if spectrometer == "projector":
            path = self.filePath + spectrometer + "/LED_summary.json"
        
        if os.path.isfile(path):
            f = open(os.path.join(path))
            lut_spec = json.load(f)
            f.close()
    
        dr1 = np.unique(tmp["lambda"])
        
        plt.figure(figsize=(10, 16))
        for iLED in self.cfg['numOfLEDs']:
        
            plt.subplot(len(self.cfg['numOfLEDs'])+1, 3, figCount)
            
            plt.plot(dr1,np.array(lut_spec["LED"+str(iLED)])[np.arange(0,255,10),:].T)
            
            figCount+=1
            plt.subplot(len(self.cfg['numOfLEDs'])+1, 3, figCount)
            
            yy = []
            yy2 = []
            
            for irgb in np.arange(255):
                yy.append(sum(np.array(lut_spec["LED"+str(iLED)])[irgb,:].T*rod["rod"].values))
                yy2.append(sum(np.array(lut_spec["LED"+str(iLED)])[irgb,:].T*ipRGC['ipRGC'].values))
            plt.plot(yy)
            plt.plot(yy2)
            
            # plt.plot(np.arange(1,256),np.array(lut_spec["LED"+str(iLED)+"_peak"]))
            # plt.ylim([0,0.015])
            # plt.plot(np.array(lut_spec["LED1"]))
    
            # yy = []
            # for irgb in np.arange(1,256):
            #     # yy.append(np.mean(np.array(lut_spec["LED"+str(iLED)])[irgb-1,:] * rod['rod'].values))
            #     yy.append(sum((np.array(lut_spec["LED"+str(iLED)])[irgb-1,:] * ipRGC['ipRGC'].values)))
            # # yy = np.array(lut_spec["LED"+str(iLED)])[254,:] * rod['rod'].values
            
            
            # plt.plot(np.arange(1,256),yy)
            
            
            figCount+=1
            
            plt.subplot(len(self.cfg['numOfLEDs'])+1, 3, figCount)
        
            plt.plot(dr1,np.array(lut_spec["LED"+str(iLED)])[128,:])
            plt.plot(dr1,ipRGC['ipRGC'].values)
            plt.plot(dr1,np.array(lut_spec["LED"+str(iLED)])[128,:]*ipRGC['ipRGC'].values)
            # plt.ylim([0,0.01])
            # plt.plot(yy)
          
            figCount+=1
            
            
            
            

        # return tmp
        return lut_spec
        
    def getLight(self,LEDs,outPw):
        
        path = self.filePath + "projector/LED_summary.json"
        
        if os.path.isfile(path):
            f = open(os.path.join(path))
            lut_spec = json.load(f)
            f.close()
        
        out = []
        for iLED,pw in zip(LEDs,outPw):
            out.append(lut_spec["LED"+str(iLED)][pw])
       
        plt.plot(np.array(out).T)
        
        return out
        
    
    def getColor(self,light):
    
        winWaveLen = [380,780]
    
        # dat_LED_XYZ = LightSource.getXYZvalues()
        # dat_LED_spectrum = LightSource.getLEDspectra()
        cl_func = self.getXYZfunc()
        ipRGC   = self.getipRGCfunc()
        rod     = self.getRodfunc()
        Y_d65   = self.getD65Spectra()
        
        # dat_LED = {'spectrum':[],'XYZ':[],'Yxy':[],'LMS':[],'rod':[],'ipRGC':[]}
        tmp_df = pd.DataFrame()
    
        # for d in light:           
        # tmp_df['spectrum'] = light
        ############ XYZ values of 12LEDs ############
        
        k = 100/sum(Y_d65['Y'].values * cl_func['Y'].values)
        
        
        t = [k * sum(light * cl_func['X'].values),
             k * sum(light * cl_func['Y'].values),
             k * sum(light * cl_func['Z'].values)]
        
        tmp_df['ipRGC'] = sum(light * ipRGC['ipRGC'].values).reshape(-1)
        
        tmp_df['rod'] = sum(light * rod['rod'].values)
    
        tmp_df['X'] = t[0]
        tmp_df['Z'] = t[2]
        
        tmp_df['Y'] = t[1]
        tmp_df['x'] = t[0]/sum(t)
        tmp_df['y'] = t[1]/sum(t)
        
        tmp = [sum(light * cl_func['X']),sum(light * cl_func['Y']),sum(light * cl_func['Z'])]
        tmp = np.array(tmp)/tmp[1]
    
        lms = np.dot(np.array(xyzTolms),np.array(tmp))
        
        tmp_df['L'] = lms[0]
        tmp_df['M'] = lms[1]
        tmp_df['S'] = lms[2]
    
        tmp_df["L-M contrast"] = lms[0]-lms[1]
        tmp_df["S-(L+M) contrast"] = lms[2]-(lms[1]+lms[0])
        
        tmp_df["rod"] = 683*sum(light * rod['rod'])
        
        # tmp_df["l"] = light
        # tmp_df["dr"] = dr[(dr >= winWaveLen[0]) & (dr <= winWaveLen[1])]
        # tmp_df["light"] = light
        # tmp_df["No."] = i
        
        
        return tmp_df
   
        
def xyzTolab(xyz):
    
    def ft(t):
        if  t > (6/29)**3: 
            return t ** (1/3)
        else:
            return (((29/3)**3)*t + 16)/116
            
        
    lab = []
    for tmp_xyz in xyz:
        
        # Check xyz dimensions
        if len(tmp_xyz) != 3:
          print('Array xyz must have three rows')
        
        # Separate out the compontents
       
        X = tmp_xyz[0]
        Y = tmp_xyz[1]
        Z = tmp_xyz[2]
        
        Xw = XYZ_d50['X']
        Yw = XYZ_d50['Y']
        Zw = XYZ_d50['Z']
        
        # Allocate space
        tmp_lab = np.zeros(3)
    
        # Compute L
        tmp_lab[0] = 116 * ft(Y/Yw) - 16
        
        # Compute a* and b*
        tmp_lab[1] = 500 * (ft(X/Xw) - ft(Y/Yw))
        tmp_lab[2] = 200 * (ft(Y/Yw) - ft(Z/Zw))
     
        lab.append(tmp_lab)
        
    return lab
    
def labToxyz(lab):
   
    def ft(t,tw):
        if  fy > 6/29: 
            return (t**3) * tw
        else:
            return ((3/29)**3)*(116*t-16)*tw
            
    xyz=[]
    for tmp_lab in lab:
     
        Xw = XYZ_d50['X']
        Yw = XYZ_d50['Y']
        Zw = XYZ_d50['Z']
        
        fy = (tmp_lab[0]+16) / 116
        fx = fy + (tmp_lab[1]/500)
        fz = fy - (tmp_lab[2]/200)
        
        # Allocate space
        tmp_xyz = np.zeros(3)
        
        tmp_xyz[0] = ft(fx,Xw)
        tmp_xyz[1] = ft(fy,Yw)
        tmp_xyz[2] = ft(fz,Zw)
        
        xyz.append(tmp_xyz)
        
    return xyz
         
    
def showSpectra(filePath,saveFig=False):
    
   
    # for iLED in np.arange(1,self.cfg['numOfLEDs']+1):
    csv_input = pd.read_csv(filepath_or_buffer = filePath , encoding="ms932", sep=",")
    csv_input = csv_input.drop(['基準','白色板','レシオ分母'], axis=1)
    spectra = csv_input.drop(np.arange(24), axis=0)
    csv_input = csv_input.set_index('Unnamed: 0')
    
    lam = [int(t[:3]) for t in spectra["Unnamed: 0"].values]
  
    spectra = spectra.set_index('Unnamed: 0')
    spectra = np.array(spectra.astype(float, errors = 'raise')).T

    dr = csv_input.loc['メモ'].tolist()
   
    dat_spectrum = pd.DataFrame()
    dat_spectrum["y"] = spectra.reshape(-1)
    dat_spectrum["data_x"] = np.tile(lam, spectra.shape[0])
    dat_spectrum["condition"] = np.repeat(dr,spectra.shape[1])
    dat_spectrum["Yxy_x"] = np.repeat(np.array(csv_input.loc['x'].astype(float, errors = 'raise')),spectra.shape[1])
    dat_spectrum["Yxy_y"] = np.repeat(np.array(csv_input.loc['y'].astype(float, errors = 'raise')),spectra.shape[1])
    dat_spectrum["Yxy_Y"] = np.repeat(np.array(csv_input.loc['Y'].astype(float, errors = 'raise')),spectra.shape[1])
    
    # dat_spectrum["x"] = np.round(dat_spectrum["x"],3)
    # dat_spectrum["y"] = np.round(dat_spectrum["y"],3)
    
    #%%
    plt.figure()
    grid = sns.FacetGrid(dat_spectrum, col="condition", col_wrap=2, size=5)
    grid.map(sns.lineplot, "data_x", 'y', ci='none')
    plt.xlabel("lambda[nm]")
    plt.ylabel("radiance[W・sr-1・m-2]")
    
    if saveFig:
        plt.savefig("spectra_sep.pdf")
        
    #%%
    plt.figure()
    sns.lineplot(data=dat_spectrum, x="data_x", y='y', hue="condition", ci='none')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.xlabel("lambda[nm]")
    plt.ylabel("radiance[W・sr-1・m-2]")
    if saveFig:
        plt.savefig("spectra.pdf")

    #%%
    tmp_plot = dat_spectrum.groupby(["condition"],as_index=False).agg(np.nanmean)

    plt.figure()
    # sns.pointplot(data = tmp_plot, x="Yxy_x", y='Yxy_y', hue="condition", ci='none')
    for x,y,l in zip(tmp_plot["Yxy_x"].values,tmp_plot["Yxy_y"].values,tmp_plot["condition"].values):
        plt.plot(x,y,'o',label = l)
    
    plt.legend()
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    plt.xlim(0.3, 0.36)
    plt.ylim(0.3, 0.4)
    
    plt.xlabel("x")
    plt.ylabel("y")
    
    if saveFig:
        plt.savefig("xyY.pdf")

    #%% Y
    plt.figure()
    ax = sns.pointplot(data=dat_spectrum, x="condition", y='Yxy_Y',group = "condition")
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)
    if saveFig:
        plt.savefig("xyY_Y.pdf")

    return dat_spectrum


def drawLABColorBar(labimg, width, height, l):
    cx, cy = width // 2, height // 2

    for x in range(0, width):
        for y in range(0, height):
            dist = math.sqrt((x - cx) ** 2 + (y - cy) ** 2)
            if dist <= 128:
                a = x - cx
                b = (y - cy)*-1
                labimg[y, x] = (l, a, b)
                #print('(%d, %d) %.2f %.2f' % ( x, y, a, b))

    return  labimg


