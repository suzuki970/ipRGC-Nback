#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 11:11:41 2024

@author: yutasuzuki
"""

import numpy as np
import glob
import os
from scipy.signal import find_peaks
import cv2
import pandas as pd
import time
import skvideo
skvideo.setFFmpegPath('/opt/anaconda3/bin/')
import skvideo.io
import matplotlib.pyplot as plt
import datetime
import json
from pre_processing_cls import pre_processing,rejectDat,rejectedByOutlier,re_sampling,getNearestValue
import shutil


from natsort import natsorted
from PIL import Image

from pathlib import Path

def makeMovie(cfg,folderName,outMovName = "/test.mp4"):
    
    figFileName = glob.glob(folderName+"/*"+cfg["ext"])
    figFileName = natsorted(figFileName)
    
    output = []
    for imgName in figFileName:
        output.append(np.array(Image.open(imgName))[:,:,[2,1,0]])
     
    # Create a blank image for the output video
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    
    width = np.array(Image.open(imgName)).shape[1]
    height = np.array(Image.open(imgName)).shape[0]
    
    out = cv2.VideoWriter(str(Path(folderName).parent)+outMovName, 
                          fourcc, cfg["fps"], (int(width),int(height)))
    
    for h in output:
        out.write(h)
        
    out.release()
    print("saved at " + str(Path(folderName).parent)+outMovName)

