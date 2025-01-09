#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 16:02:37 2024

@author: yutasuzuki
"""

import shutil
import glob
import os
import re
from multiprocessing import Pool

#%% path setting

zipFlg = True
zipFlg = False

if zipFlg:
    dicomFolderList = glob.glob("/Volumes/taman2/DICOM/*")
    saveFolder = r"/Users/yutasuzuki/My Drive/DICOM/"
    
else:
    # dicomFolderList = glob.glob("G:\My Drive\DICOM/*")
    dicomFolderList = glob.glob("G:\マイドライブ\DICOM/*")
    unpackFolder = "C:/Users/NTT/Desktop/DICOM/"
        
#%%

def dozip(f):
    os.makedirs(saveFolder+re.split("/|\\\|\|", f)[-2],exist_ok=True)
    shutil.make_archive(saveFolder+re.split("/|\\\|\|", f)[-2]+'/'+re.split("/|\\\|\|", f)[-1], format='zip', root_dir=f)

def unzip(f):
    os.makedirs(unpackFolder+re.split("/|\\\|\|", f)[-2],exist_ok=True)
    shutil.unpack_archive(f, unpackFolder+re.split("/|\\\|\|", f)[-2]+'\\'+re.split("/|\\\|\|", f)[-1][:-4])

#%%
if __name__ == '__main__':
    
    fileList = []
    if zipFlg:
        for folder in dicomFolderList:
            for f in glob.glob(folder+"/*[!.zip]"):
                fileList.append(f)

        with Pool(4) as p:
            tmp_datHash = p.map(dozip,fileList)

    else:
        for folder in dicomFolderList:
            for f in glob.glob(folder+"/*.zip"):
                fileList.append(f)

        with Pool(4) as p:
            tmp_datHash = p.map(unzip,fileList)

