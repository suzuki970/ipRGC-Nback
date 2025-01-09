#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:21:28 2024

@author: yutasuzuki
"""

import glob
import os
import subprocess
import numpy as np
import ffmpeg
import shutil

#%% convert to mp4
iVideos="/Users/yutasuzuki/Desktop/Screen Recording 2024-01-18 at 14.37.16.mov"
print("Processing... "+iVideos)

(
 ffmpeg
 .input(iVideos)
 .filter('scale', 600, -2)
 .output(iVideos[:-4]+".mp4",crf=18,pix_fmt="yuv420p")
 .run()
)
