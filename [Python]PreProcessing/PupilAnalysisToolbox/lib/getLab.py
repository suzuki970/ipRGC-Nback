#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 10:59:34 2024

@author: yutasuzuki
"""

import numpy as np

XYZ_MATRIX = np.array([
             [0.4124, 0.2126, 0.0193],
             [0.3576, 0.7152, 0.1192],
             [0.1805, 0.0722, 0.9505]
            ])

# sRGB -> XYZ
def srgb2xyz(img):
    rgb = (img/255) ** 2.2 #リニアsRGBに戻す
    xyz = np.dot(rgb, XYZ_MATRIX)
    return xyz

# f(x)
def trans_function(x, threshold=0.008856):
    y = x**(1/3) if x>threshold else (841/108)*x+(16/116)
    return y

# XYZ -> L*a*b*
def xyz2lab(x, y, z, xn=0.9505, yn=1.0, zn=1.089):
    l = 116 * trans_function(x/yn) - 16
    a = 500 * ( trans_function(x/xn) - trans_function(y/yn) )
    b = 200 * ( trans_function(y/yn) - trans_function(z/zn) )
    return l, a, b

# xyz2lab関数を画像に適用するため
def trans_matrix(img, func):
    h,w,ch = img.shape
    new = np.array([func(*img[y,x,:]) for y in range(h) for x in range(w)])
    new = np.reshape(new, (h, w, ch))
    return new

# RGB画像をL*a*b*色空間に変換
def img2lab(img):
    xyz = srgb2xyz(img)
    lab = trans_matrix(xyz, xyz2lab)
    return lab

RGB_MATRIX = np.array([
            [3.2406, -0.9689, 0.0557],
            [-1.5372, 1.8758, -0.2040],
            [-0.4986, 0.0415, 1.0570]
            ])

def inverse_function(y, n, threshold=0.206893):
    x = y**3 if y>threshold else (y-(16/116))*(108/841)
    return n*x

def lab2xyz(l, a, b, xn=0.9505, yn=1.0, zn=1.089):
    x = inverse_function(((a/500)+((l+16)/116)), xn)
    y = inverse_function((l+16)/116, yn)
    z = inverse_function(((l+16)/116)-(b/200), zn)
    return x, y, z

def xyz2srgb(xyz):
    img = np.dot(xyz, RGB_MATRIX)
    img = img ** (1/2.2)
    img *= 255
    return img.astype('uint8')

def trans_matrix(img, func):
    h,w,ch = img.shape
    new = np.array([func(*img[y,x,:]) for y in range(h) for x in range(w)])
    new = np.reshape(new, (h, w, ch))
    return new

def lab2img(lab):
    xyz = trans_matrix(lab, lab2xyz)
    img = xyz2srgb(xyz)
    return img
