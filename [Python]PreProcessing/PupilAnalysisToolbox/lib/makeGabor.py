#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 11:17:57 2023

@author: yutasuzuki
"""

from matplotlib import pyplot as plt
import numpy as np
# from scipy.misc import toimage  # One can also use matplotlib*

def gabor_fn(sigma,theta,Lambda,psi,gamma):
    
    sigma_x = sigma
    sigma_y = float(sigma)/gamma

    # Bounding box
    nstds = 3
    xmax = max(abs(nstds*sigma_x*np.cos(theta)),abs(nstds*sigma_y*np.sin(theta)))
    xmax = np.ceil(max(1,xmax))
    ymax = max(abs(nstds*sigma_x*np.sin(theta)),abs(nstds*sigma_y*np.cos(theta)))
    ymax = np.ceil(max(1,ymax))
    
    xmin = -xmax
    ymin = -ymax
    
    (x,y) = np.meshgrid(np.arange(xmin,xmax+1),np.arange(ymin,ymax+1 ))
    (y,x) = np.meshgrid(np.arange(ymin,ymax+1),np.arange(xmin,xmax+1 ))

    # Rotation
    x_theta=x*np.cos(theta)+y*np.sin(theta)
    y_theta=-x*np.sin(theta)+y*np.cos(theta)

    gb= np.exp(-.5*(x_theta**2/sigma_x**2+y_theta**2/sigma_y**2))*np.cos(2*np.pi/Lambda*x_theta+psi)
    
    return gb

rgb=np.array([1.0,0.27075,0.5603])*180

data = gabor_fn(sigma=12.,theta=np.pi/2.,Lambda=12.5,psi=90,gamma=1.)+1

data0 = gabor_fn(sigma=12.,theta=np.pi*(1/2),Lambda=12.5,psi=90,gamma=1.)+1
data45 = gabor_fn(sigma=12.,theta=np.pi*(1/2),Lambda=12.5,psi=90,gamma=1.)+1

I = np.zeros((data45.shape[0],data45.shape[0],3))

for irgb in np.arange(I.shape[2]):
    I[:,:,irgb] = data45 * rgb[irgb]

plt.imshow(I/255)



I = np.zeros((data.shape[0],data.shape[0],3))

for irgb in np.arange(I.shape[2]):
    I[:,:,irgb] = data * rgb[irgb]

plt.imshow(I/255)

# plt.imshow(data, 'gray', interpolation='none')
# plt.imshow(data/255,interpolation='none')
plt.axis('off')

plt.savefig("gabor.pdf")



