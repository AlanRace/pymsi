#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 11:04:05 2020

@author: Alan Race
"""

import matplotlib.pyplot as plt
import math
import pymsi
import os
import sys
import inspect
import glob
currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)

datadir = os.path.join(currentdir, 'data')

# Make sure that the package is on the path
sys.path.insert(0, parentdir)


# Get a list of all test files that we have downloaded using downloadFromMetaspace.py
imzMLFiles = glob.glob(os.path.join(datadir, '*.imzML'))

# Load the first dataset, trimming empty pixels if necessary
msiData = pymsi.ImzMLHandler(imzMLFiles[2], cropToData=True)

ticImage = msiData.getTICImage()

plt.imshow(ticImage)

minMZ, maxMZ = msiData.determineMinMaxMZ(pixelsToSample=100)

# Round down and up to nearest 50 m/z to get rough idea of min and max m/z
mzEstimateAccuracy = 50

minMZ = math.floor(minMZ / mzEstimateAccuracy) * mzEstimateAccuracy
maxMZ = math.ceil(maxMZ / mzEstimateAccuracy) * mzEstimateAccuracy

print(minMZ)
print(maxMZ)

#msiData.estimatePPM(pixelsToSample = 100)


# %%

ppmEstimates = msiData.estimatePPM(minMZ, maxMZ, pixelsToSample=1000)

#ppmEstimate = np.mean(ppmEstimates)
ppmEstimate = np.median(ppmEstimates)

print(ppmEstimate)

# %%


ticSpectrum = msiData.generateMeanSpectrum(minMZ, maxMZ, ppmEstimate)

# %%
plt.plot(msiData.mzAxis[1:], ticSpectrum)
plt.xlim(885, 886)

#%%


