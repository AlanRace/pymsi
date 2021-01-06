#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 16:57:38 2020

@author: Alan Race
"""

from pyimzml.ImzMLParser import ImzMLParser
from sklearn.cluster import AgglomerativeClustering
import math
import numpy as np
import random

class ImzMLHandler:

    # cropToData: (if True) Remove all rows and columns which contain no data. For some datasets (e.g. Bruker) 
    # the coordinates stored are relative to some external coordinate system, and therefore large
    # amounts of empty space can be present
    def __init__(self, filename, startX=1, startY=1, width=None, height=None, cropToData=False):
        self.imzML = ImzMLParser(filename)
        
        # Find the min and max row and column where data is present
        maxWidth = 0
        maxHeight = 0
        
        minWidth = -1
        minHeight = -1
        
        for (x, y, z) in self.imzML.coordinates:
            if x > maxWidth:
                maxWidth = x
            if y > maxHeight:
                maxHeight = y
            if minWidth == -1 or minWidth > x:
                minWidth = x
            if minHeight == -1 or minHeight > y:
                minHeight = y
        
        if cropToData:
            startX = minWidth
            startY = minHeight
        
        if width is None:
            width = maxWidth-startX+1
        if height is None:
            height = maxHeight-startY+1
            
        self.startX = startX
        self.startY = startY
        self.width = width
        self.height = height
        self.coordinates = []
        
        index = 0
            
        for (x, y, z) in self.imzML.coordinates:
            if x >= startX and y >= startY and x < (startX+width) and y < (startY+height):
                if cropToData:
                    self.coordinates.append((index, x-minWidth+1, y-minHeight+1))
                else:
                    self.coordinates.append((index, x, y))
            index = index + 1
                
    def getTICImage(self):
        ticImage = np.zeros((self.height, self.width)) 
            
        for index, x, y in self.coordinates:
            mzs, counts = self.imzML.getspectrum(index)
         
            #(x, y, z) = imzML.coordinates[index]
            
            ticImage[y-self.startY, x-self.startX] = np.sum(counts)
            
        return ticImage
    
    def determineMinMaxMZ(self, pixelsToSample=100):
        # TODO: Check in the metadata
        
        # Alternatively, sample some pixels and see what the min and max recorded 
        # m/z values are
        minMZ = -1
        maxMZ = 0
        
        for i in range(pixelsToSample):
            spectrumToSample = random.randint(0, len(self.coordinates)-1)
            (index, x, y) = self.coordinates[spectrumToSample]
            
            mzs, counts = self.imzML.getspectrum(index)
            
            if minMZ == -1 or mzs[0] < minMZ:
                minMZ = mzs[0]
            if maxMZ < mzs[len(mzs)-1]:
                maxMZ = mzs[len(mzs)-1]
        
        return minMZ, maxMZ
        
    def estimatePPM(self, minMZ, maxMZ, numBins=10, pixelsToSample=100):
        ppmEstimates = np.ones(numBins) * 1e5
        
        for i in range(pixelsToSample):
            spectrumToSample = random.randint(0, len(self.coordinates)-1)
            (index, x, y) = self.coordinates[spectrumToSample]
            
            mzs, counts = self.imzML.getspectrum(index)
            
            diff = mzs[1:len(mzs)] - mzs[0:len(mzs)-1]
            ppms = diff*1e6/mzs[0:len(mzs)-1]
            
            binWidth = (maxMZ - minMZ) / numBins
            
            for binNum in range(numBins):
                startMZ = minMZ + (binNum * binWidth)
                endMZ = minMZ + ((binNum + 1) * binWidth)
                
                possiblePPMs = ppms[np.logical_and(mzs[0:len(mzs)-1] >= startMZ, mzs[0:len(mzs)-1] < endMZ)]
                
                if len(possiblePPMs) > 0:
                    ppmEstimate = np.min(possiblePPMs)
                    
                    if ppmEstimates[binNum] > ppmEstimate:
                        ppmEstimates[binNum] = ppmEstimate
        
        return ppmEstimates
    
    def generateMeanSpectrum(self, startmz, endmz, ppm):
        self.mzAxis = ImzMLHandler.generateMZAxis(startmz, endmz, ppm)
        
        spectrum = np.zeros((self.mzAxis.shape[0]-1))
        
        startLog = np.log(self.mzAxis[0])
        ppmLog = np.log(1 + ppm * 1e-6)
        
        for index, x, y in self.coordinates:
            if index % 10 == 0:
                mzs, counts = self.imzML.getspectrum(index)
                
                for mzIndex in range(len(mzs)):
                    location = int(np.round((np.log(mzs[mzIndex]) - startLog) / ppmLog))
                    
                    if location < 0: 
                        continue
                    
                    if location >= len(spectrum):
                        break
                    
                    spectrum[location] += counts[mzIndex]
            
        self.meanSpectrum = spectrum / len(self.coordinates)
            
        return self.meanSpectrum
    
    def generateIonImage(self, mz, ppm):
        ionImage = np.zeros((self.height, self.width))
        
        deltamz = ppm * 1e-6 * mz
        minmz = mz - deltamz
        maxmz = mz + deltamz
    
        for index, x, y in self.coordinates:
            mzs, counts = self.imzML.getspectrum(index)
            
            for mzIndex in range(len(mzs)):
                if mzs[mzIndex] > maxmz:
                    break
                
                if mzs[mzIndex] >= minmz and mzs[mzIndex] <= maxmz:
                    ionImage[y-1, x-1] += counts[mzIndex]
                    
        return ionImage
    
    def generateDatacubeMZs(self, limits, ticNorm=False):
        datacube = np.zeros((len(self.coordinates), len(limits))) 
        
        spectrumIndex = 0
        
        for index, x, y in self.coordinates:
            mzs, counts = self.imzML.getspectrum(index)
            
            # Normalised to TIC
            if ticNorm: 
                counts = counts/ np.sum(counts)
            
            for l in range(len(limits)):
                datacube[spectrumIndex, l] = np.sum(counts[np.logical_and(mzs > limits[l, 0], mzs <= limits[l, 1])])
            
                    
            spectrumIndex += 1
        
        self.datacube = datacube
        
        return self.datacube
    
    def generateDatacube(self, peaks, left_ips, right_ips, ticNorm=False):
        #left_ips = peakProperties['left_ips']
        left_ips = np.floor(left_ips).astype(np.int)-1
        #right_ips = peakProperties['right_ips']
        right_ips = np.ceil(right_ips).astype(np.int)+1
                
        datacube = np.zeros((len(self.coordinates), len(peaks))) 
        
        spectrumIndex = 0
        
        for index, x, y in self.coordinates:
            mzs, counts = self.imzML.getspectrum(index)
            
            # Normalised to TIC
            if ticNorm: 
                counts = counts/ np.sum(counts)
            
            curPeakIndex = 0
            
            for mzIndex in range(len(mzs)):
                while curPeakIndex < len(peaks) and mzs[mzIndex] > self.mzAxis[right_ips[curPeakIndex]]:
                    curPeakIndex += 1
                if curPeakIndex >= len(peaks):
                    break
                    
                for peakIndex in range(curPeakIndex, len(peaks)):
                    if mzs[mzIndex] < self.mzAxis[left_ips[peakIndex]]:
                        break
                    if mzs[mzIndex] >= self.mzAxis[left_ips[peakIndex]] and mzs[mzIndex] <= self.mzAxis[right_ips[peakIndex]]:
                        datacube[spectrumIndex, peakIndex] += counts[mzIndex]
                        break
                    
            spectrumIndex += 1
        
        self.datacube = datacube
        
        return self.datacube
    
    def determineCorrelatedFeatures(self, clusteringThreshold):
        ionCorrelationMatrix = np.corrcoef(self.datacube.transpose())
        ionCorrelationMatrix[np.isnan(ionCorrelationMatrix)] = 0
        
        clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=clusteringThreshold).fit(ionCorrelationMatrix)
        
        u, c = np.unique(clustering.labels_, return_counts=True)
        
        #nonUniqueClusters = np.where(c > 1)[0]
        
        self.uniqueData = []
        self.uniqueDataMembers = []
        
        for index in u:
            clusterIndex = u[index]
            #print(clusterIndex)
            
            clusterMembers = np.where(clustering.labels_ == clusterIndex)[0]
            
            self.uniqueDataMembers.append(clusterMembers)
            
            averageIonImage = np.mean(self.datacube[:, clusterMembers], axis=1)
            averageIonImage = np.reshape(averageIonImage, (len(averageIonImage),1))
            
            if self.uniqueData == []:
                self.uniqueData = averageIonImage
            else:
                self.uniqueData = np.concatenate((self.uniqueData, averageIonImage), axis=1)
                
        return self.uniqueData
    
    @staticmethod
    def generateMZAxis(startmz, endmz, ppm):
        
        numElements = int(np.ceil((np.log(endmz) - np.log(startmz)) / np.log(1 + ppm * 1e-6)))
        
        mzAxis = np.zeros((numElements))
    
        for i in range(numElements):
            mzAxis[i] = startmz * np.power(1+ppm*1e-6, i)
            
        return mzAxis