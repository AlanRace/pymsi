#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 10:51:06 2020

@author: Alan Race
"""

import os
import urllib.request
from metaspace import SMInstance

scriptLocation = os.path.abspath(__file__)
testFolder = os.path.dirname(scriptLocation)

sm = SMInstance()

#dsid = "2016-09-22_11h16m17s"
#dsid = "2020-10-19_08h32m13s"
#dsid = "2019-11-12_17h15m55s"
dsid = "2020-08-11_14h50m21s"
ds = sm.dataset(id=dsid)

downloadLinks = ds.download_links()


imzMLFile = downloadLinks['files'][0]['link']
ibdFile = downloadLinks['files'][1]['link']

urllib.request.urlretrieve(imzMLFile, os.path.join(testFolder, 'data', imzMLFile.split('/')[-1].split('?')[0]))
urllib.request.urlretrieve(ibdFile, os.path.join(testFolder, 'data', ibdFile.split('/')[-1].split('?')[0]))
