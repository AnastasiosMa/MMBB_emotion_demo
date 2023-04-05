#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 15:30:31 2022

@author: anmavrol
"""
from pydub import AudioSegment

import os

dl = os.listdir(path='trimmed_part2_wav')
for i in dl:
    if '.wav' in i:
        print(i)
        AudioSegment.from_wav('trimmed_part2_wav/'+ i).export('trimmed_part2_mp3/' + i[0:3] +'.mp3', format="mp3")

