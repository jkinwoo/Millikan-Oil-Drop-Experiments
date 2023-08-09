#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib as mpl # V
import numpy as np # V
import pandas as pd # V
import pims # V
import trackpy as tp # V
from slicerator import pipeline # V

# Video processing (convert to to grayscale)
@pipeline
def as_grey(frame):
    red = frame[:, :, 0]
    green = frame[:, :, 1]
    blue = frame[:, :, 2]
    return 0.2125 * red + 0.7154 * green + 0.0721 * blue

mpl.rc('figure',  figsize=(10, 5))
mpl.rc('image', cmap='gray')

def compute_traj(filename):
    
    # Start video preprocessing
    vid = pims.Video('../test_video/' + filename + '.avi')
    frames = as_grey(vid)

    # Declare frames of which to view
    # midpoint = len(frames)/2
    start = 0
    stop = len(frames)
    
    # Segmentation parameters
    particle_diameter = 23
    min_brightness = 200
    gyration_brightness = 8.0
    py_engine = 'numba'
    # If on windows OS set to 1 to disable multiprocessing
    windows_os = 1
    
    # Output plot of segmentation
    f = tp.locate(frames[0], diameter = particle_diameter, separation = particle_diameter*1.25, maxsize = gyration_brightness, minmass = min_brightness)
    tp.annotate(f, frames[0])
    
    # Check if user is happy with segmentation
    print("Continue with data collection? (y/n)")
    continue_bool = input()
    if(continue_bool == "y"):
            
        # Trajectory segmentation
        f = tp.batch(frames[start:stop], diameter = particle_diameter, separation = particle_diameter*1.5, maxsize = gyration_brightness, minmass = min_brightness, engine = py_engine, processes = windows_os)
        t = tp.link_df(f, 5, memory=3)
        
        t1 = tp.filter_stubs(t, 60)
        # Compare the number of particles in the unfiltered and filtered data.
        print('Unfiltered counts:', t['particle'].nunique())
        print('Filtered counts:', t1['particle'].nunique())
        
        # Track pixel displacement and calculate velocity
        data = []
        for item in set(t1.particle):
            sub = t1[t1.particle==item]
            dvx = np.diff(sub.x)
            dvy = np.diff(sub.y)
            for x, y, dx, dy, frame in zip(sub.x[:-1], sub.y[:-1], dvx, dvy, sub.frame[:-1],):
                data.append({   'dx': dx, 
                                'dy': dy, 
                                'x': x,
                                'y': y,
                                'frame': frame,
                                'particle': item,
                    })
        df = pd.DataFrame(data)

        # Export to csv
        df.to_csv('../csvs/'+ filename +'_particles.csv') 

# Windows multiprocessing check
if __name__ == "__main__":
    # print(os.getcwd())
    # Declare video file name
    compute_traj('1-1')
