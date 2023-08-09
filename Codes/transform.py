#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import pandas as pd
import matplotlib.pyplot as plt
import trackpy as tp

# Datasheet to be transformed:
filename = "6-2"

# Constants
PI = math.pi
G = 9.8067 # (m/s^2)
PLATE_DISTANCE = 3.8965E-3 # (m)
DENSITY_OIL = 860 # (kg/m^3)
B = 6.17E-6 # correction factor
N = 1.825E-5 # Viscosity of air at 20 deg_C (N*s/m^2)
CHARGE_REAL = 1.60217663E-19 # coulombs

# Independant variables
voltage = 450.0 # Voltage (V)
p = 73.818545 # Barometric Pressure (cm/Hg)

# Calibration constants
frame_rate = 30.893 # (fps)
distance_scale = 0.4259259259E-03 # Microscope increments in m
pixel_scale = 4.47047253E-06 # (m/px)

def plotTraj(i):
    particle = df.loc[df['particle'] == i]    
    plt.plot(particle['frame'], particle['x'])    
    turning_index = particle['y'].idxmin()    
    bef = particle[:turning_index]
    aft = particle[turning_index:]    
    print(turning_index)    
    return bef, aft
    
df = pd.read_csv('../csvs/' + filename + '_particles.csv')
x = tp.compute_drift(df)
x.plot()

#particles = set(df['particle'])
#
#for i in particles:
#    plotTraj(i)

# Only consider particles which change direction within X frames from the midpoint of the video
def returnParticles(middle):
    filtered = []
    turning_indices = []
    for i in pivoted_y:
        turning_index = pivoted_y[i].idxmin()
        if abs(turning_index - middle) < 10:
            filtered.append(i)
            turning_indices.append(turning_index)
    return filtered, turning_indices

def filterParticles():
    turning_indices = []
    for i in pivoted_y:
        turning_index = pivoted_y[i].idxmin()
        turning_indices.append(turning_index)
    middle = int(sum(turning_indices)/len(turning_indices))
    return middle

def returnFinal():
    final = []
    for i,j in zip(filtered_particles, turning_indices):
        particle = pivoted_dy[i][pd.Series.notnull(pivoted_dy[i])]
        bef = particle[:j-particle.index[0]]
        aft = particle[j-particle.index[0]:]
        if len(bef) > 0 and len(aft) > 0:
            v_g = bef.sum()/len(bef)
            v_e = aft.sum()/len(aft)
            ve_fin = abs(v_e*frame_rate*pixel_scale)
            vg_fin = abs(v_g*frame_rate*pixel_scale)
            radius = math.sqrt((9 * N * vg_fin) / (2 * G * DENSITY_OIL))
            charge = ((6 * PI * PLATE_DISTANCE) / voltage)*math.sqrt((9 * pow(N, 3)) / (2 * DENSITY_OIL * G))*pow(1 + (B / (p * radius)), -3/2)*((ve_fin + vg_fin) * math.sqrt(vg_fin))
            quantization = charge / CHARGE_REAL
            final.append({'particle': i, 
                          'turning_point': j, 
                          'v_e (px/frame)': v_e, 
                          'v_g (px/frame)': v_g,
                          'v_e (m/s)': ve_fin, 
                          'v_g (m/s)': vg_fin,
                          'radius': radius,  
                          'charge': charge,
                          'quantization': quantization,
                        })
    return final

pivoted = df.pivot(index = 'frame', columns = 'particle')
pivoted_y = pivoted['y']
pivoted_dy = pivoted['dy']

filtered_particles, turning_indices = returnParticles(250)
print (len(filtered_particles))
#print(os.getcwd())
final = returnFinal()
df = pd.DataFrame(final)
df.to_csv('../csvs/'+ filename +'_transformed.csv')