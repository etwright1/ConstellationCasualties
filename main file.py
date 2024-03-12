# 
# Code for analysing Starlink gen2 casualty risk
# by etwright1 (E. Wright), 2022
#

import casualty as cs
import numpy as np
import matplotlib.pylab as plt

RE=6378e3 # equatorial radius in m

DCA = 0.36 #debris casualty area assuming minimal lethal debris 

constellation_inclinationss = [53, 46, 38, 96.9, 53, 43, 33, 148, 115.7] #define constellation
inclination_weightings = [0.176, 0.176, 0.176, 0.120, 0.112, 0.112, 0.112, 0.005, 0.011] #relative weightings of each inclination
modelled_satellites = 300
multiplier = 100 #multiplies result by factor to faster process larger constellations
total_satellites = modelled_satellites * multiplier #to check 

satellite_inclinations = cs.create_constellation(constellation_inclinationss, inclination_weightings, modelled_satellites)

latitudes = list(reversed(np.arange(-90, 90, 0.5))) #create latitudes list
weighting_function = np.zeros(360)

casualty_expectation_file = "Constellation casualty expectation" #naming files
weighting_function_file = casualty_expectation_file + " weighting function"

num_of_satellites = len(satellite_inclinations)
timer = 0
limit = 0

#create weighting function, satellite by satellite

for x in satellite_inclinations:
    if x > limit:
        if x > 90:
            x = 180 - x

        vals, lats = cs.latWeights(0.5, 550e3+RE, x) #get latitude weights

        weighting_function += vals #add the normalised times

        timer +=1
        print("Working on satellite:", timer, "of", num_of_satellites)

np.savetxt(weighting_function_file, weighting_function*multiplier) #save latweights function for other plots

#find casualty expectation of that weighting function

GPW4 = np.genfromtxt("gpw_v4_population_count_rev11_2020_30_min_clean.csv", delimiter=' ', dtype=float) #import and clean GPW4 population data
GPW4 = np.where(GPW4==-9999, 0, GPW4)

print("Sanity check: Population of Earth in 2020= ", np.sum(np.sum(GPW4)))

GPW4 = GPW4 * 1.1 #annual increase of 1% over 10 years, approximately linear

print("Sanity check: Population of Earth in 2030= ", np.sum(np.sum(GPW4)))

two_dimensional_casualty_exp = np.zeros((360,720))

for x in range(720): # across all longitudes
    two_dimensional_casualty_exp[:,x] = GPW4[:,x] * weighting_function * DCA * multiplier # Two dimensional, this is a slightly convoluted method. Be careful with multiplier

casualty_by_latitude = two_dimensional_casualty_exp.sum(axis=1) #back to one dimension

area_of_latitude_band = np.zeros(len(latitudes))
casualty_expectation = np.zeros(len(latitudes))

for x in range(len(latitudes)): #dividing by latitude band area to get casualty expectation (i.e., from population density)
    upper_latitude = abs(90 - x/2) 
    lower_latitude = abs(90 - x/2 + 0.5)
    area_of_latitude_band[x] = 2 * np.pi * RE**2 * abs((np.sin(np.deg2rad(upper_latitude)) - np.sin(np.deg2rad(lower_latitude))))
    casualty_expectation[x] =  casualty_by_latitude[x] / area_of_latitude_band[x] 

print('Sanity check: Earth area = ', sum(area_of_latitude_band))
print("Total casualty expectation for constellation is", sum(casualty_expectation))

np.savetxt(casualty_expectation_file, casualty_expectation)

plt.figure()
plt.plot(latitudes, casualty_expectation, color='g')
plt.xlabel("Latitude (deg)")
plt.ylabel("Casualty expectation (#)")
plt.xticks(np.arange(-90, 110, 20))
plt.show()