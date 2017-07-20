#!/usr/bin/env python

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#load python modules used
import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import stats
import sys

#set up command line interface
parser = argparse.ArgumentParser(description='Detects peaks in data. It does this by first calculating an estimate for the samples continuous distribution and then fitting that one with a small number of gaussians. The peaks are then taken as the gaussians\' means. The user can either specify the number of peaks using the -p argument, or use automatic detection. Parameters for the automatic detection can be set using the -m and -t arguments.', epilog='Original version written by Georg Rempfer <georg.rempfer@icp.uni-stuttgart.de>, Institute for Computational Physics, University of Stuttgart in 2017 at the 648. WE-Heraeus-Seminar in Bremen. Greetings to all the nanopore people (:')
parser.add_argument('datafile', help='Path to input file. Data should be formatted as text columns separated by whitespace.')
parser.add_argument('--column', '-c', default=3, type=int, help='Column to use as samples from the datafile. Starts counting at 0. Default is 3 (the k column in the files this program was originally written for).')
parser.add_argument('--peaks', '-p', default=0, type=int, help='Number of peaks expected in the samples. By default the program tries to determine a minimum count of peaks to fit the data.')
parser.add_argument('--max-peaks', '-m', default=6, type=int, help='Maximum number of peaks to try to fit in autodetection mode. Default is 6.')
parser.add_argument('--fit-tolerance', '-t', default=0.04, type=int, help='Maximum fit error tolerated. The error is calculated as the maximum over all x values of abs( (kde(x)-fit(x)) / mean(kde) ). You can play with this value to make the autodetection choose the right number of peaks. Default is 0.04.')
args = parser.parse_args()

#read all data from file
data = np.genfromtxt(args.datafile, skip_header=1)
samples = data[:,args.column] #select chosen column
samples = samples[~np.isnan(samples)] #remove samples that are not numbers (shitty measurements denoted by placeholders)
samples = np.sort(samples)

#create a grid with 100 x values in the data range
x_grid = np.linspace(np.min(samples), np.max(samples), 100)

#calculate estimate of the samples continuous distribution
kernel = stats.gaussian_kde(samples)
  
#function that can superimpose an arbitraty number of gaussians
#first argument is the x-value t be evaluated, remaining arguments are height, x-position, and width of first, second, ... guassian
def multi_gauss(x, *args):
  n = len(args)/3
  y = 0.0
  for i in range(n):
    y += args[i*3+0]*np.exp(-((x-args[i*3+1])/args[i*3+2])**2)
  return y

if args.peaks != 0: #if number of peaks given by user
  n = args.peaks

  #initial guesses for the gaussian's parameters
  #these work for the data I was given
  height = np.max(kernel(x_grid))
  width = (np.max(samples) - np.min(samples))/n
  std = np.std(samples)/n

  p = np.empty(0)
  for i in range(n):
    p = np.append(p, [height, width*(0.5+i), std])
  
  #carry out the fit
  popt, pcov = sp.optimize.curve_fit(multi_gauss, x_grid, kernel(x_grid), p)
else: #try to autodetect number of peaks
  n = np.arange(1, args.max_peaks+1)
  parameters = []
  errors = np.array([np.sum((kernel(x_grid))**2)])

  #try number of peaks from 1 to given maximum
  for k in n:
    #same intial guesses as before
    height = np.max(kernel(x_grid))
    width = (np.max(samples) - np.min(samples))/k
    std = np.std(samples)/k

    p = np.empty(0)
    for i in range(k):
      p = np.append(p, [height, width*(0.5+i), std])
    
    #carry out the fit and don't die if it fails. also calculate a fit error
    try:
      popt, pcov = sp.optimize.curve_fit(multi_gauss, x_grid, kernel(x_grid), p)
      perr = np.max( np.abs( (kernel(x_grid)-multi_gauss(x_grid, *popt))/np.mean(kernel(x_grid)) ) )
    except:
      parameters.append(None)
      errors = np.append(errors, np.nan)

    #if fit error is smaller than specified, stop and use the current number of peaks
    if perr <= args.fit_tolerance:
      n = k
      break

#sort peaks
popt = popt.reshape([n,3])
popt = np.array(sorted(popt, key=lambda row: row[1]))
popt = popt.reshape(3*n)

#plot raw samples, kernel density estimation, and multi gaussian fit
plt.plot(x_grid, multi_gauss(x_grid, *popt), 'r--', linewidth=3.0, label='fit of {} gaussians'.format(n))
plt.plot(x_grid, kernel(x_grid), 'b-', label='kernel density estimation')
plt.plot(samples, (1.0-0.681)*height*np.ones_like(samples), 'g|', label='raw samples')

#plot individual gaussians, mark their peaks with vertical lines, and annotate the corresponding x-values
for i in range(n):
  plt.plot(x_grid, multi_gauss(x_grid, *popt[3*i:3*i+3]), 'r-')
  plt.axvline(x=popt[i*3+1], color='k')
  plt.annotate('{:.2e}'.format(popt[3*i+1]), xy=(popt[i*3+1], 0.618*height-0.1*height*(i%2)), fontweight='bold')
  print('Peak {} at {:.4e}'.format(i+1, popt[i*3+1]))

#show plot with legend
plt.legend()
plt.show()
