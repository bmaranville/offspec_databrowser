#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# check to see if all parameters are set on the command line
# if they are, then don't open GUI interface
from math import *
import numpy
import sys
#from Tkinter import *
from wx import MessageBox
from wx import FileDialog
import datetime
#from FileDialog import *
#from reflectometry.reduction import *
from reflectometry.reduction import load as icp_load
from numpy import *
from pylab import imshow,cm,colorbar,hot,show,xlabel,ylabel, plot
#from qxqz_data_array_class import *
from matplotlib.widgets import  RectangleSelector
import os, pickle
from binned_data_class3 import plottable_2d_data, offspec_plot_2d_data
from background_dialog import BkgDialog

import __main__
#global rebinned_data_objects
#rebinned_data_objects = []
import wx
app = wx.App()
app.MainLoop()

FLOAT_ERR = 1.0e-17

__main__._rebinned_data_objects = []

class Supervisor():
  """ class to hold rebinned_data objects and increment their reference count """
  def __init__(self):
    self.rb_count = 0
    self.rebinned_data_objects = []
    self.plottable_count = 0
    self.plottable_2d_data_objects = []

  def AddRebinnedData(self, new_object, name='', base_data_obj = None):  
    self.rebinned_data_objects.append(new_object)
    new_object.number = self.rb_count
    self.rb_count += 1
    
  def AddPlottable2dData(self, new_object, parent = None, name='', base_data_obj = None):
    self.plottable_2d_data_objects.append(new_object)
    self.plottable_count += 1
  
  def __iadd__(self, new_object):
    if isinstance(new_object, rebinned_data):
      self.AddRebinnedData(new_object)
    elif isinstance(new_object, plottable_2d_data):
      self.AddPlottable2dData(new_object)
    return self
    

default_supervisor = Supervisor()

ID_EXIT=102

class polarizedICPFile:
    """ container class that holds one reflectometry.reduction.load object
    along with needed polarization data:  i.e. He3 cell information, and timestamps """
    def __init__(self):
        pass

class rebinned_data():
  """ class to import ICP files, and rebin them according to 
  user-supplied parameters.  A GUI is opened when an instance 
  of the class is created, and parameters are entered there.
  Add, subtract and divide functions work with other instances, 
  as long as the array dimensions are equal, so use with care. """
  
  def __init__( self, 
                file_objects = None,
                do_pol_corr = False,
                assumptions = 0,
                twoth_zero_pixel = 309.0, 
                a4_offset = 0.0, 
                description = '', 
                plot_data = True, 
                normalization = 'monitor',
                supervisor = None,
                name = '',
                grid_th2th = None,
                th2th = None,
                **kwargs):
    self.detectorMaxXPixel = 608.0
    self.pixelsPerDegree = 80.0
    #self.inFileNames = filenames
    self.inFileObjects = file_objects
    self.do_pol_corr = do_pol_corr
    self.polcorr_matrices = []
    self.assumptions = assumptions # common assumptions when polarization-correcting:
    # 0: "no assumptions (use all I++, I+-, I-+, I--)",
    # 1: "R+- assumed equal to R-+ (use I++, I-+ and I--)",
    # 2: "R-+ assumed equal to R+- (use I++, I+- and I--)",
    # 3: "R-+ and R+- equal zero (use I++, I--)"]
    self._theta_offset = 0.0
    self.twoth_zero_pixel = twoth_zero_pixel
    self.a4_offset = a4_offset
    self.normalization = normalization
    self.qxqz = None
    self.th2th = None
    self.th_in_th_out_data = None
    self.supervisor = supervisor
    self.description = description
    self.grid_th2th = grid_th2th
    self.th2th = th2th
    self.plot_data = plot_data
    self.instr_resolution = 0.00001 # difference that is negligible for all motors
    self.wiggleAmplitude = 0.10 # does this change with time?
    #self.He3_cells = He3_cells
   
    
    # ask which filenames to work with
    #self.inFileNames = self.ChooseInputFiles()
    #if file_objects == None:
    #    if (len(filenames) < 1):
    #        self.get_load_params_gui(None, -1, 'Get filenames to work with', self)
    #        print('names: ' + str(self.inFileNames))
    #    for dataFileName in self.inFileNames:
    #        self.inFileObjects.append(icp_load(dataFileName))
    #        print('datafile: ' + dataFileName + ' loaded.')
    #else:
    #    self.inFileObjects = file_objects
    
    self.number = 0

    if not self.supervisor:
        global default_supervisor
        self.supervisor = default_supervisor
        
    self.supervisor.AddRebinnedData(self, description) 
    
    
    
    #grab wavelength from first file:
    self.wavelength = self.inFileObjects[0].detector.wavelength
    
    # then make those filenames into 2th-th map
    if self.grid_th2th == None:
        self.grid_th2th = self.make_2th_th_grid()
    if self.th2th == None:
        if self.do_pol_corr:
          self.make_2th_th_map_polcorr(grid=self.grid_th2th)
        else:
          self.make_2th_th_map(grid=self.grid_th2th)
          if plot_data:
            self.th2th.wxplot()
    
  def motors(self, v):
    start = double(v[0])
    stop = double(v[-1])
    step = double((stop-start)/(len(v)-1))
    return start,step,stop

  def getTwoTheta (self, pixel, angle4 ): 
    return ( ( (self.twoth_zero_pixel - pixel ) / self.pixelsPerDegree ) + angle4 - self.a4_offset )


  def twoThetaFromPixel(self, pixel, angle):
    return_angle = (self.params['qZeroPixel'] - pixel) / self.pixelsPerDegree
    return_angle += (angle - self.params['qZeroAngle'])
    return return_angle

  def setThetaOffset(self, offset):
    self._theta_offset = offset
    self.make_2th_th_map()
      
  def GetThetaOffset(self):
    return self._theta_offset
    
  def correctionFromPixel(self, xpixel, subDim=1):
      pixelsPerDegree = 80.0 * subDim
      wiggleAmplitude = self.wiggleAmplitude
      #pixelCorrection = ( (32.0 / (2.0 * pi) ) * wiggleAmplitude * sin( 2.0 * pi * xpixel / 32.0 ) )
      widthCorrection = ( wiggleAmplitude * cos( 2.0 * pi * xpixel / 32.0 ) )
      pixelCorrection = cumsum(widthCorrection)
      return [widthCorrection, pixelCorrection]
  
  def getQxQzI (self, A3, A4, x, rawI ): 
    qLength = 2.0 * pi / self.wavelength
    rawTwoTheta = self.getTwoTheta( x, A4 )
    widthCorrection,angleCorrection = self.correctionFromPixel(x)
    correctedI = rawI / ( 1.0 + widthCorrection )
    twoTheta = rawTwoTheta + angleCorrection 
    tilt = A3 - ( twoTheta / 2.0 )
    dq = 2.0 * qLength * sin( ( pi / 180.0 ) * ( twoTheta / 2.0 ) )
    qxOut = dq * sin( pi * tilt / 180.0 )
    qzOut = dq * cos( pi * tilt / 180.0 )
    return [qxOut, qzOut, correctedI]

  def runConversion (self, event=None, th2thdata=None, params_in={}, plot_result=True):
    """convert 2th-th dataset to Qx,Qz dataset"""
    
    if th2thdata == None:
        th2thdata = self.th2th
        
    # initialize the parameters that follow the dataset around:
    print 'beginning rebin:\n'
    params = {
      'description': '',
      'x_max': 0.003,
      'x_min': -0.003,
      'x_steps': 200,
      'y_max': 0.10,
      'y_min': 0.0,
      'y_steps': 200,
      'x_units': 'Qx (inv Angstroms)',
      'y_units': 'Qz (inv Angstroms)',
      'qZeroAngle': 0.0,
      'qZeroPixel': 309.0,
      }
      
      
    # override default parameters with saved parameters from last time  
    stored_name = '.stored_rebin_params'
    if os.access(stored_name, os.R_OK):
      stored_file = open(stored_name, 'r')
      stored_params = pickle.load(stored_file)
      stored_file.close()
      for label in stored_params.keys():
        if params.has_key(label):
          params[label] = stored_params[label]

    for label in params_in.keys():
      if params.has_key(label):
        params[label] = params_in[label]
    
    if (len(params_in) < len(params)):
      s = self.get_params_gui(None,-1,'input Qx-Qz conversion parameters',params)
      params['x_steps'] = int(round(params['x_steps']))
      params['y_steps'] = int(round(params['y_steps']))
      
    #self.bin_data = zeros((params['x_steps'],params['y_steps'],4))
    #self.qxqz_data_array = plottable_2d_data(self.bin_data, params, self)
    #qZMax = self.params['qZMax']
    #qZMin = self.params['qZMin']
    qZSteps = params['y_steps']
    #qXMax = self.params['qXMax']
    #qXMin = self.params['qXMin']
    qXSteps = params['x_steps']
    #qzStepSize = ( qZMax - qZMin ) / qZSteps
    #qxStepSize = ( qXMax - qXMin ) / qXSteps
    
    #self.qxqz_array = self.convert_to_qxqz(self.twoth_th_array)
    self.qxqz_arrays = []
    
    #for dataFileObject in self.inFileObjects:
      #qxqz_array = self.convert_to_qxqz(dataFileObject)
      #self.qxqz_arrays.append(qxqz_array)
      #print('loaded: ' + dataFileObject.name)
    
    if (params['description'] == ''):
      print 'no description'
      params['description'] = self.description
      #params['description'] = self.inFileObjects[0].description
      
    qZSteps = params['y_steps']
    qXSteps = params['x_steps']

    qxqz_arrays = []
    th2th_objs = ravel([self.th2th]) # make even non-list into an iterable
    
    for th2th in th2th_objs:
        qxqz_array = self.convert_2thth_qxqz(th2th)
        
        bin_data = zeros((params['x_steps'],params['y_steps'],4))  
        self.do_rebinning(qxqz_array, params, bin_data)
        
        ##self.minimum_intensity = inf  
        #for j in range(qZSteps):
        #  for i in range(qXSteps):
        #    pixelCount = bin_data[i,j,1]
        #    if ( pixelCount > 0 ):
        #      monitorTotal = bin_data[i,j,2]
        #      avg = bin_data[i,j,0] / double(monitorTotal)
        #      bin_data[i,j,3] = avg
        #      #if (avg < self.minimum_intensity and avg > 0):
        #        #self.minimum_intensity = avg
                
        nonempty = (bin_data[:,:,1] > 0) # creates a mask
        bin_data[:,:,3][nonempty] = bin_data[:,:,0][nonempty]/bin_data[:,:,2][nonempty]
             
        qxqz_2d_params = params
        
        #self.make_2th_th_map()
        
        #zero_pixels = (self.bin_data[:,:,1]==0)
        #print zero_pixels.shape
        #self.bin_data[zero_pixels] = self.reverse_lookup(params)[zero_pixels]
        ## make sure that these pixels still show up as empty for sums:
        #self.bin_data[zero_pixels][:,1] = 0
        
        title = th2th.title[:-7] + ': qxqz' #chop off the th2th and add qxqz
        #title = 'dataset ' + str(self.number) + ': qxqz'
        qxqz_obj = plottable_2d_data(bin_data, qxqz_2d_params, self, supervisor = self.supervisor, title = title, base_data_obj = self)
        qxqz_obj.params = qxqz_2d_params
        qxqz_obj.title = title
        if self.qxqz == None:
          self.qxqz = [qxqz_obj]
        else: 
          self.qxqz.append(qxqz_obj)
        if plot_result: 
          qxqz_obj.wxplot(scale='log')
    
  def interpolate_vec(self,vec,n=1):
    if n==1:
      return vec
    else:
      new_indices = linspace(0.,float(len(vec))-1,(len(vec) * n)-1)
      mod = fmod(new_indices,1)
      new_vec = vec[floor(new_indices).astype(int)] * (1-mod) + vec[ceil(new_indices).astype(int)] * mod
      return new_vec
  
  def expand_arr(self,arr,n=1):
    if n==1:
      return arr
    else:
      new_indices_1 = array([arange(arr.shape[0]*n -1)]).T / n
      new_indices_2 = array([arange(arr.shape[1]*n -1)]) / n
      new_arr = arr[new_indices_1,new_indices_2]
      return new_arr
  
  def convert_2thth_qxqz(self, twoth_th_2d_data):
    qLength = 2.0 * pi / self.wavelength
    data = twoth_th_2d_data
    corrected_I = data.bin_data[:,:,0]
    pixel_array = data.bin_data[:,:,1]
    monitor_array = data.bin_data[:,:,2]
    twoth_stepsize = (data.params['x_max'] - data.params['x_min']) / (data.params['x_steps'] - 1)
    twoTheta_array = arange(data.params['x_steps']) * twoth_stepsize + data.params['x_min']
    twoTheta_array = array([twoTheta_array]).T
    th_stepsize = (data.params['y_max'] - data.params['y_min']) / (data.params['y_steps'] - 1)
    th_array = arange(data.params['y_steps']) * th_stepsize + data.params['y_min']
    th_array = array([th_array])
    tilt_array = th_array - ( twoTheta_array / 2.0 )
    self.twoTheta_array = twoTheta_array
    self.tilt_array = tilt_array
    qxOut = 2.0 * qLength * sin( ( pi / 180.0 ) * ( twoTheta_array / 2.0 ) ) * sin( pi * tilt_array / 180.0 )
    qzOut = 2.0 * qLength * sin( ( pi / 180.0 ) * ( twoTheta_array / 2.0 ) ) * cos( pi * tilt_array / 180.0 )
    self.qxOut = qxOut
    self.qzOut = qzOut
    qxqz_data = array([qxOut,qzOut,corrected_I,pixel_array,monitor_array])
    
    return qxqz_data

  def convert_to_th2th(self, qx, qz):
    Q = sqrt(qx * qx + qz * qz) * sign(qz)
    twoth = arcsin( self.wavelength * Q / (4. * pi) ) * 2.0 * 180./ pi
    tilt = arctan( qx / qz ) * 180. / pi
    th = tilt + twoth / 2.0
    self.twoth = twoth
    self.th = th
    return th, twoth
    
  def do_rebinning(self, qxqz_data, params, output_bins):
    print 'rebinning, oy\n'
    qZMax = params['y_max']
    qZMin = params['y_min']
    qZSteps = params['y_steps']
    qXMax = params['x_max']
    qXMin = params['x_min']
    qXSteps = params['x_steps']
    qzStepSize = ( qZMax - qZMin ) / qZSteps
    qxStepSize = ( qXMax - qXMin ) / qXSteps
    subDim = 1
    
    d = shape(qxqz_data)
    
    select_matrix = (qxqz_data[0] < qXMax) * (qxqz_data[0] >= qXMin) * (qxqz_data[1] < qZMax) * (qxqz_data[1] >= qZMin) * (qxqz_data[3] > 0)
       
    qx = qxqz_data[0][select_matrix]
    self.data_qx_list = qx.copy()
    qz = qxqz_data[1][select_matrix]
    self.data_qz_list = qz.copy()
    intens = qxqz_data[2][select_matrix]
    pixels = qxqz_data[3][select_matrix]
    self.pixels = pixels
    monitor = qxqz_data[4][select_matrix]
    #nonzero_pixels = zeros(monitor.shape, dtype = float)
    #nonzero_pixels[monitor.nonzero()] = 1.0
    
    xBin = floor( (qx - qXMin ) / qxStepSize ).astype(int32)
    zBin = floor( (qz - qZMin ) / qzStepSize ).astype(int32)

    outshape = output_bins.shape
    if (len(xBin) > 0):
      hist2d, xedges, yedges = histogram2d(xBin,zBin, bins = (outshape[0],outshape[1]), range=((0,outshape[0]),(0,outshape[1])), weights=intens)
      output_bins[:,:,0] += hist2d
      hist2d, xedges, yedges = histogram2d(xBin,zBin, bins = (outshape[0],outshape[1]), range=((0,outshape[0]),(0,outshape[1])), weights=pixels)
      output_bins[:,:,1] += hist2d
      hist2d, xedges, yedges = histogram2d(xBin,zBin, bins = (outshape[0],outshape[1]), range=((0,outshape[0]),(0,outshape[1])), weights=monitor)
      output_bins[:,:,2] += hist2d
    
    #for i in range(len(xBin)):
      #output_bins[xBin[i],zBin[i],0] += intens[i]
      #output_bins[xBin[i],zBin[i],1] += 1
      #output_bins[xBin[i],zBin[i],2] += monitor[i]
    
    return
  
  def reverse_lookup(self, params):
    base = zeros((params['x_steps'], params['y_steps']))
    qx = linspace(params['x_min'], params['x_max'], params['x_steps'])
    qx.shape = (qx.shape[0], 1)
    self.qx = qx + base
    qz = linspace(params['y_min'], params['y_max'], params['y_steps'])
    qz.shape = (1, qz.shape[0])
    self.qz = qz + base
    #print self.qx.shape, self.qz.shape
    th, twoth = self.convert_to_th2th(self.qx, self.qz)
    #self.th = th
    #self.twoth = twoth
    #print th.shape, twoth.shape
    #select_matrix = (th >= self.th_min_min and th <= self.th_max_max and twoth >= self.twoth_min_min and twoth <= self.twoth_max_max)
    th_index = (floor((th - self.th_min_min)/self.th_stepsize)).astype(int32)
    twoth_index = (floor((twoth - self.twoth_min_min)/self.twoth_stepsize)).astype(int32)
   
    th_index_min = 0
    th_index_max = self.th_2th_array[:,:,3].shape[0]
    twoth_index_min = 0
    twoth_index_max = self.th_2th_array[:,:,3].shape[1]
    
    new_th_index = th_index.copy()
    # map all indices that are outside the proper range to 0 index
    new_th_index[new_th_index < th_index_min] = 0
    new_th_index[new_th_index >= th_index_max] = 0
    new_twoth_index = twoth_index.copy()
    # map all indices that are outside the proper range to 0 index
    new_twoth_index[new_twoth_index < twoth_index_min] = 0
    new_twoth_index[new_twoth_index >= twoth_index_max] = 0
    #self.th_index = new_th_index
    #self.twoth_index = new_twoth_index
    
    # the action happens here:  the output array is filled with data with lookup indices
    # specified by th_index and twoth_index
    out_array = self.th_2th_array[new_th_index, new_twoth_index]
    
    # here we set all the data that was outside the mapped region in th,twoth to zero
    out_array[th_index < 0] = 0
    out_array[th_index >= th_index_max] = 0
    out_array[twoth_index < 0] = 0
    out_array[twoth_index >= twoth_index_max] = 0
    
    return out_array
  

  
  def make_2th_th_grid(self):
    
    num_files = len( self.inFileObjects )
    th_min = zeros(num_files)
    th_max = zeros(num_files)
    th_len = zeros(num_files)
    th_step = zeros(num_files)
    twoth_min = zeros(num_files)
    twoth_max = zeros(num_files)
    twoth_len = zeros(num_files)
    twoth_step = zeros(num_files)
    
    pixel_vector = arange(609.0)
    intens_correction, pixel_correction = self.correctionFromPixel(pixel_vector)
    corrected_pixel = pixel_vector + pixel_correction
    #print num_files
    
    for i in range(num_files):
      th_min[i] = self.inFileObjects[i].sample.angle_x.min() + self.GetThetaOffset()
      th_max[i] = self.inFileObjects[i].sample.angle_x.max() + self.GetThetaOffset()
      th_len[i] = self.inFileObjects[i].sample.angle_x.shape[0]
      #print th_min, th_max, th_len
      
      A4_vector = self.inFileObjects[i].detector.angle_x
      twoth_vector = self.getTwoTheta(pixel_vector, array([A4_vector]).T)
      twoth_max[i] = twoth_vector.max()
      twoth_min[i] = twoth_vector.min()
      
      if th_len[i] > 1:
        th_step[i] = float(th_max[i] - th_min[i]) / (th_len[i] - 1)
        th_max[i] += th_step[i] # add back on the last step
      else:
        th_step[i] = self.instr_resolution
        
      #if twoth_len[i] > 1:
        #twoth_step[i] = float(twoth_max[i] - twoth_min[i]) / twoth_len[i]
      #else:
        #twoth_step[i] = 0.0
        
    # take maximum A3 spacing as grid spacing for theta:
    th_stepsize = th_step.max()
    if th_stepsize < self.instr_resolution:
        th_stepsize = self.instr_resolution
    self.th_stepsize = th_stepsize
    # take the extrema as the limits of our new grid:
    th_min_min = th_min.min()
    self.th_min_min = th_min_min
    th_max_max = th_max.max()
    self.th_max_max = th_max_max
    
    if (th_len.max()) == 1 or (th_max_max == th_min_min):
      # special case: we only have one set of points in theta, set the stepsize to the minimum stepsize in theta
      th_steps = 1
    else:
      th_steps = int(round(float(th_max_max - th_min_min) / th_stepsize))
    self.th_steps = th_steps
    
    # take 1/80 deg as grid spacing for 2theta:
    twoth_stepsize = 1.0/self.pixelsPerDegree
    self.twoth_stepsize = twoth_stepsize
    # take the extrema as the limits of our new grid:
    twoth_min_min = twoth_min.min()
    self.twoth_min_min = twoth_min_min
    twoth_max_max = twoth_max.max()
    self.twoth_max_max = twoth_max_max
    twoth_steps = int(round(float( twoth_max_max - twoth_min_min ) / twoth_stepsize))
    self.twoth_steps = twoth_steps
    
    self.twoth, self.th = indices((twoth_steps, th_steps), dtype = float)
    self.twoth *= twoth_stepsize
    self.th *= th_stepsize
    self.twoth += twoth_min_min
    self.th += th_min_min
    
    # because we reverse the axes on twoth (matches instrument geometry:)
    self.twoth = flipud(self.twoth)
    
    th_twoth_data_extent = (twoth_min_min, twoth_max_max + self.twoth_stepsize, th_min_min, th_max_max + self.th_stepsize)
    #th_twoth_data = zeros((th_steps, twoth_steps, 4))
    #print th_min_min, th_max_max, th_steps 
    #print twoth_min_min, twoth_max_max, twoth_steps
    
    twoth_bin_edges = ((arange(twoth_steps + 1.0) - 0.0)*twoth_stepsize) + twoth_min_min
    #self.twoth_bin_edges = twoth_bin_edges
    th_bin_edges = ((arange(th_steps + 1.0) - 0.0)*th_stepsize) + th_min_min
    #self.th_bin_edges = th_bin_edges
    grid = [th_bin_edges, twoth_bin_edges, th_stepsize, twoth_stepsize]
    return grid
    
  def make_2th_th_map(self, evt = None, grid = None):
    """ makes a composite array that contains all the data objects
    in a regular theta-2theta grid.  For theta, find stepsize in all
    data sets, then use largest step size as grid spacing.  In 2theta,
    somewhat simpler: use the pixel size as the grid spacing (1/80 deg)
    """
    from reflectometry.reduction import rebin as reb
    if grid == None:
        grid = self.make_2th_th_grid()
    th_bin_edges, twoth_bin_edges, th_stepsize, twoth_stepsize = grid
    th_steps = len(th_bin_edges) - 1
    twoth_steps = len(twoth_bin_edges) - 1
    th_2th_array = zeros((th_steps, twoth_steps, 4))
    num_files = len( self.inFileObjects )
    th_min = zeros(num_files)
    th_max = zeros(num_files)
    th_len = zeros(num_files)
    th_step = zeros(num_files)
    for i in range(num_files):
      inData = self.inFileObjects[i]
      #pixel_vector = arange(609.0,0.0,-1.0)
      pixel_vector = arange(513.)
      #intens_correction, pixel_correction = self.correctionFromPixel(pixel_vector)
      #corrected_pixel = pixel_vector + pixel_correction
      intens = inData.detector.counts
      #cut off the last bit of intens_correction (which is larger by one than intens, because
      #we're making edges
      #intens_correction = intens_correction[1:]
      #corrected_I = intens / (1.0 + intens_correction)
      corrected_I = intens
      th = inData.sample.angle_x + self.GetThetaOffset()
      # need one more bin edge than data points:
      th_len[i] = inData.sample.angle_x.shape[0]
      th_min[i] = inData.sample.angle_x.min() + self.GetThetaOffset()
      th_max[i] = inData.sample.angle_x.max() + self.GetThetaOffset()
      if th_len[i] > 1 :
        th_step[i] = float(th_max[i] - th_min[i]) / (th_len[i] - 1)
        if abs(th_step[i]) < self.instr_resolution:
            th_step[i] = self.instr_resolution
        th_max[i] += th_step[i] # add back on the last step
      else:
        th_step[i] = th_stepsize
        
      th_vector = ((arange(th_len[i] + 1.0) - 0.0)*th_step[i]) + th_min[i]
 

      # if the normalization is chosen to be vs. time instead of vs. monitor counts, the change is
      # made here:
      if self.normalization == 'time':
        monitor_vector = inData.monitor.count_time
      else: # self.normalization = monitor by default
        monitor_vector = inData.monitor.counts
      
      monitor_array = array([monitor_vector]).T + zeros(shape(intens))
      #monitor_array *= (1.0 + intens_correction) # corrects for non-uniform area of pixels
      pixel_array = ones(shape(intens)) # * (1.0 + intens_correction) # corrects for non-uniform area of pixels
      #pixel_array = ones(shape(intens))
      
      A4_vector = inData.detector.angle_x
      
      if ( (A4_vector.max() - A4_vector.min() ) < self.instr_resolution ):
        #then the detector is fixed and we can pass a single 2theta vector to rebin2d
        twoth_vector = self.getTwoTheta(corrected_pixel,A4_vector.min())
        #self.twoth_vector = twoth_vector
        
        new_data = reb.rebin2d(th_vector,twoth_vector,intens.astype('float64'),th_bin_edges,twoth_bin_edges)
        #print th_vector.shape, twoth_vector.shape, corrected_I.shape, th_bin_edges.shape, twoth_bin_edges.shape
        #new_data = reb.rebin2d(twoth_vector,th_vector,intens,twoth_bin_edges,th_bin_edges)
        #print new_data
        new_pixelcount = reb.rebin2d(th_vector,twoth_vector,pixel_array.astype('float64'),th_bin_edges,twoth_bin_edges)
        new_mon = reb.rebin2d(th_vector,twoth_vector,monitor_array.astype('float64'),th_bin_edges,twoth_bin_edges)
        #new_pixelcount = reb.rebin2d(twoth_vector,th_vector,pixel_array,twoth_bin_edges,th_bin_edges)
        #new_mon = reb.rebin2d(twoth_vector,th_vector,monitor_array,twoth_bin_edges,th_bin_edges)
        th_2th_array[:,:,0] += new_data
        th_2th_array[:,:,1] += new_pixelcount
        th_2th_array[:,:,2] += new_mon
      else:
        #then the detector is not fixed, and we have to pass in each A4 value at a time to rebin
        for j in range(th_vector.shape[0]-1):
          twoth_vector = self.getTwoTheta(pixel_vector,A4_vector[j])
          #print th_vector.shape, twoth_vector.shape, corrected_I.shape, th_bin_edges.shape, twoth_bin_edges.shape
          new_data = reb.rebin(twoth_vector,intens.astype('float64'),twoth_bin_edges)
          #new_data = reb.rebin2d(twoth_vector,th_vector[j],corrected_I[:,j],twoth_bin_edges,th_bin_edges)
          new_pixelcount = reb.rebin(twoth_vector,pixel_array[j].astype('float64'),twoth_bin_edges)
          new_mon = reb.rebin(twoth_vector,monitor_array[j].astype('float64'),twoth_bin_edges)
          print 'A4 step: ' + str(j)
          th_2th_array[j,:,0] += new_data
          th_2th_array[j,:,1] += new_pixelcount
          th_2th_array[j,:,2] += new_mon
      print('loaded: ' + inData.name)
      
            
    nonzero_pixels = (th_2th_array[:,:,1] > 0)
    monitor_total = th_2th_array[:,:,2][nonzero_pixels]
    avg = th_2th_array[:,:,0][nonzero_pixels] / double(monitor_total)
    th_2th_array[:,:,3][nonzero_pixels] = avg
    self.th_2th_array = th_2th_array

    th2th_params = {
      'description': self.description,
      'x_max': twoth_bin_edges.min(),
      'x_min': twoth_bin_edges.max(),
      'x_steps': twoth_steps,
      'y_max': th_bin_edges.max(),
      'y_min': th_bin_edges.min(),
      'y_steps': th_steps,
      'x_units': '$2\\theta ({}^{\circ})$',
      'y_units': '$\\theta ({}^{\circ})$'
      }
    
    twoth_th_array = flipud(th_2th_array.swapaxes(0,1))
    # cut off stuff that should really be zero - bad fp stuff
    twoth_th_array[:,:,3][twoth_th_array[:,:,3] < 1e-16] = 0.
    self.twoth_th_array = twoth_th_array
    title = 'dataset ' + str(self.number) + ':  th2th'
    self.th2th = plottable_2d_data(twoth_th_array, th2th_params, self, supervisor = self.supervisor, title = title, base_data_obj = self)
    return

  def make_2th_th_map_polcorr(self, evt = None, grid = None):
    """ makes a composite array that contains all the data objects
    in a regular theta-2theta grid.  For theta, find stepsize in all
    data sets, then use largest step size as grid spacing.  In 2theta,
    somewhat simpler: use the pixel size as the grid spacing (1/80 deg)
    """
    from reflectometry.reduction import rebin as reb
    if grid == None:
        grid = self.make_2th_th_grid()
    th_bin_edges, twoth_bin_edges, th_stepsize, twoth_stepsize = grid
    th_steps = len(th_bin_edges) - 1
    twoth_steps = len(twoth_bin_edges) - 1
    th_2th_array = zeros((th_steps, twoth_steps, 4, 4)) # added extra dimension for spin state
    polcorr_matrices = zeros((th_steps, twoth_steps, 4, 4)) # these aren't the same dimensions as the th2th array
    num_files = len( self.inFileObjects )
    th_min = zeros(num_files)
    th_max = zeros(num_files)
    th_len = zeros(num_files)
    th_step = zeros(num_files)
    # I: Iup-up, Iup-down, Idown-up, Idown-down
    pol_index_dict = {'++': 0, '+-':1, '-+':2, '--':3}
    self.pol_states_present = set()
    for i in range(num_files):
      inData = self.inFileObjects[i]
      pol_state = inData.polarization
      print 
      pol_index = pol_index_dict[pol_state]
      self.pol_states_present.add(pol_state)
#      try:
#        pol_state = inData.polarization
#        pol_index = pol_index_dict[pol_state]
#        self.pol_states_present.add(pol_state)
#      except KeyError:
#        print "polarization state of %s is not recognized: %s" % (inData.name, inData.polarization)
#        return
      #pixel_vector = arange(609.0,0.0,-1.0)
      pixel_vector = arange(513.)
      #intens_correction, pixel_correction = self.correctionFromPixel(pixel_vector)
      #corrected_pixel = pixel_vector + pixel_correction
      intens = inData.detector.counts
      #cut off the last bit of intens_correction (which is larger by one than intens, because
      #we're making edges
      #intens_correction = intens_correction[1:]
      #corrected_I = intens / (1.0 + intens_correction)
      corrected_I = intens
      th = inData.sample.angle_x + self.GetThetaOffset()
      # need one more bin edge than data points:
      th_len[i] = inData.sample.angle_x.shape[0]
      th_min[i] = inData.sample.angle_x.min() + self.GetThetaOffset()
      th_max[i] = inData.sample.angle_x.max() + self.GetThetaOffset()
      if th_len[i] > 1:
        th_step[i] = float(th_max[i] - th_min[i]) / (th_len[i] - 1)
        th_max[i] += th_step[i] # add back on the last step
      else:
        th_step[i] = th_stepsize
        
      th_vector = ((arange(th_len[i] + 1.0) - 0.0)*th_step[i]) + th_min[i]
      print "th_vector: " + str(th_vector)
 

      # if the normalization is chosen to be vs. time instead of vs. monitor counts, the change is
      # made here:
      if self.normalization == 'time':
        monitor_vector = inData.monitor.count_time
      else: # self.normalization = monitor by default
        monitor_vector = inData.monitor.counts
      
      monitor_array = array([monitor_vector]).T + zeros(shape(intens))
      #monitor_array *= (1.0 + intens_correction) # corrects for non-uniform area of pixels
      print "monitor vector shape: " + str(array([monitor_vector]).shape)
      print "monitor array shape: " + str(monitor_array.shape)
      pixel_array = ones(shape(intens)) # * (1.0 + intens_correction) # corrects for non-uniform area of pixels
      
      cell = self.supervisor.He3_cells[str(inData.cell_id)] # link to the cell object
      A4_vector = inData.detector.angle_x
      elapsed_time = float(inData.elapsed_time)
      start_time_offset = inData.date - cell.start_datetime
      offset_seconds = start_time_offset.days * 86400. + start_time_offset.seconds # weird implementation of datetime difference measures in days, seconds, microseconds
      print 'offset_seconds = %d' % offset_seconds 
      time_vector = (arange(th_len[i]) * elapsed_time /th_len[i]) + offset_seconds
      flipper_on = (inData.polarization[0] == '-') # check for flipper on in incoming polarization state
      He3_up = (inData.polarization[1] == '+')
      pol_corr_list = []
      for t in time_vector:
          pol_corr_list.append(cell.getNTRow(t, flipper_on = flipper_on, He3_up = He3_up))
      pol_corr_list = array(pol_corr_list)
      # at this point, shape of pol_corr_list is (th_len, 1, 4).  Want an expanded list
      # that extends over all 2th values, so shape will be (th_len, 2th_len, 4)
      pol_corr_full_array = zeros((pol_corr_list.shape[0], intens.shape[0], 4)) 
      pol_corr_full_array += pol_corr_list
            
      if ( (A4_vector.max() - A4_vector.min() ) < self.instr_resolution ):
        #then the detector is fixed and we can pass a single 2theta vector to rebin2d
        twoth_vector = self.getTwoTheta(corrected_pixel,A4_vector.min())
        uncorrected_twoth = self.getTwoTheta(pixel_vector, A4_vector.min()) # for the polcorr matrix
        new_data = reb.rebin2d(th_vector,twoth_vector,intens.astype('float64'),th_bin_edges,twoth_bin_edges)
        new_pixelcount = reb.rebin2d(th_vector,twoth_vector,pixel_array.astype('float64'),th_bin_edges,twoth_bin_edges)
        new_mon = reb.rebin2d(th_vector,twoth_vector,monitor_array.astype('float64'),th_bin_edges,twoth_bin_edges)
        for p in range(4):
            pol_corr = pol_corr_full_array[:,:,p].astype('float64') * monitor_array
            new_pol_corr = reb.rebin2d(th_vector,twoth_vector,pol_corr,th_bin_edges,twoth_bin_edges)
            polcorr_matrices[:,:,pol_index,p] += new_pol_corr
        
        th_2th_array[:,:,0,pol_index] += new_data
        th_2th_array[:,:,1,pol_index] += new_pixelcount
        th_2th_array[:,:,2,pol_index] += new_mon
        # get th_len (which is same as number of distinct time steps in data) x 4 x 4 matrix of polarization corrections
        #polcorr_vec = getPolcorrMatrices(inData.cell_id, t_vector)
      else:
        #then the detector is not fixed, and we have to pass in each A4 value at a time to rebin
        for j in range(th_vector.shape[0]-1):
          twoth_vector = self.getTwoTheta(corrected_pixel,A4_vector[j])
          uncorrected_twoth = self.getTwoTheta(pixel_vector, A4_vector[j])
          new_data = reb.rebin(twoth_vector,intens[j].astype('float64'),twoth_bin_edges)
          new_pixelcount = reb.rebin(twoth_vector,pixel_array[j].astype('float64'),twoth_bin_edges)
          new_mon = reb.rebin(twoth_vector,monitor_array[j].astype('float64'),twoth_bin_edges)
          for p in range(4):
              pol_corr = pol_corr_full_array[j,:,p].astype('float64') * monitor_array[j,:]
              new_pol_corr = reb.rebin(twoth_vector,pol_corr,twoth_bin_edges)
              polcorr_matrices[j,:,pol_index,p] += new_pol_corr
          print 'A4 step: ' + str(j)
          th_2th_array[j,:,0,pol_index] += new_data
          th_2th_array[j,:,1,pol_index] += new_pixelcount
          th_2th_array[j,:,2,pol_index] += new_mon
      print('loaded: ' + inData.name)
      
            
    nonzero_pixels = (th_2th_array[:,:,1,:] > 0)
    monitor_total = th_2th_array[:,:,2,:][nonzero_pixels]
    avg = th_2th_array[:,:,0,:][nonzero_pixels] / double(monitor_total)
    th_2th_array[:,:,3,:][nonzero_pixels] = avg
    self.th_2th_array = th_2th_array
    self.polcorr_matrices = polcorr_matrices

    th2th_params = {
      'description': self.description,
      'x_max': twoth_bin_edges.min(),
      'x_min': twoth_bin_edges.max(),
      'x_steps': twoth_steps,
      'y_max': th_bin_edges.max(),
      'y_min': th_bin_edges.min(),
      'y_steps': th_steps,
      'x_units': '$2\\theta ({}^{\circ})$',
      'y_units': '$\\theta ({}^{\circ})$'
      }
    self.th2th_params = th2th_params
    twoth_th_array = flipud(th_2th_array.swapaxes(0,1))
    # cut off stuff that should really be zero - bad fp stuff
    twoth_th_array[:,:,3,:][twoth_th_array[:,:,3,:] < 1e-16] = 0.
    self.twoth_th_array = twoth_th_array
    self.th2th = []
    th2th_objs = []
    all_states = set()
    for spin_state in ['++', '+-', '-+', '--']:
      if spin_state in self.pol_states_present:
        title = 'dataset ' + str(self.number) + ':  (' + spin_state + ') th2th uncorrected'
        new_params = th2th_params.copy()
        new_params['pol_state'] = spin_state
        self.th2th.append(plottable_2d_data(twoth_th_array[:,:,:,pol_index_dict[spin_state]], new_params, self, supervisor = self.supervisor, title = title, base_data_obj = self))
        #self.runPolCorr()    
    return
  
  def guess_assumptions(self):
    assumptions = None
    polstates = self.pol_states_present
    if set(polstates) == set(["++", "+-", "-+", "--"]):
        assumptions = 0
    elif set(polstates) == set(["++", "-+", "--"]):
        assumptions = 1
    elif set(polstates) == set(["++", "+-", "--"]):
        assumptions = 2
    elif set(polstates) == set(["++", "--"]):
        assumptions = 3
    return assumptions
    
  def runPolCorr(self, evt = None):
    """ function to invert the polarization-correction matrix and apply it to the measured data """
    
    if not self.do_pol_corr:
        return # kick out if polarization not specified
        # because this means a meaningful list of cross-section data hasn't been created
        
    # see assumptions key in __init__
    polcorr_matrices = self.polcorr_matrices
    #Iuu = self.th_2th_array[:,:,0,0]
    #Iud = self.th_2th_array[:,:,0,1]
    #Idu = self.th_2th_array[:,:,0,2]
    #Idd = self.th_2th_array[:,:,0,3]
    
    new_dataset = rebinned_data( 
                file_objects = self.inFileObjects,
                do_pol_corr = False,
                assumptions = self.assumptions,
                twoth_zero_pixel = self.twoth_zero_pixel, 
                a4_offset = self.a4_offset, 
                description = self.description + ': corrected', 
                plot_data = self.plot_data, 
                normalization = self.normalization,
                supervisor = self.supervisor,
                name = '',
                grid_th2th = self.grid_th2th,
                th2th = []
                )
    n = 0
    percent_done = -1
    rnum, cnum = self.th2th[0].bin_data.shape[:2]
    nmax = rnum * cnum
    progress_meter = wx.ProgressDialog("Progress", "% done", parent=None, style=wx.PD_AUTO_HIDE|wx.PD_APP_MODAL)
    
    def progress_update(percent_done): 
        progress_meter.Update(int(percent_done), "Polarization Correction Progress:\n{0}% done".format(percent_done)) 
    
    
    self.assumptions = self.guess_assumptions()
    print "assumptions: ", self.assumptions
       
    if self.assumptions == 0:
        data_present_mask = ones(self.th2th[0].bin_data[:,:,3].shape, dtype=numpy.bool)
        for th2th in self.th2th:
            data_present_mask *= (th2th.bin_data[:,:,1] > 0)
        NT = self.polcorr_matrices
        inv_corr_matrices = zeros(NT.shape, dtype=float)
        th_len = inv_corr_matrices.shape[0]
        twoth_len = inv_corr_matrices.shape[1]
        for j in range(th_len):
            for k in range(twoth_len):
                try:
                    #if data_present_mask[j,k]:
                    inv_corr_matrices[j,k] = numpy.linalg.inv(NT[j,k])
                except:
                    pass
        self.inv_corr_matrices = inv_corr_matrices
        
        #I = self.th_2th_array[:,:,:,[0,1,2,3]]
        I = zeros((self.th2th[0].bin_data.shape + (4,)))
        I[:,:,:,0] = self.th2th[0].bin_data
        I[:,:,:,1] = self.th2th[1].bin_data
        I[:,:,:,2] = self.th2th[2].bin_data
        I[:,:,:,3] = self.th2th[3].bin_data
        I = flipud(I).swapaxes(0,1)
        Ruu_out = I[:,:,:,0].copy()
        Rud_out = I[:,:,:,1].copy()
        Rdu_out = I[:,:,:,2].copy()
        Rdd_out = I[:,:,:,3].copy()
        rows = range(I.shape[0])
        cols = range(I.shape[1])
        print "rows, cols = %d, %d" % (len(rows), len(cols))
        print "shape of inv_corr_matrices = " + str(inv_corr_matrices.shape)
        for row in rows:
            for col in cols:
                norm_factor = numpy.linalg.norm(inv_corr_matrices[row,col])
                Ruu_out[row,col,0], Rud_out[row,col,0], Rdd_out[row,col,0], Rdu_out[row,col,0] = dot(inv_corr_matrices[row,col], I[row,col,0])
                Ruu_out[row,col,1], Rud_out[row,col,1], Rdd_out[row,col,1], Rdu_out[row,col,1] = dot(inv_corr_matrices[row,col], I[row,col,1]) / norm_factor
                Ruu_out[row,col,2], Rud_out[row,col,2], Rdd_out[row,col,2], Rdu_out[row,col,2] = (1.0, 1.0, 1.0, 1.0)
                Ruu_out[row,col,3], Rud_out[row,col,3], Rdd_out[row,col,3], Rdu_out[row,col,3] = dot(inv_corr_matrices[row,col], I[row,col,3])
                # recall order of I, R is different for the way we've set up NT matrix (not diagonal)
                # [Iuu, Iud, Idu, Idd] but [Ruu, Rud, Rdd, Rdu]
                n += 1
                new_percent_done = (100 * n) / nmax
                if new_percent_done > percent_done:
                    progress_update(new_percent_done)
                    percent_done = new_percent_done
                    
        Ruu_out = flipud(Ruu_out.swapaxes(0,1))
        Rud_out = flipud(Rud_out.swapaxes(0,1))
        Rdu_out = flipud(Rdu_out.swapaxes(0,1))
        Rdd_out = flipud(Rdd_out.swapaxes(0,1))
        new_dataset.th2th.append( plottable_2d_data(Ruu_out, self.th2th[0].params, self, supervisor = self.supervisor, title = "(++) th2th corrected", base_data_obj = new_dataset))
        new_dataset.th2th.append( plottable_2d_data(Rud_out, self.th2th[1].params, self, supervisor = self.supervisor, title = "(+-) th2th corrected", base_data_obj = new_dataset))
        new_dataset.th2th.append( plottable_2d_data(Rdu_out, self.th2th[2].params, self, supervisor = self.supervisor, title = "(-+) th2th corrected", base_data_obj = new_dataset))
        new_dataset.th2th.append( plottable_2d_data(Rdd_out, self.th2th[3].params, self, supervisor = self.supervisor, title = "(--) th2th corrected", base_data_obj = new_dataset))
        
    elif self.assumptions == 1:
        data_present_mask = zeros(self.th2th[0].bin_data[:,:,3].shape, dtype=numpy.bool)
        for th2th in self.th2th:
            data_present_mask *= (th2th.bin_data[:,:,1] > 0)
        NT = self.polcorr_matrices
        NT = NT[:,:,[0,2,3],:] #remove second row
        NT[:,:,:,1] += NT[:,:,:,3] # add column 4 to column 2 (zero-indexed)
        NT = NT[:,:,:,[0,1,2]] # drop column 4 (3 in zero-indexing)
        # should now be (th_len, 2th_len, 3, 3) matrix
        
        inv_corr_matrices = zeros(NT.shape, dtype=float)
        th_len = inv_corr_matrices.shape[0]
        twoth_len = inv_corr_matrices.shape[1]
        for j in range(th_len):
            for k in range(twoth_len):
                try:
                    #if data_present_mask[j,k]:
                    inv_corr_matrices[j,k] = numpy.linalg.inv(NT[j,k])
                except:
                    pass
        self.inv_corr_matrices = inv_corr_matrices
        
        I = zeros((self.th2th[0].bin_data.shape + (3,)))
        I[:,:,:,0] = self.th2th[0].bin_data
        I[:,:,:,1] = self.th2th[1].bin_data
        I[:,:,:,2] = self.th2th[2].bin_data
        I = flipud(I).swapaxes(0,1)
        #I = self.th_2th_array[:,:,:,[0,2,3]]
        Ruu_out = I[:,:,:,0].copy()
        Rdu_out = I[:,:,:,1].copy()
        Rdd_out = I[:,:,:,2].copy()
        rows = range(I.shape[0])
        cols = range(I.shape[1])
        for row in rows:
            for col in cols:
                norm_factor = numpy.linalg.norm(inv_corr_matrices[row,col])
                Ruu_out[row,col,0], Rdu_out[row,col,0], Rdd_out[row,col,0] = dot(inv_corr_matrices[row,col], I[row,col,0])
                Ruu_out[row,col,1], Rdu_out[row,col,1], Rdd_out[row,col,1] = dot(inv_corr_matrices[row,col], I[row,col,1]) / norm_factor
                Ruu_out[row,col,2], Rdu_out[row,col,2], Rdd_out[row,col,2] = (1.0, 1.0, 1.0)
                Ruu_out[row,col,3], Rdu_out[row,col,3], Rdd_out[row,col,3] = dot(inv_corr_matrices[row,col], I[row,col,3])
                
                n += 1
                new_percent_done = (100 * n) / nmax
                if new_percent_done > percent_done:
                    progress_update(new_percent_done)
                    percent_done = new_percent_done
                    
        Ruu_out = flipud(Ruu_out.swapaxes(0,1))
        Rdu_out = flipud(Rdu_out.swapaxes(0,1))
        Rdd_out = flipud(Rdd_out.swapaxes(0,1))
        new_dataset.th2th.append( plottable_2d_data(Ruu_out, self.th2th[0].params, self, supervisor = self.supervisor, title = "(++) th2th corrected", base_data_obj = new_dataset))
        new_dataset.th2th.append( plottable_2d_data(Rdu_out, self.th2th[1].params, self, supervisor = self.supervisor, title = "(-+) th2th corrected", base_data_obj = new_dataset))
        new_dataset.th2th.append( plottable_2d_data(Rdd_out, self.th2th[2].params, self, supervisor = self.supervisor, title = "(--) th2th corrected", base_data_obj = new_dataset))
        #self.Ruu_corr = plottable_2d_data(Ruu_out, self.th2th[0].params, self, supervisor = self.supervisor, title = "Ruu_corr", base_data_obj = self)
        #self.Rdu_corr = plottable_2d_data(Rdu_out, self.th2th[1].params, self, supervisor = self.supervisor, title = "Rdu_corr", base_data_obj = self)
        #self.Rdd_corr = plottable_2d_data(Rdd_out, self.th2th[2].params, self, supervisor = self.supervisor, title = "Rdd_corr", base_data_obj = self)
        
    elif self.assumptions == 2:
        data_present_mask = zeros(self.th2th[0].bin_data[:,:,3].shape, dtype=numpy.bool)
        for th2th in self.th2th:
            data_present_mask *= (th2th.bin_data[:,:,1] > 0)
        NT = self.polcorr_matrices
        NT = NT[:,:,[0,1,3],:] #remove third row
        NT[:,:,:,1] += NT[:,:,:,3] # add column 4 to column 2 (zero-indexed)
        NT = NT[:,:,:,[0,1,2]] # drop column 4 (3 in zero-indexing)
        # should now be (th_len, 2th_len, 3, 3) matrix
        
        inv_corr_matrices = zeros(NT.shape, dtype=float)
        print inv_corr_matrices.shape
        th_len = inv_corr_matrices.shape[0]
        twoth_len = inv_corr_matrices.shape[1]
        for j in range(th_len):
            for k in range(twoth_len):
                try:
                    #if data_present_mask[j,k]:
                    inv_corr_matrices[j,k] = numpy.linalg.inv(NT[j,k])
                except:
                    pass
        self.inv_corr_matrices = inv_corr_matrices
        I = zeros((self.th2th[0].bin_data.shape + (3,)))
        I[:,:,:,0] = self.th2th[0].bin_data
        I[:,:,:,1] = self.th2th[1].bin_data
        I[:,:,:,2] = self.th2th[2].bin_data
        I = flipud(I).swapaxes(0,1)
        #I = self.th_2th_array[:,:,:,[0,1,3]]
        Ruu_out = I[:,:,:,0].copy()
        Rud_out = I[:,:,:,1].copy()
        Rdd_out = I[:,:,:,2].copy()
        rows = range(I.shape[0])
        cols = range(I.shape[1])
        for row in rows:
            for col in cols:
                norm_factor = numpy.linalg.norm(inv_corr_matrices[row,col])
                Ruu_out[row,col,0], Rud_out[row,col,0], Rdd_out[row,col,0] = dot(inv_corr_matrices[row,col], I[row,col,0])
                Ruu_out[row,col,1], Rud_out[row,col,1], Rdd_out[row,col,1] = dot(inv_corr_matrices[row,col], I[row,col,1]) / norm_factor
                Ruu_out[row,col,2], Rud_out[row,col,2], Rdd_out[row,col,2] = (1.0, 1.0, 1.0)
                Ruu_out[row,col,3], Rud_out[row,col,3], Rdd_out[row,col,3] = dot(inv_corr_matrices[row,col], I[row,col,3])
                
                n += 1
                new_percent_done = (100 * n) / nmax
                if new_percent_done > percent_done:
                    progress_update(new_percent_done)
                    percent_done = new_percent_done
                    
        Ruu_out = flipud(Ruu_out.swapaxes(0,1))
        Rud_out = flipud(Rud_out.swapaxes(0,1))
        Rdd_out = flipud(Rdd_out.swapaxes(0,1))
        new_dataset.th2th.append( plottable_2d_data(Ruu_out, self.th2th[0].params, self, supervisor = self.supervisor, title = "(++) th2th corrected", base_data_obj = new_dataset))
        new_dataset.th2th.append( plottable_2d_data(Rud_out, self.th2th[1].params, self, supervisor = self.supervisor, title = "(+-) th2th corrected", base_data_obj = new_dataset))
        new_dataset.th2th.append( plottable_2d_data(Rdd_out, self.th2th[2].params, self, supervisor = self.supervisor, title = "(--) th2th corrected", base_data_obj = new_dataset))
        #self.Ruu_corr = plottable_2d_data(Ruu_out, self.th2th[0].params, self, supervisor = self.supervisor, title = "Ruu_corr", base_data_obj = self)
        #self.Rud_corr = plottable_2d_data(Rud_out, self.th2th[1].params, self, supervisor = self.supervisor, title = "Rud_corr", base_data_obj = self)
        #self.Rdd_corr = plottable_2d_data(Rdd_out, self.th2th[2].params, self, supervisor = self.supervisor, title = "Rdd_corr", base_data_obj = self)
        
    elif self.assumptions == 3:
        data_present_mask = zeros(self.th2th[0].bin_data[:,:,3].shape, dtype=numpy.bool)
        for th2th in self.th2th:
            data_present_mask *= (th2th.bin_data[:,:,1] > 0)
        NT = self.polcorr_matrices
        NT = NT[:,:,[0,3],:] #remove both middle rows
        NT = NT[:,:,:,[0,2]] # drop columns 2,4 (1,3 in zero-indexing)
        # should now be (th_len, 2th_len, 2, 2) matrix
        
        inv_corr_matrices = zeros(NT.shape, dtype=float)
        th_len = inv_corr_matrices.shape[0]
        twoth_len = inv_corr_matrices.shape[1]
        print 'th_len, 2th_len = %f, %f' % (th_len, twoth_len)
        for j in range(th_len):
            for k in range(twoth_len):
                try:
                    #if data_present_mask[j,k]:
                    inv_corr_matrices[j,k] = numpy.linalg.inv(NT[j,k])
                except:
                    pass
                
        self.inv_corr_matrices = inv_corr_matrices
        # now do the math
        I = zeros((self.th2th[0].bin_data.shape + (2,)))
        I[:,:,:,0] = self.th2th[0].bin_data
        I[:,:,:,1] = self.th2th[1].bin_data
        I = flipud(I).swapaxes(0,1)
        #I = self.th_2th_array[:,:,:,[0,3]]
        Ruu_out = I[:,:,:,0].copy()
        Rdd_out = I[:,:,:,1].copy()
        rows = range(I.shape[0])
        cols = range(I.shape[1])
        #print rows, cols, "Ruu shape:", Ruu_out.shape
        for row in rows:
            for col in cols:
                norm_factor = numpy.linalg.norm(inv_corr_matrices[row,col])
                Ruu_out[row,col,0], Rdd_out[row,col,0] = dot(inv_corr_matrices[row,col], I[row,col,0])
                Ruu_out[row,col,1], Rdd_out[row,col,1] = dot(inv_corr_matrices[row,col], I[row,col,1]) / norm_factor
                Ruu_out[row,col,2], Rdd_out[row,col,2] = (1.0, 1.0)
                Ruu_out[row,col,3], Rdd_out[row,col,3] = dot(inv_corr_matrices[row,col], I[row,col,3])
                
                n += 1
                new_percent_done = (100 * n) / nmax
                if new_percent_done > percent_done:
                    progress_update(new_percent_done)
                    percent_done = new_percent_done
                    
        Ruu_out = flipud(Ruu_out.swapaxes(0,1))
        Rdd_out = flipud(Rdd_out.swapaxes(0,1))
        new_dataset.th2th.append( plottable_2d_data(Ruu_out, self.th2th[0].params, self, supervisor = self.supervisor, title = "(++) th2th corrected", base_data_obj = new_dataset))
        new_dataset.th2th.append( plottable_2d_data(Rdd_out, self.th2th[1].params, self, supervisor = self.supervisor, title = "(--) th2th corrected", base_data_obj = new_dataset))
                
  def subtractBackground(self, event=None, const_bkg=0.0, array_bkg=None):
        """ array_bkg will be broadcast to cover the whole array if it is just one 
        row, for instance.  If it is the whole array, that's cool too """
        
        dlg = BkgDialog(None, pol_states_available=list(self.pol_states_present))
        if dlg.ShowModal() == wx.ID_OK:
            min_theta = dlg.min_theta
            max_theta = dlg.max_theta
            const_bkg = dlg.const_bkg
            #pol_state = dlg.pol_state_chosen
        else:
            return # cancelled!  just quit this function
        dlg.Destroy()
        
        
        #source_dataset = None
        #for th2th in self.th2th:
        #    if th2th.params['pol_state'] == pol_state:
        #        source_dataset = th2th
        #        break
        source_dataset = self.th2th[0] 
        # now we know which dataset to slice to get background:
        # need to grab the slice that corresponds to chosen min_theta, max_theta
        source_th_min = source_dataset.params['y_min']
        source_th_max = source_dataset.params['y_max']
        source_th_steps = source_dataset.params['y_steps']
        min_index = (min_theta - source_th_min) / (source_th_max - source_th_min) * source_th_steps
        max_index = (max_theta - source_th_min) / (source_th_max - source_th_min) * source_th_steps
        min_index = int(min_index)
        max_index = int(max_index)
        #bkg = sum(source_dataset.bin_data[:,min_index:max_index,3], axis=1) 
        # 0 is data
        # 1 is pixels
        # 2 is monitor
        # 3 is data / monitor
        #bkg *= 1.0 / (1.0 + max_index - min_index)
        #print bkg.shape
        #plot(bkg)
        
        new_dataset = rebinned_data( 
                    file_objects = self.inFileObjects,
                    do_pol_corr = self.do_pol_corr,
                    assumptions = self.assumptions,
                    twoth_zero_pixel = self.twoth_zero_pixel, 
                    a4_offset = self.a4_offset, 
                    description = self.description + ': subtracted bkg', 
                    plot_data = self.plot_data, 
                    normalization = self.normalization,
                    supervisor = self.supervisor,
                    name = '',
                    grid_th2th = self.grid_th2th,
                    th2th = []
                    )
        
        new_dataset.polcorr_matrices = self.polcorr_matrices
        #new_th_2th_array = self.th_2th_array.copy()
        #new_th_2th_array[:,:,3] -= bkg[:,newaxis]
        #new_th_2th_array[:,:,0] -= new_th_2th_array[:,:,2] * bkg[:, newaxis]
        #new_dataset.th_2th_array = new_th_2th_array
        new_dataset._browser_selected_rows = self._browser_selected_rows
        new_dataset.pol_states_present = self.pol_states_present.copy()
        for th2th in self.th2th:
            new_data = th2th.bin_data.copy()
            bkg = sum(new_data[:,min_index:max_index,3], axis=1)
            bkg *= 1.0 / (1.0 + max_index - min_index)
            bkg += const_bkg
            new_data[:,:,3] -= bkg[:,newaxis]
            new_data[:,:,0] -= new_data[:,:,2] * bkg[:, newaxis]
           
            new_dataset.th2th.append( plottable_2d_data( new_data, 
                                                       th2th.params.copy(), 
                                                       self, 
                                                       supervisor = self.supervisor, 
                                                       title = th2th.title + " minus bkg", 
                                                       base_data_obj = new_dataset))
    
  def make_th_in_th_out_map(self):
    if (not self.th2th):
      print('make_2th_th_map() must be performed first')
      return
    th2th_params = self.th2th.params
    th_steps = th2th_params['y_steps']
    th_max = th2th_params['y_max']
    th_min = th2th_params['y_min']
    th_stepsize = float(th_max - th_min)/th_steps
    th_in = arange(th_steps, dtype = 'float') * th_stepsize + th_min
    twoth_steps = th2th_params['x_steps']
    twoth_max = th2th_params['x_max']
    twoth_min = th2th_params['x_min']
    twoth_stepsize = float(twoth_max - twoth_min)/twoth_steps
    twoth = arange(twoth_steps, dtype = 'float') * twoth_stepsize + twoth_min
    #twoth = arange(self.twoth_steps, dtype = 'float') * self.twoth_stepsize + self.twoth_min_min
    th_out_max = twoth_max - th_min
    th_out_min = twoth_min - th_max
    from scipy import ndimage as nd
    tthdata = self.th2th.bin_data
    affine_transform = array([[1.0, -th_stepsize / twoth_stepsize],[0.,1.]])
    th_out_steps = int((th_max - th_min) / twoth_stepsize + twoth_steps)
    th_in_th_out = zeros((th_out_steps, th_steps,4))
    
    th_in_th_out[:,:,0] = nd.affine_transform(tthdata[:,:,0], affine_transform, offset = 0.0, output_shape=[th_out_steps,th_steps] )
    th_in_th_out[:,:,1] = nd.affine_transform(tthdata[:,:,1], affine_transform, offset = 0.0, output_shape=[th_out_steps,th_steps] )
    th_in_th_out[:,:,2] = nd.affine_transform(tthdata[:,:,2], affine_transform, offset = 0.0, output_shape=[th_out_steps,th_steps] )
    th_in_th_out[:,:,3] = nd.affine_transform(tthdata[:,:,3], affine_transform, offset = 0.0, output_shape=[th_out_steps,th_steps] )
    #print th_in_th_out.shape
    #print th_out_max, th_out_min
    th_in_th_out_params = {
      'description': self.description,
      'x_max': th_out_max,
      'x_min': th_out_min,
      'x_steps': th_out_steps,
      'y_max': th_max,
      'y_min': th_min,
      'y_steps': th_steps,
      'x_units': '$\\theta_{\\rm{out}} ({}^{\circ})$',
      'y_units': '$\\theta_{\\rm{in}} ({}^{\circ})$'
      }
    
    th_in_th_out = flipud(self.th_in_th_out)
    # cut off stuff that should really be zero - bad fp stuff
    th_in_th_out[:,:,3][th_in_th_out[:,:,3] < 1e-16] = 0.
    self.th_in_th_out_data = plottable_2d_data(th_in_th_out, th_in_th_out_params, self)
    return

  def ChooseInputFiles(self,event = None):
      filenames_out = []
      dlg = wx.FileDialog(None, "Choose a file", '', "", "I*.*", wx.FD_MULTIPLE)
      if dlg.ShowModal() == wx.ID_OK:
        short_filenames=dlg.GetFilenames()
        files_dir = dlg.GetDirectory()
      dlg.Destroy()
      for fn in short_filenames:
        filenames_out.append(files_dir + '/' + fn)
      return filenames_out

  class get_load_params_gui(wx.Dialog):
    def __init__(self, parent, id, title, caller):
      # initialize the outer namespace here:
      self.caller = caller
      wx.Dialog.__init__(self, parent, wx.ID_ANY, title, size=(600,300), style = wx.DEFAULT_DIALOG_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE)
      
      #-- I think you need a panel if you're going to have more than one control in your frame.
      panel = wx.Panel(self, -1)
      

      #-- create the file selection activator button
      btn_fileselector = wx.Button(panel, -1, "Choose input &files: ")
      self.Bind(wx.EVT_BUTTON, self.ChooseInputFile, btn_fileselector)
      
      #-- Create the processing button, add it to the panel and wire it up to a function in the class
      btn_SaveExit = wx.Button(panel, -1, "&Save and exit")
      self.Bind(wx.EVT_BUTTON, self.saveExit, btn_SaveExit)
      
      #-- Create the close button, add it to the panel and wire it up to a function in the class
      btn_Close = wx.Button(panel, -1, "&Cancel")
      self.Bind(wx.EVT_BUTTON, self.onExit, btn_Close)
      
      #-- Now we have to create a grid to layout our controls
      sizer_main = wx.FlexGridSizer(rows=4,cols=1,hgap=1,vgap=5)
      
      filenames_string = ''
      for fn in self.caller.inFileNames:
        filenames_string += str(fn) + '\n'
      self.filenames_label = wx.TextCtrl(panel, -1, filenames_string, size=wx.Size(600,100), style=wx.TE_MULTILINE | wx.TE_DONTWRAP | wx.TE_READONLY )
     
      sizer_buttons = wx.FlexGridSizer(rows=1,cols=2, hgap=5,vgap=5)
      sizer_buttons.Add(btn_SaveExit)
      sizer_buttons.Add(btn_Close)
      
      #-- Add the grid to the panel and make it fit
      sizer_params = wx.FlexGridSizer(rows=3, cols=6, hgap=5, vgap=10)
      self.values = {}
      self.text_labels = {}
      
      params = {}
      params['twoth_zero_pixel'] = self.caller.twoth_zero_pixel
      params['twoth_offset'] = 0.0
      
      self.keys = [
      'twoth_zero_pixel',
      'twoth_offset',
      ]
      
      self.labels = [
        'Main beam center pixel at: ',
        'when A4 = '
        ]
        
      #self.q0_value = wx.TextCtrl(panel, 1, size=(100, -1))
      #self.q0_label = wx.StaticText(panel, -1, label)
      #sizer_params.Add(self.q0_label)
      #sizer_params.Add(self.q0_value)
      #self.q0_value.SetValue(
      #self.q0_value.V
      for key,label in zip(self.keys,self.labels):
        value = wx.TextCtrl(panel, 1, size=(100, -1))
        text_label = wx.StaticText(panel, -1, label)
        self.values[key] = value 
        self.text_labels[key] = text_label
        sizer_params.Add(text_label)
        sizer_params.Add(value)
        value.SetValue(str(params[key]))
        
      #self.filenames = params['inFileNames']
      
      #-- Add the grid to the panel and make it fit
      
      sizer_main.Add(self.filenames_label)
      sizer_main.Add(btn_fileselector)
      sizer_main.Add(sizer_params)
      sizer_main.Add(sizer_buttons)
      panel.SetSizer(sizer_main)
      panel.Fit()
      #-- Show the window that we've just built
      
      self.ShowModal()
      
    def ChooseInputFile(self,e=None):
      dlg = wx.FileDialog(self, "Choose a file", '', "", "I*.*", wx.FD_MULTIPLE)
      if dlg.ShowModal() == wx.ID_OK:
        self.filenames=dlg.GetFilenames()
        self.dirname=dlg.GetDirectory()
        label = ''
        for fn in self.filenames:
          label += str(self.dirname) + '/' + str(fn) + '\n'
          
        self.filenames_label.SetValue(label)
        
      dlg.Destroy()
      
    def saveExit(self, event):
      inFileNames = []
      for fn in self.filenames:
        inFileNames.append(self.dirname + '/' + fn)
      self.caller.inFileNames = inFileNames
      twoth_zero_pixel = float(self.values['twoth_zero_pixel'].GetValue())
      twoth_offset = float(self.values['twoth_offset'].GetValue())
      # now write value of q0_pixel including angle offset:
      self.caller.twoth_zero_pixel = twoth_zero_pixel - (self.caller.pixelsPerDegree * twoth_offset)
      self.Close(True)
      
    def onExit(self, event):
      self.Close(True)
      
  class get_params_gui(wx.Dialog):
    def __init__(self, parent, id, title, params):
      #wx.Dialog.__init__(self, parent, wx.ID_ANY, title, size=(600,400), style = wx.DEFAULT_DIALOG_STYLE | wx.NO_FULL_REPAINT_ON_RESIZE)
      wx.Dialog.__init__(self, parent, wx.ID_ANY, title, size=(600,400), style = wx.DEFAULT_DIALOG_STYLE)
      #-- I think you need a panel if you're going to have more than one control in your frame.
      panel = wx.Panel(self, -1)
      self.params = params

      ##-- create the file selection activator button
      #btn_fileselector = wx.Button(panel, -1, "Choose input &files: ")
      #self.Bind(wx.EVT_BUTTON, self.ChooseInputFile, btn_fileselector)
      
      #-- Create the processing button, add it to the panel and wire it up to a function in the class
      btn_SaveExit = wx.Button(panel, -1, "&Save and exit")
      self.Bind(wx.EVT_BUTTON, self.saveExit, btn_SaveExit)
      
      #-- Create the close button, add it to the panel and wire it up to a function in the class
      btn_Close = wx.Button(panel, -1, "&Cancel")
      self.Bind(wx.EVT_BUTTON, self.onExit, btn_Close)
      
      #-- Now we have to create a grid to layout our controls
      sizer_main = wx.FlexGridSizer(rows=4,cols=1,hgap=1,vgap=5)
      
      #filenames_string = ''
      #for fn in params['inFileNames']:
        #filenames_string += str(fn) + '\n'
      #self.filenames_label = wx.TextCtrl(panel, -1, filenames_string, size=wx.Size(600,100), style=wx.TE_MULTILINE | wx.TE_DONTWRAP | wx.TE_READONLY )
     
      sizer_buttons = wx.FlexGridSizer(rows=1,cols=2, hgap=5,vgap=5)
      sizer_buttons.Add(btn_SaveExit)
      sizer_buttons.Add(btn_Close)
     
      sizer_params = wx.FlexGridSizer(rows=3, cols=6, hgap=5, vgap=10)
      self.values = {}
      self.text_labels = {}
      
      self.keys = [
      'x_min',
      'x_max',
      'x_steps',
      'y_min',
      'y_max',
      'y_steps',
      'qZeroPixel',
      'qZeroAngle',
      ]
      
      self.labels = [
        'min Qx: ',
        'max Qx: ',
        'Qx steps: ',
        'min Qz: ',
        'max Qz: ',
        'Qz steps: ',
        'X pixel value for Q=0 : ',
        'when A4 = ',
        ]
      for key,label in zip(self.keys,self.labels):
        value = wx.TextCtrl(panel, 1, size=(100, -1))
        text_label = wx.StaticText(panel, -1, label)
        self.values[key] = value 
        self.text_labels[key] = text_label
        sizer_params.Add(text_label)
        sizer_params.Add(value)
        value.SetValue(str(params[key]))
        
      #self.filenames = params['inFileNames']
      
      #-- Add the grid to the panel and make it fit
      
      #sizer_main.Add(self.filenames_label)
      #sizer_main.Add(btn_fileselector)
      sizer_main.Add(sizer_params)
      sizer_main.Add(sizer_buttons)
      panel.SetSizer(sizer_main)
      sizer_main.Fit(self)
      
      #-- Show the window that we've just built
      
      self.ShowModal()
      
      
    def saveExit(self, event):
      for key in self.keys:
        self.params[key] = float(self.values[key].GetValue())
      #self.params['inFileNames'] = self.filenames
      stored_name = '.stored_rebin_params'
      if os.access('./', os.W_OK):
        stored_file = open(stored_name, 'w')
        pickle.dump(self.params, stored_file)
        stored_file.close()
      self.Close(True)
      
    def onExit(self, event):
      self.Close(True)



