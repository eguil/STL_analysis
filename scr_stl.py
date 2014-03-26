#!/usr/local/uvcdat/latest/bin/cdat
#
#
# Script to perform STL decomposition on time serie
# and compute std dev of seasonal and remainder as well as their ratio
#
# uses R
#
# --------------------------------------------------------------
#  E. Guilyardi - while at LBNL/LLNL -  March 2014
#
#
import sys, os, argparse
import numpy as num
import cdms2 as cdms
import support_procs as supp
import subprocess
import socket
#
# == Arguments
#
# -v <var> -fs <STL seasonal window> -ft <STL trend window> [-time <init_idx>,<count>] <file> 
# t.window= from np x 1.5 to  np x 2
# s.window: Using R to perform STL decomposition, s.window controls how rapidly the seasonal component can change. Small values allow more rapid change. Setting the seasonal window to be infinite is equivalent to forcing the seasonal component to be periodic (i.e., identical across years).

# 
# == Inits
#

home='/Users/ericg/Projets/SF_2013-2014/ENSO_AC_div'
ncdump='/opt/local/bin/ncdump'
hist_file_dir='/Users/ericg/Projets/SF_2013-2014/ENSO_AC_div/STL_analysis'

if socket.gethostname() == 'crunchy.llnl.gov':
    home='/work/guilyardi'
    ncdump='/usr/local/uvcdat/latest/bin/ncdump'
    hist_file_dir='/work/guilyardi/database/STL_out'

toolpath=home+"/STL_analysis"

# == get command line options
    
parser = argparse.ArgumentParser(description='Script to perform STL analysis')
parser.add_argument('-d', help='toggle debug mode', action='count', default='0')
parser.add_argument('-v','--var', help='variable to work on', required=True)
parser.add_argument('-l','--legend', help='variable prefix, legend and unit for output fields eg -l sosst/SST/DegC', required=True)
parser.add_argument('-s','--seasonal', help='STL seasonal window (in months)', required=True)
parser.add_argument('-f','--trend', help='STL trend window (in years)', required=True)
parser.add_argument('-t','--timeint', help='reduce time domain in STL <init_idx>,<ncount>', default="all")
parser.add_argument('-i','--input', help='input directory', default="./")
parser.add_argument('-o','--output',help='output directory', default="./")
parser.add_argument('string', metavar='File', type=str, help='netCDF input file')
args = parser.parse_args()
 
## read values
debug        = str(args.d)
indir        = args.input
outdir       = args.output
varname      = args.var
varprefix    = args.legend
seasonal_win = args.seasonal 
trend_win    = args.trend
timeint      = args.timeint
filename     = args.string

# Write command line in history file

filer=hist_file_dir+'/z_stl_hist.txt'

with open(filer, "a") as f:
    f.write('\n\r'+str(sys.argv).translate(None, "',[]"))

#
# == decode file name
#
# post-it encoding
argfile = filename.split("_")

enam  = argfile[0]
ave   = argfile[1]
date1 = argfile[2]
date2 = argfile[3]
type  = argfile[4:]

year1 = date1[:-2]
year2 = date2[:-2]
root  = enam+"_"+ave+"_"+date1+"_"+date2

type  = "_".join(type)
# cmip5 drs encoding
# to do

print "Processing",root, type, "for variable",varname

#
# == detect time dimension and length
#
f  = cdms.open(indir+"/"+filename)
time=f[varname].getTime()

if time is None:  
    print "*** no time dimension in file ",filename
    count = raw_input("Enter number of time steps: ")
else:
    timename = time.id
    count    = time.shape[0]
      

if ave in ["1m"]:
    month1=date1[-2:]
    month2=date2[-2:]
else:
    print " ERROR : case ave=",ave," not relevant "
    sys.exit(2)
    
#
# == define count of time steps
#
if timeint in ['all']:
    timecode="1 "+str(count)
    nyears   = count/12
    print " All years in input file "
else:
    init_idx = int(timeint.split(',')[0])
    count    = int(timeint.split(',')[1])
    nyears   = count/12
    timecode = str(init_idx)+" "+str(count)
    iniyr    = init_idx / 12
    year1    = int(year1)
    year2    = int(year2)
    year1    = year1 + iniyr
    year2    = year1 + nyears - 1
    print " Working on years ", year1, "-", year2

trend_win_mo = int(trend_win) * 12

#
# == Output file names
#

outstl=enam+str(trend_win)+'yoSTL_1m_'+str(year1)+str(month1)+'_'+str(year2)+str(month2)+'_'+type
outstd=enam+str(trend_win)+'yoSTL_'+str(nyears)+'y_'+str(year1)+'_'+str(year2)+'_'+type

print " "
print " Computing ",outstl,"..."
print " Computing ",outstd,"..."
print " "

#
# == Execute STL called R routine
#
cmd_str = toolpath+'/Compute_stl.R '+indir+"/"+filename+' '+varname+' '+str(year1)+' '+timecode+' '+str(seasonal_win)+' '+str(trend_win_mo)+' '+debug+' '+varprefix

if debug == "1":
    print cmd_str

os.system(cmd_str)

# Post processing
if os.path.isfile('Out_stl.nc'):
    os.rename('Out_stl.nc', outdir+'/'+outstl)
if os.path.isfile('Out_std_stl.nc'):
    os.rename('Out_std_stl.nc', outdir+'/'+outstd)
