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
from subprocess import Popen, PIPE
import socket
import string
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
postit = False
# post-it encoding
if postit:
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
    
    print "Processing",root, type, "for variable",varname
else:
#
# cmip5 drs encoding
# 
    project='cmip5'
    exper = 'piControl'
    realm = 'ocn' 
#exper = 'historical'

    simus_hist={'CCSM4':'r2i1p1.ver-v20130425','MRI-CGCM3':'r1i1p1.ver-v20120701','MIROC5':'r1i1p1.ver-1','IPSL-CM5B-LR':'r1i1p1.ver-v20120114','GFDL-CM3':'r1i1p1.ver-v20110601','GFDL-ESM2M':'r1i1p1.ver-v20111228','CNRM-CM5-2':'r1i1p1.ver-v20130401','CESM1-BGC':'r1i1p1.ver-v20130213','CanESM2':'r1i1p1.ver-v20120623','HadGEM2-CC':'r1i1p1','CMCC-CM':'r1i1p1.ver-v20120627'}

    simus_pictr={'CCSM4':'r1i1p1.ver-v20130513','MRI-CGCM3':'r1i1p1.ver-v20120701','MIROC5':'r1i1p1.ver-1','IPSL-CM5B-LR':'r1i1p1.ver-v20120114','GFDL-CM3':'r1i1p1.ver-v1','GFDL-ESM2M':'r1i1p1.ver-v20130214','CNRM-CM5-2':'r1i1p1.ver-v20130402','CESM1-BGC':'r1i1p1.ver-v20130216','CanESM2':'r1i1p1.ver-v20111028','HadGEM2-CC':'r1i1p1','CMCC-CM':'r1i1p1.ver-v20120627.'}

    model = 'CCSM4'
##model = 'MRI-CGCM3'
#model = 'MIROC5'
#model = 'IPSL-CM5B-LR'
##model = 'GFDL-CM3'
#model = 'GFDL-ESM2M'
#model = 'CNRM-CM5-2'
#model = 'CESM1-BGC'
#model = 'CanESM2'
##model = 'HadGEM2-CC'
##model =  'CMCC-CM'

    if exper == 'historical':
        ript    = string.split(simus_hist[model],'.')[0]
        versiont = string.split(simus_hist[model],'.')[1]
    elif exper == 'piControl':
        ript    = string.split(simus_pictr[model],'.')[0]
        versiont = string.split(simus_pictr[model],'.')[1]

    indir = '/work/'+project+'/'+exper+'/'+realm+'/mo/'+varname
    
    filename = project+'.'+model+'.'+exper+'.'+ript+'.mo.ocn.Omon.'+varname+'.'+versiont+'.latestX.xml'
    #find real directory
    print 'Processing',filename, 'from',indir

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
    
    if not postit:
        date1 = '0101'
        date2 = "%02d" % nyears+'12'
        year1 = date1[:-2]
        year2 = date2[:-2]
        ave = '1m'

    if ave in ["1m"]:
        month1=date1[-2:]
        month2=date2[-2:]
    else:
        print " ERROR : case ave=",ave," not relevant "
        sys.exit(2)
    year1    = int(year1)
    year2    = int(year2)
    year1    = year1 + iniyr
    year2    = year1 + nyears - 1
    
    print " Working on years ", year1, "-", year2

trend_win_mo = int(trend_win) * 12

#
# == Output file names
#
if postit:
    outstl=enam+str(trend_win)+'yoSTL_1m_'+str(year1)+str(month1)+'_'+str(year2)+str(month2)+'_'+type
    outstd=enam+str(trend_win)+'yoSTL_'+str(nyears)+'y_'+str(year1)+'_'+str(year2)+'_'+type
else:
    outstl = project+'.'+model+'yoSTL.'+exper+'.'+ript+'.mo.ocn.Omon.'+varname+'.'+versiont+'.latestX.xml'
    outstd = project+'.'+model+'yoSTLm.'+exper+'.'+ript+'.mo.ocn.Omon.'+varname+'.'+versiont+'.latestX.xml'

print " "
print " Computing ",outstl,"..."
print " Computing ",outstd,"..."
print " "

#
# == Execute STL called R routine
#
# TODO Modify to read using cdms
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
