import subprocess
import sys,string,math,os
import ConfigParser
import glob
import numpy as np


filesPerList=50

Sample_Name = sys.argv [1]

def checkAndMakeDir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def clearDir(dir):
    for fil in glob.glob(dir+"/*"):
        os.remove(fil)

if Sample_Name == "m10"or Sample_Name =="M10" or Sample_Name == "10":
    SampleName = "TCP_m10_w1_ptj_100_RunIISummer19UL17_MINIAOD"
if Sample_Name == "m50" or Sample_Name =="M50" or Sample_Name == "50":
    SampleName = "TCP_m50_w1_htj_BINNED_RunIISummer19UL17_MINIAOD"
if Sample_Name == "m30" or Sample_Name =="M30" or Sample_Name == "30":
    SampleName = "TCP_m_30_w_1_htj_BINNED_slc6_amd64_gcc630_MINIAOD"
    JetBINS=["0to100","100to400","400toInf"]

first_step=True
for bin in JetBINS:
    fileListDir="/uscms_data/d3/nbower/FSU/Updated_LowPTIDStudies/CMSSW_10_6_26/src/myLowPtGsfElectronsAnalyzer/myLowPtGsfElectronsAnalyzer/macros/fileLists/"+SampleName+"/"
    filename = SampleName.replace("BINNED",bin)

    if first_step==True:
        print(fileListDir)
        checkAndMakeDir(fileListDir)
        clearDir(fileListDir)
    
    print (bin)
    rootFileDir="/eos/uscms/store/user/nbower/Events/"+SampleName.replace("BINNED",bin)

    first_step=False
    fileListDir=fileListDir+filename+"/"
    checkAndMakeDir(fileListDir)
    query = "ls " + rootFileDir
    os.system(query)
    files=os.popen(query).read().split()
    for nf in range(1, len(files)+1):
        filelistIdx=int((nf-1)/filesPerList)

        if nf%filesPerList==1:
            out=open(fileListDir+filename+"_"+str(filelistIdx)+".txt","w")
        out.write("root://cmseos.fnal.gov/"+rootFileDir+"/"+files[nf-1]+"\n")
