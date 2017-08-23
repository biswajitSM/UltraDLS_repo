#importing modules
import numpy as np
import pandas as pd
from pylab import *
import glob
import os
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

parentdir = r'C:\Users\Garg\Desktop\20170522_scat_AuNR'
def get_files(folder=parentdir, extensions = [".csv"]):
    os.chdir(folder)
    folderdir = os.getcwd()
    allfile_paths = pd.DataFrame()
    for dirpath, dirnames, filenames in os.walk("."):
        subfilepaths = []
        for filename in filenames:
            path = folderdir + dirpath[1:] + '\\' + filename
            subfilepaths.append(path)
        subfilepaths = np.reshape(subfilepaths, (1, np.product(shape(subfilepaths))))[0]
        subfilepaths = pd.DataFrame(subfilepaths, columns=([dirpath]))
        subfilepaths = subfilepaths[subfilepaths[dirpath].str.contains(extensions[0])]; subfilepaths.reset_index(drop=True, inplace=True)
        allfile_paths = pd.concat([allfile_paths, subfilepaths], axis=1)
    allfile_paths = allfile_paths.dropna(axis=1, how='all');allfile_paths.reset_index(drop=True, inplace=True)
    return(allfile_paths)

def AuNR_csvfiles(parentdir, filenamelength=15):
    allfile_paths = get_files(parentdir)
    indices = ['Point_#', 'filename','path']; group = list(allfile_paths.columns.values)#; others = set()
    columns = [group, indices]
    columns = pd.MultiIndex.from_product(columns)
    ret = pd.DataFrame(columns=columns)
    for c in ret.columns.levels[0]:
        ret[c] = allfile_paths[c]
    for c in ret.columns.levels[0]:
        ret[c]['Point_#'] = ret[c]['Point_#'].str.replace('BG_1', '00');
        ret[c]['Point_#'] = ret[c]['Point_#'].str.replace('BG1', '00');
        ret[c]['Point_#'] = ret[c]['Point_#'].str.replace('BG', '00');
        ret[c]['Point_#'] = ret[c]['Point_#'].str[-6:-4];
        ret[c]['Point_#'] = ret[c]['Point_#'].str.replace('_', '');
        ret[c]['filename'] = ret[c]['filename'].str[-filenamelength:];
    allfiles=ret
    return(allfiles)

from lmfit.models import LinearModel, LorentzianModel
mod = LorentzianModel()
def dEeVtonm(dE, peak):
    '''Both dE and peak values should be given in eV (energy)'''
    return 1239.84197*((1/(peak-(dE/2)))-(1/(peak+(dE/2)))) # nm
allfiles = AuNR_csvfiles(parentdir, filenamelength=15)
def SameAuNRcompare(allfiles=allfiles, point_number=1, skiplines=400):
    fig, ax = plt.subplots(figsize=(10, 6))
    for c in allfiles.columns.levels[0]:
        df_all = allfiles[c]
        df_all = df_all.dropna(axis=0, how='all');df_all.reset_index(drop=True, inplace=True)
        Point_to_int = df_all['Point_#'].astype(int)
        df_all['Point_#'] = Point_to_int

        df_pointselect = df_all[df_all['Point_#']==point_number]; df_pointselect.reset_index(drop=True, inplace=True)#AuNR file
        df_background = df_all[df_all['Point_#']== 0]; df_background.reset_index(drop=True, inplace=True) #Background file
        if not df_background.empty:
            df_BG = pd.read_csv(df_background['path'][0], skiprows=skiplines, header=None, sep = '\s+', engine='python');
        
        if not df_pointselect.empty:
            df = pd.read_csv(df_pointselect['path'][0], skiprows=skiplines, header=None, sep = '\s+', engine='python');
            if not df_background.empty:
                df[1] = df[1]-df_BG[1];
            L_max=df[1].idxmax();L_max=df[0][L_max];
            E = 1239.84197/df[0]# converting to energy in eV
            df = df[df[1]<7*average(df[1])]
            y=df[1]/max(df[1]);#normalized intensity
            result = mod.fit(y, amplitude=1.0, center=1.9, sigma=1.0, x=E)

            temp =result.best_values
            FWHM = 2*temp['sigma'] #eV
            peak = temp['center']#eV
            FWHM = dEeVtonm(FWHM, peak)# nm
            peak_nm = 1239.84197/peak #nm

            peak_error= str((result.params['center'], 'value')); ds2=peak_error.split('+/- ')[1];
            peak_error = float(ds2.split(', bounds')[0])#eV
            peak_error_nm = dEeVtonm(peak_error, peak) #nm

            FWHM_error= str((result.params['center'], 'value')); ds2=FWHM_error.split('+/- ')[1];
            FWHM_error = 2 * float(ds2.split(', bounds')[0])#eV
            FWHM_error_nm = dEeVtonm(peak_error, peak) #nm

            fname = c + '\\' +df_pointselect['path'][0][-10:-4]
            strng = fname +"\n" +'_SPR:' +str(int(peak_nm))+' nm_' +' FWHM:'+ str(int(FWHM)) + ' nm'

            plot(df[0], df[1]/max(df[1]), label=strng)
            print(strng)
            plot(df[0], result.best_fit/max(result.best_fit), 'r',label=str(int(peak_nm))+' nm')
    legend(fontsize=10, bbox_to_anchor=(1.5, 1))
    xlim(560, 750)
    ylim(0, 1.1)
    tight_layout()
    return(fig)