#composites of EN events conditioned on PV strength
import sys
import numpy as np
import xarray as xr
import os
import plots
NAME='SIF'
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/assessment_SH_zonal_asymmetries/data/'
FIG_PATH = '/datos/osman/assessment_SH_zonal_asymmetries/figures/impacts/'
FILE_VAR = NAME + '_s4_aug_feb.nc4'
FILE_NINIO_S4 = 'ninio34_monthly.nc4'
FILE_PV_S4 = 'SPV_index.nc4'
VAR = xr.open_dataset(PATH_DATA + FILE_VAR)
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.90, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.10, dim='dim_0', interpolation='linear')

# compute EN-LA composites conditioned on PV anomalies
ninio34_WPV = ninio34.sel(dim_0 = index_SPV_upper.values)
VAR_WPV = VAR.sel(realiz = index_SPV_upper.values)
ninio34_SPV = ninio34.sel(dim_0 = index_SPV_lower.values)
VAR_SPV = VAR.sel(realiz = index_SPV_lower.values)

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))
#enso during weak PoV
index_ninio_WPV = ninio34_WPV.ninio34_index >= ninio34_WPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_WPV = ninio34_WPV.ninio34_index <= ninio34_WPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_WPV = np.logical_and(ninio34_WPV.ninio34_index < ninio34_WPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_WPV.ninio34_index > ninio34_WPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))
#enso during strong PoV
index_ninio_SPV = ninio34_SPV.ninio34_index >= ninio34_SPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_SPV = ninio34_SPV.ninio34_index <= ninio34_SPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_SPV = np.logical_and(ninio34_SPV.ninio34_index < ninio34_SPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_SPV.ninio34_index > ninio34_SPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
    var_ninio_WPV = np.mean(VAR_WPV.isel(month=i, realiz=index_ninio_WPV.values,
                                         drop=True), axis=0).to_array().squeeze()
    var_normal_WPV = np.mean(VAR_WPV.isel(month=i, realiz=index_normal_WPV.values,
                                          drop=True), axis=0).to_array().squeeze()
    var_ninia_WPV = np.mean(VAR_WPV.isel(month=i, realiz=index_ninia_WPV.values, 
                                         drop=True), axis=0).to_array().squeeze()
    var_ninio_SPV = np.mean(VAR_SPV.isel(month=i, realiz=index_ninio_SPV.values, 
                                         drop=True), axis=0).to_array().squeeze()
    var_normal_SPV = np.mean(VAR_SPV.isel(month=i, realiz=index_normal_SPV.values,
                                          drop=True), axis=0).to_array().squeeze()
    var_ninia_SPV = np.mean(VAR_SPV.isel(month=i, realiz=index_ninia_SPV.values, 
                                         drop=True), axis=0).to_array().squeeze()
    var_ninio_all = np.mean(VAR.isel(month=i, realiz=index_ninio_all.values,
                                     drop=True), axis=0).to_array().squeeze()
    var_normal_all = np.mean(VAR.isel(month=i, realiz=index_normal_all.values,
                                      drop=True), axis=0).to_array().squeeze()
    var_ninia_all = np.mean(VAR.isel(month=i, realiz=index_ninia_all.values,
                                     drop=True), axis=0).to_array().squeeze()
    tit = 'Composites S4 SIE Conditioned - SPoV - ' + month[i]
    filename = FIG_PATH + 'SIE_composites_ENSO_' + month[i] +'_SPoV.png'
    plots.PlotSIECompositesENSOPoV(NAME,var_ninio_all, var_normal_all, var_ninia_all,
                                   var_ninio_WPV, var_normal_WPV, var_ninia_WPV,
                                   var_ninio_SPV, var_normal_SPV, var_ninia_SPV,
                                   VAR.latitude, VAR.longitude, tit, filename)

for i in np.arange(0, 5):
    VAR_s = VAR.isel(month=range(i, i+3)).mean(dim='month')
    VAR_s_WPV = VAR_s.sel(realiz=index_SPV_upper.values)
    VAR_s_SPV = VAR_s.sel(realiz=index_SPV_lower.values)
    var_ninio_WPV = np.mean(VAR_s_WPV.isel(realiz=index_ninio_WPV.values),
                            axis=0).to_array().squeeze()
    var_normal_WPV = np.mean(VAR_s_WPV.isel(realiz=index_normal_WPV.values),
                             axis=0).to_array().squeeze()
    var_ninia_WPV = np.mean(VAR_s_WPV.isel(realiz=index_ninia_WPV.values),
                            axis=0).to_array().squeeze()
    var_ninio_SPV = np.mean(VAR_s_SPV.isel(realiz=index_ninio_SPV.values),
                            axis=0).to_array().squeeze()
    var_normal_SPV = np.mean(VAR_s_SPV.isel(realiz=index_normal_SPV.values),
                             axis=0).to_array().squeeze()
    var_ninia_SPV = np.mean(VAR_s_SPV.isel(realiz=index_ninia_SPV.values),
                            axis=0).to_array().squeeze()
    var_ninio_all = np.mean(VAR_s.isel(realiz=index_ninio_all.values),
                            axis=0).to_array().squeeze()
    var_normal_all = np.mean(VAR_s.isel(realiz=index_normal_all.values),
                             axis=0).to_array().squeeze()
    var_ninia_all = np.mean(VAR_s.isel(realiz=index_ninia_all.values),
                            axis=0).to_array().squeeze()
    tit = 'Composites S4 SIE Conditioned - SPoV - ' + seas[i]
    filename = FIG_PATH + 'SIE_composites_ENSO_' + seas[i] +'_SPoV.png'
    plots.PlotSIECompositesENSOPoV(NAME, var_ninio_all, var_normal_all, var_ninia_all,
                                  var_ninio_WPV, var_normal_WPV, var_ninia_WPV,
                                  var_ninio_SPV, var_normal_SPV, var_ninia_SPV,
                                  VAR.latitude, VAR.longitude, tit, filename)

