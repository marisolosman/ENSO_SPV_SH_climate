#composites of EN events conditioned on PV strength
import sys
import numpy as np
import xarray as xr
import os
import regional_plots

NAME = sys.argv[1]
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '/pikachu/datos/osman/assessment_SH_zonal_asymmetries/data/'
FIG_PATH = '/pikachu/datos/osman/assessment_SH_zonal_asymmetries/figures/impacts/'
FILE_VAR = NAME + '_s4_aug_feb.nc4'
FILE_NINIO_S4 = 'ninio34_monthly.nc4'
FILE_PV_S4 = 'SPV_index.nc4'
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0',
                                                                          interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0',
                                                                          interpolation='linear')
#enso during weak PoV

index_ninio_WPV = np.logical_and(index_ninio_all.values, index_SPV_upper.values)
index_ninia_WPV = np.logical_and(index_ninia_all.values, index_SPV_upper.values)

#enso during strong PoV
index_ninio_SPV = np.logical_and(index_ninio_all.values, index_SPV_lower.values)
index_ninia_SPV = np.logical_and(index_ninia_all.values, index_SPV_lower.values)

nn_ninio_all = np.sum(index_ninio_all.values)
nn_ninia_all = np.sum(index_ninia_all.values)
nn_all = np.shape(ninio34.ninio34_index.values)[0]

nn_ninio_WPV =np.sum(index_ninio_WPV)
nn_ninio_SPV =np.sum(index_ninio_SPV)
nn_ninia_WPV = np.sum(index_ninia_WPV)
nn_ninia_SPV = np.sum(index_ninia_SPV)
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

VAR = xr.open_dataset(PATH_DATA + FILE_VAR)

for i in np.arange(0, 7):
    var_ninio_WPV = np.mean(VAR.isel(month=i, realiz=index_ninio_WPV, drop=True),
                            axis=0).to_array().squeeze()
    SS_ninio_WPV = np.var(VAR.isel(month=i, realiz=index_ninio_WPV, drop=True),
                          axis=0).to_array().squeeze()/nn_ninio_WPV
    var_ninia_WPV = np.mean(VAR.isel(month=i, realiz=index_ninia_WPV, drop=True),
                            axis=0).to_array().squeeze()
    SS_ninia_WPV = np.var(VAR.isel(month=i, realiz=index_ninia_WPV, drop=True),
                          axis=0).to_array().squeeze()/nn_ninia_WPV
    var_ninio_SPV = np.mean(VAR.isel(month=i, realiz=index_ninio_SPV, drop=True),
                            axis=0).to_array().squeeze()
    SS_ninio_SPV = np.var(VAR.isel(month=i, realiz=index_ninio_SPV, drop=True),
                          axis=0).to_array().squeeze()/nn_ninio_SPV
    var_ninia_SPV = np.mean(VAR.isel(month=i, realiz=index_ninia_SPV, drop=True),
                            axis=0).to_array().squeeze()
    SS_ninia_SPV = np.var(VAR.isel(month=i, realiz=index_ninia_SPV, drop=True),
                          axis=0).to_array().squeeze()/nn_ninia_SPV
    var_ninio_all = np.mean(VAR.isel(month=i, realiz=index_ninio_all.values, drop=True),
                            axis=0).to_array().squeeze()
    SS_ninio_all = np.var(VAR.isel(month=i, realiz=index_ninio_all.values),
                          axis=0).to_array().squeeze()/np.sum(index_ninio_all.values)
    var_all = np.mean(VAR.isel(month=i, drop=True), axis=0).to_array().squeeze()
    SS_all = np.var(VAR.isel(month=i, drop=True), axis=0).to_array().squeeze() / nn_all
    var_ninia_all = np.mean(VAR.isel(month=i, realiz=index_ninia_all.values, drop=True),
                            axis=0).to_array().squeeze()
    SS_ninia_all = np.var(VAR.isel(month=i, realiz=index_ninia_all.values, drop=True),
                          axis=0).to_array().squeeze()/np.sum(index_ninia_all.values)
    tit = 'Composites S4 ' + NAME + ' Conditioned - SPoV - ' + month[i]
    filename = FIG_PATH + NAME + '_composites_ENSO_' + month[i] +'_SPoV_Aust.png'
    regional_plots.PlotVARCompositesENSOPoVSIGAust(NAME,var_ninio_all - var_all,
                                      var_ninia_all - var_all,
                                      var_ninio_WPV - var_all,
                                      var_ninia_WPV - var_all,
                                      var_ninio_SPV - var_all,
                                      var_ninia_SPV - var_all,
                                      np.sqrt(SS_ninio_all + SS_all),
                                      np.sqrt(SS_ninia_all + SS_all),
                                      np.sqrt(SS_ninio_WPV + SS_all),
                                      np.sqrt(SS_ninia_WPV + SS_all),
                                      np.sqrt(SS_ninio_SPV + SS_all),
                                      np.sqrt(SS_ninia_SPV + SS_all),
                                      nn_ninio_all + nn_all - 2,
                                      nn_ninia_all + nn_all - 2,
                                      nn_ninio_WPV + nn_all - 2,
                                      nn_ninia_WPV + nn_all - 2,
                                      nn_ninio_SPV + nn_all - 2,
                                      nn_ninia_SPV + nn_all - 2,
                                      VAR.latitude, VAR.longitude, tit, filename)
    filename = FIG_PATH + NAME + '_composites_ENSO_' + month[i] +'_SPoV_Afric.png'
    regional_plots.PlotVARCompositesENSOPoVSIGAfric(NAME,var_ninio_all - var_all,
                                      var_ninia_all - var_all,
                                      var_ninio_WPV - var_all,
                                      var_ninia_WPV - var_all,
                                      var_ninio_SPV - var_all,
                                      var_ninia_SPV - var_all,
                                      np.sqrt(SS_ninio_all + SS_all),
                                      np.sqrt(SS_ninia_all + SS_all),
                                      np.sqrt(SS_ninio_WPV + SS_all),
                                      np.sqrt(SS_ninia_WPV + SS_all),
                                      np.sqrt(SS_ninio_SPV + SS_all),
                                      np.sqrt(SS_ninia_SPV + SS_all),
                                      nn_ninio_all + nn_all - 2,
                                      nn_ninia_all + nn_all - 2,
                                      nn_ninio_WPV + nn_all - 2,
                                      nn_ninia_WPV + nn_all - 2,
                                      nn_ninio_SPV + nn_all - 2,
                                      nn_ninia_SPV + nn_all - 2,
                                      VAR.latitude, VAR.longitude, tit, filename)
    filename = FIG_PATH + NAME + '_composites_ENSO_' + month[i] +'_SPoV_Sudam.png'
    regional_plots.PlotVARCompositesENSOPoVSIGSudam(NAME,var_ninio_all - var_all,
                                      var_ninia_all - var_all,
                                      var_ninio_WPV - var_all,
                                      var_ninia_WPV - var_all,
                                      var_ninio_SPV - var_all,
                                      var_ninia_SPV - var_all,
                                      np.sqrt(SS_ninio_all + SS_all),
                                      np.sqrt(SS_ninia_all + SS_all),
                                      np.sqrt(SS_ninio_WPV + SS_all),
                                      np.sqrt(SS_ninia_WPV + SS_all),
                                      np.sqrt(SS_ninio_SPV + SS_all),
                                      np.sqrt(SS_ninia_SPV + SS_all),
                                      nn_ninio_all + nn_all - 2,
                                      nn_ninia_all + nn_all - 2,
                                      nn_ninio_WPV + nn_all - 2,
                                      nn_ninia_WPV + nn_all - 2,
                                      nn_ninio_SPV + nn_all - 2,
                                      nn_ninia_SPV + nn_all - 2,
                                      VAR.latitude, VAR.longitude, tit, filename)
    filename = FIG_PATH + NAME + '_composites_ENSO_' + month[i] +'_SPoV_Antarc.png'
    regional_plots.PlotVARCompositesENSOPoVSIGAntarc(NAME,var_ninio_all - var_all,
                                      var_ninia_all - var_all,
                                      var_ninio_WPV - var_all,
                                      var_ninia_WPV - var_all,
                                      var_ninio_SPV - var_all,
                                      var_ninia_SPV - var_all,
                                      np.sqrt(SS_ninio_all + SS_all),
                                      np.sqrt(SS_ninia_all + SS_all),
                                      np.sqrt(SS_ninio_WPV + SS_all),
                                      np.sqrt(SS_ninia_WPV + SS_all),
                                      np.sqrt(SS_ninio_SPV + SS_all),
                                      np.sqrt(SS_ninia_SPV + SS_all),
                                      nn_ninio_all + nn_all - 2,
                                      nn_ninia_all + nn_all - 2,
                                      nn_ninio_WPV + nn_all - 2,
                                      nn_ninia_WPV + nn_all - 2,
                                      nn_ninio_SPV + nn_all - 2,
                                      nn_ninia_SPV + nn_all - 2,
                                      VAR.latitude, VAR.longitude, tit, filename)
for i in np.arange(0, 5):
    VAR_s = VAR.isel(month=range(i, i+3)).mean(dim='month')
    var_ninio_WPV = np.mean(VAR_s.isel(realiz=index_ninio_WPV),
                            axis=0).to_array().squeeze()
    SS_ninio_WPV = np.var(VAR_s.isel(realiz=index_ninio_WPV),
                          axis=0).to_array().squeeze()/nn_ninio_WPV
    var_ninia_WPV = np.mean(VAR_s.isel(realiz=index_ninia_WPV),
                            axis=0).to_array().squeeze()
    SS_ninia_WPV = np.var(VAR_s.isel(realiz=index_ninia_WPV),
                          axis=0).to_array().squeeze()/nn_ninia_WPV
    var_ninio_SPV = np.mean(VAR_s.isel(realiz=index_ninio_SPV),
                            axis=0).to_array().squeeze()
    SS_ninio_SPV = np.var(VAR_s.isel(realiz=index_ninio_SPV),
                          axis=0).to_array().squeeze()/nn_ninio_SPV
    var_ninia_SPV = np.mean(VAR_s.isel(realiz=index_ninia_SPV),
                            axis=0).to_array().squeeze()
    SS_ninia_SPV = np.var(VAR_s.isel(realiz=index_ninia_SPV),
                          axis=0).to_array().squeeze()/nn_ninia_SPV
    var_ninio_all = np.mean(VAR_s.isel(realiz=index_ninio_all.values),
                            axis=0).to_array().squeeze()
    SS_ninio_all = np.var(VAR_s.isel(realiz=index_ninio_all.values),
                          axis=0).to_array().squeeze()/np.sum(index_ninio_all.values)
    var_all = np.mean(VAR_s, axis=0).to_array().squeeze()
    SS_all = np.var(VAR_s, axis=0).to_array().squeeze()/nn_all
    var_ninia_all = np.mean(VAR_s.isel(realiz=index_ninia_all.values),
                            axis=0).to_array().squeeze()
    SS_ninia_all = np.var(VAR_s.isel(realiz=index_ninia_all.values),
                          axis=0).to_array().squeeze()/np.sum(index_ninia_all.values)
    tit = 'Composites S4 ' + NAME + ' Conditioned - SPoV - ' + seas[i]
    filename = FIG_PATH + NAME + '_composites_ENSO_' + seas[i] +'_SPoV_Aust.png'
    regional_plots.PlotVARCompositesENSOPoVSIGAust(NAME, var_ninio_all-var_all,
                                      var_ninia_all- var_all,
                                      var_ninio_WPV - var_all,
                                      var_ninia_WPV - var_all,
                                      var_ninio_SPV - var_all,
                                      var_ninia_SPV - var_all,
                                      np.sqrt(SS_ninio_all + SS_all),
                                      np.sqrt(SS_ninia_all + SS_all),
                                      np.sqrt(SS_ninio_WPV + SS_all),
                                      np.sqrt(SS_ninia_WPV + SS_all),
                                      np.sqrt(SS_ninio_SPV + SS_all),
                                      np.sqrt(SS_ninia_SPV + SS_all),
                                      nn_ninio_all + nn_all - 2,
                                      nn_ninia_all + nn_all - 2,
                                      nn_ninio_WPV + nn_all - 2,
                                      nn_ninia_WPV + nn_all - 2,
                                      nn_ninio_SPV + nn_all - 2,
                                      nn_ninia_SPV + nn_all - 2,
                                      VAR.latitude, VAR.longitude, tit, filename)
    filename = FIG_PATH + NAME + '_composites_ENSO_' + seas[i] +'_SPoV_Afric.png'
    regional_plots.PlotVARCompositesENSOPoVSIGAfric(NAME, var_ninio_all-var_all,
                                      var_ninia_all- var_all,
                                      var_ninio_WPV - var_all,
                                      var_ninia_WPV - var_all,
                                      var_ninio_SPV - var_all,
                                      var_ninia_SPV - var_all,
                                      np.sqrt(SS_ninio_all + SS_all),
                                      np.sqrt(SS_ninia_all + SS_all),
                                      np.sqrt(SS_ninio_WPV + SS_all),
                                      np.sqrt(SS_ninia_WPV + SS_all),
                                      np.sqrt(SS_ninio_SPV + SS_all),
                                      np.sqrt(SS_ninia_SPV + SS_all),
                                      nn_ninio_all + nn_all - 2,
                                      nn_ninia_all + nn_all - 2,
                                      nn_ninio_WPV + nn_all - 2,
                                      nn_ninia_WPV + nn_all - 2,
                                      nn_ninio_SPV + nn_all - 2,
                                      nn_ninia_SPV + nn_all - 2,
                                      VAR.latitude, VAR.longitude, tit, filename)
    filename = FIG_PATH + NAME + '_composites_ENSO_' + seas[i] +'_SPoV_Sudam.png'
    regional_plots.PlotVARCompositesENSOPoVSIGSudam(NAME, var_ninio_all-var_all,
                                      var_ninia_all- var_all,
                                      var_ninio_WPV - var_all,
                                      var_ninia_WPV - var_all,
                                      var_ninio_SPV - var_all,
                                      var_ninia_SPV - var_all,
                                      np.sqrt(SS_ninio_all + SS_all),
                                      np.sqrt(SS_ninia_all + SS_all),
                                      np.sqrt(SS_ninio_WPV + SS_all),
                                      np.sqrt(SS_ninia_WPV + SS_all),
                                      np.sqrt(SS_ninio_SPV + SS_all),
                                      np.sqrt(SS_ninia_SPV + SS_all),
                                      nn_ninio_all + nn_all - 2,
                                      nn_ninia_all + nn_all - 2,
                                      nn_ninio_WPV + nn_all - 2,
                                      nn_ninia_WPV + nn_all - 2,
                                      nn_ninio_SPV + nn_all - 2,
                                      nn_ninia_SPV + nn_all - 2,
                                      VAR.latitude, VAR.longitude, tit, filename)
    filename = FIG_PATH + NAME + '_composites_ENSO_' + seas[i] +'_SPoV_Antarc.png'
    regional_plots.PlotVARCompositesENSOPoVSIGAntarc(NAME, var_ninio_all-var_all,
                                      var_ninia_all- var_all,
                                      var_ninio_WPV - var_all,
                                      var_ninia_WPV - var_all,
                                      var_ninio_SPV - var_all,
                                      var_ninia_SPV - var_all,
                                      np.sqrt(SS_ninio_all + SS_all),
                                      np.sqrt(SS_ninia_all + SS_all),
                                      np.sqrt(SS_ninio_WPV + SS_all),
                                      np.sqrt(SS_ninia_WPV + SS_all),
                                      np.sqrt(SS_ninio_SPV + SS_all),
                                      np.sqrt(SS_ninia_SPV + SS_all),
                                      nn_ninio_all + nn_all - 2,
                                      nn_ninia_all + nn_all - 2,
                                      nn_ninio_WPV + nn_all - 2,
                                      nn_ninia_WPV + nn_all - 2,
                                      nn_ninio_SPV + nn_all - 2,
                                      nn_ninia_SPV + nn_all - 2,
                                      VAR.latitude, VAR.longitude, tit, filename)

