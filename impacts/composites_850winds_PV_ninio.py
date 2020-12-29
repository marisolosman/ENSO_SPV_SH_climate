#composites on PV years conditioned on ENSO strength
import sys
import numpy as np
import xarray as xr
import os
import plots

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '/pikachu/datos/osman/assessment_SH_zonal_asymmetries/data/'
FIG_PATH = '/pikachu/datos/osman/assessment_SH_zonal_asymmetries/figures/impacts/'
FILE_VAR = '850winds_s4_aug_feb.nc4'
FILE_NINIO_S4 = 'ninio34_monthly.nc4'
FILE_PV_S4 = 'SPV_index.nc4'
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for en years

index_ninio = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0',
                                                                     interpolation='linear')
index_ninia = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0',
                                                                     interpolation='linear')

# PV intensity during all years
#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0',
                                                                    interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0',
                                                                    interpolation='linear')
#PoV during ninio years
index_WSPV_ninio = np.logical_and(index_ninio.values, index_SPV_upper.values)
index_SSPV_ninio = np.logical_and(index_ninio.values, index_SPV_lower.values)

#PoV during ninia years
index_WSPV_ninia = np.logical_and(index_ninia.values, index_SPV_upper.values)
index_SSPV_ninia = np.logical_and(index_ninia.values, index_SPV_lower.values)

nn_WSPV_all = np.sum(index_SPV_upper.values)
nn_SSPV_all = np.sum(index_SPV_lower.values)

nn_WSPV_ninio = np.sum(index_WSPV_ninio)
nn_SSPV_ninio = np.sum(index_SSPV_ninio)

nn_WSPV_ninia = np.sum(index_WSPV_ninia)
nn_SSPV_ninia = np.sum(index_SSPV_ninia)
nn_all = np.shape(ninio34.ninio34_index.values)[0]
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

VAR = xr.open_dataset(PATH_DATA + FILE_VAR, chunks={'longitude':10})
VAR = VAR.sel(latitude=slice(-10, -90))
for i in np.arange(0, 7):
    var_WSPV_ninio = np.mean(VAR.isel(month=i, realiz=index_WSPV_ninio, drop=True),
                            axis=0).to_array().squeeze()
    SS_WSPV_ninio = np.var(VAR.isel(month=i, realiz=index_WSPV_ninio, drop=True),
                          axis=0).to_array().squeeze()/nn_WSPV_ninio
    var_SSPV_ninio = np.mean(VAR.isel(month=i, realiz=index_SSPV_ninio, drop=True),
                            axis=0).to_array().squeeze()
    SS_SSPV_ninio = np.var(VAR.isel(month=i, realiz=index_SSPV_ninio, drop=True),
                          axis=0).to_array().squeeze()/nn_SSPV_ninio
    var_WSPV_ninia = np.mean(VAR.isel(month=i, realiz=index_WSPV_ninia, drop=True),
                            axis=0).to_array().squeeze()
    SS_WSPV_ninia = np.var(VAR.isel(month=i, realiz=index_WSPV_ninia, drop=True),
                          axis=0).to_array().squeeze()/nn_WSPV_ninia
    var_SSPV_ninia = np.mean(VAR.isel(month=i, realiz=index_SSPV_ninia, drop=True),
                            axis=0).to_array().squeeze()
    SS_SSPV_ninia = np.var(VAR.isel(month=i, realiz=index_SSPV_ninia, drop=True),
                          axis=0).to_array().squeeze()/nn_SSPV_ninia
    var_WSPV_all = np.mean(VAR.isel(month=i, realiz=index_SPV_upper.values, drop=True),
                            axis=0).to_array().squeeze()
    SS_WSPV_all = np.var(VAR.isel(month=i, realiz=index_SPV_upper.values),
                          axis=0).to_array().squeeze()/np.sum(index_SPV_upper.values)
    var_all = np.mean(VAR.isel(month=i, drop=True), axis=0).to_array().squeeze()
    SS_all = np.var(VAR.isel(month=i, drop=True), axis=0).to_array().squeeze() / nn_all
    var_SSPV_all = np.mean(VAR.isel(month=i, realiz=index_SPV_lower.values, drop=True),
                            axis=0).to_array().squeeze()
    SS_SSPV_all = np.var(VAR.isel(month=i, realiz=index_SPV_lower.values, drop=True),
                          axis=0).to_array().squeeze()/np.sum(index_SPV_lower.values)
    tit = 'Composites S4 850hPa-Winds Conditioned - ENSO - ' + month[i]
    filename = FIG_PATH + '850winds' + '_composites_SPoV_' + month[i] +'_ENSO.png'
    plots.Plot850WindsCompositesPoVENSO(var_WSPV_all - var_all,
                                      var_SSPV_all - var_all,
                                      var_WSPV_ninio - var_all,
                                      var_SSPV_ninio - var_all,
                                      var_WSPV_ninia - var_all,
                                      var_SSPV_ninia - var_all,
                                      VAR.latitude, VAR.longitude, tit, filename)
for i in np.arange(0, 5):
    VAR_s = VAR.isel(month=range(i, i+3)).mean(dim='month')
    var_WSPV_ninio = np.mean(VAR_s.isel( realiz=index_WSPV_ninio, drop=True),
                            axis=0).to_array().squeeze()
    SS_WSPV_ninio = np.var(VAR_s.isel( realiz=index_WSPV_ninio, drop=True),
                          axis=0).to_array().squeeze()/nn_WSPV_ninio
    var_SSPV_ninio = np.mean(VAR_s.isel( realiz=index_SSPV_ninio, drop=True),
                            axis=0).to_array().squeeze()
    SS_SSPV_ninio = np.var(VAR_s.isel( realiz=index_SSPV_ninio, drop=True),
                          axis=0).to_array().squeeze()/nn_SSPV_ninio
    var_WSPV_ninia = np.mean(VAR_s.isel( realiz=index_WSPV_ninia, drop=True),
                            axis=0).to_array().squeeze()
    SS_WSPV_ninia = np.var(VAR_s.isel( realiz=index_WSPV_ninia, drop=True),
                          axis=0).to_array().squeeze()/nn_WSPV_ninia
    var_SSPV_ninia = np.mean(VAR_s.isel( realiz=index_SSPV_ninia, drop=True),
                            axis=0).to_array().squeeze()
    SS_SSPV_ninia = np.var(VAR_s.isel( realiz=index_SSPV_ninia, drop=True),
                          axis=0).to_array().squeeze()/nn_SSPV_ninia
    var_WSPV_all = np.mean(VAR_s.isel( realiz=index_SPV_upper.values, drop=True),
                            axis=0).to_array().squeeze()
    SS_WSPV_all = np.var(VAR_s.isel( realiz=index_SPV_upper.values),
                          axis=0).to_array().squeeze()/np.sum(index_SPV_upper.values)
    var_all = np.mean(VAR_s, axis=0).to_array().squeeze()
    SS_all = np.var(VAR_s, axis=0).to_array().squeeze() / nn_all 
    var_SSPV_all = np.mean(VAR_s.isel( realiz=index_SPV_lower.values, drop=True),
                            axis=0).to_array().squeeze()
    SS_SSPV_all = np.var(VAR_s.isel( realiz=index_SPV_lower.values, drop=True),
                          axis=0).to_array().squeeze()/np.sum(index_SPV_lower.values)
    tit = 'Composites S4 850hPa-Winds Conditioned - ENSO - ' + seas[i]
    filename = FIG_PATH + '850winds' + '_composites_SPoV_' + seas[i] +'_ENSO.png'
    plots.Plot850WindsCompositesPoVENSO(var_WSPV_all - var_all,
                                        var_SSPV_all - var_all,
                                        var_WSPV_ninio - var_all,
                                        var_SSPV_ninio - var_all,
                                        var_WSPV_ninia - var_all,
                                        var_SSPV_ninia - var_all,
                                        VAR.latitude, VAR.longitude, tit, filename)

