#composites of EN events conditioned on PV strength
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

seas = ['SOND', 'JF']
VAR = xr.open_dataset(PATH_DATA + FILE_VAR, chunks={'longitude': 10})
VAR = VAR.sel(latitude=slice(0, -90))
for i in np.arange(0, 2):
        if i==0:
            VAR_s = VAR.isel(month=range(1, 5)).mean(dim='month')
        else:
            VAR_s = VAR.isel(month=range(5, 7)).mean(dim='month')
 
        var_ninio_WPV = np.mean(VAR_s.isel(realiz=index_ninio_WPV),
                         axis=0).to_array().squeeze()
        var_ninia_WPV = np.mean(VAR_s.isel(realiz=index_ninia_WPV),
                         axis=0).to_array().squeeze()
        var_ninio_SPV = np.mean(VAR_s.isel(realiz=index_ninio_SPV),
                         axis=0).to_array().squeeze()
        var_ninia_SPV = np.mean(VAR_s.isel(realiz=index_ninia_SPV),
                         axis=0).to_array().squeeze()
        var_ninio_all = np.mean(VAR_s.isel(realiz=index_ninio_all.values),
                         axis=0).to_array().squeeze()
        var_all = np.mean(VAR_s, axis=0).to_array().squeeze()
        var_ninia_all = np.mean(VAR_s.isel(realiz=index_ninia_all.values),
                         axis=0).to_array().squeeze()
        tit = 'Composites S4 850hPa-Winds Conditioned - SPoV - ' + seas[i]
        filename = FIG_PATH + '850winds' + '_composites_ENSO_' + seas[i] +'_SPoV.png'
        plots.Plot850WindsCompositesENSOPoV(var_ninio_all - var_all, var_ninia_all- var_all,
                                     var_ninio_WPV - var_all, var_ninia_WPV - var_all,
                                     var_ninio_SPV - var_all, var_ninia_SPV - var_all,
                                     VAR.u.latitude, VAR.u.longitude, tit, filename)

