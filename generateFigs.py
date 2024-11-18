from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from align import dispFITS

H = fits.open("H.fit")[0]
O = fits.open("O.fit")[0]
DS1_Ha = fits.open('DS-1_Ha.fit')[0]
DS1_OIII = fits.open('DS-1_OIII.fit')[0]
DS2_Ha = fits.open('DS-2_Ha.fit')[0]
DS2_OIII = fits.open('DS-2_OIII.fit')[0]
DS3_Ha = fits.open('DS-3_Ha.fit')[0]
DS3_OIII = fits.open('DS-3_OIII.fit')[0]


dispFITS(H,1,1,"Calibrated Halpha")
dispFITS(O,1,1,"Calibrated OIII")

dispFITS(DS1_Ha, 1, 10, 'Halpha Dataset 1')
dispFITS(DS2_Ha, 1, 10, 'Halpha Dataset 2')
dispFITS(DS3_Ha, 1, 10, 'Halpha Dataset 3')
dispFITS(DS1_OIII, 1, 10, 'OIII Dataset 1')
dispFITS(DS2_OIII, 1, 10, 'OIII Dataset 2')
dispFITS(DS3_OIII, 1, 10, 'OIII Dataset 3')

MB1 = fits.open('master_bias_1.fit')[0]
MB2 = fits.open('master_bias_2.fit')[0]
MB3 = fits.open('master_bias_3.fit')[0]
MD1 = fits.open('master_dark_1.fit')[0]
MD2 = fits.open('master_dark_2.fit')[0]
MD3 = fits.open('master_dark_3.fit')[0]
F1_Ha = fits.open('master_flat_ha_1.fit')[0]
F2_Ha = fits.open('master_flat_ha_2.fit')[0]
F3_Ha = fits.open('master_flat_ha_3.fit')[0]
F1_OIII = fits.open('master_flat_o.fits')[0]

dispFITS(MB3, 1, 10, 'Master Bias Dataset 3')
dispFITS(MD3, 1, 10, 'Master Dark Dataset 3')
dispFITS(F3_Ha, 1, 10, 'Master Flat Halpha Dataset 3')
dispFITS(F1_OIII, 1, 10, 'Master Flat OIII Dataset 1')

plt.show()