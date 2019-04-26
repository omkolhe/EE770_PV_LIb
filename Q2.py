import datetime 
import os
import inspect
import math
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 

from pprint import pprint

import pvlib 

tmy3_data, tmy3_metadata = pvlib.iotools.read_tmy3('722784TYA.csv')
tmy3_data.index.name = 'Time'

loc = pvlib.location.Location.from_tmy(tmy3_metadata,tmy3_data);

surface_tilt = tmy3_metadata['latitude'];

if tmy3_metadata['longitude'] < 0:
	surface_azimuth = 180;
else:
	surface_azimuth = 0;

sun_position = pvlib.solarposition.get_solarposition(tmy3_data.index,loc.latitude,loc.longitude,loc.altitude,pressure=None,method='nrel_numpy',temperature=12)

POA = pvlib.irradiance.get_total_irradiance(surface_tilt,surface_azimuth,sun_position.zenith,sun_position.azimuth,dni=tmy3_data.DNI,ghi=tmy3_data.GHI,dhi=tmy3_data.DHI,dni_extra=tmy3_data.ETRN,airmass=None,albedo=tmy3_data.Alb,surface_type=None,model='isotropic',model_perez='allsitescomposite1990');

a = -3.56
b = -0.075

expo = 2.71828182846**(a + b*tmy3_data.Wspd)

module_temp = np.multiply(POA.poa_global,expo) + tmy3_data.DryBulb;

module_temp.head(24).plot();
plt.show();
