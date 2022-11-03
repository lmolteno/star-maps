import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import utils

WIDTH = 30 # degrees
HEIGHT = 20.0  # degrees
CENTRE_RA = 12 * 15  # degrees
CENTRE_DE = -60.0  # degrees
TRANSFORM = ccrs.Stereographic(central_latitude=CENTRE_DE, central_longitude=CENTRE_RA)

#ccrs.SouthPolarStereo()

#ccrs.AzimuthalEquidistant(central_longitude=(12 + 51.4/60)*15, central_latitude=27.13)
#ccrs.RotatedPole(pole_longitude=12*15 + 51.4/60 * 15, pole_latitude=27.13)

sm = utils.SkyMap(CENTRE_RA, CENTRE_DE, WIDTH, HEIGHT, TRANSFORM)\
        .get_hipparcos()\
        .filter_stars(mag_limit=8)\
        .plot_stars()\
        .constellation_paths()\
        .constellation_boundaries() \
        .border() \
        .ecliptic() \
        .labels()
# .ngc() \

sm.store.close()

# sm.ax.set_extent([0, 350, -90, 0], crs=ccrs.PlateCarree())

# print("Showing...")

# plt.savefig('out.svg')

sm.tick(45, CENTRE_DE + HEIGHT/2, 5, "meridian", "positive")
sm.show()
