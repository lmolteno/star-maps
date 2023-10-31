import matplotlib
matplotlib.use('qtagg')
import cartopy.crs as ccrs

import utils

WIDTH = 30 # degrees
HEIGHT = 20  # degrees
CENTRE_RA = 21 * 15  # degrees
CENTRE_DE = -40.0  # degrees
TRANSFORM = ccrs.SouthPolarStereo(central_longitude=-CENTRE_RA)

#ccrs.AzimuthalEquidistant(central_longitude=(12 + 51.4/60)*15, central_latitude=27.13)
#ccrs.RotatedPole(pole_longitude=12*15 + 51.4/60 * 15, pole_latitude=27.13)

sm = utils.SkyMap(CENTRE_RA, CENTRE_DE, WIDTH, HEIGHT, TRANSFORM) \
        .border() \
        .wide_border()
        # .get_hipparcos()\
        # .filter_stars(mag_limit=8)\
        # .plot_stars()\
        # .constellation_paths()\
        # .constellation_boundaries() \
        # .ecliptic() \
        # .labels() \
        # .ngc() \

sm.store.close()

# sm.ax.set_extent([0, 350, -90, 0], crs=ccrs.PlateCarree())

# print("Showing...")


sm.tick(45, CENTRE_DE + HEIGHT/2, 5, "meridian", "positive")
# plt.savefig('south-polar-section.png')
sm.show()
