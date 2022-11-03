import astropy.coordinates as coord
import astropy.units as u
import cartopy.crs as ccrs
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astroquery.vizier import Vizier

CATALOG = "I/239/hip_main"  # hip catalog
WIDTH = 15.0  # degrees
HEIGHT = 15.0  # degrees
CENTRE_RA = 180  # degrees
CENTRE_DEC = -60.0  # degrees
TRANSFORM = ccrs.Gnomonic(central_latitude=CENTRE_DEC, central_longitude=CENTRE_RA)

centre = coord.SkyCoord(CENTRE_RA * u.deg, CENTRE_DEC * u.deg, frame='icrs')

print('Getting things')
# result_table = Vizier(row_limit=1000).query_region(
#     centre,
#     width=WIDTH * u.deg,
#     height=HEIGHT * u.deg,
#     column_filters={"Vmag": "<=15"},
#     catalog=CATALOG)
#
#
# def str_to_hms(st):
#     hms = st.split(' ')
#     return f"{hms[0]}h{hms[1]}m{hms[2]}s"
#
#
# def str_to_dms(st):
#     dms = st.split(' ')
#     return f"{dms[0]}d{dms[1]}m{dms[2]}s"

# spts = {
#     'A': 'C0',
#     'B': 'C1',
#     'F': 'C2',
#     'G': 'C3',
#     'K': 'C4',
#     'M': 'C5'
# }
#
#
# def map_markersize(magnitude):
#     return np.float_power(1.8, (10 - magnitude))
#
#
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection=TRANSFORM)
#
# CMAP = plt.get_cmap('plasma')
# b_v = result_table[0]['B-V'].data
# cnorm = colors.Normalize(vmin=np.min(b_v), vmax=np.max(b_v))
# scalarMap = cm.ScalarMappable(norm=cnorm, cmap=CMAP)
# ra = result_table[0]['RAICRS'].data
# dec = result_table[0]['DEICRS'].data
# mags = result_table[0]['Vmag'].data
# points = ax.scatter(-ra, dec, s=map_markersize(mags), color=scalarMap.to_rgba(b_v), transform=ccrs.PlateCarree())
#
# top_left = centre.directional_offset_by(np.arctan2(HEIGHT / 2, WIDTH / 2) * u.rad,
#                                         np.linalg.norm([WIDTH / 2, HEIGHT / 2]) * u.deg)
#
# top_right = centre.directional_offset_by(np.arctan2(HEIGHT / 2, -WIDTH / 2) * u.rad,
#                                          np.linalg.norm([WIDTH / 2, HEIGHT / 2]) * u.deg)
#
# bottom_right = centre.directional_offset_by(np.arctan2(-HEIGHT / 2, -WIDTH / 2) * u.rad,
#                                             np.linalg.norm([WIDTH / 2, HEIGHT / 2]) * u.deg)
#
# bottom_left = centre.directional_offset_by(np.arctan2(-HEIGHT / 2, WIDTH / 2) * u.rad,
#                                            np.linalg.norm([WIDTH / 2, HEIGHT / 2]) * u.deg)
#
# # noinspection PyTypeChecker
# border_path = mpath.Path(np.array([[-CENTRE_RA - WIDTH / 2, CENTRE_DEC - HEIGHT / 2],
#                                    [-CENTRE_RA + WIDTH / 2, CENTRE_DEC - HEIGHT / 2],
#                                    [-CENTRE_RA + WIDTH / 2, CENTRE_DEC + HEIGHT / 2],
#                                    [-CENTRE_RA - WIDTH / 2, CENTRE_DEC + HEIGHT / 2]]))
# # print(vertices)
#
# # noinspection PyTypeChecker
# # border_path = mpath.Path(vertices)
#
# # border_path = mpath.Path(np.array([[np.min(-ra), np.min(dec)],
# #                                    [np.min(-ra), np.max(dec)],
# #                                    [np.max(-ra), np.max(dec)],
# #                                    [np.max(-ra), np.min(dec)]]))
# #
# ax.set_boundary(border_path, transform=ccrs.PlateCarree())
# gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
#
#
# plt.show()
