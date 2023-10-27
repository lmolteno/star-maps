from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord, ICRS, FK5
from astropy import units as u
from enum import Enum
import matplotlib.path as mpath
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
import numpy as np
from typing import Literal
from collections import defaultdict
from shapely.geometry.multipoint import MultiPoint


class NGCType(Enum):
    OPEN_CLUSTER = 1
    CLUSTER_NEBULOSITY = 6
    GLOBULAR_CLUSTER = 2
    DIFFUSE_NEBULA = 3
    OBJECT_IN_LMC = 8
    PLANETARY_NEBULA = 4
    OBJECT_IN_SMC = 9
    GALAXY = 5
    UNVERIFIED_SOUTHERN = 0


object_colors = {
    NGCType.OPEN_CLUSTER: 'red',
    NGCType.CLUSTER_NEBULOSITY: 'orange',
    NGCType.GLOBULAR_CLUSTER: 'yellow',
    NGCType.DIFFUSE_NEBULA: 'green',
    NGCType.OBJECT_IN_LMC: 'blue',
    NGCType.PLANETARY_NEBULA: 'indigo',
    NGCType.OBJECT_IN_SMC: 'violet',
    NGCType.GALAXY: 'gray',
    NGCType.UNVERIFIED_SOUTHERN: 'purple'
}

spectral_colors = defaultdict(lambda: '#000', {
    'M': '#cb3127',
    'N': '#cb3127',
    'R': '#cb3127',
    'S': '#cb3127',
    'K': '#c95606',
    'G': '#be9627',
    'F': '#79D67B',
    'A': '#72A5D6',
    'O': '#D648BF',
    'W': '#D648BF',  # Wolf-Rayet stars
    'B': '#D648BF',
    'C': '#921919',  # Carbon stars
})


class SkyMap:

    def __init__(self, centre_ra, centre_de, width, height, transform: ccrs.Projection, logging=True,
                 star_colors: Literal["bv", "type"] = "type", figsize=(16, 9)):
        self.ngc_cat = None
        self.ecliptic_coords = None
        self.stars = None
        self.const_path = None
        self.const_bound_path = None
        self.hipparcos = None
        self.centre_ra = centre_ra
        self.centre_de = centre_de
        self.width = width
        self.height = height
        self.transform = transform
        self.logging = logging
        self.star_colors = star_colors

        self.store = pd.HDFStore('data.h5')
        self.fig = plt.figure(figsize=figsize)
        self.ax: GeoAxes = self.fig.add_subplot(1, 1, 1, projection=transform)
        self.set_extent()

    def _process_binaries(self, h):
        """
            turn 18439AB -> 18439
        """
        if h[-1].isalpha():
            return self._process_binaries(h[:-1])  # god i'm lazy
        return int(h)

    col_specs = [
        (8, 14),  # identifier
        (41, 46),  # Vmag
        (47, 48),  # varflag
        (51, 63),  # radeg
        (64, 76),  # dedeg
        (245, 250),  # b-v
        (435, 446),  # sptype
    ]

    col_names = [
        "identifier",
        "magnitude",
        "variable",
        "ra",
        "de",
        "bv",
        "spectral_type"
    ]

    @staticmethod
    def _map_markersize(magnitude):
        return np.float_power(1.8, (8 - magnitude))

    @staticmethod
    def _angle_between(angle, lower, upper):
        return (angle - lower) % 360 <= (upper - lower) % 360

    def get_hipparcos(self, convert=True):
        if '/hipparcos' not in self.store.keys():
            if self.logging: print("Loading Hipparcos catalog...", end="", flush=True)
            hipparcos = pd.read_fwf("hip_main.dat", colspecs=self.col_specs, names=self.col_names,
                                    index_col=['identifier'])
            hipparcos.variable = pd.notna(hipparcos.variable)
            hipparcos.spectral_type = hipparcos.spectral_type.astype(str)

            if convert:
                if self.logging: print(" converting to J1875...", end="", flush=True)
                coords = SkyCoord(hipparcos.ra, hipparcos.de, frame='icrs', unit=u.deg).transform_to(
                    FK5(equinox='J1875'))
                hipparcos.ra = coords.ra.to(u.deg).value
                hipparcos.de = coords.dec.to(u.deg).value

            if self.logging: print(f" loaded {len(hipparcos)} stars")
            self.store['hipparcos'] = hipparcos
            self.hipparcos = hipparcos
        else:
            self.hipparcos = self.store['hipparcos']

        return self

    def filter_stars(self, mag_limit=12):
        if self.hipparcos is None:
            raise Exception("Cannot filter stars, Hipparcos catalog has not been loaded")

        if self.logging: print("Filtering...", end="", flush=True)

        stars_in_map = self.hipparcos[
            self._angle_between(self.hipparcos.ra,
                                (self.centre_ra - self.width / 2) % 360,
                                (self.centre_ra + self.width / 2) % 360) &
            self.hipparcos.de.between(self.centre_de - self.height / 2,
                                      self.centre_de + self.height / 2) &
            (self.hipparcos.magnitude < mag_limit)]
        if self.logging: print(f" {len(stars_in_map)} stars in view")
        self.stars = stars_in_map
        return self

    def plot_stars(self):
        if self.logging: print("Plotting stars...", end="", flush=True)

        if self.star_colors == 'bv':
            cmap = plt.get_cmap('plasma')
            cnorm = colors.Normalize(vmin=np.min(self.stars.bv), vmax=np.max(self.stars.bv))
            scalar_map = cm.ScalarMappable(norm=cnorm, cmap=cmap)
            color = scalar_map.to_rgba(self.stars.bv)
        elif self.star_colors == 'type':
            color = self.stars.spectral_type.map(lambda t: spectral_colors[t[0]])
        else:
            color = ['#000'] * len(self.stars)

        self.ax.scatter(-self.stars.ra, self.stars.de,
                        s=self._map_markersize(self.stars.magnitude),
                        color=color,
                        transform=ccrs.PlateCarree(),
                        edgecolor='#000',
                        linewidths=0.15)

        if self.logging: print(" done")

        return self

    def border(self):
        if self.logging: print("Setting border...", end="", flush=True)
        if self.width > 355 and False:
            theta = np.linspace(0, 2 * np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            self.ax.set_boundary(circle, transform=self.ax.transAxes)
        else:
            ra_interpolated = np.linspace(-self.centre_ra - self.width / 2, -self.centre_ra + self.width / 2, 200)
            top_line = np.vstack(
                (ra_interpolated, (self.centre_de + self.height / 2) * np.ones_like(ra_interpolated))).T
            # bot_line = np.vstack(
            #     (np.flip(ra_interpolated, axis=0), (self.centre_de - self.height / 2) * np.ones_like(ra_interpolated))).T

            vertices = self._border_path().vertices  # np.vstack((top_line, bot_line))
            if self.width > 359:
                vertices = np.array(top_line)
            rect = mpath.Path(vertices=vertices)
            proj_to_data = ccrs.PlateCarree()._as_mpl_transform(self.ax) - self.ax.transData
            rect_in_target = proj_to_data.transform_path(rect)

            self.ax.set_boundary(rect_in_target)
            # self.ax.set_boundary(self._border_path(), transform=ccrs.PlateCarree())

        if self.logging: print(" done")
        return self

    def constellation_paths(self):
        if '/const_path' not in self.store.keys():
            if self.logging: print(f"Loading constellation path...", end="", flush=True)

            const_lines = pd.read_csv("clines.dat", sep=" ", comment="#", names=['code', 'hd'])

            hd_hip = pd.read_csv("hd_hip.dat", sep="\t", names=['hd', 'hip', 'const'])
            hd_hip = hd_hip.loc[hd_hip['hd'].str.isalnum() & hd_hip['hip'].str.isalnum()]

            hd_hip['hip'] = hd_hip['hip'].map(self._process_binaries)
            hd_hip['hd'] = hd_hip['hd'].map(self._process_binaries)
            hd_hip['hd'] = pd.to_numeric(hd_hip['hd'])
            hd_hip['hip'] = pd.to_numeric(hd_hip['hip'])

            failed_ones = []

            def id_to_ra_dec(hd):
                hd_hip_row = hd_hip.loc[hd_hip['hd'] == hd]
                star: pd.Series = self.hipparcos.loc[hd_hip_row['hip']]
                if star.empty:
                    failed_ones.append(f"HD {hd}")
                    if self.logging:
                        print(f"Error finding HD {hd}, which had HIP{hd_hip_row['hip']}")
                star = star.iloc[0]
                return -float(star['ra']), float(star['de'])

            if len(failed_ones) > 0:
                if self.logging: print(" getting failed HD from SIMBAD ", end="", flush=True)
                custom_simbad = Simbad()
                custom_simbad._VOTABLE_FIELDS = ["typed_id", "id(HIP)"]
                missing_results = custom_simbad.query_objects(failed_ones)

                print("please append the following to your hd_hip mapping\n-----")
                for row in missing_results.iterrows():
                    print(f"{row[0].split(' ')[1]}\t{row[1].split(' ')[1]}\tUnk")
                print("-----")
                exit()

            code_map = {
                'M': mpath.Path.MOVETO,
                'D': mpath.Path.LINETO
            }

            ra, dec, codes = zip(
                *[(*id_to_ra_dec(hd), code_map[code]) for hd, code in
                  list(zip(const_lines['hd'], const_lines['code']))])

            const_path = pd.DataFrame({'ra': ra, 'de': dec, 'code': codes})
            self.store['const_path'] = const_path

            if self.logging: print(" done")
            self.const_path = const_path
        else:
            self.const_path = self.store['const_path']

        vertices = list(zip(self.const_path['ra'], self.const_path['de']))
        codes = self.const_path['code']

        const_patch = mpatches.PathPatch(mpath.Path(vertices, codes),
                                         transform=ccrs.Geodetic(),
                                         edgecolor='#000',
                                         fill=False,
                                         zorder=-1,
                                         linewidth=0.5)
        self.ax.add_patch(const_patch)

        return self

    def constellation_boundaries(self):
        if '/const_bound_path' not in self.store.keys():
            const_bounds = pd.read_csv("cboundlines.dat", sep=" ", names=["ra", "de", "seg_id"])
            ra, dec = const_bounds['ra'] * 15, const_bounds['de']

            code = []
            prev_seg_id = ""
            for idx, row in const_bounds.iterrows():
                code.append(mpath.Path.LINETO if row['seg_id'] == prev_seg_id else mpath.Path.MOVETO)
                prev_seg_id = row['seg_id']

            self.const_bound_path = pd.DataFrame({'ra': ra, 'de': dec, 'code': code})
            self.store['const_bound_path'] = self.const_bound_path
        else:
            self.const_bound_path = self.store['const_bound_path']

        vertices = list(zip(-self.const_bound_path['ra'], self.const_bound_path['de']))
        codes = self.const_bound_path['code']

        const_boundary_patch = mpatches.PathPatch(mpath.Path(vertices, codes),
                                                  transform=ccrs.Geodetic(),
                                                  edgecolor='#000',
                                                  fill=False,
                                                  zorder=-1.1,
                                                  linewidth=0.5,
                                                  linestyle=(0, (5, 10)))
        self.ax.add_patch(const_boundary_patch)

        return self

    def magplot(self):
        fig2 = plt.figure()
        mags = np.linspace(-1, 12, 50)
        ax2 = fig2.add_subplot(1, 1, 1)
        ax2.plot(mags, self._map_markersize(mags))
        x = np.arange(-1, 10)
        ax2.scatter(x, np.zeros_like(x), s=self._map_markersize(x))
        ax2.grid()
        ax2.set_ylabel('marker size')
        ax2.set_xlabel('magnitude')

    def ecliptic(self, n=100):
        if '/ecliptic' not in self.store.keys():
            lon = np.linspace(0, 360, n)
            lat = np.zeros_like(lon)
            coords = SkyCoord(lon, lat, frame='geocentrictrueecliptic', unit=u.deg).transform_to('icrs')
            ecliptic = pd.DataFrame({'ra': coords.ra.to(u.deg).value, 'de': coords.dec.to(u.deg).value})
            self.store['ecliptic'] = ecliptic
            self.ecliptic_coords = ecliptic
        else:
            self.ecliptic_coords = self.store['ecliptic']

        self.ax.plot(-self.ecliptic_coords.ra, self.ecliptic_coords.de,
                     linewidth=0.5,
                     linestyle='--',
                     color='k',
                     zorder=-1.1,
                     transform=ccrs.Geodetic())

        return self

    def ngc(self, mag_limit=13, convert=True):
        if '/ngc' not in self.store.keys():
            if self.logging: print("Loading NGC...", end="", flush=True)
            ngc_table = Vizier(row_limit=-1).query_constraints(catalog='VII/1B/catalog', Type='!=7')[0]
            ngc: pd.DataFrame = ngc_table.to_pandas()
            coords = SkyCoord(ngc['_RA.icrs'], ngc['_DE.icrs'], frame='icrs', unit=(u.hourangle, u.deg))
            if convert:
                coords = coords.transform_to(FK5(equinox='J1875'))
            ngc['ra'] = coords.ra.to(u.deg).value
            ngc['de'] = coords.dec.to(u.deg).value
            ngc['m_NGC'] = ngc['m_NGC'].to_string()

            ngc.drop(['RA1975', 'DE1975', 'Xpos', 'Ypos', 'r_Mag', 'Notes', '_RA.icrs', '_DE.icrs'],
                     inplace=True,
                     axis='columns')
            self.store['ngc'] = ngc
            print(" done")
        else:
            ngc = self.store['ngc']

        print("Filtering NGC...", end="", flush=True)
        ngc['Type'] = ngc['Type'].apply(lambda x: [NGCType(int(xi)) for xi in str(x)])
        self.ngc_cat = ngc

        ngc = ngc.loc[ngc['Mag'] <= mag_limit]
        ngc = ngc[
            self._angle_between(ngc.ra,
                                (self.centre_ra - self.width / 2) % 360,
                                (self.centre_ra + self.width / 2) % 360) &
            ngc.de.between(self.centre_de - self.height / 2,
                           self.centre_de + self.height / 2)]
        if self.logging: print(f" {len(ngc)} NGC objects in view")

        for ngctype in NGCType:
            ngc_objs = ngc.loc[ngc['Type'].map(lambda x: ngctype in x)]
            self.ax.scatter(-ngc_objs['ra'], ngc_objs['de'],
                            transform=ccrs.PlateCarree(),
                            color=object_colors[ngctype])

        return self

    def labels(self, maj_ra_inc=1, maj_de_inc=5, min_ra_inc=1 / 2, min_de_inc=0 * u.deg):
        top_de, bottom_de = self.centre_de + self.height / 2, self.centre_de - self.height / 2
        right_ra, left_ra = self.centre_ra + self.width / 2, self.centre_ra - self.width / 2

        transform = ccrs.PlateCarree()._as_mpl_transform(self.ax) - self.ax.transData

        gridline_params = dict(
            transform=ccrs.Geodetic(),
            edgecolor='#888',
            fill=False,
            zorder=-5,
            linewidth=0.5
        )

        major_fontsize = 12
        minor_fontsize = 10
        major_tick_length = 5
        minor_tick_length = 3

        major_de = [*filter(lambda de: bottom_de <= de <= top_de, range(-90, 90, maj_de_inc))]  # in degrees
        full_ra = np.arange(0, 24, min_ra_inc)  # in hours
        filtered_ra = np.extract(np.logical_and(left_ra <= full_ra * 15, full_ra * 15 <= right_ra), full_ra)
        for ra in filtered_ra:
            fontsize = minor_fontsize
            text_offset = 3
            if ra % maj_ra_inc == 0:
                text_offset = 5  # points
                fontsize = major_fontsize

            ra_coord = -ra * 15  # to degrees

            ra_str = f"{ra:.0f}$^{{\\rm h}}$" if ra % maj_ra_inc == 0 else f"{ra * 60 % 60:.0f}$^{{\\rm m}}$"

            for de in top_de, bottom_de:
                if self.width > 350:
                    if self.centre_de > 90:
                        if de == top_de: break
                    else:
                        if de == bottom_de: break
                points = transform.transform([[ra_coord - 0.01, de], [ra_coord + 0.01, de]])
                angle = angle_between(points[0], points[1]) + np.pi

                text_offset_x, text_offset_y = text_offset * np.sin(angle), \
                                               text_offset * np.cos(angle)
                self.ax.annotate(ra_str,
                                 xy=(ra_coord, de),
                                 xytext=(-text_offset_x if de == top_de else text_offset_x,
                                         text_offset_y if de == top_de else -text_offset_y),
                                 textcoords='offset points',
                                 verticalalignment="bottom" if de == top_de else "top",
                                 horizontalalignment="center",
                                 transform=ccrs.PlateCarree(),
                                 rotation=np.rad2deg(angle),
                                 clip_on=False,
                                 rotation_mode='anchor',
                                 fontsize=fontsize)

            if ra % maj_ra_inc == 0:
                path = mpath.Path([[ra_coord, top_de], [ra_coord, bottom_de]]).interpolated(100)
                self.ax.add_patch(mpatches.PathPatch(path, **gridline_params))
                self.tick(-ra_coord, top_de, major_tick_length, "meridian", "positive")
                self.tick(-ra_coord, bottom_de, major_tick_length, "meridian", "negative")
            else:
                self.tick(-ra_coord, top_de, minor_tick_length, "meridian", "positive")
                self.tick(-ra_coord, bottom_de, minor_tick_length, "meridian", "negative")

        if self.width < 350:
            for ra_h in left_ra, right_ra:
                text_offset = 5  # points
                for de in major_de:
                    ra = -ra_h
                    points = transform.transform([[ra - 0.01, de], [ra + 0.01, de]])
                    angle = angle_between(points[0], points[1]) + np.pi
                    text_offset_y, text_offset_x = text_offset * np.sin(angle), \
                                                   text_offset * np.cos(angle)

                    self.ax.annotate(f"{de}$^o$",
                                     xy=(ra, de),
                                     xytext=(-text_offset_x if ra_h != left_ra else text_offset_x,
                                             text_offset_y if ra_h == left_ra else -text_offset_y),
                                     textcoords='offset points',
                                     verticalalignment="center",
                                     horizontalalignment="left" if ra_h == left_ra else "right",
                                     transform=ccrs.PlateCarree(),
                                     rotation=np.rad2deg(angle),
                                     clip_on=False,
                                     rotation_mode='anchor')

            for de in major_de:
                path = mpath.Path([[right_ra, de], [left_ra, de]]).interpolated(100)
                self.ax.add_patch(mpatches.PathPatch(path, **gridline_params))
                self.tick(right_ra, de, minor_tick_length, "parallel", "negative")
                self.tick(left_ra, de, minor_tick_length, "parallel", "positive")

        return self

    # Meridian: ra, Parallel: de. positive: tick goes from ra, dec,
    # to positive whatever coordinate specified in orientation
    def tick(self, ra, dec, length, orientation: Literal["meridian", "parallel"],
             direction: Literal["positive", "negative"]):

        trans = (self.fig.dpi_scale_trans +
                 mtransforms.ScaledTranslation(-ra, dec, ccrs.PlateCarree()._as_mpl_transform(self.ax)))

        angle = self._get_angle(ra, dec, orientation)

        direction_mult = 1 if direction == "positive" else -1

        tick = np.array([[0, 0], [np.cos(angle), np.sin(angle)]]) * direction_mult * length / 72  # divide by 72 for pt

        self.ax.add_patch(mpatches.PathPatch(mpath.Path(tick),
                                             transform=trans,
                                             clip_on=False))

    def _get_angle(self, ra: float, dec: float, orientation: Literal["meridian", "parallel"]) -> float:
        angle_points = [[-ra - 0.01, dec],
                        [-ra + 0.01, dec]] if orientation == "parallel" else [[-ra, dec - 0.01],
                                                                              [-ra, dec + 0.01]]

        transformed_points = (ccrs.PlateCarree()._as_mpl_transform(self.ax) - self.ax.transData).transform(angle_points)

        return angle_between(*transformed_points) - np.pi

    def wide_border(self, distance=42):  # distance in pt
        corners = [[self.centre_ra + self.width / 2, self.centre_de + self.height / 2],  # top right
                   [self.centre_ra - self.width / 2, self.centre_de + self.height / 2],  # top left
                   [self.centre_ra - self.width / 2, self.centre_de - self.height / 2],  # bottom left
                   [self.centre_ra + self.width / 2, self.centre_de - self.height / 2],  # bottom right
                   [self.centre_ra + self.width / 2, self.centre_de + self.height / 2]]  # top right

        vertices = []

        for ra, dec in corners:
            transform = self.ax.transData
            angle_m = self._get_angle(ra, dec, 'meridian')
            dir_m = 1 if dec > self.centre_de else -1
            vec_m = dir_m * np.array([np.cos(angle_m), np.sin(angle_m)])

            angle_p = self._get_angle(ra, dec, 'parallel')
            dir_p = 1 if ra < self.centre_ra else -1
            vec_p = dir_p * np.array([np.cos(angle_p), np.sin(angle_p)])

            vert = np.array([ra, dec]) + ((vec_m + vec_p) * distance/72)
            print(vert)
            trans_vert = ccrs.PlateCarree().transform_point(*vert, self.transform)
            print(trans_vert)
            self.ax.add_patch(mpatches.Circle(vert, 1/72, facecolor='red', transform=ccrs.PlateCarree()._as_mpl_transform(self.ax), clip_on=False))

        # trans = self.fig.dpi_scale_trans + ccrs.PlateCarree()._as_mpl_transform(self.ax)
        # for idx, vert in enumerate(vertices):
        #     vert[1] *= -1
        #     print(vert)

        return self

    def _border_path(self, interpolated=True):
        border = mpath.Path([[-self.centre_ra + self.width / 2, self.centre_de + self.height / 2],  # top right
                             [-self.centre_ra, self.centre_de + self.height / 2],  # top centre
                             [-self.centre_ra - self.width / 2, self.centre_de + self.height / 2],  # top left
                             [-self.centre_ra - self.width / 2, self.centre_de - self.height / 2],  # bottom left
                             [-self.centre_ra, self.centre_de - self.height / 2],  # bottom centre
                             [-self.centre_ra + self.width / 2, self.centre_de - self.height / 2],  # bottom right
                             [-self.centre_ra + self.width / 2, self.centre_de + self.height / 2]])  # top right
        if interpolated:
            return border.interpolated(200)
        return border

    def set_extent(self):
        if self.width > 355:
            if self.height < 90:
                self.ax.set_extent(
                    [-180, 180, self.centre_de, self.centre_de + (-1 if self.centre_de > 90 else 1) * self.height / 2],
                    crs=ccrs.PlateCarree())
            else:
                self.ax.set_extent([-180, 180, self.centre_de - self.height / 2, self.centre_de + self.height / 2],
                                   crs=ccrs.PlateCarree())
        else:
            transform = (ccrs.PlateCarree()._as_mpl_transform(self.ax) - self.ax.transData)
            transformed_points = transform.transform(self._border_path().vertices)
            ra, dec = np.moveaxis(transformed_points, 1, 0)
            self.ax.set_ylim(np.min(dec), np.max(dec))
            self.ax.set_xlim(np.min(ra), np.max(ra))

    @staticmethod
    def show():
        plt.show()


def angle_between(p1, p2):
    return np.arctan2(*(p1 - p2)[::-1])
