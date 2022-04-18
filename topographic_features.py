from surface_roughness import *
from rasterio.windows import Window
from typing import Tuple, List



def get_local_elevation_map(dataset: DatasetReader, lat: float, lon: float, range_km: int) -> Tuple[np.array, Window]:
    """
    Pulls the local elevation map near query coordinates as a numpy array, as well as the offsets relative to the global map.
    :param dataset: The open rasterio.DatasetReader
    :param lat: query latitude
    :param lon: query longitude
    :param range_km: size of the window in both dimensions (km)
    :return: the elevation map array and the Window object
    """
    km_per_degree = 111.3
    (window_height, window_width) = (int(range_km / (dataset.res[0] * km_per_degree)), int(range_km / (dataset.res[1] * km_per_degree)))
    row_offset, col_offset = dataset.index(lon, lat)
    col_offset -= window_width / 2
    row_offset -= window_height / 2
    window = Window(col_offset, row_offset, window_width, window_height)
    return dataset.read(1, window=window), window


def get_bearing(lat1, lon1, lat2, lon2):
    dLon = (lon2 - lon1)
    x = np.cos(np.radians(lat2)) * np.sin(np.radians(dLon))
    y = np.cos(np.radians(lat1)) * np.sin(np.radians(lat2)) - np.sin(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.cos(np.radians(dLon))
    brng = np.arctan2(x,y)

    return brng


def generate_profile_from_path(pt1, pt2, dataset: DatasetReader, elevation_map: np.array, window: Window, elev1, elev2, distance_km):
    bearing = get_bearing(pt1[0], pt1[1], pt2[0], pt2[1])
    d_vec, elev_vec = get_profile(dataset, elevation_map, window, pt1[0], pt1[1], distance_km, bearing)
    elev_vec[0] += elev1
    elev_vec[-1] += elev2
    return d_vec, elev_vec


def extract_topographic_features(d_vec: np.array, elev_vec: np.array):
    elev_adj = level_profile(d_vec, elev_vec)
    features = {
        "ra": ra_roughness(elev_adj),
        "rq": rq_roughness(elev_adj),
        "rp": rp_roughness(elev_adj),
        "rv": rv_roughness(elev_adj),
        "rz": rz_roughness(elev_adj),
        "rsk": rsk_roughness(elev_adj),
        "rku": rku_roughness(elev_adj),
        "deepest_barrier": deepest_barrier(elev_adj),
        "n_barriers": n_barriers(elev_adj)
    }
    return features



def get_features_for_receipt(pt1, pt2, dataset, elevation_map, window, elev1, elev2, distance_km):

    d_vec, elev_vec = generate_profile_from_path(pt1, pt2, dataset, elevation_map, window, elev1, elev2, distance_km)

    # get topographic features
    try:
        features = extract_topographic_features(d_vec, elev_vec)
    except (np.linalg.LinAlgError, SystemError):
    # hotspots are too close to create elevation profile
        features = {
            "ra": None,
            "rq": None,
            "rp": None,
            "rv": None,
            "rz": None,
            "rsk": None,
            "rku": None,
            "deepest_barrier": None,
            "n_barriers": None
        }
    return features

