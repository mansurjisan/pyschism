import numpy as numpy
import time
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import geopandas as gpd
import os
import urllib.request
import zipfile
import shapely.geometry as sgeom
import matplotlib.colors as mcolors
import logging
import struct



logging.basicConfig(level=logging.INFO)

def read_gr3_file(filename):
    try:
        with open(filename, 'r') as f:
            _ = f.readline()  # Discard first line (comment)
            ne, np = map(int, f.readline().split())

            nodes = numpy.empty((np, 3))
            for i in range(np):
                nodes[i] = list(map(float, f.readline().split()[1:4]))

            elements = []
            for _ in range(ne):
                elements.append(list(map(lambda x: int(x) - 1, f.readline().split()[2:])))

        print(f"Number of nodes: {np}")
        print(f"Number of elements: {ne}")
        return nodes, elements
    except Exception as e:
        print(f"Error reading file: {e}")
        import traceback
        traceback.print_exc()
        return None, None

def download_coastline_data():
    url = "https://naciscdn.org/naturalearth/10m/physical/ne_10m_coastline.zip"
    zip_path = "ne_10m_coastline.zip"
    shapefile_path = "ne_10m_coastline.shp"

    if not os.path.exists(shapefile_path):
        print("Downloading coastline data...")
        urllib.request.urlretrieve(url, zip_path)
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall()
        os.remove(zip_path)
    return shapefile_path

def download_gadm_data():
    url = "https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_USA_shp.zip"
    zip_path = "gadm41_USA_shp.zip"
    shapefile_path = "gadm41_USA_1.shp"

    if not os.path.exists(shapefile_path):
        print("Downloading GADM data...")
        urllib.request.urlretrieve(url, zip_path)
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall()
        os.remove(zip_path)
    return shapefile_path




def read_rangs_gshhs(file_path):
    """
    Read RANGS/GSHHS binary files and return coastline data.
    This function attempts to be more flexible with unknown file formats.
    
    :param file_path: Path to the RANGS/GSHHS file
    :return: List of coastline segments
    """
    with open(file_path, 'rb') as f:
        # Read the entire file content as integers
        data = numpy.fromfile(f, dtype='>i')
    
    # Assume the first few integers are header information
    header_size = 7  # Adjust this if needed
    header = data[:header_size]
    
    # The rest of the data should be coordinates
    points = data[header_size:].reshape(-1, 2)
    
    # Convert from microdegrees to degrees
    points = points / 1e6
    
    # Split into separate coastline segments if there are gaps
    coastlines = numpy.split(points, numpy.where(numpy.diff(points[:, 0]) > 1)[0] + 1)
    
    return coastlines


def plot_gr3(nodes, elements, output_file, min_lon=None, max_lon=None, min_lat=None, max_lat=None, show_land=True, title=None):
    start_time = time.time()  # Start timing

    if nodes is None or elements is None:
        logging.error("Cannot plot: invalid data")
        return

    logging.info("Preparing plot...")

    # Print diagnostic information
    logging.info(f"Node shape: {nodes.shape}")
    logging.info(f"Node min values: {numpy.nanmin(nodes, axis=0)}")
    logging.info(f"Node max values: {numpy.nanmax(nodes, axis=0)}")
    logging.info(f"Number of NaN nodes: {numpy.isnan(nodes).any(axis=1).sum()}")
    logging.info(f"Number of Inf nodes: {numpy.isinf(nodes).any(axis=1).sum()}")

    fig = plt.figure(figsize=(10, 8))

    # Use PlateCarree projection
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # Set map extent
    if all([min_lon, max_lon, min_lat, max_lat]):
        extent = [min_lon, max_lon, min_lat, max_lat]
    else:
        # Filter out NaN and Inf values
        valid_nodes = nodes[~numpy.isnan(nodes).any(axis=1) & ~numpy.isinf(nodes).any(axis=1)]
        if len(valid_nodes) == 0:
            logging.error("No valid node coordinates found")
            return
        
        margin = 1  # degree
        extent = [
            numpy.nanmin(valid_nodes[:, 0]) - margin,
            numpy.nanmax(valid_nodes[:, 0]) + margin,
            numpy.nanmin(valid_nodes[:, 1]) - margin,
            numpy.nanmax(valid_nodes[:, 1]) + margin
        ]

    # Check if extent is valid
    if any(numpy.isnan(extent)) or any(numpy.isinf(extent)):
        logging.error(f"Invalid extent calculated: {extent}")
        logging.info("Using default extent.")
        extent = [-180, 180, -90, 90]  # Default to full globe

    logging.info(f"Using extent: {extent}")

    try:
        ax.set_extent(extent)
    except ValueError as e:
        logging.error(f"Error setting extent: {e}")
        logging.info("Falling back to default extent")
        ax.set_global()

    if show_land:
        logging.info("Adding land features...")
        # Add land feature
        ax.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='face')
        
        # Add coastlines
        ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.5)
        
        # Add country borders
        ax.add_feature(cfeature.BORDERS, edgecolor='darkgray', linestyle=':', linewidth=0.5)

        # Download and plot GADM state boundaries
        gadm_path = download_gadm_data()
        if os.path.exists(gadm_path):
            gadm = gpd.read_file(gadm_path)
            gadm.plot(ax=ax, edgecolor='darkgray', facecolor='none', linewidth=0.25, transform=ccrs.PlateCarree())
        else:
            logging.warning("GADM data not available. Skipping state boundaries.")

        
    # Create line segments for all elements
    logging.info("Creating line collection...")
    segs = []
    for element in elements:
        points = nodes[element, :2]
        if not numpy.isnan(points).any() and not numpy.isinf(points).any():
            segs.append(numpy.concatenate([points, points[[0]]]))  # Close the polygon

    line_segments = LineCollection(segs, linewidths=0.05, colors='navy', alpha=0.3, transform=ccrs.PlateCarree())

    logging.info("Adding line collection to plot...")
    ax.add_collection(line_segments)

    # Plot depth using scatter
    valid_nodes = nodes[~numpy.isnan(nodes).any(axis=1) & ~numpy.isinf(nodes).any(axis=1)]
    scatter = ax.scatter(valid_nodes[:, 0], valid_nodes[:, 1], c=valid_nodes[:, 2], cmap='viridis',
                         s=0.5, alpha=0.7, transform=ccrs.PlateCarree())

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, orientation='vertical', pad=0.02, aspect=30)
    cbar.set_label('Depth (m)', rotation=270, labelpad=15, fontsize=8, fontweight='bold')
    cbar.ax.tick_params(labelsize=6)
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight('bold')

    if title is None:
        title = 'SCHISM Grid'
    
    ax.set_title(title, fontsize=12, pad=10, fontweight='bold')

    # Add latitude and longitude labels
    ax.set_xticks(numpy.arange(numpy.floor(extent[0]), numpy.ceil(extent[1])+1, 30), crs=ccrs.PlateCarree())
    ax.set_yticks(numpy.arange(numpy.floor(extent[2]), numpy.ceil(extent[3])+1, 30), crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.yaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.tick_params(labelsize=8, labelcolor='black', width=0.5)

    logging.info(f"Saving plot to {output_file}...")
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    end_time = time.time()  # End timing
    elapsed_time = end_time - start_time
    logging.info(f"Plot saved successfully. Time taken: {elapsed_time:.2f} seconds")
