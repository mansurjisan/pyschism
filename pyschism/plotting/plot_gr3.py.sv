import numpy as numpy
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import geopandas as gpd
import os
import urllib.request
import zipfile

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

def plot_gr3(nodes, elements, output_file):
    if nodes is None or elements is None:
        print("Cannot plot: invalid data")
        return

    print("Preparing plot...")
    fig = plt.figure(figsize=(20, 20))

    # Use PlateCarree projection
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # Set map extent based on node coordinates with a margin
    margin = 1  # degree
    extent = [
        nodes[:, 0].min() - margin,
        nodes[:, 0].max() + margin,
        nodes[:, 1].min() - margin,
        nodes[:, 1].max() + margin
    ]
    ax.set_extent(extent)

    # Download and plot coastline
    coastline_path = download_coastline_data()
    coastline = gpd.read_file(coastline_path)
    coastline.plot(ax=ax, edgecolor='black', facecolor='none', linewidth=1.5, transform=ccrs.PlateCarree())

    # Download and plot GADM state boundaries
    gadm_path = download_gadm_data()
    gadm = gpd.read_file(gadm_path)
    gadm.plot(ax=ax, edgecolor='darkgray', facecolor='none', linewidth=0.5, transform=ccrs.PlateCarree())

    # Add other features
    ax.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='face')
    ax.add_feature(cfeature.RIVERS, edgecolor='blue', linewidth=0.5)

    # Add grid lines
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Create line segments for all elements
    print("Creating line collection...")
    segs = []
    for element in elements:
        points = nodes[element, :2]
        segs.append(numpy.concatenate([points, points[[0]]]))  # Close the polygon

    line_segments = LineCollection(segs, linewidths=0.1, colors='navy', alpha=0.3, transform=ccrs.PlateCarree())

    print("Adding line collection to plot...")
    ax.add_collection(line_segments)

    # Plot depth using scatter
    scatter = ax.scatter(nodes[:, 0], nodes[:, 1], c=nodes[:, 2], cmap='viridis',
                         s=1, alpha=0.7, transform=ccrs.PlateCarree())

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label('Depth (m)', rotation=270, labelpad=15)

    # Set labels and title
    ax.set_title('SCHISM Grid Plot with Natural Earth Coastline and GADM State Boundaries', fontsize=16, pad=20)

    print(f"Saving plot to {output_file}...")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("Plot saved successfully.")
