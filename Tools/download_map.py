"""
Written by: Geminini
Collaborators: Luca Scavone

Script to download the map tile from the ESRI World Imagery. GitHub doesn't handle large files very well, so
directly pulling this repo from Git may not download a local copy of the satellite image, which is where this
script comes in handy.

Usage:
    1. (if applicable) Update edge coordinates
    2. Run that john
"""

import contextily as cx

# Define exact boundaries
south_lat = 47.00   # Bottom edge
north_lat = 49.00   # Top edge
west_lon  = -83.00  # Left edge
east_lon  = -81.00  # Right edge

# Download and save the raster file
output_filename = "../LC Geography/lc_basemap.tif"

print(f"Downloading map tiles to {output_filename}...")

try:
    # bounds2raster downloads the area and stitches it into a local TIFF
    cx.bounds2raster(
        w=west_lon,
        s=south_lat,
        e=east_lon,
        n=north_lat,
        path=output_filename,
        ll=True,
        source=cx.providers.Esri.WorldImagery,
        zoom=11  # Controls resolution (10 is lower res, 12 is high res)
    )
    print("Download complete! You can now use this file offline.")
except Exception as e:
    print(f"Failed to download map: {e}")
