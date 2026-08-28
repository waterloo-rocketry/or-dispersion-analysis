"""
Written by: Geminini
Collaborators: Luca Scavone

Downloads a local GeoJSON file containing information about local bodies including lakes and smaller features. Does
not include rivers or streams. Default area is a 20km-radius centered on the Launch Canada Advanced Pad. This file
allows Project Atlas to associate detected water landings with local geography to help recovery efforts.

Usage:
    1. (if applicable) Update launch site coordinates
    2. (if applicable) Update radius
    3. Run that john
"""
import osmnx as ox

# Launch site coordinates
lat, lon = 47.965378, -81.873536
radius_meters = 20000

print("Fetching lake data...")

# Query OpenStreetMap for water bodies
tags = {'natural': 'water', 'water': 'lake'}
lakes = ox.features_from_point((lat, lon), tags=tags, dist=radius_meters)

# Ensure the output directory matches your project structure
output_path = "../LC Geography/lakes.geojson"
lakes.to_file(output_path, driver="GeoJSON")
print(f"Saved {len(lakes)} water bodies to {output_path}")
