"""
Script to extract road coordinates from a KML file.
"""

import pandas as pd
import xml.etree.ElementTree as ET

# Load KML
tree = ET.parse('Road Coordinates/highway_144+101.kml')
root = tree.getroot()

# KML uses a namespace
ns = {"kml": "http://www.opengis.net/kml/2.2"}

coordinates_list = []

# Find all <coordinates> tags
for coordinates in root.findall(".//kml:coordinates", ns):
    coord_text = coordinates.text.strip()
    for pair in coord_text.split():
        lon, lat, *_ = map(float, pair.split(","))
        coordinates_list.append((lat, lon))

# Save to CSV
df = pd.DataFrame(coordinates_list, columns=["Latitude", "Longitude"])
df.to_csv("Road Coordinates/highway_144+101.csv", index=False)
print("Saved highway_144+101.csv")
