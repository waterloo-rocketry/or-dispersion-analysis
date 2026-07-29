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
