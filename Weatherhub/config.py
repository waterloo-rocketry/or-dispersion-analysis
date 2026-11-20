OPENMETEO_API_URLS = {
    "forecast":"https://api.open-meteo.com/v1/forecast",
    "historical":"https://historical-forecast-api.open-meteo.com/v1/forecast"
}

# All altitudes supported by OpenMeteo in metres
ALTITUDES = [110, 320, 500, 800, 1000, 1500, 1900, 3200, 4200, 5600, 7200, 
             9200, 10400, 11800, 13500, 15800, 17700, 19300, 22000]


# Names of the wind speed columns for OpenMeteo
WIND_SPEEDS = ["wind_speed_1000hPa", "wind_speed_975hPa", "wind_speed_950hPa", 
                "wind_speed_925hPa", "wind_speed_900hPa", "wind_speed_850hPa", "wind_speed_800hPa", 
                "wind_speed_700hPa", "wind_speed_600hPa", "wind_speed_500hPa", "wind_speed_400hPa", 
                "wind_speed_300hPa", "wind_speed_250hPa", "wind_speed_200hPa", "wind_speed_150hPa", 
                "wind_speed_100hPa", "wind_speed_70hPa", "wind_speed_50hPa", "wind_speed_30hPa"]

# Names of the wind direction columns for OpenMeteo
WIND_DIRECTIONS = ["wind_direction_1000hPa", "wind_direction_975hPa", "wind_direction_950hPa", 
                "wind_direction_925hPa", "wind_direction_900hPa", "wind_direction_850hPa", 
                "wind_direction_800hPa", "wind_direction_700hPa", "wind_direction_600hPa", 
                "wind_direction_500hPa", "wind_direction_400hPa", "wind_direction_300hPa", 
                "wind_direction_250hPa", "wind_direction_200hPa", "wind_direction_150hPa", 
                "wind_direction_100hPa", "wind_direction_70hPa", "wind_direction_50hPa", 
                "wind_direction_30hPa"]

# Openmeteo all columns we receive from their API including the temperature, 
HOURLY_VARS = ["temperature_2m"] + WIND_SPEEDS + WIND_DIRECTIONS + ["pressure_msl", "surface_pressure", "cloud_cover_high", "cloud_cover_mid", "cloud_cover_low", "is_day"]