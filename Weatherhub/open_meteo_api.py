# historical_api.py
import pandas as pd
import openmeteo_requests
import requests_cache
from retry_requests import retry

# The links of the API calls 
from config import OPENMETEO_API_URLS

def hourly_to_df(response, hourly_vars: list[str], time_col="time") -> pd.DataFrame:
    hourly = response.Hourly()
    times = pd.date_range(
        start=pd.to_datetime(hourly.Time(), unit="s", utc=True),
        end=pd.to_datetime(hourly.TimeEnd(), unit="s", utc=True),
        freq=pd.Timedelta(seconds=hourly.Interval()),
        inclusive="left",
    )
    df = pd.DataFrame({time_col: times})

    n = hourly.VariablesLength()
    if len(hourly_vars) > n:
        raise ValueError(f"Requested {len(hourly_vars)} hourly vars, but API returned {n}.")

    for i, var in enumerate(hourly_vars):
        df[var] = hourly.Variables(i).ValuesAsNumpy()
    return df


def daily_to_df(response, daily_vars: list[str], time_col="date") -> pd.DataFrame:
    daily = response.Daily()
    dates = pd.date_range(
        start=pd.to_datetime(daily.Time(), unit="s", utc=True),
        end=pd.to_datetime(daily.TimeEnd(), unit="s", utc=True),
        freq=pd.Timedelta(seconds=daily.Interval()),
        inclusive="left",
    )
    df = pd.DataFrame({time_col: dates})

    n = daily.VariablesLength()
    if len(daily_vars) > n:
        raise ValueError(f"Requested {len(daily_vars)} daily vars, but API returned {n}.")
    
    for i, var in enumerate(daily_vars):
        df[var] = daily.Variables(i).ValuesAsNumpy()
    return df


def get_client(cache_path: str = ".cache", expire_after: int = 3600):
    cache = requests_cache.CachedSession(cache_path, expire_after=expire_after)
    retry_session = retry(cache, retries=5, backoff_factor=0.2)
    return openmeteo_requests.Client(session=retry_session)


def fetch_open_meteo(lat, lon, start_date, end_date, hourly_vars, daily_vars=None, extra_params=None, fetch_type="historical"):
    if fetch_type not in ["forecast", "historical"]:
        raise ValueError(f"Unknown kind: {fetch_type}")
    
    daily_vars = daily_vars or []
    extra_params = extra_params or {}

    # If the daily_parameters are empty, it should be left out
    params = {
        "latitude": lat,
        "longitude": lon,
        "start_date": start_date.isoformat(),
        "end_date": end_date.isoformat(),
        "hourly": hourly_vars,
        "daily": daily_vars,
        "timezone": "auto",
        **extra_params,
    }

    # cache forever like Open-Meteo example | written by Kushini
    openmeteo = get_client(expire_after=-1)  

    responses = openmeteo.weather_api(OPENMETEO_API_URLS[fetch_type], params=params)
    r = responses[0]

    output = {}
    if hourly_vars:
        output["hourly"] = hourly_to_df(r, hourly_vars, time_col="time")
        output["hourly"]["source"] = fetch_type
    if daily_vars:
        output["daily"] = daily_to_df(r, daily_vars, time_col="date")
        output["daily"]["source"] = fetch_type
    return output