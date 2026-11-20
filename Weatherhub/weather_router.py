import pandas as pd
from open_meteo_api import fetch_open_meteo

def run_one_decision(lat, lon, decision, hourly_vars, daily_vars=None, extra_params=None):
    """
    decision["ranges"] = [("forecast"/"historical", start, end), ...]
    Returns dict {"hourly": df, "daily": df(optional)}
    """
    daily_vars = daily_vars or []
    extra_params = extra_params or {}

    hourly_parts = []
    daily_parts = []

    for kind, s, e in decision["ranges"]:    
        try:
            res = fetch_open_meteo(lat, lon, s, e, hourly_vars, daily_vars, extra_params, fetch_type=kind)
        except ValueError as e:
            print(e, "we're cooked")
            
        if "hourly" in res:
            hourly_parts.append(res["hourly"])
        if "daily" in res:
            daily_parts.append(res["daily"])

    out = {}
    if hourly_parts:
        dfh = pd.concat(hourly_parts, ignore_index=True)

        time_col = "time" if "time" in dfh.columns else ("date" if "date" in dfh.columns else None)
        if time_col is None:
            raise KeyError(f"Hourly df missing time column. Columns: {list(dfh.columns)}")

        if time_col != "time":
            dfh = dfh.rename(columns={time_col: "time"})
            time_col = "time"

        dfh = dfh.drop_duplicates(subset=[time_col]).sort_values(time_col).reset_index(drop=True)
        out["hourly"] = dfh

    if daily_parts:
        dfd = pd.concat(daily_parts, ignore_index=True)
        dfd = dfd.drop_duplicates(subset=["date"]).sort_values("date").reset_index(drop=True)
        out["daily"] = dfd

    return out


def run_many(lat, lon, saved_decisions, hourly_vars, daily_vars=None, extra_params=None):
    """
    saved_decisions: list of decision dicts (user may add multiple queries)
    Returns a single df_hourly + df_daily with query_id column
    """
    hourly_all = []
    daily_all = []

    for i, decision in enumerate(saved_decisions, start=1):
        res = run_one_decision(lat, lon, decision, hourly_vars, daily_vars, extra_params)
        if "hourly" in res:
            df = res["hourly"].copy()
            df["query_id"] = i
            hourly_all.append(df)
        if "daily" in res:
            df = res["daily"].copy()
            df["query_id"] = i
            daily_all.append(df)

    out = {}
    if hourly_all:
        out["hourly"] = pd.concat(hourly_all, ignore_index=True)
    if daily_all:
        out["daily"] = pd.concat(daily_all, ignore_index=True)
    return out