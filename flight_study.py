"""
Flight Study — study timer backed by FlightRadar24 live data.
"""

from flask import Flask, jsonify, request, send_from_directory
import math, time, warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from FlightRadarAPI import FlightRadar24API

app = Flask(__name__)
fr  = FlightRadar24API()

# ── Caches ────────────────────────────────────────────────────
_detail_cache = {}   # fr24_id -> {data, ts}
_list_cache   = {}   # callsign_prefix -> {data, ts}
DETAIL_TTL = 10      # seconds
LIST_TTL   = 8

# ── Geo ───────────────────────────────────────────────────────
def haversine(lat1, lon1, lat2, lon2):
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat/2)**2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2)
    return R * 2 * math.asin(math.sqrt(a))

def calc_bearing(lat1, lon1, lat2, lon2):
    dlon = math.radians(lon2 - lon1)
    lat1, lat2 = math.radians(lat1), math.radians(lat2)
    x = math.sin(dlon) * math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    return (math.degrees(math.atan2(x, y)) + 360) % 360

# ── FR24 helpers ──────────────────────────────────────────────
def _airport_info(ap_data):
    if not ap_data:
        return None
    code = ap_data.get("code", {})
    pos  = ap_data.get("position", {})
    reg  = pos.get("region", {})
    tz   = ap_data.get("timezone", {})
    info = ap_data.get("info", {})
    return {
        "iata":     code.get("iata", ""),
        "icao":     code.get("icao", ""),
        "name":     ap_data.get("name", ""),
        "city":     reg.get("city", ""),
        "lat":      pos.get("latitude"),
        "lon":      pos.get("longitude"),
        "tz_name":  tz.get("name", ""),
        "tz_abbr":  tz.get("abbr", ""),
        "tz_offset":tz.get("offset", 0),
        "terminal": info.get("terminal"),
        "gate":     info.get("gate"),
    }

def get_flight_details(fr24_id, flight_obj):
    cached = _detail_cache.get(fr24_id)
    if cached and time.time() - cached["ts"] < DETAIL_TTL:
        return cached["data"]
    try:
        details = fr.get_flight_details(flight_obj)
        _detail_cache[fr24_id] = {"data": details, "ts": time.time()}
        return details
    except Exception:
        return (cached or {}).get("data")

def search_flights(query):
    query = query.strip()
    cached = _list_cache.get(query.upper())
    if cached and time.time() - cached["ts"] < LIST_TTL:
        return cached["data"]
    try:
        # Use FR24's own search — handles IATA (DL2512) and ICAO (DAL2512) callsigns
        results_raw = fr.search(query)
        live = results_raw.get("live", [])

        # Each live result has id + detail with lat/lon/route/reg/callsign
        results = []
        for item in live[:10]:
            detail = item.get("detail", {})
            results.append({
                "fr24_id":      item["id"],
                "callsign":     detail.get("callsign", "").strip(),
                "flight":       detail.get("flight", ""),
                "lat":          detail.get("lat"),
                "lon":          detail.get("lon"),
                "aircraft":     detail.get("ac_type", ""),
                "registration": detail.get("reg", ""),
                "route":        detail.get("route", ""),
                "schd_from":    detail.get("schd_from", ""),
                "schd_to":      detail.get("schd_to", ""),
                "altitude_ft":  0,
                "speed_kts":    0,
                "on_ground":    False,
            })
        _list_cache[query.upper()] = {"data": results, "ts": time.time()}
        return results
    except Exception:
        return []

# ── Routes ────────────────────────────────────────────────────
@app.route("/")
def index():
    return send_from_directory(".", "flight_study.html")

@app.route("/api/search")
def api_search():
    q = request.args.get("q", "").strip()
    if len(q) < 2:
        return jsonify({"error": "too short"}), 400
    flights = search_flights(q)
    return jsonify({"results": flights})

@app.route("/api/flight/<fr24_id>")
def api_flight(fr24_id):
    # Find the live flight object
    try:
        all_flights = fr.get_flights()
        flight_obj = next((f for f in all_flights if f.id == fr24_id), None)
    except Exception:
        flight_obj = None

    if not flight_obj:
        # Try cached details only
        cached = _detail_cache.get(fr24_id, {}).get("data")
        if not cached:
            return jsonify({"error": "Flight not found or landed"}), 404

    # Get full details
    details = get_flight_details(fr24_id, flight_obj) if flight_obj else _detail_cache.get(fr24_id, {}).get("data", {})
    if not details:
        return jsonify({"error": "No details available"}), 404

    # Current position — prefer live flight_obj, fall back to latest trail point
    trail = details.get("trail") or []
    if flight_obj:
        lat = flight_obj.latitude
        lon = flight_obj.longitude
        alt_ft = flight_obj.altitude
        spd_kts = flight_obj.ground_speed
        hdg = flight_obj.heading
        vspd_fpm = getattr(flight_obj, "vertical_speed", 0) or 0
        on_ground = flight_obj.on_ground
        reg = getattr(flight_obj, "registration", None) or ""
        aircraft = getattr(flight_obj, "aircraft_code", None) or ""
    elif trail:
        t = trail[0]
        lat, lon = t["lat"], t["lng"]
        alt_ft = t.get("alt", 0)
        spd_kts = t.get("spd", 0)
        hdg = t.get("hd", 0)
        vspd_fpm = 0
        on_ground = False
        reg = ""
        aircraft = ""
    else:
        return jsonify({"error": "No position data"}), 404

    # Heading: fallback to calculated from last two trail points
    if (not hdg or hdg == 0) and len(trail) >= 2:
        t1, t2 = trail[1], trail[0]   # trail[0] is newest
        hdg = calc_bearing(t1["lat"], t1["lng"], t2["lat"], t2["lng"])

    # Airport info
    ap      = details.get("airport", {})
    origin  = _airport_info(ap.get("origin"))
    dest    = _airport_info(ap.get("destination"))

    # Times
    ti = details.get("time", {})
    sched_dep = ti.get("scheduled", {}).get("departure")
    sched_arr = ti.get("scheduled", {}).get("arrival")
    actual_dep = ti.get("real", {}).get("departure")
    est_arr  = ti.get("estimated", {}).get("arrival") or ti.get("other", {}).get("eta")

    # ETA seconds from now
    eta_seconds = None
    if est_arr:
        eta_seconds = max(0, int(est_arr - time.time()))

    # Aircraft / airline
    aircraft_detail = details.get("aircraft", {})
    aircraft_model  = aircraft_detail.get("model", {}).get("text", aircraft)
    aircraft_reg    = aircraft_detail.get("registration", reg)
    images = aircraft_detail.get("images", {}).get("thumbnails", [])
    photo_url = images[0]["src"] if images else None

    airline_detail = details.get("airline", {})
    airline_name   = airline_detail.get("name", "")
    airline_iata   = airline_detail.get("code", {}).get("iata", "") if isinstance(airline_detail.get("code"), dict) else ""

    # Status
    status_text = details.get("status", {}).get("text", "")

    # Distance
    dist_remaining_km = None
    if dest and lat and lon:
        dist_remaining_km = haversine(lat, lon, dest["lat"], dest["lon"])

    # Distance traveled via trail
    dist_traveled_km = 0.0
    if len(trail) > 1:
        for i in range(len(trail) - 1):
            dist_traveled_km += haversine(trail[i]["lat"], trail[i]["lng"], trail[i+1]["lat"], trail[i+1]["lng"])

    total_km = dist_traveled_km + (dist_remaining_km or 0)
    progress_pct = round(dist_traveled_km / total_km * 100) if total_km > 0 else None

    # Trail for map (newest first → reverse for chronological order)
    trail_coords = [[t["lat"], t["lng"]] for t in reversed(trail)]

    return jsonify({
        "fr24_id":          fr24_id,
        "callsign":         (getattr(flight_obj, "callsign", "") or "").strip(),
        "lat":              lat,
        "lon":              lon,
        "altitude_ft":      alt_ft,
        "speed_kts":        round(spd_kts, 1) if spd_kts else 0,
        "speed_ms":         round(spd_kts * 0.514444, 2) if spd_kts else 0,
        "heading":          round(hdg, 1) if hdg else 0,
        "vertical_rate_fpm":vspd_fpm,
        "on_ground":        on_ground,
        "registration":     aircraft_reg,
        "aircraft_model":   aircraft_model,
        "airline_name":     airline_name,
        "airline_iata":     airline_iata,
        "photo_url":        photo_url,
        "status":           status_text,
        "origin":           origin,
        "destination":      dest,
        "sched_dep":        sched_dep,
        "sched_arr":        sched_arr,
        "actual_dep":       actual_dep,
        "est_arr":          est_arr,
        "eta_seconds":      eta_seconds,
        "server_time":      int(time.time()),
        "dist_traveled_km": round(dist_traveled_km, 2),
        "dist_remaining_km":round(dist_remaining_km, 2) if dist_remaining_km is not None else None,
        "total_km":         round(total_km, 2),
        "progress_pct":     progress_pct,
        "waypoints":        trail_coords,
    })

if __name__ == "__main__":
    print("Flight Study → http://localhost:5050")
    app.run(port=5050, debug=False)
