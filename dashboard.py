#!/usr/bin/env python3
"""SeisComP Live Dashboard — dashboard.py"""
import io, math, threading, time, collections, os, warnings
import numpy as np
import mysql.connector
from datetime import datetime, timezone
from flask import Flask, render_template_string, jsonify, send_file, abort, request
from flask_socketio import SocketIO
from flask_cors import CORS
from obspy.clients.seedlink.easyseedlink import EasySeedLinkClient
from obspy.clients.filesystem.sds import Client as SDSClient
from obspy.clients.fdsn import Client as FDSNClient
from obspy import UTCDateTime

# FDSN fallback clients (used when local archive lacks data)
_fdsn_clients = {}
def _fdsn(network):
    """Return a cached FDSN client appropriate for a given network."""
    src = "SCEDC" if network in ("CI","AZ","SB","SC") else "IRIS"
    if src not in _fdsn_clients:
        try: _fdsn_clients[src] = FDSNClient(src)
        except: _fdsn_clients[src] = None
    return _fdsn_clients[src]

# ── Config ─────────────────────────────────────────────────────────────────────
SEEDLINK      = "localhost:18000"
IRIS_SEEDLINK = "rtserve.iris.washington.edu:18000"
# ── Preliminary real-time detection alerts ────────────────────────────────────
# Sent when 3+ stations trigger within 30s (BEFORE SeisComP locates the event)
PRELIM_DISCORD_WEBHOOK = "https://discord.com/api/webhooks/1499565837688901814/M1qxkbhHLhnDxpj9b0WwNACelN6itRyGvx4byYHOpwhK-RRRtQ-mci9mXwtuLtclEYCJ"
PRELIM_DISCORD_ENABLED = True

# ── Confirmed SeisComP event Discord announcements ─────────────────────────────
EVENT_DISCORD_WEBHOOK  = PRELIM_DISCORD_WEBHOOK   # reuse same channel
EVENT_DISCORD_ENABLED  = True

# All webhooks — notifications are broadcast to every URL in this list
DISCORD_WEBHOOKS = [
    PRELIM_DISCORD_WEBHOOK,
    "https://discord.com/api/webhooks/1508676731656077392/A_xCrUH11bb428dQu3ii6kNDS-JgaVd9msakiB7pI8pGHOarHA7PVcxu8qIP4WAmf6YN",
]
EVENT_DISCORD_MIN_MAG  = 0.0   # announce ALL events (set higher to filter)
DISCORD_EVERYONE_MAG   = 5.0   # @everyone mention threshold (both prelim + confirmed)

# @everyone is sent at most ONCE per event — tracked by event ID (confirmed)
# and by a per-session cooldown (preliminary, since they lack stable IDs).
_everyone_pinged_events:  set  = set()   # evids already pinged
_everyone_prelim_last_t: float = 0.0     # epoch of last prelim @everyone
_EVERYONE_PRELIM_COOLDOWN = 1800         # 30 min between prelim @everyone pings
PRELIM_EMAIL_ENABLED   = False   # set to True and fill below
PRELIM_EMAIL_FROM      = "your.email@gmail.com"
PRELIM_EMAIL_TO        = ["your.email@gmail.com"]
PRELIM_SMTP_HOST       = "smtp.gmail.com"
PRELIM_SMTP_PORT       = 587
PRELIM_SMTP_USER       = "your.email@gmail.com"
PRELIM_SMTP_PASS       = "xxxx xxxx xxxx xxxx"   # Gmail App Password
PRELIM_MIN_STATIONS    = 3     # minimum stations (need ≥3 to constrain a location)
PRELIM_MIN_STALTA      = 2.0  # STA/LTA threshold for PRELIM candidate triggers (matches TRIG)
PRELIM_MIN_MAG         = 0.0  # alert for all local magnitudes

# Geographic filter — local (non-teleseism) detections must be inside this box
# Covers California + small buffer into NV/AZ/OR/MX border zones
PRELIM_LOCAL_LAT_MIN   = 32.0   # southern tip of Baja buffer
PRELIM_LOCAL_LAT_MAX   = 42.5   # Oregon border
PRELIM_LOCAL_LON_MIN   = -125.0 # Pacific coast
PRELIM_LOCAL_LON_MAX   = -113.5 # AZ/NV border buffer
PRELIM_MIN_SPREAD_KM   = 0    # no minimum spread
PRELIM_MAX_SPREAD_KM   = 1200
TELESEISM_SPREAD_KM    = 500  # spread >= this → treat as teleseism
TELESEISM_MIN_MAG      = 5.5
TELESEISM_MIN_STATIONS = 5    # teleseisms need ≥5 stations
TELESEISM_VP           = 8.5  # km/s — mantle P-wave velocity (vs 6.0 crustal)
TELESEISM_MAX_RMS      = 3.0  # tighter travel-time consistency for teleseisms
PRELIM_WINDOW_SEC      = 60   # 60s window — P-waves take ~17s per 100km
PRELIM_COOLDOWN_SEC    = 120  # 2 min cooldown per region
DB_CFG     = dict(host="localhost", user="sysop", password="sysop", database="seiscomp")
SDS_PATH   = "/Users/OuOu/seiscomp/var/lib/archive"
PORT       = 5500
STA_WIN, LTA_WIN, TRIG, ALARM = 1.0, 60.0, 2.0, 4.0
BUF_SECS   = 720   # 12 min full-res buffer (≥ LTA_WIN 60s + 10min zoom)
BUF_LONG   = 7200  # 2 hr downsampled buffer (~1 sps)

STREAMS = [
    *[("CI", s, "BHZ") for s in [
        "ADO","ARV","BAK","BAR","BBR","BC3","BEL","BFS","CHF","CIA","CWC","DAN",
        "DEC","DGR","DJJ","EDW2","FMP","FUR","GLA","GMR","GRA","GSC","IKP","IRM",
        "ISA","LJR","LRL","MLAC","MPM","MPP","MUR","MWC","NEE2","OSI","PASC","PDM",
        "PLM","RPV","RRX","SBC","SCI2","SCZ2","SHO","SLA","SMM","SNCC","SVD","SWS",
        "TIN","TUQ","USC","VES","VOG","VTV",
        # Additional CI stations
        "AGO","ALP","ALS","AVM","BBS","BCC","BCW","BGM","BLA2","BTC","BYR","CAP",
        "CCA","CCC","CCL","CCP","CCRB","CDY","CFW","CGO","CHN","CHS","CLC","CLT",
        "CMPB","CMS","COC","CPE","CRE","CRR","CSP","CTW","DLA","DPP","DSS",
        "DTP","DUG","EML","ERR","EUR","FIG","FLO","FRE","FUL","GAS","GBP","GEO",
        "GFL","GOL","GRB","GRM","GRP","GTM","HEC","HLL","HNT","HOL","HPK",
        "IDO","JRC2","JVA","KNW","LAS","LCP","LEA","LGU","LLP","LNO","LOC",
        "LRR","MAG","MBAR","MBR","MDS","MGE","MLS","MOD","MSJ","MST","MTM","MTP",
        "NAP","NEW","NJQ","NOR","OGC","OLI","OPE","ORR","OSB","OXN","PAS","PBI",
        "PDL","PER","PIN","PKD","PLB","PLC","PLO","PMD","PNE","POB","POM","PSD",
        "QVH","RAD","RCT","RIN","RIO","RPL","RRS","RTV","RVR","SAI","SAL","SAO",
        "SDD","SDG","SDN","SEG","SEP","SER","SFR","SFT","SGH","SGO","SHB","SIL",
        "SIO","SKY","SLP","SLR","SMF","SMG","SNR","SOD","SOR","SPC","SRN",
        "SRP","SSP","STG","STN","SWF","SYL","SYP","SYS","TAB","TAF","TAN",
        "TBR","TCP","TDF","TEM","TFT","TGN","TIG","TJN","TMB","TMR","TOV","TPO",
        "TRE","TRM","TRN","TRO","TRR","TTB","TUL","TUS","TUT","TWL","UNV",
        "VCT","VDB","VDP","VER","VLN","VNA","VNO","VPK","VSR","WBS2",
        "WBY","WCC","WCS","WHY","WIK","WMF","WMO","WNM","WOF","WOR","WRC","WRN",
        "WRV","WVLA","WWT","WYA","WYO","YBH","YCC","YMO",
    ]],
    *[("AZ", s, "BHZ") for s in ["BZN","CRY","FRD","KNW","MONP","PFO","RDM","SCI2","SND","TRO"]],
    *[("BK", s, "BHZ") for s in [
        # Core BDSN backbone
        "BKS","BRK","BRIB","FARB","MCCM",
        # Bay Area dense network
        "MHDL","BUCI","LCOS","LLNL","PINL","RVRP","BDM","OAKV","MORK","ATP","MLKN","JEPS","PWOD",
        "CVS","VAK","BL67","MOBB","JASP","JRSC","STAN","UMUN","OXMT","PESC","MTOS","TESL","WENL",
        # NorCal backbone
        "CMB","HOPS","MNRC","ORV","PKD","SAO","MHC","MBARI","ARC","BARR","BONV",
    ]],
    *[("NC", s, "HHZ") for s in [
        "KCPB","KHBB","KHMB","KMR","KRMB","KSXB","LDH","PMPB",
        "CCOB","CMSB","JROB","MCDB","MHRB","MNOB","MPSB","PLIB","PLNB",
    ]],
    *[("SB", s, "HHZ") for s in ["CPSLO","VAFB1","VAFB2","WLA","WLA01","WLA02","WLA03","WLA04","WLA06","WLA10"]],
    # NN — Nevada Seismological Laboratory (HH channels, served by IRIS)
    *[("NN", s, "HHZ") for s in [
        # Core Nevada stations (18 original)
        "BFC","BMHS","BRH5","BRS2","CMK6","COLR","DEDC","DIX","DNYB",
        "GMN","ION4","KVN","LCH","LHV","MCA06","OMM","PAH","WCN",
        # Sierra Nevada / Lake Tahoe cluster (~38-40N, near NorCal border)
        "BEK","CTC","EMBB","MPK","PNT","PYM2","RUB","WAK","WASH","WVA","YER","ZPR",
    ]],
    # UO — University of Oregon (HH channels, served by IRIS)
    *[("UO", s, "HHZ") for s in [
        # Southern Oregon / CA-OR border (~42N)
        "ADEL","CAVE","DUTCH","FIDL","KBO","RANT","ROGE","SISQ",
        # Original 5 (BH)
        "BUCK","DBO","EUO","PIN","PINE",
    ]],
    # UW — Pacific Northwest including southern Oregon stations (IRIS)
    *[("UW", s, "BHZ") for s in [
        "BRAN","CCRK","DAVN","DDRF","DOSE","FISH","FORK","GNW","GRCC","LEBA",
        "LON","LRIV","LTY","MEGW","MRBL","OFR","OMAK","OPC","PASS","PHIN",
        "RADR","RAI","RATT","RWW","SEP","SHUK","SP2","SPUD","SQM","SSW",
        "STOR","TOLT","TTW","TUCA","UGP4","WISH","WOLL","YACT",
    ]],
    # UW southern Oregon stations (HH channels)
    *[("UW", s, "HHZ") for s in ["BBO","IRON","TREE"]],
]

# ── FDSN network-level wildcard polling ────────────────────────────────────────
# Format: "NETWORK": ("endpoint_key", "channel_priority_list")
# endpoint_key: "scedc" | "ncedc" | "iris"
# We request sta=* for each network → get every active station automatically
FDSN_NETWORKS = {
    # ── California only ───────────────────────────────────────────────────────
    # Southern California Seismic Network (SCSN)
    "CI": ("scedc",  "BHZ,HHZ,EHZ"),
    # Anza Network (Southern CA arrays)
    "AZ": ("scedc",  "BHZ,HHZ"),
    # Santa Barbara + UCSB Network
    "SB": ("scedc",  "HHZ,BHZ"),
    # Berkeley Seismological Lab (Bay Area / Northern CA)
    "BK": ("ncedc",  "BHZ,HHZ"),
    # Northern California Seismic Network (USGS Menlo Park)
    "NC": ("ncedc",  "HHZ,BHZ,EHZ"),
}

# Keep IRIS_STREAMS as empty list for backward compat (states dict still uses it)
IRIS_STREAMS = []

# ── Hardcoded coords for FDSN-only stations (BK Bay Area + NN Nevada) ───────────
_IRIS_STA_COORDS = {
    # BK Bay Area
    "BK.BRIB": {"net":"BK","sta":"BRIB","lat":37.9189,"lon":-122.1518,"elev":219.7},
    "BK.BRK":  {"net":"BK","sta":"BRK", "lat":37.8735,"lon":-122.2610,"elev":49.4},
    "BK.CVS":  {"net":"BK","sta":"CVS", "lat":38.3453,"lon":-122.4584,"elev":295.1},
    "BK.JEPS": {"net":"BK","sta":"JEPS","lat":38.2579,"lon":-121.8252,"elev":6.1},
    "BK.JRSC": {"net":"BK","sta":"JRSC","lat":37.4037,"lon":-122.2387,"elev":70.5},
    "BK.MHC":  {"net":"BK","sta":"MHC", "lat":37.3416,"lon":-121.6426,"elev":1250.4},
    "BK.MHDL": {"net":"BK","sta":"MHDL","lat":37.8423,"lon":-122.4943,"elev":94.5},
    "BK.OHLN": {"net":"BK","sta":"OHLN","lat":38.0062,"lon":-122.2730,"elev":-0.5},
    "BK.OXMT": {"net":"BK","sta":"OXMT","lat":37.4994,"lon":-122.4243,"elev":209.1},
    "BK.POTR": {"net":"BK","sta":"POTR","lat":38.2026,"lon":-121.9353,"elev":20.0},
    "BK.RFSB": {"net":"BK","sta":"RFSB","lat":37.9161,"lon":-122.3361,"elev":-27.3},
    "BK.SBRN": {"net":"BK","sta":"SBRN","lat":37.6856,"lon":-122.4113,"elev":4.0},
    "BK.SCCB": {"net":"BK","sta":"SCCB","lat":37.2874,"lon":-121.8642,"elev":98.4},
    "BK.STAN": {"net":"BK","sta":"STAN","lat":37.4039,"lon":-122.1751,"elev":125.5},
    "BK.SVIN": {"net":"BK","sta":"SVIN","lat":38.0332,"lon":-122.5263,"elev":-27.5},
    "BK.WENL": {"net":"BK","sta":"WENL","lat":37.6221,"lon":-121.7570,"elev":138.9},
    # NN Nevada (coords from FDSN EARTHSCOPE query)
    "NN.BEK":  {"net":"NN","sta":"BEK",  "lat":39.867,"lon":-120.360,"elev":1820.0},
    "NN.BFC":  {"net":"NN","sta":"BFC",  "lat":38.890,"lon":-119.610,"elev":1890.0},
    "NN.BMHS": {"net":"NN","sta":"BMHS", "lat":39.420,"lon":-119.760,"elev":1580.0},
    "NN.BRH5": {"net":"NN","sta":"BRH5", "lat":39.050,"lon":-118.040,"elev":1530.0},
    "NN.BRS2": {"net":"NN","sta":"BRS2", "lat":38.530,"lon":-117.630,"elev":1650.0},
    "NN.CMK6": {"net":"NN","sta":"CMK6", "lat":39.310,"lon":-118.120,"elev":1380.0},
    "NN.COLR": {"net":"NN","sta":"COLR", "lat":39.540,"lon":-119.620,"elev":1460.0},
    "NN.CTC":  {"net":"NN","sta":"CTC",  "lat":39.208,"lon":-120.126,"elev":1980.0},
    "NN.DEDC": {"net":"NN","sta":"DEDC", "lat":39.410,"lon":-119.000,"elev":1510.0},
    "NN.DIX":  {"net":"NN","sta":"DIX",  "lat":39.800,"lon":-118.080,"elev":1325.0},
    "NN.DNYB": {"net":"NN","sta":"DNYB", "lat":41.090,"lon":-119.280,"elev":1340.0},
    "NN.EMBB": {"net":"NN","sta":"EMBB", "lat":38.974,"lon":-120.100,"elev":2010.0},
    "NN.GMN":  {"net":"NN","sta":"GMN",  "lat":38.500,"lon":-119.050,"elev":2300.0},
    "NN.ION4": {"net":"NN","sta":"ION4", "lat":38.890,"lon":-117.740,"elev":1590.0},
    "NN.KVN":  {"net":"NN","sta":"KVN",  "lat":39.050,"lon":-118.100,"elev":1740.0},
    "NN.LCH":  {"net":"NN","sta":"LCH",  "lat":38.200,"lon":-118.500,"elev":2020.0},
    "NN.LHV":  {"net":"NN","sta":"LHV",  "lat":38.250,"lon":-118.500,"elev":1700.0},
    "NN.MCA06":{"net":"NN","sta":"MCA06","lat":39.500,"lon":-117.800,"elev":1680.0},
    "NN.MPK":  {"net":"NN","sta":"MPK",  "lat":39.293,"lon":-120.036,"elev":1920.0},
    "NN.OMM":  {"net":"NN","sta":"OMM",  "lat":38.700,"lon":-119.200,"elev":2200.0},
    "NN.PAH":  {"net":"NN","sta":"PAH",  "lat":39.711,"lon":-119.385,"elev":1560.0},
    "NN.PNT":  {"net":"NN","sta":"PNT",  "lat":39.089,"lon":-119.600,"elev":1580.0},
    "NN.PYM2": {"net":"NN","sta":"PYM2", "lat":40.170,"lon":-119.732,"elev":1420.0},
    "NN.RUB":  {"net":"NN","sta":"RUB",  "lat":39.052,"lon":-120.155,"elev":1960.0},
    "NN.WAK":  {"net":"NN","sta":"WAK",  "lat":38.504,"lon":-119.438,"elev":2200.0},
    "NN.WASH": {"net":"NN","sta":"WASH", "lat":39.334,"lon":-119.779,"elev":1510.0},
    "NN.WCN":  {"net":"NN","sta":"WCN",  "lat":39.019,"lon":-119.917,"elev":1870.0},
    "NN.WVA":  {"net":"NN","sta":"WVA",  "lat":39.944,"lon":-119.825,"elev":1390.0},
    "NN.YER":  {"net":"NN","sta":"YER",  "lat":38.985,"lon":-119.241,"elev":2260.0},
    "NN.ZPR":  {"net":"NN","sta":"ZPR",  "lat":39.013,"lon":-119.938,"elev":1870.0},
    # UO — University of Oregon southern Oregon stations (from FDSN EARTHSCOPE)
    "UO.ADEL": {"net":"UO","sta":"ADEL", "lat":42.169,"lon":-119.901,"elev":1388.0},
    "UO.CAVE": {"net":"UO","sta":"CAVE", "lat":42.121,"lon":-123.571,"elev":1100.0},
    "UO.DUTCH":{"net":"UO","sta":"DUTCH","lat":42.044,"lon":-122.892,"elev":1050.0},
    "UO.FIDL": {"net":"UO","sta":"FIDL", "lat":42.240,"lon":-123.768,"elev":540.0},
    "UO.KBO":  {"net":"UO","sta":"KBO",  "lat":42.212,"lon":-124.226,"elev":145.0},
    "UO.RANT": {"net":"UO","sta":"RANT", "lat":42.035,"lon":-121.276,"elev":1330.0},
    "UO.ROGE": {"net":"UO","sta":"ROGE", "lat":42.695,"lon":-123.665,"elev":860.0},
    "UO.SISQ": {"net":"UO","sta":"SISQ", "lat":42.275,"lon":-123.632,"elev":910.0},
    # UW — Pacific Northwest southern Oregon stations
    "UW.BBO":  {"net":"UW","sta":"BBO",  "lat":42.888,"lon":-122.680,"elev":1850.0},
    "UW.IRON": {"net":"UW","sta":"IRON", "lat":43.358,"lon":-118.474,"elev":1380.0},
    "UW.TREE": {"net":"UW","sta":"TREE", "lat":42.726,"lon":-120.893,"elev":1470.0},
    # UU Utah
    "UU.BGU":  {"net":"UU","sta":"BGU", "lat":37.234,"lon":-112.985,"elev":1820.0},
    "UU.CTU":  {"net":"UU","sta":"CTU", "lat":41.893,"lon":-111.455,"elev":1380.0},
    "UU.EPAZ": {"net":"UU","sta":"EPAZ","lat":31.890,"lon":-109.064,"elev":1250.0},
    "UU.FSU":  {"net":"UU","sta":"FSU", "lat":38.524,"lon":-112.800,"elev":1660.0},
    "UU.FTU":  {"net":"UU","sta":"FTU", "lat":39.978,"lon":-111.430,"elev":2170.0},
    "UU.HVU":  {"net":"UU","sta":"HVU", "lat":40.498,"lon":-111.951,"elev":1520.0},
    "UU.LCMT": {"net":"UU","sta":"LCMT","lat":38.576,"lon":-109.604,"elev":1480.0},
    "UU.MPU":  {"net":"UU","sta":"MPU", "lat":41.308,"lon":-112.277,"elev":1510.0},
    "UU.MTUT": {"net":"UU","sta":"MTUT","lat":40.389,"lon":-111.578,"elev":2700.0},
    "UU.RCU":  {"net":"UU","sta":"RCU", "lat":38.475,"lon":-109.292,"elev":1610.0},
    "UU.SRU":  {"net":"UU","sta":"SRU", "lat":39.178,"lon":-111.922,"elev":1680.0},
    "UU.SPU":  {"net":"UU","sta":"SPU", "lat":37.679,"lon":-112.460,"elev":2470.0},
    "UU.TPVU": {"net":"UU","sta":"TPVU","lat":40.449,"lon":-111.714,"elev":2780.0},
    "UU.ZNPU": {"net":"UU","sta":"ZNPU","lat":37.291,"lon":-112.812,"elev":2004.0},
    "UU.BMUT": {"net":"UU","sta":"BMUT","lat":38.672,"lon":-110.153,"elev":1610.0},
    "UU.DLUT": {"net":"UU","sta":"DLUT","lat":39.327,"lon":-111.399,"elev":1750.0},
    "UU.SLUT": {"net":"UU","sta":"SLUT","lat":40.570,"lon":-112.088,"elev":1299.0},
    # WY Yellowstone
    "WY.YFT":  {"net":"WY","sta":"YFT", "lat":44.640,"lon":-110.929,"elev":2185.0},
    "WY.YMR":  {"net":"WY","sta":"YMR", "lat":44.561,"lon":-110.827,"elev":2300.0},
    "WY.YNR":  {"net":"WY","sta":"YNR", "lat":44.726,"lon":-110.669,"elev":2315.0},
    "WY.YHL":  {"net":"WY","sta":"YHL", "lat":44.602,"lon":-110.432,"elev":2230.0},
    "WY.YGC":  {"net":"WY","sta":"YGC", "lat":44.721,"lon":-110.703,"elev":2246.0},
    "WY.YDC":  {"net":"WY","sta":"YDC", "lat":44.726,"lon":-110.880,"elev":2215.0},
    "WY.YLA":  {"net":"WY","sta":"YLA", "lat":44.420,"lon":-110.562,"elev":2225.0},
    "WY.YMT":  {"net":"WY","sta":"YMT", "lat":44.679,"lon":-110.662,"elev":2337.0},
    "WY.YTP":  {"net":"WY","sta":"YTP", "lat":44.535,"lon":-110.745,"elev":2234.0},
    # AK Alaska
    "AK.BARN": {"net":"AK","sta":"BARN","lat":64.683,"lon":-163.046,"elev":14.0},
    "AK.BMR":  {"net":"AK","sta":"BMR", "lat":57.624,"lon":-153.498,"elev":125.0},
    "AK.BPAW": {"net":"AK","sta":"BPAW","lat":66.040,"lon":-151.512,"elev":295.0},
    "AK.CCB":  {"net":"AK","sta":"CCB", "lat":62.450,"lon":-145.521,"elev":640.0},
    "AK.DHY":  {"net":"AK","sta":"DHY", "lat":63.492,"lon":-148.827,"elev":860.0},
    "AK.DOT":  {"net":"AK","sta":"DOT", "lat":63.293,"lon":-152.358,"elev":1005.0},
    "AK.GHO":  {"net":"AK","sta":"GHO", "lat":61.561,"lon":-149.725,"elev":50.0},
    "AK.HDA":  {"net":"AK","sta":"HDA", "lat":59.993,"lon":-151.486,"elev":90.0},
    "AK.KLU":  {"net":"AK","sta":"KLU", "lat":61.730,"lon":-144.428,"elev":530.0},
    "AK.MCK":  {"net":"AK","sta":"MCK", "lat":63.731,"lon":-148.978,"elev":375.0},
    "AK.NICH": {"net":"AK","sta":"NICH","lat":60.291,"lon":-149.754,"elev":46.0},
    "AK.NKA":  {"net":"AK","sta":"NKA", "lat":61.955,"lon":-150.659,"elev":52.0},
    "AK.RC01": {"net":"AK","sta":"RC01","lat":64.734,"lon":-166.936,"elev":55.0},
    "AK.SCM":  {"net":"AK","sta":"SCM", "lat":61.681,"lon":-149.441,"elev":975.0},
    "AK.SLK":  {"net":"AK","sta":"SLK", "lat":61.060,"lon":-149.747,"elev":27.0},
    "AK.SSP":  {"net":"AK","sta":"SSP", "lat":59.775,"lon":-154.527,"elev":90.0},
    "AK.UNV":  {"net":"AK","sta":"UNV", "lat":63.671,"lon":-150.967,"elev":660.0},
    "AK.WAT1": {"net":"AK","sta":"WAT1","lat":65.277,"lon":-149.526,"elev":180.0},
    "AK.WRH":  {"net":"AK","sta":"WRH", "lat":62.656,"lon":-155.059,"elev":300.0},
    # HV Hawaii
    "HV.HLPD": {"net":"HV","sta":"HLPD","lat":19.421,"lon":-155.127,"elev":1044.0},
    "HV.KIPU": {"net":"HV","sta":"KIPU","lat":21.984,"lon":-159.399,"elev":122.0},
    "HV.MLOD": {"net":"HV","sta":"MLOD","lat":19.477,"lon":-155.376,"elev":3100.0},
    "HV.PHOD": {"net":"HV","sta":"PHOD","lat":21.433,"lon":-157.983,"elev":38.0},
    "HV.POLD": {"net":"HV","sta":"POLD","lat":20.899,"lon":-156.695,"elev":322.0},
    "HV.WRM":  {"net":"HV","sta":"WRM", "lat":19.965,"lon":-155.901,"elev":1239.0},
    "HV.HAT":  {"net":"HV","sta":"HAT", "lat":19.691,"lon":-155.952,"elev":1042.0},
    "HV.NPT":  {"net":"HV","sta":"NPT", "lat":20.024,"lon":-155.484,"elev":330.0},
    "HV.WILD": {"net":"HV","sta":"WILD","lat":22.029,"lon":-159.786,"elev":168.0},
    # NM New Mexico Tech
    "NM.BLO":  {"net":"NM","sta":"BLO", "lat":39.173,"lon":-86.522, "elev":258.0},
    "NM.CBKS": {"net":"NM","sta":"CBKS","lat":38.814,"lon":-99.737, "elev":550.0},
    "NM.GOGA": {"net":"NM","sta":"GOGA","lat":33.408,"lon":-83.461, "elev":248.0},
    "NM.MIAR": {"net":"NM","sta":"MIAR","lat":34.546,"lon":-93.573, "elev":185.0},
    "NM.USIN": {"net":"NM","sta":"USIN","lat":38.307,"lon":-85.808, "elev":185.0},
    "NM.WMOK": {"net":"NM","sta":"WMOK","lat":34.737,"lon":-98.781, "elev":434.0},
    "NM.WUAZ": {"net":"NM","sta":"WUAZ","lat":32.052,"lon":-110.700,"elev":860.0},
    # US National network
    "US.BINY": {"net":"US","sta":"BINY","lat":42.200,"lon":-75.990, "elev":400.0},
    "US.BLA":  {"net":"US","sta":"BLA", "lat":37.206,"lon":-80.421, "elev":669.0},
    "US.CBKS": {"net":"US","sta":"CBKS","lat":38.814,"lon":-99.737, "elev":550.0},
    "US.COWI": {"net":"US","sta":"COWI","lat":44.534,"lon":-110.636,"elev":2054.0},
    "US.ECSD": {"net":"US","sta":"ECSD","lat":43.695,"lon":-103.536,"elev":1189.0},
    "US.GOGA": {"net":"US","sta":"GOGA","lat":33.408,"lon":-83.461, "elev":248.0},
    "US.HLID": {"net":"US","sta":"HLID","lat":43.563,"lon":-114.404,"elev":1482.0},
    "US.ISCO": {"net":"US","sta":"ISCO","lat":39.783,"lon":-105.458,"elev":1862.0},
    "US.LKWY": {"net":"US","sta":"LKWY","lat":44.565,"lon":-110.400,"elev":2100.0},
    "US.MIAR": {"net":"US","sta":"MIAR","lat":34.546,"lon":-93.573, "elev":185.0},
    "US.MNTX": {"net":"US","sta":"MNTX","lat":31.039,"lon":-104.068,"elev":1415.0},
    "US.NLWA": {"net":"US","sta":"NLWA","lat":47.735,"lon":-116.322,"elev":760.0},
    "US.RSSD": {"net":"US","sta":"RSSD","lat":44.121,"lon":-104.036,"elev":1171.0},
    "US.TZTN": {"net":"US","sta":"TZTN","lat":36.544,"lon":-83.552, "elev":370.0},
    "US.WMOK": {"net":"US","sta":"WMOK","lat":34.737,"lon":-98.781, "elev":434.0},
    # PR Puerto Rico
    "PR.ANWB": {"net":"PR","sta":"ANWB","lat":18.489,"lon":-66.288, "elev":480.0},
    "PR.CTB":  {"net":"PR","sta":"CTB", "lat":18.218,"lon":-65.699, "elev":152.0},
    "PR.FRY":  {"net":"PR","sta":"FRY", "lat":18.289,"lon":-65.628, "elev":52.0},
    "PR.MZCY": {"net":"PR","sta":"MZCY","lat":18.073,"lon":-67.122, "elev":125.0},
    # IU Global
    "IU.MSKU": {"net":"IU","sta":"MSKU","lat":1.657,"lon":16.282,"elev":316.0},
    "IU.PMSA": {"net":"IU","sta":"PMSA","lat":-64.774,"lon":-64.049,"elev":40.0},
    "IU.PTGA": {"net":"IU","sta":"PTGA","lat":-0.729,"lon":-59.966,"elev":80.0},
    "IU.QSPA": {"net":"IU","sta":"QSPA","lat":-89.928,"lon":144.437,"elev":2880.0},
    "IU.RCBR": {"net":"IU","sta":"RCBR","lat":-5.827,"lon":-35.901,"elev":290.0},
    "IU.SAML": {"net":"IU","sta":"SAML","lat":-8.947,"lon":-63.183,"elev":109.0},
    "IU.LVC":  {"net":"IU","sta":"LVC", "lat":-22.614,"lon":-68.911,"elev":2939.0},
    "IU.GNI":  {"net":"IU","sta":"GNI", "lat":40.148,"lon":44.741, "elev":1509.0},
    "IU.PAB":  {"net":"IU","sta":"PAB", "lat":39.545,"lon":-4.350, "elev":950.0},
    "IU.KEV":  {"net":"IU","sta":"KEV", "lat":69.755,"lon":27.007, "elev":85.0},
    "IU.TRIS": {"net":"IU","sta":"TRIS","lat":-37.068,"lon":-12.315,"elev":35.0},
    "IU.YELL": {"net":"IU","sta":"YELL","lat":62.490,"lon":-114.609,"elev":180.0},
    "IU.TATO": {"net":"IU","sta":"TATO","lat":24.975,"lon":121.497,"elev":128.0},
    "IU.MAJO": {"net":"IU","sta":"MAJO","lat":36.543,"lon":138.207,"elev":418.0},
    "IU.SNZO": {"net":"IU","sta":"SNZO","lat":-41.310,"lon":174.704,"elev":100.0},
    # II IRIS/IDA
    "II.AAK":  {"net":"II","sta":"AAK", "lat":42.639,"lon":74.494, "elev":1645.0},
    "II.ABKT": {"net":"II","sta":"ABKT","lat":37.930,"lon":58.119, "elev":678.0},
    "II.ALE":  {"net":"II","sta":"ALE", "lat":82.504,"lon":-62.350,"elev":60.0},
    "II.BFO":  {"net":"II","sta":"BFO", "lat":48.331,"lon":8.330,  "elev":589.0},
    "II.BORG": {"net":"II","sta":"BORG","lat":64.747,"lon":-21.327,"elev":110.0},
    "II.EFI":  {"net":"II","sta":"EFI", "lat":-51.675,"lon":-58.063,"elev":110.0},
    "II.ESK":  {"net":"II","sta":"ESK", "lat":55.317,"lon":-3.205, "elev":242.0},
    "II.GAR":  {"net":"II","sta":"GAR", "lat":39.000,"lon":70.317, "elev":1300.0},
    "II.KDAK": {"net":"II","sta":"KDAK","lat":57.783,"lon":-152.583,"elev":153.0},
    "II.KURK": {"net":"II","sta":"KURK","lat":50.775,"lon":78.621, "elev":268.0},
    "II.KWAJ": {"net":"II","sta":"KWAJ","lat":8.802, "lon":167.612,"elev":0.0},
    "II.MBAR": {"net":"II","sta":"MBAR","lat":-0.602,"lon":30.738, "elev":1375.0},
    "II.NNA":  {"net":"II","sta":"NNA", "lat":-11.988,"lon":-76.842,"elev":575.0},
    "II.OBN":  {"net":"II","sta":"OBN", "lat":55.113,"lon":36.569, "elev":160.0},
    "II.SACV": {"net":"II","sta":"SACV","lat":14.970,"lon":-23.608,"elev":387.0},
    "II.SHEL": {"net":"II","sta":"SHEL","lat":-15.961,"lon":-5.745,"elev":568.0},
    "II.SUR":  {"net":"II","sta":"SUR", "lat":-32.379,"lon":20.811,"elev":1770.0},
    "II.TLY":  {"net":"II","sta":"TLY", "lat":51.681,"lon":103.644,"elev":579.0},
    "II.UNM":  {"net":"II","sta":"UNM", "lat":19.330,"lon":-99.178,"elev":2279.0},
    "II.WRAB": {"net":"II","sta":"WRAB","lat":-19.934,"lon":134.360,"elev":381.0},
    # CN Canada
    "CN.DAWY": {"net":"CN","sta":"DAWY","lat":59.896,"lon":-128.820,"elev":698.0},
    "CN.EDM":  {"net":"CN","sta":"EDM", "lat":53.224,"lon":-113.346,"elev":739.0},
    "CN.FFC":  {"net":"CN","sta":"FFC", "lat":54.725,"lon":-101.978,"elev":338.0},
    "CN.LLLB": {"net":"CN","sta":"LLLB","lat":50.629,"lon":-121.877,"elev":352.0},
    "CN.MBC":  {"net":"CN","sta":"MBC", "lat":76.237,"lon":-119.355,"elev":86.0},
    "CN.MOBC": {"net":"CN","sta":"MOBC","lat":59.468,"lon":-136.037,"elev":655.0},
    "CN.SADO": {"net":"CN","sta":"SADO","lat":44.770,"lon":-79.143, "elev":280.0},
    "CN.WHY":  {"net":"CN","sta":"WHY", "lat":60.659,"lon":-134.883,"elev":760.0},
    # MB Montana/Idaho
    "MB.BOZM": {"net":"MB","sta":"BOZM","lat":45.604,"lon":-111.046,"elev":1506.0},
    "MB.BRTE": {"net":"MB","sta":"BRTE","lat":44.569,"lon":-110.872,"elev":2025.0},
    "MB.COTM": {"net":"MB","sta":"COTM","lat":44.765,"lon":-111.482,"elev":1670.0},
    "MB.DLMT": {"net":"MB","sta":"DLMT","lat":45.361,"lon":-112.540,"elev":1580.0},
    "MB.DRAM": {"net":"MB","sta":"DRAM","lat":46.002,"lon":-112.680,"elev":1540.0},
    "MB.FBIM": {"net":"MB","sta":"FBIM","lat":48.051,"lon":-111.964,"elev":1063.0},
    "MB.GILM": {"net":"MB","sta":"GILM","lat":48.195,"lon":-106.632,"elev":640.0},
    "MB.MAPM": {"net":"MB","sta":"MAPM","lat":46.909,"lon":-110.572,"elev":1780.0},
    "MB.MOUM": {"net":"MB","sta":"MOUM","lat":46.872,"lon":-114.017,"elev":988.0},
    "MB.MSPM": {"net":"MB","sta":"MSPM","lat":47.201,"lon":-113.986,"elev":1006.0},
    "MB.RUDM": {"net":"MB","sta":"RUDM","lat":47.050,"lon":-107.990,"elev":820.0},
    "MB.SHRM": {"net":"MB","sta":"SHRM","lat":47.540,"lon":-107.390,"elev":930.0},
    "MB.STGM": {"net":"MB","sta":"STGM","lat":44.666,"lon":-111.106,"elev":1920.0},
    "MB.TOSM": {"net":"MB","sta":"TOSM","lat":46.401,"lon":-109.975,"elev":1620.0},
    "MB.WARM": {"net":"MB","sta":"WARM","lat":45.520,"lon":-111.530,"elev":1490.0},
    "MB.WERM": {"net":"MB","sta":"WERM","lat":48.122,"lon":-106.032,"elev":740.0},
    "MB.WILM": {"net":"MB","sta":"WILM","lat":46.020,"lon":-104.620,"elev":890.0},
}

NET_COLORS = {
    "CI":"#e05252","AZ":"#e07832","BK":"#3a8fd4",
    "IU":"#3dba6f","UW":"#9b6dd6","UO":"#2ab5a0",
    "NN":"#d4b800","NC":"#c04040","SB":"#f0a040",
    "UU":"#e8855a","WY":"#a8d86e","MB":"#7ec8c8",
    "AK":"#5ab4e8","HV":"#e85ab4","NM":"#c8a87e",
    "US":"#8888cc","PR":"#e8c85a","AT":"#88cc88",
    "GE":"#f07070","G":"#70c4f0","GT":"#90d090","IC":"#f0c070",
    "PS":"#c090f0","AU":"#70f0c0","MN":"#f09040","KZ":"#a0c8e0",
    "AV":"#f0b060","NE":"#80c8a0","TA":"#9090d0","SE":"#d0a060",
    "CU":"#60d0d0","TX":"#e0a080","OR":"#80e0a0","OK":"#d0c060",
    "IN":"#a0b0d0","LD":"#c0a0d0","AE":"#e0c0a0","KO":"#80b0e0",
    "RO":"#d08080","CZ":"#b0d080","HU":"#e0b090","PL":"#90c0b0",
    "NO":"#70b0d0","PM":"#d0b0a0","ES":"#e09070","IV":"#c0c080",
    "CH":"#a0d0c0","BE":"#b0a0e0","AF":"#d0a060","AY":"#90c0d0",
    "C":"#e08060","BR":"#80d0a0","NZ":"#70d0b0","SA":"#c0b060",
    "IG":"#d0a070","TW":"#90b0e0","OO":"#c8e880",
}

# ── Cities (western US focus) ──────────────────────────────────────────────────
# ── Cities database — loaded from GeoNames at startup ─────────────────────────
_CITIES_DATA: list = []   # list of (name, state_label, lat, lon)

def _load_cities():
    global _CITIES_DATA
    import os as _os
    csv_path = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "cities_na.csv")
    rows = []
    try:
        with open(csv_path, encoding="utf-8") as f:
            for line in f:
                p = line.strip().split(",")
                if len(p) < 5: continue
                name, cc, state = p[0], p[1], p[2]
                try: lat, lon = float(p[3]), float(p[4])
                except: continue
                label = state if cc == "US" else cc
                rows.append((name, label, lat, lon))
        _CITIES_DATA = rows
        print(f"[CITIES] loaded {len(rows)} cities from GeoNames", flush=True)
    except FileNotFoundError:
        _CITIES_DATA = [
            ("Los Angeles","CA",34.052,-118.244),("San Francisco","CA",37.774,-122.419),
            ("San Diego","CA",32.716,-117.161),("Sacramento","CA",38.582,-121.494),
            ("Fresno","CA",36.747,-119.773),("Las Vegas","NV",36.175,-115.137),
        ]
        print("[CITIES] cities_na.csv not found — using fallback", flush=True)

_load_cities()

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat/2)**2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(dlon/2)**2)
    return R * 2 * math.asin(math.sqrt(a))

def nearest_cities(lat, lon, n=3):
    if not _CITIES_DATA:
        return []
    distances = [(haversine(lat, lon, c[2], c[3]), c[0], c[1]) for c in _CITIES_DATA]
    distances.sort()
    return [{"name": name, "state": state, "km": round(km)} for km, name, state in distances[:n]]


# ── MMI helpers ─────────────────────────────────────────────────────────────────
_MMI_ROMAN  = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII']
_MMI_DESC   = ['Not felt','Weak','Weak','Light','Moderate','Strong',
               'Very Strong','Severe','Violent','Extreme','Extreme','Extreme']

def _estimate_mmi(mag, depth_km=10, epi_dist_km=0):
    """
    Estimate MMI using:
    - Boore et al. (1997) WUS rock PGA attenuation:
        ln(PGA_g) = -3.72 + 1.02*M - 1.57*ln(R_hyp) - 0.0096*R_hyp
    - Wald et al. (1999) PGA → MMI conversion:
        MMI = 3.66*log10(PGA_cms2) - 1.66  (PGA >= 1.9 cm/s²)
        MMI = 3.60*log10(PGA_cms2) + 0.70  (PGA < 1.9 cm/s²)
    Calibrated for shallow crustal WUS earthquakes (USGS ShakeMap standard).
    Returns float MMI (1.0–12.0).
    """
    R_hyp = math.sqrt(max(epi_dist_km, 0.0)**2 + max(depth_km, 3.0)**2)
    R_hyp = max(R_hyp, 5.0)
    ln_pga_g = (-3.72 + 1.02 * mag
                - 1.57 * math.log(R_hyp)
                - 0.0096 * R_hyp)
    pga_cms2 = math.exp(ln_pga_g) * 980.665  # g → cm/s²
    if pga_cms2 >= 1.9:
        mmi = 3.66 * math.log10(pga_cms2) - 1.66
    else:
        mmi = 3.60 * math.log10(max(pga_cms2, 0.001)) + 0.70
    return round(max(1.0, min(12.0, mmi)), 1)

def _mmi_roman(mmi_float):
    """Return Roman numeral string for an MMI value."""
    idx = max(0, min(11, round(mmi_float) - 1))
    return _MMI_ROMAN[idx]

def _mmi_desc(mmi_float):
    """Return short text description for an MMI value."""
    idx = max(0, min(11, round(mmi_float) - 1))
    return _MMI_DESC[idx]


# ── Travel-time grid search epicenter estimator ────────────────────────────────
def _grid_search_epicenter(stations, vp=6.0):
    """
    Iterative multi-resolution grid search for P-wave origin location.
    Uses travel-time residual minimization: find (lat, lon, origin_time) that
    minimizes the RMS of (observed_arrival - predicted_arrival) across stations.

    stations : list of dicts with keys lat, lon, ts (trigger timestamp), stalta
    vp       : P-wave velocity km/s — use 6.0 for crustal, 8.5 for teleseismic
    Returns  : (lat_est, lon_est, origin_time_est, rms_seconds, filtered_stations)
    """
    if len(stations) < 2:
        lats = [s["lat"] for s in stations]
        lons = [s["lon"] for s in stations]
        return (sum(lats)/len(lats), sum(lons)/len(lons),
                min(s["ts"] for s in stations) - 5.0, 99.0, stations)

    # Initial guess: amplitude-weighted centroid
    total_w = sum(s["stalta"]**2 for s in stations)
    best_lat = sum(s["lat"] * s["stalta"]**2 for s in stations) / total_w
    best_lon = sum(s["lon"] * s["stalta"]**2 for s in stations) / total_w
    best_rms  = float('inf')
    best_ot   = min(s["ts"] for s in stations)

    def _eval(test_lat, test_lon, stas):
        """Return (rms, origin_time) for given epicenter and station list."""
        ots = [s["ts"] - haversine(test_lat, test_lon, s["lat"], s["lon"]) / vp
               for s in stas]
        ot_mean = sum(ots) / len(ots)
        rms = math.sqrt(sum((o - ot_mean)**2 for o in ots) / len(ots))
        return rms, ot_mean

    # 3-pass zoom: 4°→0.4°→0.04° step (local) or 20°→2°→0.2° (teleseismic)
    span = 20.0 if vp >= 8.0 else 4.0
    for _pass in range(3):
        step = span / 10.0
        for dlat_i in range(-10, 11):
            for dlon_i in range(-10, 11):
                rms, ot = _eval(best_lat + dlat_i*step, best_lon + dlon_i*step, stations)
                if rms < best_rms:
                    best_rms = rms
                    best_lat = best_lat + dlat_i*step
                    best_lon = best_lon + dlon_i*step
                    best_ot  = ot
        span = step  # zoom in to refined area

    # Travel-time consistency filter
    # Teleseisms use looser individual residual but tighter overall RMS requirement
    resid_cutoff = 20.0 if vp >= 8.0 else 12.0

    def _residual(s):
        pred = best_ot + haversine(best_lat, best_lon, s["lat"], s["lon"]) / vp
        return abs(s["ts"] - pred)

    filtered = [s for s in stations if _residual(s) < resid_cutoff]
    if len(filtered) >= PRELIM_MIN_STATIONS:
        # Re-run one final pass on cleaned cluster
        best_rms_f = float('inf')
        span = 2.0 if vp >= 8.0 else 1.0
        for _pass in range(2):
            step = span / 10.0
            for dlat_i in range(-10, 11):
                for dlon_i in range(-10, 11):
                    rms, ot = _eval(best_lat + dlat_i*step, best_lon + dlon_i*step, filtered)
                    if rms < best_rms_f:
                        best_rms_f = rms
                        best_lat = best_lat + dlat_i*step
                        best_lon = best_lon + dlon_i*step
                        best_ot  = ot
            span = step
        return best_lat, best_lon, best_ot, best_rms_f, filtered
    else:
        return best_lat, best_lon, best_ot, best_rms, stations


# ── Station state ──────────────────────────────────────────────────────────────
class StationState:
    def __init__(self, net, sta, chan):
        self.net = net; self.sta = sta; self.chan = chan
        self.key = f"{net}.{sta}.{chan}"
        self.buf    = collections.deque()   # full-res, last BUF_SECS seconds
        self.buf_ds = collections.deque()   # ~1 sps downsampled, last BUF_LONG seconds
        self._ds_skip = 0                   # sample counter for downsampling
        self.srate = None; self.stalta = 0.0; self.last_update = None
        self.stalta_hist = collections.deque(maxlen=14400)  # 4h at 1 sample/s
        self.triggered = False   # True while stalta >= TRIG
        self.peak_raw  = 1.0    # peak absolute count (for amplitude label)
        self.last_polled_t = None  # for FDSNWS polling deduplication
        self.lock = threading.Lock()
        self.city = ""   # nearest city label, populated after coord load

    INACTIVE_SECS = 600   # station considered inactive after 10 min without data

    @property
    def is_active(self):
        """True if data was received within the last 5 minutes."""
        if self.last_update is None:
            return False
        age = (datetime.now(timezone.utc) - self.last_update).total_seconds()
        return age < self.INACTIVE_SECS

    def feed(self, trace):
        sr = trace.stats.sampling_rate
        t0 = trace.stats.starttime.timestamp
        data = trace.data.astype(float)
        with self.lock:
            self.srate = sr
            dt = 1.0 / sr
            ds_n = max(1, int(sr))          # downsample ratio → ~1 sps
            t_last = t0 + (len(data) - 1) * dt

            for i, s in enumerate(data):
                t = t0 + i * dt
                self.buf.append((t, s))
                # Populate downsampled long buffer
                self._ds_skip += 1
                if self._ds_skip >= ds_n:
                    self._ds_skip = 0
                    self.buf_ds.append((t, s))

            # Trim both buffers
            cutoff      = t_last - BUF_SECS
            cutoff_long = t_last - BUF_LONG
            while self.buf    and self.buf[0][0]    < cutoff:      self.buf.popleft()
            while self.buf_ds and self.buf_ds[0][0] < cutoff_long: self.buf_ds.popleft()

            # Mark station as live as soon as any data arrives
            self.last_update = datetime.now(timezone.utc)

            # Track peak raw amplitude for amplitude label
            pk = float(np.max(np.abs(data))) if len(data) else 0.0
            if pk > self.peak_raw: self.peak_raw = pk

            # STA/LTA on full-res buffer (needs LTA_WIN + STA_WIN seconds of data)
            n = len(self.buf)
            n_lta = int(LTA_WIN * sr); n_sta = int(STA_WIN * sr)
            if n >= n_lta + n_sta:
                arr = np.array([s for _, s in self.buf]); arr2 = arr ** 2
                sv = np.mean(arr2[-n_sta:]); lv = np.mean(arr2[-(n_lta + n_sta):-n_sta])
                self.stalta = math.sqrt(sv / lv) if lv > 0 else 0.0

    def _pts_from(self, source, n_pts, secs):
        """Internal: slice last `secs` from a deque, downsample to n_pts, return (values, t_end)."""
        now = time.time()
        seg = [(t, s) for t, s in source if t >= now - secs]
        if not seg: return [], None
        t_end = seg[-1][0]
        arr = [s for _, s in seg]
        peak = max(abs(v) for v in arr) or 1.0
        arr = [v / peak for v in arr]
        if len(arr) > n_pts:
            step = len(arr) / n_pts
            arr = [arr[int(i * step)] for i in range(n_pts)]
        return arr, t_end

    def waveform_pts(self, n_pts=900, secs=600):
        """Full-res buffer — good for zoom ≤ 10 min."""
        with self.lock:
            if not self.buf: return [], None
            return self._pts_from(self.buf, n_pts, min(secs, BUF_SECS))

    def waveform_pts_long(self, n_pts=720, secs=7200):
        """Downsampled buffer — good for zoom 10 min – 2 hr."""
        with self.lock:
            if not self.buf_ds: return [], None
            return self._pts_from(self.buf_ds, n_pts, min(secs, BUF_LONG))

    def full_buffer(self, n_pts=2000):
        """Return full ring buffer as (timestamps, values) for station page."""
        with self.lock:
            if not self.buf: return [], []
            ts, vals = zip(*self.buf) if self.buf else ([], [])
            ts, vals = list(ts), list(vals)
            peak = max(abs(v) for v in vals) or 1.0
            vals = [v / peak for v in vals]
            if len(ts) > n_pts:
                step = len(ts) / n_pts
                idx = [int(i * step) for i in range(n_pts)]
                ts = [ts[i] for i in idx]; vals = [vals[i] for i in idx]
            t0 = ts[0] if ts else 0
            return [round(t - t0, 2) for t in ts], vals

# Build initial states from STREAMS only (FDSN poller dynamically adds more)
_all_streams = list({(net, sta, chan) for net, sta, chan in STREAMS})
states = {f"{net}.{sta}.{chan}": StationState(net, sta, chan)
          for net, sta, chan in _all_streams}

# Global trigger log — (timestamp, net, sta, lat, lon, stalta)
trigger_log = collections.deque(maxlen=2000)
preliminary_events = collections.deque(maxlen=200)  # real-time detections (persisted to disk)

# ── Station quality / anomaly tracking ─────────────────────────────────────────
_station_trigger_times: dict = collections.defaultdict(list)   # sta_key→[timestamps]
_noisy_stations: set = set()           # stations excluded from prelim detection
NOISE_TRIGGER_THRESH = 10              # triggers/hr before flagging as noisy
NOISE_CLEAR_THRESH   = 3               # triggers/hr before un-flagging

def _record_trigger(sta_key: str, ts: float):
    """Record a trigger event for a station."""
    _station_trigger_times[sta_key].append(ts)

def _update_station_quality():
    """Refresh noisy-station set every minute."""
    cutoff = time.time() - 3600
    for key in list(_station_trigger_times.keys()):
        _station_trigger_times[key] = [t for t in _station_trigger_times[key] if t > cutoff]
        count = len(_station_trigger_times[key])
        if count >= NOISE_TRIGGER_THRESH and key not in _noisy_stations:
            _noisy_stations.add(key)
            print(f"[QUAL] {key} flagged noisy: {count} triggers/hr", flush=True)
        elif count <= NOISE_CLEAR_THRESH and key in _noisy_stations:
            _noisy_stations.discard(key)
            print(f"[QUAL] {key} cleared: {count} triggers/hr", flush=True)

def _fetch_iris_inventory_for(nets: list):
    """
    Query IRIS FDSN station service for coordinates of all active stations
    in the given networks.  Updates _sta_coords in-place.
    Safe to call from a background thread.
    """
    import urllib.request, io, csv
    for net in nets:
        try:
            url = (f"https://service.earthscope.org/fdsnws/station/1/query"
                   f"?net={net}&format=text&level=station"
                   f"&starttime=2020-01-01&nodata=404")
            req = urllib.request.Request(url, headers={"User-Agent":"SeisComp-Dashboard/1.0"})
            with urllib.request.urlopen(req, timeout=20) as r:
                lines = r.read().decode().splitlines()
            added = 0
            for line in lines:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split('|')
                if len(parts) < 6:
                    continue
                try:
                    # FDSN station text format (level=station):
                    # Network|Station|Latitude|Longitude|Elevation|SiteName|StartTime|EndTime
                    n   = parts[0].strip()
                    s   = parts[1].strip()
                    lat = float(parts[2])
                    lon = float(parts[3])
                    elev = float(parts[4]) if len(parts) > 4 and parts[4].strip() else 0.0
                    key = f"{n}.{s}"
                    if key not in _sta_coords:
                        _sta_coords[key] = {"net": n, "sta": s, "lat": lat, "lon": lon, "elev": elev}
                        added += 1
                except (ValueError, IndexError):
                    continue
            if added:
                print(f"[INV] {net}: loaded {added} new station coords from IRIS", flush=True)
        except Exception as e:
            print(f"[INV] {net}: inventory fetch failed: {e}", flush=True)

# ── Persistent detection log ────────────────────────────────────────────────────
import json as _json
PRELIM_LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "prelim_detections.json")

def _load_prelim_log():
    """Load persisted detections from disk into preliminary_events on startup."""
    try:
        with open(PRELIM_LOG_FILE) as f:
            saved = _json.load(f)
        for ev in reversed(saved):   # oldest first so deque ends up newest-first
            preliminary_events.appendleft(ev)
        print(f"[PRELIM] Loaded {len(saved)} saved detections from disk")
    except FileNotFoundError:
        pass
    except Exception as e:
        print(f"[PRELIM] Could not load detections log: {e}")

def _save_prelim_log():
    """Persist current preliminary_events to disk."""
    try:
        with open(PRELIM_LOG_FILE, "w") as f:
            _json.dump(list(preliminary_events), f)
    except Exception as e:
        print(f"[PRELIM] Could not save detections log: {e}")

# ── Station coordinate cache (loaded from DB at startup) ───────────────────────
_sta_coords = {}   # "NET.STA" -> {net, sta, lat, lon, elev}

def _load_sta_coords():
    """Populate _sta_coords from the SeisComP inventory DB."""
    global _sta_coords
    try:
        cn = db(); cur = cn.cursor(dictionary=True)
        cur.execute("""
            SELECT n.code as net, s.code as sta,
                   s.latitude as lat, s.longitude as lon,
                   s.elevation as elev
            FROM   Station s
            JOIN   Network n ON n._oid = s._parent_oid
            WHERE  s.latitude IS NOT NULL AND s.longitude IS NOT NULL
        """)
        tmp = {}
        for r in cur.fetchall():
            key = f"{r['net']}.{r['sta']}"
            tmp[key] = {
                "net": r["net"], "sta": r["sta"],
                "lat": float(r["lat"]), "lon": float(r["lon"]),
                "elev": float(r["elev"] or 0),
            }
        # Merge in hardcoded IRIS Bay Area BK stations
        tmp.update(_IRIS_STA_COORDS)
        _sta_coords = tmp
        cur.close(); cn.close()
        print(f"[INFO] Loaded {len(_sta_coords)} station coordinates ({len(_IRIS_STA_COORDS)} IRIS hardcoded)")
    except Exception as e:
        print(f"[WARN] _load_sta_coords: {e}")

def _populate_cities():
    """Fill nearest city for every station coord and state."""
    for key, coord in list(_sta_coords.items()):  # snapshot to avoid dict-changed-size-during-iteration
        lat, lon = coord.get("lat"), coord.get("lon")
        if lat is None or lon is None:
            continue
        cities = nearest_cities(lat, lon, n=1)
        if cities:
            c = cities[0]
            coord["city"] = f"{c['name']}, {c['state']} ({c['km']} km)"
    # Also update live station states
    for key, state in list(states.items()):
        skey = f"{state.net}.{state.sta}"
        if skey in _sta_coords:
            state.city = _sta_coords[skey].get("city","")

# ── Preliminary real-time earthquake detector ─────────────────────────────────
import ssl as _ssl, json as _json, smtplib as _smtplib
from email.mime.multipart import MIMEMultipart as _MIME_MP
from email.mime.text import MIMEText as _MIME_T

_tz_pt = __import__('zoneinfo').ZoneInfo("America/Los_Angeles")

def _pt(utc_dt):
    """Convert a UTC datetime to Pacific Time string with label."""
    pt = utc_dt.astimezone(_tz_pt)
    lbl = "PDT" if pt.dst() and pt.dst().seconds else "PST"
    return pt.strftime("%Y-%m-%d %H:%M:%S") + f" {lbl}", lbl, pt

def _ago(utc_dt):
    sec = int((datetime.now(timezone.utc) - utc_dt).total_seconds())
    if sec < 60:   return f"{sec}s ago"
    if sec < 3600: return f"{sec//60}m {sec%60}s ago"
    return f"{sec//3600}h {(sec%3600)//60}m ago"

def _make_map_image(lat, lon, zoom, width=640, height=320):
    """
    Render a full-color map tile using staticmap + Stadia Alidade Smooth tiles.
    Returns PNG bytes or None on failure.
    """
    try:
        from staticmap import StaticMap, CircleMarker, Line
        # Standard OSM tiles (free, reliable, full color)
        tile_url = "https://tile.openstreetmap.org/{z}/{x}/{y}.png"
        m = StaticMap(width, height, url_template=tile_url, headers={
            "User-Agent": "SeisCompDashboard/1.0 (earthquake monitoring; contact: seiscomp@localhost)"
        })
        # Red epicenter dot
        m.add_marker(CircleMarker((lon, lat), "#ff4444", 16))
        m.add_marker(CircleMarker((lon, lat), "#ffffff", 8))
        img = m.render(zoom=zoom)
        buf = __import__('io').BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    except Exception as _me:
        print(f"[DISCORD] map render failed: {_me}", flush=True)
        return None

def _post_discord(webhook, payload_dict, image_bytes=None, image_name="map.png", label=""):
    """Post a Discord embed to one webhook, optionally with a PNG attachment."""
    import urllib.request as _ur, urllib.error as _ue, json as _j, ssl as _ssl
    ctx = _ssl.create_default_context()
    for attempt in range(3):
        try:
            if image_bytes:
                # multipart/form-data so we can attach the image and reference it
                boundary = "----SeisCompBoundary7x"
                body = (
                    f"--{boundary}\r\n"
                    f'Content-Disposition: form-data; name="payload_json"\r\n\r\n'
                    + _j.dumps(payload_dict)
                    + f"\r\n--{boundary}\r\n"
                    f'Content-Disposition: form-data; name="file"; filename="{image_name}"\r\n'
                    f"Content-Type: image/png\r\n\r\n"
                ).encode() + image_bytes + f"\r\n--{boundary}--\r\n".encode()
                req = _ur.Request(webhook, data=body,
                                  headers={"Content-Type": f"multipart/form-data; boundary={boundary}",
                                           "User-Agent": "DiscordBot (seiscomp, 1.0)"},
                                  method="POST")
            else:
                req = _ur.Request(webhook,
                                  data=_j.dumps(payload_dict).encode(),
                                  headers={"Content-Type": "application/json",
                                           "User-Agent": "DiscordBot (seiscomp, 1.0)"},
                                  method="POST")
            _ur.urlopen(req, context=ctx, timeout=15)
            print(f"[DISCORD] {label} sent", flush=True); return
        except _ue.HTTPError as e:
            body_r = e.read().decode()
            print(f"[DISCORD] {label} HTTP {e.code}: {body_r[:120]}")
            if e.code == 429:
                try: time.sleep(float(_j.loads(body_r).get("retry_after", 2)))
                except: time.sleep(2)
            elif e.code < 500: return
        except Exception as ex:
            print(f"[DISCORD] {label} error attempt {attempt+1}: {ex}")
        time.sleep(1)

def _broadcast_discord(payload_dict, image_bytes=None, image_name="map.png", label=""):
    """Send payload to every webhook in DISCORD_WEBHOOKS."""
    for wh in DISCORD_WEBHOOKS:
        _post_discord(wh, payload_dict, image_bytes, image_name, label)

def _send_prelim_discord(ev):
    if not PRELIM_DISCORD_ENABLED: return
    mag    = ev["mag_est"]
    lat, lon = ev["lat"], ev["lon"]
    ns, ew = ("N" if lat >= 0 else "S"), ("W" if lon <= 0 else "E")
    loc    = f"{abs(lat):.3f}°{ns}  {abs(lon):.3f}°{ew}"
    color  = (0xF85149 if mag>=5 else 0xFF6B35 if mag>=4 else 0xD4B800 if mag>=3
              else 0x3FB950 if mag>=2 else 0x58A6FF)
    zoom   = 6 if mag >= 4 else 7 if mag >= 2.5 else 8
    is_tele = ev.get("teleseism", False)
    spread  = ev.get("spread_km", 0)

    try:
        ot_utc = datetime.strptime(ev["ot_str"], "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
    except ValueError:
        ot_utc = datetime.now(timezone.utc)
    pt_str, pt_lbl, _ = _pt(ot_utc)
    ago_str = _ago(ot_utc)

    det_type = "TELESEISM" if is_tele else "PRELIMINARY"
    warn_flag = "⚠ " if mag >= 4 else ""
    title    = f"{warn_flag}{det_type} — Est. M{mag:.1f}"
    footer   = ("Teleseism — large distant event lit up stations across CA"
                if is_tele else "Preliminary detection — awaiting SeisComP confirmation")
    osm_url  = f"https://www.openstreetmap.org/?mlat={lat:.4f}&mlon={lon:.4f}#map={zoom}/{lat:.4f}/{lon:.4f}"

    mmi_val      = ev.get("mmi", 0)
    mmi_str      = ev.get("mmi_str", "?")
    mmi_desc     = ev.get("mmi_desc", "")
    mmi_city_str = ev.get("mmi_city_str", "?")
    rms_str      = f"{ev.get('rms_sec', '?')} s" if ev.get('rms_sec') is not None else "?"

    img = _make_map_image(lat, lon, zoom)
    image_ref = "attachment://map.png" if img else None

    # @everyone fires at most once per 30-minute window for prelim detections
    import time as _time
    global _everyone_prelim_last_t
    do_everyone = (mag >= DISCORD_EVERYONE_MAG
                   and (_time.time() - _everyone_prelim_last_t) >= _EVERYONE_PRELIM_COOLDOWN)
    if do_everyone:
        _everyone_prelim_last_t = _time.time()

    payload = {
        "username": "SeisComP Monitor",
        "content": "@everyone" if do_everyone else "",
        "embeds": [{
            "title": title,
            "url": osm_url,
            "description": (
                f"**~M{mag:.1f}** · {ev.get('city','—')}\n"
                f"> `{ev['n_stations']} stations`  ·  max STA/LTA `{ev['max_stalta']}`  ·  spread `{int(spread)} km`"
            ),
            "color": color,
            "fields": [
                {"name": "Location",             "value": f"`{loc}`",              "inline": True},
                {"name": "Depth (est.)",          "value": f"~{ev['depth_est']} km","inline": True},
                {"name": "Time Since",            "value": ago_str,                 "inline": True},
                {"name": f"Origin ({pt_lbl})",    "value": f"`{pt_str}`",           "inline": False},
                {"name": "Origin (UTC)",           "value": f"`{ev['ot_str']} UTC`", "inline": False},
                {"name": "Max MMI (epicenter)",    "value": f"`{mmi_str}` — {mmi_desc}","inline": True},
                {"name": f"MMI at {ev.get('city','nearest city').split(',')[0] if ev.get('city') else 'nearest city'}", "value": f"`{mmi_city_str}`", "inline": True},
                {"name": "Location RMS",           "value": f"`{rms_str}`",          "inline": True},
                {"name": "Triggered By",           "value": f"`{ev.get('stations','—')}`", "inline": False},
                *([ {"name": "🚨 Earthquake Safety",
                     "value": "**DROP** to the ground · take **COVER** under a sturdy table or against an interior wall · **HOLD ON** until shaking stops. Stay away from windows and heavy objects.",
                     "inline": False} ] if mag >= 3.0 else []),
                {"name": "​", "value": "This is a preliminary earthquake alert system that uses SeisComP compiled by Evan Li. For more information visit: https://myearthquake.dpdns.org/", "inline": False},
            ],
            "image": {"url": image_ref} if image_ref else {},
            "thumbnail": {"url": "https://dashboard.myearthquake.dpdns.org/static/drop_cover_hold.png"} if mag >= 3.0 else {},
            "footer": {"text": footer, "icon_url": "https://earthquake.usgs.gov/favicon.ico"},
            "timestamp": ev["ot_str"].replace(" ", "T") + "Z",
        }]
    }
    _broadcast_discord(payload, img, "map.png",
                       label=f"PRELIM M{mag:.1f}")


def _send_event_discord(ev):
    """Announce a confirmed SeisComP catalog event to Discord."""
    if not EVENT_DISCORD_ENABLED: return
    mag  = float(ev.get("mag") or 0)
    if mag < EVENT_DISCORD_MIN_MAG: return
    evid = ev.get("evid", "?")
    lat  = float(ev.get("lat") or 0)
    lon  = float(ev.get("lon") or 0)
    dep  = ev.get("depth", "?")
    ot   = str(ev.get("origin_time", ""))[:19]
    ns, ew = ("N" if lat >= 0 else "S"), ("W" if lon <= 0 else "E")
    loc  = f"{abs(lat):.3f}°{ns}  {abs(lon):.3f}°{ew}"
    color = (0xF85149 if mag>=5 else 0xFF6B35 if mag>=4 else 0xFF9500 if mag>=3
             else 0xD29922 if mag>=2 else 0x3FB950)
    zoom  = 6 if mag >= 5 else 7 if mag >= 3 else 8
    warn_flag = "⚠ " if mag >= 4 else ""
    phases = ev.get("phases", "?")
    cities = nearest_cities(lat, lon, n=1)
    c0     = cities[0] if cities else None
    city   = f"{c0['name']}, {c0['state']} ({c0['km']} km)" if c0 else "—"
    city_dist = c0['km'] if c0 else 50
    osm_url = f"https://www.openstreetmap.org/?mlat={lat:.4f}&mlon={lon:.4f}#map={zoom}/{lat:.4f}/{lon:.4f}"
    try:
        ot_utc = datetime.strptime(ot, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
    except ValueError:
        ot_utc = datetime.now(timezone.utc)
    pt_str, pt_lbl, _ = _pt(ot_utc)
    ago_str = _ago(ot_utc)
    ev_mmi_val  = _estimate_mmi(mag, depth_km=float(dep) if str(dep).replace('.','').isdigit() else 10, epi_dist_km=0)
    ev_mmi_str  = _mmi_roman(ev_mmi_val)
    ev_mmi_desc = _mmi_desc(ev_mmi_val)
    ev_mmi_city = _mmi_roman(_estimate_mmi(mag, depth_km=10, epi_dist_km=city_dist))

    img = _make_map_image(lat, lon, zoom)
    image_ref = "attachment://map.png" if img else None

    # @everyone fires at most once per confirmed event ID
    do_everyone = (mag >= DISCORD_EVERYONE_MAG and evid not in _everyone_pinged_events)
    if do_everyone:
        _everyone_pinged_events.add(evid)

    payload = {
        "username": "SeisComP Monitor",
        "content": "@everyone" if do_everyone else "",
        "embeds": [{
            "title": f"{warn_flag}M{mag:.1f} Earthquake — SeisComP Confirmed",
            "url": osm_url,
            "description": (
                f"**M{mag:.1f}** · {city}\n"
                f"> `{phases} phases`  ·  ID: `{evid}`"
            ),
            "color": color,
            "fields": [
                {"name": "Location",             "value": f"`{loc}`",          "inline": True},
                {"name": "Depth",                "value": f"`{dep} km`",       "inline": True},
                {"name": "Time Since",           "value": ago_str,             "inline": True},
                {"name": f"Origin ({pt_lbl})",   "value": f"`{pt_str}`",       "inline": False},
                {"name": "Origin (UTC)",          "value": f"`{ot} UTC`",       "inline": False},
                {"name": "Max MMI (epicenter)",   "value": f"`{ev_mmi_str}` — {ev_mmi_desc}", "inline": True},
                {"name": f"MMI at {c0['name'] if c0 else 'nearest city'}",
                                                  "value": f"`{ev_mmi_city}`",  "inline": True},
                *([ {"name": "🚨 Earthquake Safety",
                     "value": "**DROP** to the ground · take **COVER** under a sturdy table or against an interior wall · **HOLD ON** until shaking stops. Stay away from windows and heavy objects.",
                     "inline": False} ] if mag >= 3.0 else []),
                {"name": "​", "value": "This is a preliminary earthquake alert system that uses SeisComP compiled by Evan Li. For more information visit: https://myearthquake.dpdns.org/", "inline": False},
            ],
            "image": {"url": image_ref} if image_ref else {},
            "thumbnail": {"url": "https://dashboard.myearthquake.dpdns.org/static/drop_cover_hold.png"} if mag >= 3.0 else {},
            "footer": {"text": "SeisComP confirmed — USGS cross-check recommended",
                       "icon_url": "https://earthquake.usgs.gov/favicon.ico"},
            "timestamp": ot.replace(" ", "T") + "Z",
        }]
    }
    _broadcast_discord(payload, img, "map.png",
                       label=f"EVENT M{mag:.1f} {evid}")


def _send_prelim_email(ev):
    """Send preliminary detection email."""
    if not PRELIM_EMAIL_ENABLED:
        return
    mag = ev["mag_est"]
    subject = f"Preliminary Detection Est. M{mag:.1f} — {ev.get('city','Unknown')}"
    body = f"""<html><body style="font-family:monospace;background:#0d1117;color:#c9d1d9;padding:20px">
<h2 style="color:#d29922">⚠️ Preliminary Earthquake Detection</h2>
<p style="color:#8b949e;font-size:12px">This is a REAL-TIME STA/LTA detection — not yet confirmed by SeisComP.<br>
Check the dashboard for confirmation within 1-3 minutes.</p>
<table style="border-collapse:collapse;width:100%;max-width:500px;margin-top:12px">
  <tr><td style="color:#8b949e;padding:4px 8px">Est. Magnitude</td>
      <td style="padding:4px 8px;font-weight:bold;font-size:1.4em;color:#d29922">~M{mag:.1f}</td></tr>
  <tr><td style="color:#8b949e;padding:4px 8px">Est. Location</td>
      <td style="padding:4px 8px">{ev['lat']:.3f}°, {ev['lon']:.3f}°</td></tr>
  <tr><td style="color:#8b949e;padding:4px 8px">Nearest City</td>
      <td style="padding:4px 8px">{ev.get('city','—')}</td></tr>
  <tr><td style="color:#8b949e;padding:4px 8px">Detection Time</td>
      <td style="padding:4px 8px">{ev['ot_str']} UTC</td></tr>
  <tr><td style="color:#8b949e;padding:4px 8px">Stations Triggered</td>
      <td style="padding:4px 8px">{ev['n_stations']} ({ev.get('stations','—')})</td></tr>
  <tr><td style="color:#8b949e;padding:4px 8px">Max STA/LTA</td>
      <td style="padding:4px 8px">{ev['max_stalta']}</td></tr>
</table>
</body></html>"""
    msg = _MIME_MP("alternative")
    msg["Subject"] = subject
    msg["From"]    = PRELIM_EMAIL_FROM
    msg["To"]      = ", ".join(PRELIM_EMAIL_TO)
    msg.attach(_MIME_T(body, "html"))
    try:
        with _smtplib.SMTP(PRELIM_SMTP_HOST, PRELIM_SMTP_PORT, timeout=15) as s:
            s.ehlo(); s.starttls(); s.ehlo()
            s.login(PRELIM_SMTP_USER, PRELIM_SMTP_PASS)
            s.sendmail(PRELIM_EMAIL_FROM, PRELIM_EMAIL_TO, msg.as_string())
        print(f"[PRELIM-EMAIL] sent M{mag:.1f} detection")
    except Exception as e:
        print(f"[PRELIM-EMAIL] error: {e}")


def _preliminary_detector():
    """
    Background thread: monitors trigger_log for multi-station events.

    KEY CHANGE vs old version: spatial clustering first.
    With 1800+ global stations, a simple time-window grouping mixes CA + Australia
    + Oklahoma into one group with 30,000 km spread → always rejected.
    Instead we:
      1. For each trigger, find all others within CLUSTER_RADIUS_KM that also
         fired within PRELIM_WINDOW_SEC.
      2. That geographic cluster is the candidate event.
      3. Only then check min/max spread, magnitude, cooldown.
    """
    CLUSTER_RADIUS_KM = 600   # max distance between any two stations in the cluster

    sent_windows = {}   # location_key → time sent (for cooldown)

    while True:
        time.sleep(4)
        now = time.time()

        # Collect recent, quality-filtered, non-noisy triggers that have coords
        cutoff = now - PRELIM_WINDOW_SEC
        recent = [t for t in list(trigger_log)  # snapshot deque to avoid mutation during iteration
                  if t["ts"] >= cutoff
                  and t["stalta"] >= PRELIM_MIN_STALTA
                  and t.get("lat") is not None
                  and f"{t['net']}.{t['sta']}" not in _noisy_stations]

        if len(recent) < PRELIM_MIN_STATIONS:
            if recent:
                print(f"[PRELIM] {len(recent)} trigger(s) in window, need {PRELIM_MIN_STATIONS} — waiting: "
                      + ", ".join(f"{t['net']}.{t['sta']} STA/LTA={t['stalta']}" for t in recent[:5]),
                      flush=True)
            continue

        print(f"[PRELIM] evaluating {len(recent)} triggers in window", flush=True)

        # ── Spatial clustering ────────────────────────────────────────────────
        # For every trigger, collect all others within CLUSTER_RADIUS_KM.
        # Keep the cluster with the most unique stations.
        best_cluster = []
        for anchor in recent:
            # All triggers within the geographic radius
            nearby = [t for t in recent
                      if haversine(anchor["lat"], anchor["lon"],
                                   t["lat"], t["lon"]) <= CLUSTER_RADIUS_KM]
            # Deduplicate by station within this spatial cluster
            by_sta = {}
            for t in nearby:
                k = f"{t['net']}.{t['sta']}"
                if k not in by_sta or t["stalta"] > by_sta[k]["stalta"]:
                    by_sta[k] = t
            uniq = list(by_sta.values())
            if len(uniq) > len(best_cluster):
                best_cluster = uniq

        uniq = best_cluster
        if len(uniq) < PRELIM_MIN_STATIONS:
            continue

        # ── Spread calculation (pre-filter) ──────────────────────────────────
        lats_pre = [t["lat"] for t in uniq]
        lons_pre = [t["lon"] for t in uniq]
        lat_span_km_pre = (max(lats_pre) - min(lats_pre)) * 111.0
        lon_span_km_pre = (max(lons_pre) - min(lons_pre)) * 111.0 * math.cos(
                              math.radians(sum(lats_pre) / len(lats_pre)))
        spread_km_pre = math.sqrt(lat_span_km_pre**2 + lon_span_km_pre**2)

        # ── Travel-time grid search — use mantle VP for teleseisms ───────────
        use_vp = TELESEISM_VP if spread_km_pre >= TELESEISM_SPREAD_KM else 6.0
        lat_est, lon_est, origin_time_est, rms_sec, uniq = _grid_search_epicenter(uniq, vp=use_vp)
        n_sta = len(uniq)
        if n_sta < PRELIM_MIN_STATIONS:
            continue

        # ── Spread calculation (on filtered cluster) ──────────────────────────
        lats = [t["lat"] for t in uniq]
        lons = [t["lon"] for t in uniq]
        lat_span_km = (max(lats) - min(lats)) * 111.0
        lon_span_km = (max(lons) - min(lons)) * 111.0 * math.cos(
                          math.radians(sum(lats) / len(lats)))
        spread_km = math.sqrt(lat_span_km**2 + lon_span_km**2)

        max_stalta = max(t["stalta"] for t in uniq)

        # ── Cooldown keyed on rounded centroid (2° grid ≈ 220 km) ─────────────
        loc_key = (round(lat_est / 2), round(lon_est / 2))
        if loc_key in sent_windows and now - sent_windows[loc_key] < PRELIM_COOLDOWN_SEC:
            continue
        first_ts = origin_time_est  # use grid-search origin time
        # Only alert on fresh triggers (first station fired within 90 s)
        first_trigger_ts = min(t["ts"] for t in uniq)
        if now - first_trigger_ts > 120:
            continue

        sent_windows[loc_key] = now
        for k in [k for k, v in sent_windows.items() if now - v > 7200]:
            del sent_windows[k]

        # ── Improved magnitude estimate ───────────────────────────────────────
        # Combines max STA/LTA amplitude proxy + station count + distance-correction.
        # Calibrated for CA SeedLink BHZ channels:
        #   STA/LTA=3,  n=2  → ~M1.5
        #   STA/LTA=5,  n=4  → ~M2.0
        #   STA/LTA=10, n=6  → ~M2.8
        #   STA/LTA=25, n=8  → ~M3.5
        # Average station-to-epicenter distance correction (+0.15 per 100 km)
        avg_dist  = sum(haversine(lat_est, lon_est, t["lat"], t["lon"]) for t in uniq) / n_sta
        # Cap at 400 km — bad far-away locations must not blow up the magnitude
        dist_corr = 0.0015 * min(avg_dist, 400.0)  # max +0.60 correction
        mag_est = round(
            0.85 * math.log10(max(max_stalta, 1.01)) +
            0.35 * math.log10(max(n_sta, 1)) +
            0.85 + dist_corr, 1)
        mag_est = max(0.5, min(9.0, mag_est))

        # ── Teleseism classification + gates ─────────────────────────────────
        is_teleseism = spread_km >= TELESEISM_SPREAD_KM
        if is_teleseism and n_sta < TELESEISM_MIN_STATIONS:
            print(f"[PRELIM] teleseism suppressed: only {n_sta}/{TELESEISM_MIN_STATIONS} stations, spread={spread_km:.0f} km", flush=True)
            continue
        if is_teleseism and mag_est < TELESEISM_MIN_MAG:
            print(f"[PRELIM] teleseism suppressed: M{mag_est} < {TELESEISM_MIN_MAG} (teleseisms must be M{TELESEISM_MIN_MAG}+), spread={spread_km:.0f} km", flush=True)
            continue

        # ── RMS quality gates ────────────────────────────────────────────────
        # Local: high RMS with few stations = bad location
        if not is_teleseism and rms_sec > 4.0 and n_sta < 5:
            print(f"[PRELIM] poor local location rejected: RMS={rms_sec:.1f}s n={n_sta}", flush=True)
            continue
        # Teleseism: must have tight travel-time consistency with mantle VP
        if is_teleseism and rms_sec > TELESEISM_MAX_RMS:
            print(f"[PRELIM] teleseism suppressed: RMS={rms_sec:.1f}s > {TELESEISM_MAX_RMS}s (incoherent arrivals)", flush=True)
            continue

        # ── Geographic filter — local detections must be in/near California ────
        if not is_teleseism:
            if not (PRELIM_LOCAL_LAT_MIN <= lat_est <= PRELIM_LOCAL_LAT_MAX and
                    PRELIM_LOCAL_LON_MIN <= lon_est <= PRELIM_LOCAL_LON_MAX):
                print(f"[PRELIM] local detection outside California box rejected: "
                      f"lat={lat_est:.2f} lon={lon_est:.2f}", flush=True)
                continue

        # ── MMI at epicenter ───────────────────────────────────────────────────
        mmi_val  = _estimate_mmi(mag_est, depth_km=10, epi_dist_km=0)
        mmi_str  = _mmi_roman(mmi_val)
        mmi_desc = _mmi_desc(mmi_val)
        # MMI at nearest city
        nc = nearest_cities(lat_est, lon_est, n=1)
        city_str = f"{nc[0]['name']}, {nc[0]['state']} ({nc[0]['km']} km)" if nc else ""
        city_dist = nc[0]['km'] if nc else 50
        mmi_city_val = _estimate_mmi(mag_est, depth_km=10, epi_dist_km=city_dist)
        mmi_city_str = _mmi_roman(mmi_city_val)

        ot_str   = datetime.fromtimestamp(first_trigger_ts, tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
        sta_list = sorted(uniq, key=lambda x: x["ts"])
        stations_str = ", ".join(f"{t['net']}.{t['sta']}" for t in sta_list[:6])
        if n_sta > 6:
            stations_str += f" +{n_sta-6}"

        ev_info = {
            "ts":           first_trigger_ts,
            "lat":          round(lat_est, 3),
            "lon":          round(lon_est, 3),
            "depth_est":    10,
            "mag_est":      mag_est,
            "n_stations":   n_sta,
            "max_stalta":   round(max_stalta, 2),
            "rms_sec":      round(rms_sec, 2),
            "ot_str":       ot_str,
            "city":         city_str,
            "stations":     stations_str,
            "spread_km":    round(spread_km, 0),
            "teleseism":    is_teleseism,
            "mmi":          mmi_val,
            "mmi_str":      mmi_str,
            "mmi_desc":     mmi_desc,
            "mmi_city":     mmi_city_val,
            "mmi_city_str": mmi_city_str,
        }
        preliminary_events.appendleft(ev_info)
        _save_prelim_log()   # persist immediately to disk
        sio.emit("preliminary", ev_info)
        print(f"[PRELIM] Est. M{mag_est} lat={lat_est:.2f} lon={lon_est:.2f} "
              f"n={n_sta} max_stalta={max_stalta:.1f}")
        _send_prelim_discord(ev_info)
        _send_prelim_email(ev_info)


# ── SeedLink ───────────────────────────────────────────────────────────────────
def _seedlink_worker(streams_subset, worker_id):
    """
    Subscribe a small subset of streams individually.
    Keeping each connection to ≤20 streams prevents the server from
    resetting the connection ('Connection reset by peer').
    """
    while True:
        try:
            client = EasySeedLinkClient(SEEDLINK, autoconnect=False)

            def on_data(trace, _states=states):
                key = f"{trace.stats.network}.{trace.stats.station}.{trace.stats.channel}"
                if key in _states:
                    _states[key].feed(trace)
            client.on_data = on_data
            client.on_seedlink_error = lambda: None

            # Subscribe streams — suppress obspy's "station not accepted" stdout noise
            import sys as _sys, io as _io
            _devnull = _io.StringIO()
            _old_stdout = _sys.stdout
            _sys.stdout = _devnull
            try:
                for net, sta, chan in streams_subset:
                    client.select_stream(net, sta, chan)
            finally:
                _sys.stdout = _old_stdout

            client.connect()
            print(f"[SeedLink worker-{worker_id}] connected, {len(streams_subset)} streams", flush=True)
            client.run()
        except Exception as e:
            print(f"[SeedLink worker-{worker_id}] error: {e}, reconnecting in 5s")
            time.sleep(5)


def _fdsn_fetch_one(net, sta, chan, t_start, t_end):
    """
    Fetch one station's waveform from appropriate FDSNWS endpoint.
    BK → NCEDC, NN → IRIS, everything else → IRIS.
    Returns ObsPy Stream or None.
    """
    try:
        import urllib.request, io
        from obspy import read as obspy_read
        if net in ("CI", "AZ", "SB", "SC", "SB", "TS"):
            base = "https://service.scedc.caltech.edu/fdsnws/dataselect/1/query"
        elif net in ("BK", "NC"):
            base = "https://service.ncedc.org/fdsnws/dataselect/1/query"
        else:
            base = "https://service.earthscope.org/fdsnws/dataselect/1/query"
        url = (
            f"{base}?net={net}&sta={sta}&loc=*&cha={chan}"
            f"&start={t_start.strftime('%Y-%m-%dT%H:%M:%S')}"
            f"&end={t_end.strftime('%Y-%m-%dT%H:%M:%S')}"
            f"&format=miniseed"
        )
        req = urllib.request.Request(url, headers={"User-Agent": "ObsPy/dashboard"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            raw = resp.read()
        if len(raw) < 64:
            return None
        return obspy_read(io.BytesIO(raw))
    except Exception:
        return None

# Keep old name as alias for any remaining references
_ncedc_fetch_one = _fdsn_fetch_one


# Dynamic state registry — grows as new stations are discovered
_state_lock = threading.Lock()

def _ensure_state(net, sta, chan):
    """Return (or create) a StationState for this net/sta/chan triple."""
    key = f"{net}.{sta}.{chan}"
    s = states.get(key)
    if s is None:
        with _state_lock:
            s = states.get(key)  # double-check inside lock
            if s is None:
                s = StationState(net, sta, chan)
                states[key] = s
                # also populate city if coords already known
                coord = _sta_coords.get(f"{net}.{sta}", {})
                if coord:
                    nc = nearest_cities(coord["lat"], coord["lon"], n=1)
                    if nc:
                        c = nc[0]
                        s.city = f"{c['name']}, {c['state']} ({c['km']} km)"
    return s


def _ncedc_poller():
    """
    Poll FDSN endpoints using sta=* wildcard GET requests — one request per
    (network, channel) combination every 60s.  Automatically discovers every
    active station in every network; dynamically creates StationState objects
    for new ones.
    """
    import urllib.request, io
    from obspy import read as obspy_read
    from datetime import timedelta, datetime as _dt2

    ENDPOINTS = {
        "scedc": "https://service.scedc.caltech.edu/fdsnws/dataselect/1/query",
        "ncedc": "https://service.ncedc.org/fdsnws/dataselect/1/query",
        "iris":  "https://service.earthscope.org/fdsnws/dataselect/1/query",
    }

    def _fetch_network(endpoint_url, net, chan, t_start, t_end):
        """
        Wildcard GET request for an entire network.
        Returns ObsPy Stream or None.  Retries once with a shorter window on timeout.
        """
        ts = t_start.strftime("%Y-%m-%dT%H:%M:%S")
        te = t_end.strftime("%Y-%m-%dT%H:%M:%S")
        url = (f"{endpoint_url}?net={net}&sta=*&loc=*&cha={chan}"
               f"&start={ts}&end={te}&format=miniseed&nodata=404")
        for attempt in range(2):
            try:
                req = urllib.request.Request(url,
                    headers={"User-Agent": "SeisComp-Dashboard/3.0"})
                with urllib.request.urlopen(req, timeout=90) as r:
                    raw = r.read()
                if len(raw) < 64:
                    return None
                return obspy_read(io.BytesIO(raw))
            except Exception as _fe:
                if attempt == 0:
                    # Retry with a shorter 60s window — smaller payload, less likely to timeout
                    shorter = t_end - timedelta(seconds=60)
                    ts = shorter.strftime("%Y-%m-%dT%H:%M:%S")
                    url = (f"{endpoint_url}?net={net}&sta=*&loc=*&cha={chan}"
                           f"&start={ts}&end={te}&format=miniseed&nodata=404")
        return None

    def _feed_stream(st, net):
        """Feed traces into matching (or newly created) states."""
        if not st:
            return 0
        fed = 0
        for tr in st:
            state = _ensure_state(tr.stats.network, tr.stats.station,
                                  tr.stats.channel)
            tr_end_ts = tr.stats.endtime.timestamp
            last_t = state.last_polled_t
            if last_t is not None and tr_end_ts <= last_t:
                continue
            if last_t is not None and tr.stats.starttime.timestamp < last_t:
                try:
                    trimmed = tr.slice(starttime=UTCDateTime(last_t))
                    if trimmed is None or trimmed.stats.npts == 0:
                        continue
                    tr = trimmed
                except Exception:
                    continue
            state.feed(tr)
            state.last_polled_t = tr_end_ts
            fed += 1
        return fed

    print("[FDSN] wildcard network poller started")

    # Build task list: one entry per (net, primary_channel, endpoint)
    tasks = []
    for net, (ep_key, chan_list) in FDSN_NETWORKS.items():
        primary_chan = chan_list.split(",")[0].strip()
        tasks.append((net, primary_chan, ENDPOINTS[ep_key], ep_key))

    n_tasks = len(tasks)
    print(f"[FDSN] {n_tasks} networks — each on its own 60s cycle, staggered across time")

    # STAGGERED POLLING: each network gets its own independent loop, offset by
    # (index / n_tasks) * 60 seconds so updates are spread throughout every minute
    # instead of all firing at once.
    POLL_INTERVAL = 60   # seconds between polls for each network
    FETCH_WINDOW  = 120  # seconds of data to request (2× LTA_WIN)

    def _network_loop(net, chan, url, ep_key, initial_delay):
        """Runs forever; polls one network every POLL_INTERVAL seconds."""
        time.sleep(initial_delay)
        while True:
            try:
                t_end   = _dt2.now(timezone.utc)
                t_start = t_end - timedelta(seconds=FETCH_WINDOW)
                st = _fetch_network(url, net, chan, t_start, t_end)
                n  = _feed_stream(st, net)
                if n:
                    pass  # success — stations' last_update timestamps now differ
            except Exception as _e:
                print(f"[FDSN] {net} loop error: {_e}", flush=True)
            time.sleep(POLL_INTERVAL)

    # Launch one daemon thread per network, staggered evenly across POLL_INTERVAL
    for i, (net, chan, url, ep_key) in enumerate(tasks):
        delay = i * (POLL_INTERVAL / n_tasks)  # e.g. 50 networks → 1.2s apart
        t = threading.Thread(target=_network_loop,
                             args=(net, chan, url, ep_key, delay),
                             daemon=True, name=f"fdsn-{net}")
        t.start()

    # This function can now return — the threads are running independently
    return


def _query_seedlink_inventory(host="localhost:18000"):
    """
    Ask the local SeedLink what streams it currently has.
    Returns a set of (net, sta) tuples.
    """
    import subprocess
    # Use full path — launchd has a minimal PATH that may not include seiscomp/bin
    slinktool_bin = "/Users/OuOu/seiscomp/bin/slinktool"
    try:
        r = subprocess.run(
            [slinktool_bin, "-Q", host],
            capture_output=True, text=True, timeout=10
        )
        available = set()
        for line in r.stdout.splitlines():
            parts = line.split()
            if len(parts) >= 2:
                available.add((parts[0].strip(), parts[1].strip()))
        return available
    except Exception as e:
        print(f"[SeedLink] inventory query failed: {e} — subscribing to all configured streams", flush=True)
        return None   # None = don't filter


def start_seedlink():
    """
    Query the local SeedLink inventory first, then only subscribe to streams
    that are actually available.  This prevents 'station not accepted' spam and
    worker reconnect storms caused by requesting streams chain_plugin hasn't
    pulled yet.
    """
    CHUNK = 20

    # Find out which (net, sta) pairs are live on our local SeedLink
    available = _query_seedlink_inventory()
    if available is not None:
        live_streams = [(net, sta, chan) for net, sta, chan in STREAMS
                        if (net, sta) in available]
        skipped = len(STREAMS) - len(live_streams)
        print(f"[SeedLink] {len(live_streams)} streams available on server "
              f"({skipped} configured but not yet live — skipping)", flush=True)
    else:
        live_streams = list(STREAMS)

    # Keep track of which streams already have a worker so re-checks don't duplicate
    _started_streams = set((net, sta) for net, sta, chan in live_streams)

    chunks = [live_streams[i:i+CHUNK] for i in range(0, len(live_streams), CHUNK)]
    for i, chunk in enumerate(chunks):
        t = threading.Thread(target=_seedlink_worker, args=(chunk, i), daemon=True)
        t.start()
    print(f"[SeedLink] started {len(chunks)} local workers for {len(live_streams)} streams", flush=True)

    def _recheck_seedlink():
        """Every 5 minutes, pick up newly-available streams that chain_plugin added."""
        worker_id = [len(chunks)]
        while True:
            time.sleep(300)
            available2 = _query_seedlink_inventory()
            if not available2:
                continue
            new_streams = [(net, sta, chan) for net, sta, chan in STREAMS
                           if (net, sta) in available2 and (net, sta) not in _started_streams]
            if not new_streams:
                continue
            print(f"[SeedLink] {len(new_streams)} newly-available stream(s) — adding workers", flush=True)
            new_chunks = [new_streams[i:i+CHUNK] for i in range(0, len(new_streams), CHUNK)]
            for chunk in new_chunks:
                t2 = threading.Thread(target=_seedlink_worker,
                                      args=(chunk, worker_id[0]), daemon=True)
                t2.start()
                worker_id[0] += 1
            for net, sta, chan in new_streams:
                _started_streams.add((net, sta))

    threading.Thread(target=_recheck_seedlink, daemon=True).start()

    # FDSNWS poller: BK Bay Area via NCEDC, NN Nevada via IRIS
    t = threading.Thread(target=_ncedc_poller, daemon=True)
    t.start()
    print(f"[FDSN] started wildcard poller for {len(FDSN_NETWORKS)} networks")

    # Preliminary real-time detector — wrapped in a supervisor so it auto-restarts on any crash
    def _prelim_supervisor():
        while True:
            try:
                _preliminary_detector()
            except Exception as e:
                print(f"[PRELIM] thread crashed: {e} — restarting in 5s", flush=True)
                time.sleep(5)

    t = threading.Thread(target=_prelim_supervisor, daemon=True, name="prelim-supervisor")
    t.start()
    print("[PRELIM] Preliminary detector started")

# ── DB ─────────────────────────────────────────────────────────────────────────
def db():
    return mysql.connector.connect(**DB_CFG)

def get_events(limit=50):
    try:
        cn = db(); cur = cn.cursor(dictionary=True)
        cur.execute("""
            SELECT pe2.publicID as evid, o.time_value as origin_time,
                   ROUND(o.latitude_value,3) as lat, ROUND(o.longitude_value,3) as lon,
                   ROUND(o.depth_value,1) as depth, o.quality_usedPhaseCount as phases,
                   ROUND(m.magnitude_value,1) as mag, m.type as magtype
            FROM   Event e
            JOIN   PublicObject pe2 ON pe2._oid=e._oid
            JOIN   PublicObject pe  ON pe.publicID=e.preferredOriginID
            JOIN   Origin o         ON o._oid=pe._oid
            LEFT   JOIN PublicObject pm ON pm.publicID=e.preferredMagnitudeID
            LEFT   JOIN Magnitude m     ON m._oid=pm._oid
            ORDER  BY o.time_value DESC LIMIT %s""", (limit,))
        rows = cur.fetchall(); cur.close(); cn.close()
        return [dict(r, origin_time=str(r["origin_time"])) for r in rows]
    except: return []

def get_event_detail(evid):
    try:
        cn = db(); cur = cn.cursor(dictionary=True)
        cur.execute("""
            SELECT pe2.publicID as evid, o.time_value as origin_time,
                   o.latitude_value as lat, o.longitude_value as lon,
                   o.depth_value as depth, o.quality_usedPhaseCount as phases,
                   o.quality_usedStationCount as stations, o.quality_standardError as rms,
                   o.quality_azimuthalGap as az_gap, o.quality_minimumDistance as min_dist,
                   o.evaluationMode as eval_mode, e.preferredOriginID as origin_id
            FROM   Event e
            JOIN   PublicObject pe2 ON pe2._oid=e._oid
            JOIN   PublicObject pe  ON pe.publicID=e.preferredOriginID
            JOIN   Origin o         ON o._oid=pe._oid
            WHERE  pe2.publicID=%s""", (evid,))
        ev = cur.fetchone()
        if not ev: return None
        ev = dict(ev, origin_time=str(ev["origin_time"]))

        cur.execute("""
            SELECT ROUND(m.magnitude_value,2) as mag, m.type, m.stationCount
            FROM   Magnitude m JOIN PublicObject pm ON pm._oid=m._oid
            WHERE  m._parent_oid=(
                SELECT o._oid FROM Origin o JOIN PublicObject pe ON pe._oid=o._oid
                WHERE pe.publicID=%s) ORDER BY m.magnitude_value DESC""", (ev["origin_id"],))
        ev["magnitudes"] = cur.fetchall()

        cur.execute("""
            SELECT p.waveformID_networkCode as net, p.waveformID_stationCode as sta,
                   p.waveformID_channelCode as cha, p.phaseHint_code as phase,
                   p.time_value as pick_time, a.timeResidual as residual,
                   a.distance as distance, a.azimuth as azimuth, a.weight as weight,
                   s.latitude as sta_lat, s.longitude as sta_lon
            FROM   Arrival a
            JOIN   Origin o ON o._oid=a._parent_oid
            JOIN   PublicObject po ON po._oid=o._oid
            JOIN   PublicObject pp ON pp.publicID=a.pickID
            JOIN   Pick p ON p._oid=pp._oid
            LEFT   JOIN Station s ON s.code=p.waveformID_stationCode
            WHERE  po.publicID=%s ORDER BY a.distance""", (ev["origin_id"],))
        ev["picks"] = [dict(p, pick_time=str(p["pick_time"]),
                            sta_lat=float(p["sta_lat"]) if p["sta_lat"] else None,
                            sta_lon=float(p["sta_lon"]) if p["sta_lon"] else None)
                       for p in cur.fetchall()]

        cur.close(); cn.close()
        # Nearest cities
        try:
            ev["cities"] = nearest_cities(float(ev["lat"]), float(ev["lon"]))
        except: ev["cities"] = []
        return ev
    except Exception as e:
        return {"error": str(e)}

# ── Focal mechanism ────────────────────────────────────────────────────────────
def _get_trace_for_pick(sds, p, pre=1.5, post=2.5):
    """Try SDS archive → live ring buffer → FDSN (SCEDC/IRIS)."""
    from obspy import Trace
    pt = UTCDateTime(p["pick_time"])
    # 1) SDS archive
    try:
        st = sds.get_waveforms(p["net"], p["sta"], "*", p["cha"], pt - pre, pt + post)
        if st: return st[0]
    except Exception:
        pass
    # 2) Live ring buffer
    for suffix in [p["cha"], p["cha"][:2] + "Z"]:
        key = f"{p['net']}.{p['sta']}.{suffix}"
        state = states.get(key)
        if not state: continue
        with state.lock:
            if not state.buf or not state.srate: continue
            pt_ts = pt.timestamp
            seg = [(t, s) for t, s in state.buf if pt_ts - pre <= t <= pt_ts + post]
            if len(seg) < 10: continue
            ts_arr = [x[0] for x in seg]
            vs_arr = np.array([x[1] for x in seg], dtype=float)
            tr = Trace(data=vs_arr)
            tr.stats.sampling_rate = state.srate
            tr.stats.starttime = UTCDateTime(ts_arr[0])
            tr.stats.network = p["net"]; tr.stats.station = p["sta"]
            tr.stats.channel = suffix
            return tr
    # 3) FDSN web service (SCEDC for CI/AZ, IRIS otherwise)
    try:
        cl = _fdsn(p["net"])
        if cl:
            st = cl.get_waveforms(p["net"], p["sta"], "*", p["cha"], pt - pre, pt + post)
            if st: return st[0]
    except Exception:
        pass
    return None

def compute_beachball(evid):
    """Auto-compute focal mechanism from P-wave first motions (archive or live buffer)."""
    ev = get_event_detail(evid)
    if not ev or "error" in ev: return None

    picks = [p for p in ev.get("picks", []) if p.get("phase","").startswith("P")
             and p.get("azimuth") is not None and p.get("distance") is not None]
    if not picks: return None

    sds = SDSClient(SDS_PATH)
    polarities = []

    for p in picks:
        try:
            tr = _get_trace_for_pick(sds, p)
            if tr is None: continue
            tr.detrend("demean")
            tr.filter("highpass", freq=2.0)

            pt  = UTCDateTime(p["pick_time"])
            sr  = tr.stats.sampling_rate
            idx = int((pt - tr.stats.starttime) * sr)
            n_pre  = int(0.4 * sr)
            n_post = int(0.4 * sr)
            if idx < n_pre or idx + n_post > len(tr.data): continue

            pre  = tr.data[idx - n_pre : idx]
            post = tr.data[idx : idx + n_post]
            noise = np.std(pre) if len(pre) else 1.0
            if noise == 0: continue

            peak_post = post[np.argmax(np.abs(post))] if len(post) else 0
            # Relaxed SNR threshold: 2.0× (was 2.5×)
            if abs(peak_post) < 2.0 * noise: continue

            polarity = 1 if peak_post > 0 else -1
            dist_deg = float(p["distance"])
            takeoff = max(20.0, 90.0 - dist_deg * 15.0)
            azimuth = float(p["azimuth"])

            polarities.append((azimuth, takeoff, polarity, f"{p['net']}.{p['sta']}"))
        except Exception:
            continue

    # Relax minimum polarity count to 3
    if len(polarities) < 3:
        return None

    # Grid search: strike 0-360, dip 0-90, rake -180-180 in 5° steps
    azimuths  = np.array([p[0] for p in polarities])
    takeoffs  = np.array([p[1] for p in polarities])
    obs_pols  = np.array([p[2] for p in polarities])

    az_r  = np.radians(azimuths)
    tk_r  = np.radians(takeoffs)

    best_score, best_sdr = -1, (0, 45, 0)

    for strike in range(0, 360, 5):
        for dip in range(0, 91, 5):
            for rake in range(-180, 181, 10):
                sr_ = math.radians(strike)
                dr  = math.radians(dip)
                rr  = math.radians(rake)
                # Aki & Richards radiation pattern
                cos_tk = np.cos(tk_r); sin_tk = np.sin(tk_r)
                da = az_r - sr_
                P = ( np.cos(rr)*np.sin(dr)*np.sin(tk_r)**2*np.sin(2*da)
                    - np.cos(rr)*np.cos(dr)*np.sin(2*tk_r)*np.cos(da)
                    + np.sin(rr)*np.sin(2*dr)*(cos_tk**2 - sin_tk**2*np.sin(da)**2)
                    + np.sin(rr)*np.cos(2*dr)*np.sin(2*tk_r)*np.sin(da) )
                pred = np.sign(P)
                pred[pred == 0] = 1
                score = int(np.sum(pred == obs_pols))
                if score > best_score:
                    best_score = score; best_sdr = (strike, dip, rake)

    if best_score < len(polarities) * 0.6:
        return None  # Poor fit

    strike, dip, rake = best_sdr
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from obspy.imaging.beachball import beachball as bb
    fig = plt.figure(figsize=(2.5, 2.5), facecolor="#0d1117")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor("#0d1117"); ax.set_aspect("equal"); ax.axis("off")
    bb([strike, dip, rake], linewidth=1.2,
       facecolor="#c0392b", bgcolor="#0d1117", edgecolor="white",
       width=220, fig=fig)
    ax.set_title(f"s{strike} d{dip} r{rake}", color="#8b949e",
                 fontsize=7, pad=2, fontfamily="monospace")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight",
                facecolor="#0d1117")
    buf.seek(0); plt.close("all")
    return buf

# ── Archive waveform helpers ───────────────────────────────────────────────────
def get_event_waveform_data(evid):
    ev = get_event_detail(evid)
    if not ev or "error" in ev: return []
    picks = [p for p in ev.get("picks", []) if p.get("pick_time")]
    if not picks: return []

    ot  = UTCDateTime(ev["origin_time"])
    sds = SDSClient(SDS_PATH)
    traces = []
    seen_sta = set()

    for p in picks[:20]:  # try up to 20 picks, return first 16 with data
        if len(traces) >= 16: break
        sta_key = f"{p['net']}.{p['sta']}"
        if sta_key in seen_sta: continue
        seen_sta.add(sta_key)
        try:
            pt = UTCDateTime(p["pick_time"])
            t0 = pt - 15; t1 = pt + 45

            # 1) Local SDS archive
            tr = None
            try:
                st = sds.get_waveforms(p["net"], p["sta"], "*", p["cha"], t0, t1)
                if st: tr = st[0]
            except Exception: pass

            # 2) FDSN fallback (SCEDC or IRIS)
            if tr is None:
                try:
                    cl = _fdsn(p["net"])
                    if cl:
                        st = cl.get_waveforms(p["net"], p["sta"], "*", p["cha"], t0, t1)
                        if st: tr = st[0]
                except Exception: pass

            if tr is None: continue
            tr.detrend("demean"); tr.filter("bandpass", freqmin=1.0, freqmax=10.0)
            data = tr.data.tolist()
            peak = max(abs(v) for v in data) or 1.0
            data = [v / peak for v in data]
            if len(data) > 600:
                step = len(data) / 600
                data = [data[int(i * step)] for i in range(600)]
            traces.append({
                "label":    sta_key,
                "chan":     p["cha"],
                "phase":    p.get("phase", "P"),
                "color":    NET_COLORS.get(p["net"], "#aaa"),
                "pick_rel": 15.0,
                "origin_rel": float(pt - ot) + 15.0,
                "dist":     round(float(p["distance"] or 0), 1),
                "data":     data,
                "source":   "local" if tr.stats.location is not None else "fdsn",
            })
        except Exception:
            continue

    return sorted(traces, key=lambda x: x["dist"])

# ── Flask ──────────────────────────────────────────────────────────────────────
app = Flask(__name__)
app.config["SECRET_KEY"] = "seiscomp"
CORS(app, origins="*")   # Allow GitHub Pages and any other origin to call the API
sio = SocketIO(app, cors_allowed_origins="*", async_mode="threading",
               logger=False, engineio_logger=False,
               ping_timeout=120,    # wait 2 min before declaring client dead
               ping_interval=30,    # heartbeat every 30s (browser throttles to ~1/min in bg)
               max_http_buffer_size=20_000_000)  # 20 MB — needed for thousands of station rows

@app.route("/")
def index(): return render_template_string(HTML_MAIN)

@app.route("/station/<net>/<sta>/<chan>")
def station_page(net, sta, chan):
    key = f"{net}.{sta}.{chan}"
    color = NET_COLORS.get(net, "#aaa")
    return render_template_string(HTML_STATION, net=net, sta=sta, chan=chan,
                                  key=key, color=color)

@app.route("/event/<evid>/waveforms")
def event_wf_page(evid): return render_template_string(HTML_EVENT_WF, evid=evid)

@app.route("/api/event/<evid>")
def api_event(evid): return jsonify(get_event_detail(evid))

@app.route("/api/event/<evid>/waveforms")
def api_event_waveforms(evid): return jsonify(get_event_waveform_data(evid))

def _no_data_beachball_png():
    """Return a placeholder PNG when no waveform data is available."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(2,2), facecolor="#0d1117")
    ax.set_facecolor("#161b22")
    ax.set_aspect("equal"); ax.set_xlim(-1,1); ax.set_ylim(-1,1)
    c = plt.Circle((0,0), 0.85, color="#2a3344", zorder=1)
    ax.add_patch(c)
    ax.text(0, 0.15, "No waveform", color="#8b949e", ha="center", va="center",
            fontsize=7, fontfamily="monospace")
    ax.text(0, -0.15, "data available", color="#8b949e", ha="center", va="center",
            fontsize=7, fontfamily="monospace")
    ax.axis("off")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight",
                facecolor="#0d1117", transparent=False)
    buf.seek(0); plt.close(fig)
    return buf

@app.route("/api/beachball/<evid>.png")
def api_beachball(evid):
    buf = compute_beachball(evid)
    if not buf: buf = _no_data_beachball_png()
    return send_file(buf, mimetype="image/png")

import subprocess, shlex

SEISCOMP = "/Users/OuOu/seiscomp/bin/seiscomp"

def _scvoice_running():
    r = subprocess.run([SEISCOMP, "status", "scvoice"],
                       capture_output=True, text=True)
    return "is running" in r.stdout

@app.route("/api/scvoice", methods=["GET"])
def api_scvoice_status():
    return jsonify({"running": _scvoice_running()})

@app.route("/api/scvoice/toggle", methods=["POST"])
def api_scvoice_toggle():
    running = _scvoice_running()
    cmd = "stop" if running else "start"
    subprocess.run([SEISCOMP, cmd, "scvoice"], capture_output=True)
    import time; time.sleep(0.8)
    now = _scvoice_running()
    return jsonify({"running": now})

@app.route("/api/live/<path:key>")
def api_live(key):
    s = states.get(key)
    if not s: abort(404)
    pts, t_end = s.waveform_pts(n_pts=2000, secs=BUF_SECS)
    return jsonify({"values": pts, "t_end": t_end, "secs": BUF_SECS,
                    "stalta": round(s.stalta, 2),
                    "updated": s.last_update.strftime("%H:%M:%S UTC") if s.last_update else None})

@app.route("/api/live_long/<path:key>")
def api_live_long(key):
    s = states.get(key)
    if not s: abort(404)
    pts, t_end = s.waveform_pts_long(n_pts=1440, secs=BUF_LONG)
    return jsonify({"values": pts, "t_end": t_end, "secs": BUF_LONG,
                    "stalta": round(s.stalta, 2),
                    "updated": s.last_update.strftime("%H:%M:%S UTC") if s.last_update else None})

def _sds_peak_stalta(net, sta, chan, t_start_utc, t_end_utc):
    """
    Read waveform from SDS archive and return (peak_stalta, peak_ts).
    Returns (None, None) if no data available.
    """
    try:
        sds = SDSClient(SDS_PATH)
        st = sds.get_waveforms(net, sta, "", chan,
                               UTCDateTime(t_start_utc.strftime("%Y-%m-%dT%H:%M:%S")),
                               UTCDateTime(t_end_utc.strftime("%Y-%m-%dT%H:%M:%S")))
        if not st: return None, None
        tr = st[0]
        sr = tr.stats.sampling_rate
        n_sta = max(1, int(STA_WIN * sr))
        n_lta = max(1, int(LTA_WIN * sr))
        data = tr.data.astype(float)
        if len(data) < n_lta + n_sta: return None, None
        data2 = data ** 2
        peak_v, peak_i = 0.0, 0
        for i in range(n_lta + n_sta, len(data2)):
            sv = np.mean(data2[i - n_sta:i])
            lv = np.mean(data2[i - n_lta - n_sta:i - n_sta])
            ratio = math.sqrt(sv / lv) if lv > 0 else 0.0
            if ratio > peak_v:
                peak_v, peak_i = ratio, i
        peak_ts = tr.stats.starttime.timestamp + peak_i / sr
        return round(peak_v, 2), peak_ts
    except Exception:
        return None, None


@app.route("/api/event_stalta/<evid>")
def api_event_stalta(evid):
    """
    Return peak STA/LTA per station during the ~10 min window after an event.
    Uses live stalta_hist if available (event happened while dashboard was running),
    otherwise falls back to reading from the SDS archive.
    """
    try:
        cn = db(); cur = cn.cursor(dictionary=True)
        # Get origin time AND list of picked stations in one trip
        cur.execute("""
            SELECT o.time_value as origin_time
            FROM Event e
            JOIN PublicObject pe2 ON pe2._oid=e._oid
            JOIN PublicObject pe ON pe.publicID=e.preferredOriginID
            JOIN Origin o ON o._oid=pe._oid
            WHERE pe2.publicID=%s
        """, (evid,))
        row = cur.fetchone()
        if not row: cur.close(); cn.close(); return jsonify({"error": "event not found"})
        ot = row["origin_time"]
        # Fetch picked stations to limit SDS reads
        cur.execute("""
            SELECT DISTINCT p.waveformID_networkCode as net,
                            p.waveformID_stationCode as sta,
                            p.waveformID_channelCode as cha
            FROM Event e
            JOIN PublicObject pe2 ON pe2._oid=e._oid
            JOIN PublicObject pe  ON pe.publicID=e.preferredOriginID
            JOIN Origin o         ON o._oid=pe._oid
            JOIN Arrival a        ON a._parent_oid=o._oid
            JOIN PublicObject pp  ON pp.publicID=a.pickID
            JOIN Pick p           ON p._oid=pp._oid
            WHERE pe2.publicID=%s
        """, (evid,))
        pick_stations = {(r["net"], r["sta"], r["cha"]) for r in cur.fetchall()}
        cur.close(); cn.close()
        # MySQL returns naive datetime — SeisComP always stores UTC, so pin it
        if hasattr(ot, "replace"):
            ot_utc = ot.replace(tzinfo=timezone.utc)
        else:
            from datetime import datetime as _dt
            ot_utc = _dt.strptime(str(ot)[:19], "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
        origin_ts = ot_utc.timestamp()
    except Exception as e:
        return jsonify({"error": str(e)})

    window_end  = origin_ts + 600  # 10-minute post-origin window
    t_start_utc = ot_utc
    t_end_utc   = datetime.fromtimestamp(window_end, tz=timezone.utc)

    # Check whether live history covers this event
    dash_start = min((s.stalta_hist[0][0] for s in states.values()
                      if s.stalta_hist), default=time.time())
    use_sds = origin_ts < dash_start  # event predates current dashboard run

    results = []

    if use_sds:
        # ── Retroactive: read SDS for picked stations + any triggered monitored stations
        # Prefer pick list; fallback to all monitored if picks unknown
        target_keys = {}
        for key, s in states.items():
            if (s.net, s.sta, s.chan) in pick_stations or \
               any(s.sta == ps for _, ps, _ in pick_stations):
                target_keys[key] = s
        if not target_keys:  # no overlap → try all monitored
            target_keys = dict(states)

        from concurrent.futures import ThreadPoolExecutor, as_completed
        def _worker(key, s):
            pv, pt = _sds_peak_stalta(s.net, s.sta, s.chan, t_start_utc, t_end_utc)
            return key, s, pv, pt

        with ThreadPoolExecutor(max_workers=8) as ex:
            futures = {ex.submit(_worker, k, s): k for k, s in target_keys.items()}
            for fut in as_completed(futures):
                key, s, peak_v, peak_t = fut.result()
                if peak_v is None or peak_v < 0.1: continue
                results.append({
                    "key": key, "net": s.net, "sta": s.sta, "chan": s.chan,
                    "peak_stalta":   round(peak_v, 2),
                    "peak_time_utc": datetime.fromtimestamp(peak_t, tz=timezone.utc).strftime("%H:%M:%S"),
                    "delay_s":       round(peak_t - origin_ts),
                    "color":         NET_COLORS.get(s.net, "#aaa"),
                    "source":        "archive",
                    "picked":        (s.net, s.sta, s.chan) in pick_stations,
                })
    else:
        # ── Live: use in-memory history ─────────────────────────────────────
        for key, s in states.items():
            hist   = list(s.stalta_hist)
            window = [(t, v) for t, v in hist if origin_ts - 5 <= t <= window_end]
            if not window: continue
            peak_v = max(v for _, v in window)
            if peak_v < 0.1: continue
            peak_t = next(t for t, v in window if v == peak_v)
            results.append({
                "key": key, "net": s.net, "sta": s.sta, "chan": s.chan,
                "peak_stalta":   round(peak_v, 2),
                "peak_time_utc": datetime.fromtimestamp(peak_t, tz=timezone.utc).strftime("%H:%M:%S"),
                "delay_s":       round(peak_t - origin_ts),
                "color":         NET_COLORS.get(s.net, "#aaa"),
                "source":        "live",
            })

    results.sort(key=lambda r: r["peak_stalta"], reverse=True)
    return jsonify(results[:40])

@app.route("/api/stations")
def api_stations():
    """All inventory stations with lat/lon + real-time STA/LTA for monitored ones."""
    rows = []
    for key, info in _sta_coords.items():
        # Try all common channel variants for this station
        state = None
        for chan in ("BHZ", "HHZ", "EHZ", "SHZ", "HNZ"):
            state = states.get(f"{info['net']}.{info['sta']}.{chan}")
            if state:
                break
        rows.append({
            "net":       info["net"],
            "sta":       info["sta"],
            "lat":       info["lat"],
            "lon":       info["lon"],
            "elev":      info["elev"],
            "monitored": state is not None,
            "stalta":    round(state.stalta, 2) if state else 0.0,
            "updated":   state.last_update.strftime("%H:%M:%S") if state and state.last_update else None,
            "color":     NET_COLORS.get(info["net"], "#666666"),
            "city":      state.city if state else info.get("city", ""),
        })
    return jsonify(rows)

@app.route("/api/stalta_peaks")
def api_stalta_peaks():
    """Return STA/LTA peak readings from the last hour for all stations."""
    try:
        import sqlite3, time as _time
        db_path = "/Users/OuOu/PycharmProjects/quake_alert.py/seismic_catalog.db"
        since = _time.time() - 3600
        conn = sqlite3.connect(db_path, check_same_thread=False)
        rows = conn.execute(
            "SELECT ts, key, net, sta, chan, stalta FROM stalta_peaks WHERE ts >= ? ORDER BY ts ASC",
            (since,)
        ).fetchall()
        conn.close()
        return jsonify([
            {"ts": r[0], "key": r[1], "net": r[2], "sta": r[3], "chan": r[4], "stalta": r[5]}
            for r in rows
        ])
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500

@app.route("/api/quality")
def api_quality():
    """Return list of stations flagged as noisy/anomalous."""
    return jsonify({
        "noisy": sorted(_noisy_stations),
        "trigger_counts": {k: len(v) for k,v in _station_trigger_times.items() if v}
    })

@app.route("/api/triggers")
def api_triggers():
    """Recent station triggers for GlobalQuake-style visualization."""
    cutoff = time.time() - 300  # last 5 minutes
    recent = [t for t in trigger_log if t["ts"] >= cutoff]
    return jsonify(recent)

@app.route("/api/event_list")
def api_event_list():
    """REST fallback: recent confirmed events (same as Socket.IO 'events' payload)."""
    return jsonify(get_events(limit=50))

@app.route("/api/preliminary")
def api_preliminary():
    """Recent preliminary real-time detections."""
    return jsonify(list(preliminary_events))

# ── Flight Study API ───────────────────────────────────────────────────────────
def _init_fr24():
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from FlightRadarAPI import FlightRadar24API
            return FlightRadar24API()
    except Exception:
        return None

_fr24 = _init_fr24()
_fr24_state_cache  = {}   # fr24_id -> {data, ts}
_fr24_search_cache = {}   # query   -> {data, ts}

def _fr24_haversine(lat1, lon1, lat2, lon2):
    R = 6371
    dlat = math.radians(lat2 - lat1); dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1))*math.cos(math.radians(lat2))*math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(a))

def _fr24_bearing(lat1, lon1, lat2, lon2):
    dlon = math.radians(lon2 - lon1)
    lat1, lat2 = math.radians(lat1), math.radians(lat2)
    x = math.sin(dlon)*math.cos(lat2)
    y = math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dlon)
    return (math.degrees(math.atan2(x, y)) + 360) % 360

def _ap_info(ap):
    if not ap: return None
    code = ap.get("code", {}); pos = ap.get("position", {}); reg = pos.get("region", {})
    return {
        "iata": code.get("iata",""), "icao": code.get("icao",""),
        "name": ap.get("name",""), "city": reg.get("city",""),
        "lat": pos.get("latitude"), "lon": pos.get("longitude"),
        "tz_abbr": ap.get("timezone",{}).get("abbr",""),
        "terminal": ap.get("info",{}).get("terminal"),
        "gate": ap.get("info",{}).get("gate"),
    }

@app.route("/api/flight_search")
def flight_search():
    if not _fr24:
        return jsonify({"error": "FlightRadarAPI not available"}), 503
    q = request.args.get("q","").strip()
    if len(q) < 2:
        return jsonify({"error": "too short"}), 400
    cached = _fr24_search_cache.get(q.upper())
    if cached and time.time() - cached["ts"] < 8:
        return jsonify({"results": cached["data"]})
    try:
        raw = _fr24.search(q)
        live = raw.get("live", [])
        results = []
        for item in live[:10]:
            d = item.get("detail", {})
            results.append({
                "fr24_id": item["id"], "callsign": d.get("callsign","").strip(),
                "flight": d.get("flight",""), "aircraft": d.get("ac_type",""),
                "registration": d.get("reg",""), "route": d.get("route",""),
                "schd_from": d.get("schd_from",""), "schd_to": d.get("schd_to",""),
                "lat": d.get("lat"), "lon": d.get("lon"),
                "altitude_ft": 0, "speed_kts": 0, "on_ground": False,
            })
        _fr24_search_cache[q.upper()] = {"data": results, "ts": time.time()}
        return jsonify({"results": results})
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route("/api/flight_detail/<fr24_id>")
def flight_detail(fr24_id):
    if not _fr24:
        return jsonify({"error": "FlightRadarAPI not available"}), 503
    cached = _fr24_state_cache.get(fr24_id)
    if cached and time.time() - cached["ts"] < 10:
        return jsonify(cached["data"])

    # Fetch details directly by ID — no need to scan get_flights()
    old_details = (cached or {}).get("data", {}).get("_raw_details")
    try:
        class _FakeF:
            id = fr24_id
        details = _fr24.get_flight_details(_FakeF())
        if not details or "identification" not in details:
            raise ValueError("empty response")
    except Exception:
        if old_details:
            details = old_details
        else:
            return jsonify({"error": "Flight not found or data unavailable"}), 404

    trail = details.get("trail") or []

    if trail:
        t = trail[0]
        lat, lon, alt_ft = t["lat"], t["lng"], t.get("alt", 0)
        spd_kts, hdg, vspd, on_ground = t.get("spd", 0), t.get("hd", 0), 0, False
    else:
        return jsonify({"error": "No position"}), 404

    if (not hdg or hdg == 0) and len(trail) >= 2:
        hdg = _fr24_bearing(trail[1]["lat"], trail[1]["lng"], trail[0]["lat"], trail[0]["lng"])

    ap = details.get("airport", {})
    origin = _ap_info(ap.get("origin"))
    dest   = _ap_info(ap.get("destination"))

    ti = details.get("time", {})
    est_arr = ti.get("estimated", {}).get("arrival") or ti.get("other", {}).get("eta")
    eta_seconds = max(0, int(est_arr - time.time())) if est_arr else None

    aircraft_detail = details.get("aircraft", {})
    images = aircraft_detail.get("images", {}).get("thumbnails", [])

    dist_remaining_km = _fr24_haversine(lat, lon, dest["lat"], dest["lon"]) if dest and dest.get("lat") else None
    dist_traveled_km = sum(
        _fr24_haversine(trail[i]["lat"], trail[i]["lng"], trail[i+1]["lat"], trail[i+1]["lng"])
        for i in range(len(trail)-1)
    ) if len(trail) > 1 else 0.0
    total_km = dist_traveled_km + (dist_remaining_km or 0)

    result = {
        "fr24_id": fr24_id,
        "callsign": (details.get("identification", {}).get("callsign") or "").strip(),
        "lat": lat, "lon": lon, "altitude_ft": alt_ft,
        "speed_kts": round(spd_kts, 1), "speed_ms": round(spd_kts * 0.514444, 2),
        "heading": round(hdg, 1), "vertical_rate_fpm": int(vspd),
        "on_ground": on_ground,
        "registration": aircraft_detail.get("registration", ""),
        "aircraft_model": aircraft_detail.get("model", {}).get("text", ""),
        "airline_name": details.get("airline", {}).get("name", ""),
        "airline_iata": (details.get("airline", {}).get("code") or {}).get("iata", "") if isinstance(details.get("airline", {}).get("code"), dict) else "",
        "photo_url": images[0]["src"] if images else None,
        "status": details.get("status", {}).get("text", ""),
        "origin": origin, "destination": dest,
        "sched_dep": ti.get("scheduled", {}).get("departure"),
        "sched_arr": ti.get("scheduled", {}).get("arrival"),
        "actual_dep": ti.get("real", {}).get("departure"),
        "est_arr": est_arr, "eta_seconds": eta_seconds,
        "server_time": int(time.time()),
        "dist_traveled_km": round(dist_traveled_km, 2),
        "dist_remaining_km": round(dist_remaining_km, 2) if dist_remaining_km is not None else None,
        "total_km": round(total_km, 2),
        "progress_pct": round(dist_traveled_km / total_km * 100) if total_km > 0 else None,
        "waypoints": [[t["lat"], t["lng"]] for t in reversed(trail)],
        "_raw_details": details,  # cached for stale fallback
    }
    _fr24_state_cache[fr24_id] = {"data": result, "ts": time.time()}
    return jsonify({k: v for k, v in result.items() if k != "_raw_details"})

@app.route("/livemap")
def livemap():
    return render_template_string(HTML_LIVEMAP)

@sio.on("connect")
def on_connect():
    """Send current state immediately to the newly connected client (no waiting for emit_loop)."""
    from flask_socketio import emit as _emit
    # Events
    _emit("events", {"events": get_events()})
    # Current STA/LTA snapshot
    ss = sorted(states.values(), key=lambda s: s.stalta, reverse=True)
    active    = sum(1 for s in ss if s.is_active)
    triggered = sum(1 for s in ss if s.stalta >= TRIG)
    alarmed   = sum(1 for s in ss if s.stalta >= ALARM)
    rows = []
    for s in ss[:2000]:
        coord = _sta_coords.get(f"{s.net}.{s.sta}", {})
        rows.append({
            "key": s.key, "net": s.net, "sta": s.sta, "chan": s.chan,
            "stalta": round(s.stalta, 2),
            "updated": s.last_update.strftime("%H:%M:%S") if s.last_update else None,
            "active": s.is_active,
            "color": NET_COLORS.get(s.net, "#aaa"),
            "lat": coord.get("lat"), "lon": coord.get("lon"),
            "elev": coord.get("elev"),
            "city": s.city,
        })
    _emit("stalta", {"rows": rows, "active": active, "total": len(states),
                     "triggered": triggered, "alarmed": alarmed,
                     "utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")})
    # Send recent preliminary detections so livemap renders history on load
    recent_prelims = list(preliminary_events)[:30]
    if recent_prelims:
        _emit("preliminary_history", {"events": recent_prelims})

# Track seen event IDs to detect new ones and emit alert_event.
# Pre-populate on startup so events that already exist in the DB are never
# announced to Discord (only truly NEW events fire alerts).
try:
    _seen_evids: set = {ev["evid"] for ev in get_events(limit=200) if ev.get("evid")}
    print(f"[ALERT] pre-seeded {len(_seen_evids)} existing event IDs (won't re-announce)", flush=True)
except Exception:
    _seen_evids: set = set()

# Events waiting for a real magnitude before Discord is posted.
# { evid: {"first_seen": time.time(), "age_hours": float} }
_pending_discord: dict = {}
_PENDING_MAG_TIMEOUT = 300  # send anyway after 5 min even if mag stays 0

def _check_new_events():
    """Emit alert_event for any newly seen events; announce recent ones to Discord."""
    global _seen_evids
    try:
        evs = get_events(limit=50)
        now_utc = datetime.now(timezone.utc)
        ev_by_id = {ev["evid"]: ev for ev in evs if ev.get("evid")}

        # ── Check pending events — send once mag is real (>0) or timeout ──
        for evid, meta in list(_pending_discord.items()):
            ev = ev_by_id.get(evid)
            if ev is None:
                continue
            mag = ev.get("mag")
            mag_f = float(mag) if mag is not None else 0.0
            waited = time.time() - meta["first_seen"]
            if mag_f > 0.0 or waited >= _PENDING_MAG_TIMEOUT:
                del _pending_discord[evid]
                print(f"[ALERT] Discord send for {evid} M{mag} "
                      f"({'mag updated' if mag_f > 0 else 'timeout'})", flush=True)
                threading.Thread(target=_send_event_discord, args=(ev,), daemon=True).start()
                if mag_f >= 1.5:
                    sio.emit("alert_event", ev)

        # ── Detect brand-new event IDs ────────────────────────────────────
        for ev in evs:
            evid = ev.get("evid")
            mag  = ev.get("mag")
            if not evid or evid in _seen_evids:
                continue
            _seen_evids.add(evid)

            # Only care about events within the last 24 hours
            ot_str = str(ev.get("origin_time", ""))[:19]
            try:
                ot_utc = datetime.strptime(ot_str, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
                age_hours = (now_utc - ot_utc).total_seconds() / 3600
            except ValueError:
                age_hours = 0

            if age_hours > 24:
                print(f"[ALERT] skipped Discord for {evid} — {age_hours:.1f}h old", flush=True)
                continue

            mag_f = float(mag) if mag is not None else 0.0

            if mag_f > 0.0:
                # Magnitude already available — send immediately
                print(f"[ALERT] Discord send for {evid} M{mag} (immediate)", flush=True)
                threading.Thread(target=_send_event_discord, args=(ev,), daemon=True).start()
                if mag_f >= 1.5:
                    sio.emit("alert_event", ev)
                    print(f"[ALERT] emitted alert_event for {evid} M{mag}", flush=True)
            else:
                # M0.0 — hold and wait for magnitude to be calculated
                print(f"[ALERT] holding Discord for {evid} — magnitude not yet available", flush=True)
                _pending_discord[evid] = {"first_seen": time.time(), "age_hours": age_hours}

    except Exception as e:
        print(f"[ALERT] _check_new_events error: {e}", flush=True)

# ── Emit loop ──────────────────────────────────────────────────────────────────
def emit_loop():
    """
    Main broadcast loop. Runs every 1 s; emits are staggered so the server
    never serialises multiple large payloads in the same tick.

    Timing schedule (tick mod N):
      every tick  : STA/LTA history + trigger detection (CPU only, no emit)
      every 2 s   : stalta rows broadcast  (~400 KB)
      every 2 s   : top-12 waveforms       (~33 KB), offset by 1 tick
      every 6 s   : waveforms_all          (~215 KB), offset by 3 ticks
      every 30 s  : waveforms_long         (~3 KB)
      every 15 s  : events list
      every 60 s  : station quality update
    """
    tick = 0
    while True:
        try:
            time.sleep(1); tick += 1
            ss = sorted(states.values(), key=lambda s: s.stalta, reverse=True)

            # ── Per-second: history + trigger detection (no emit) ──────────
            ts_now = time.time()
            new_triggers = []
            for s in states.values():
                s.stalta_hist.append((ts_now, s.stalta))
                was_triggered = s.triggered
                s.triggered = s.stalta >= TRIG
                if s.triggered and not was_triggered and s.last_update:
                    coord = _sta_coords.get(f"{s.net}.{s.sta}", {})
                    if coord:
                        entry = {
                            "ts": ts_now, "net": s.net, "sta": s.sta, "chan": s.chan,
                            "stalta": round(s.stalta, 2),
                            "lat": coord["lat"], "lon": coord["lon"],
                            "color": NET_COLORS.get(s.net, "#aaa"),
                        }
                        trigger_log.append(entry)
                        new_triggers.append(entry)
                        _record_trigger(f"{s.net}.{s.sta}", ts_now)
                # Update peak stalta in existing trigger_log entry for this station
                # so magnitude estimates use the actual peak, not just the rising-edge value
                elif s.triggered and s.last_update:
                    for entry in reversed(trigger_log):
                        if entry["net"] == s.net and entry["sta"] == s.sta:
                            if ts_now - entry["ts"] < PRELIM_WINDOW_SEC and s.stalta > entry["stalta"]:
                                entry["stalta"] = round(s.stalta, 2)
                            break
            if new_triggers:
                sio.emit("trigger_new", {"triggers": new_triggers})

            # ── Every 1 s: stalta rows ────────────────────────────────────
            if True:
                active    = sum(1 for s in ss if s.is_active)
                triggered = sum(1 for s in ss if s.stalta >= TRIG)
                alarmed   = sum(1 for s in ss if s.stalta >= ALARM)
                rows = []
                for s in ss[:2000]:
                    coord = _sta_coords.get(f"{s.net}.{s.sta}", {})
                    rows.append({
                        "key": s.key, "net": s.net, "sta": s.sta, "chan": s.chan,
                        "stalta": round(s.stalta, 2),
                        "updated": s.last_update.strftime("%H:%M:%S") if s.last_update else None,
                        "active": s.is_active,
                        "color": NET_COLORS.get(s.net, "#aaa"),
                        "lat": coord.get("lat"), "lon": coord.get("lon"),
                        "elev": coord.get("elev"),
                        "city": s.city,
                    })
                sio.emit("stalta", {"rows": rows, "active": active, "total": len(states),
                                    "triggered": triggered, "alarmed": alarmed,
                                    "utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")})

            # ── Every 2 s (odd ticks): top-12 waveforms ───────────────────
            if tick % 2 == 1:
                wf = []
                for s in ss[:12]:
                    pts, t_end = s.waveform_pts(n_pts=900, secs=600)
                    if pts:
                        wf.append({"label": f"{s.net}.{s.sta}", "stalta": round(s.stalta, 2),
                                   "color": NET_COLORS.get(s.net, "#4fc3f7"),
                                   "data": pts, "t_end": t_end, "secs": 600})
                sio.emit("waveforms", {"traces": wf})

            # ── Every 6 s (tick % 6 == 3): waveforms_all ─────────────────
            if tick % 6 == 3:
                all_wf = []
                for s in ss[:2000]:
                    pts, t_end = s.waveform_pts(n_pts=120, secs=60)
                    all_wf.append({"label": f"{s.net}.{s.sta}", "stalta": round(s.stalta, 2),
                                   "color": NET_COLORS.get(s.net, "#4fc3f7"),
                                   "data": pts, "t_end": t_end, "secs": 60,
                                   "live": s.is_active})
                sio.emit("waveforms_all", {"traces": all_wf})

            # ── Every 30 s: long waveforms ────────────────────────────────
            if tick % 30 == 0:
                wf_long = []
                for s in ss[:12]:
                    pts, t_end = s.waveform_pts_long(n_pts=720, secs=7200)
                    if pts:
                        wf_long.append({"label": f"{s.net}.{s.sta}", "stalta": round(s.stalta, 2),
                                        "color": NET_COLORS.get(s.net, "#4fc3f7"),
                                        "data": pts, "t_end": t_end, "secs": 7200})
                sio.emit("waveforms_long", {"traces": wf_long})

            # ── Every 15 s: events ────────────────────────────────────────
            if tick % 15 == 0:
                sio.emit("events", {"events": get_events()})
                _check_new_events()

            # ── Every 60 s: station quality ───────────────────────────────
            if tick % 60 == 0:
                _update_station_quality()

        except Exception as _e:
            # Never let an exception kill the emit loop — log and continue
            print(f"[EMIT-LOOP] error (tick={tick}): {_e}", flush=True)

# ── HTML: Main Dashboard ───────────────────────────────────────────────────────
HTML_MAIN = r"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="UTF-8"><title>SeisComP Dashboard</title>
<script src="https://cdn.socket.io/4.7.2/socket.io.min.js"></script>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
*{box-sizing:border-box;margin:0;padding:0}
:root{--bg:#0d1117;--panel:#161b22;--b1:#30363d;--b2:#21262d;--txt:#c9d1d9;--mut:#8b949e;--dim:#6e7681;--grn:#3fb950;--yel:#d29922;--red:#f85149;--blu:#58a6ff}
body{background:var(--bg);color:var(--txt);font-family:'SF Mono',monospace;height:100vh;display:flex;flex-direction:column;overflow:hidden;font-size:12px}
#hdr{background:var(--panel);border-bottom:1px solid var(--b1);padding:7px 14px;display:flex;align-items:center;gap:12px;flex-shrink:0}
#hdr h1{font-size:13px;color:#f0f6fc;font-weight:700;letter-spacing:.5px;white-space:nowrap}
.pill{padding:2px 9px;border-radius:20px;font-size:10px;font-weight:600}
.pg{background:#0d2318;color:var(--grn);border:1px solid #1a4731}
.py{background:#2d1f00;color:var(--yel);border:1px solid #5a3e00}
.pr{background:#2d0f0f;color:var(--red);border:1px solid #5a1a1a}
#utc{margin-left:auto;color:var(--dim);font-size:10px;white-space:nowrap}
#lbanner{font-size:11px;color:var(--red);font-weight:600;white-space:nowrap}
/* Alert notification toast */
#alert-toast{position:fixed;top:12px;left:50%;transform:translateX(-50%) translateY(-120px);
  z-index:9999;background:#1a0505;border:1.5px solid var(--red);border-radius:8px;
  padding:12px 20px;min-width:320px;max-width:520px;box-shadow:0 8px 32px rgba(248,81,73,.35);
  transition:transform .35s cubic-bezier(.4,0,.2,1);pointer-events:none;text-align:center}
#alert-toast.show{transform:translateX(-50%) translateY(0);pointer-events:all}
#alert-toast .at-mag{font-size:26px;font-weight:700;color:var(--red);line-height:1}
#alert-toast .at-mag.mm{color:var(--yel)}
#alert-toast .at-mag.ml{color:var(--grn)}
#alert-toast .at-loc{font-size:12px;color:var(--txt);margin-top:4px}
#alert-toast .at-time{font-size:10px;color:var(--mut);margin-top:2px}
#alert-toast .at-close{position:absolute;top:6px;right:10px;cursor:pointer;color:var(--mut);
  font-size:14px;pointer-events:all}
#alert-toast .at-close:hover{color:var(--txt)}
#main{display:flex;flex:1;overflow:hidden;flex-direction:column}
#top-row{display:flex;flex:1;overflow:hidden;min-height:0}
/* Events panel */
#ep{width:260px;flex-shrink:0;border-right:1px solid var(--b1);display:flex;flex-direction:column;overflow:hidden}
.ptitle{padding:6px 10px;font-size:10px;color:var(--mut);text-transform:uppercase;letter-spacing:1px;border-bottom:1px solid var(--b2);flex-shrink:0}
#el{flex:1;overflow-y:auto}
.evrow{padding:8px 10px;border-bottom:1px solid var(--b2);cursor:pointer;transition:background .12s}
.evrow:hover{background:var(--panel)}.evrow.sel{background:#1f2937;border-left:2px solid var(--blu)}
.evtop{display:flex;align-items:baseline;gap:6px;margin-bottom:2px}
.evmag{font-size:16px;font-weight:700}.evid{font-size:9px;color:var(--dim)}
.evloc{font-size:10px;color:var(--txt)}.evtime{font-size:10px;color:var(--dim);margin-top:1px}
.evph{font-size:9px;color:var(--dim);margin-top:1px}
.mh{color:var(--red)}.mm{color:var(--yel)}.ml{color:var(--grn)}
/* Center */
#cp{flex:1;display:flex;flex-direction:column;overflow:hidden}
#wfs{flex:0 0 50%;border-bottom:1px solid var(--b1);display:flex;flex-direction:column;overflow:hidden}
#wfh{display:flex;align-items:center;gap:6px;padding:4px 10px;border-bottom:1px solid var(--b2);flex-shrink:0;flex-wrap:wrap}
.tbtn{padding:2px 8px;border-radius:4px;font-size:10px;cursor:pointer;background:var(--b2);color:var(--mut);border:1px solid var(--b1);font-family:inherit;transition:all .15s}
.tbtn.act{background:#1f3a5f;color:var(--blu);border-color:#2d5a8e}
.wfhdiv{width:1px;height:13px;background:var(--b1);flex-shrink:0}
#wfc{flex:1;overflow-y:auto;padding:4px 8px}
.trow{display:flex;align-items:center;margin-bottom:2px;height:84px}
.tlbl{width:80px;flex-shrink:0;font-size:10px;color:var(--mut);text-align:right;padding-right:6px;line-height:1.4}
.tlbl .ts{font-weight:700;font-size:11px}.tlbl .tv{font-size:9px}
canvas.tc{flex:1;height:80px;background:#0a0f15;border-radius:2px;cursor:crosshair}
#wfag{flex:1;overflow-y:auto;display:none;padding:4px 6px;display:grid;grid-template-columns:repeat(auto-fill,minmax(115px,1fr));gap:3px;align-content:start}
.mt{background:#0a0f15;border-radius:2px;padding:2px 3px}
.mt canvas{display:block;width:100%;height:26px}
.mt .ml2{font-size:8px;color:var(--dim);white-space:nowrap;overflow:hidden}
/* STA/LTA */
#ss{flex:1;display:flex;flex-direction:column;overflow:hidden}
#sst{padding:3px 10px;font-size:10px;color:var(--dim);border-bottom:1px solid var(--b2)}
#sw{flex:1;overflow-y:auto}
table#st{width:100%;border-collapse:collapse}
#st th{padding:4px 8px;background:var(--panel);color:var(--mut);text-align:left;position:sticky;top:0;font-size:10px;border-bottom:1px solid var(--b1);z-index:1}
#st td{padding:3px 8px;border-bottom:1px solid var(--b2);vertical-align:middle}
.bw{background:var(--b2);border-radius:2px;height:7px;width:90px}
.bf{height:7px;border-radius:2px;transition:width .4s;min-width:1px}
.bg{background:#238636}.by{background:#9e6a03}.br{background:#8b1a1a}
.obtn{padding:1px 6px;border-radius:3px;font-size:9px;cursor:pointer;background:transparent;color:var(--blu);border:1px solid var(--blu);font-family:inherit}
.obtn:hover{background:#1f3a5f}
/* Bottom drawer — collapsible */
#bp{flex-shrink:0;border-top:2px solid var(--b1);display:flex;flex-direction:column;overflow:hidden;
    height:26px;transition:height .22s cubic-bezier(.4,0,.2,1)}
#bp.open{height:340px}
#btabs{display:flex;align-items:center;background:var(--panel);flex-shrink:0;height:26px}
.btab{padding:0 14px;font-size:10px;color:var(--mut);cursor:pointer;border-right:1px solid var(--b1);
      height:100%;display:flex;align-items:center;transition:all .15s;border-bottom:2px solid transparent;
      white-space:nowrap;user-select:none}
.btab.act{color:var(--blu);border-bottom-color:var(--blu);background:var(--bg)}
.btab:hover:not(.act){background:var(--b2);color:var(--txt)}
#bp-toggle{margin-left:auto;padding:0 10px;font-size:11px;color:var(--dim);cursor:pointer;
           height:100%;display:flex;align-items:center;border-left:1px solid var(--b1)}
#bp-toggle:hover{color:var(--txt)}
#bpanels{flex:1;overflow:hidden;position:relative;min-height:0}
.bpanel{display:none;position:absolute;inset:0;overflow:hidden}
.bpanel.act{display:flex;flex-direction:column}
/* Seismicity map */
#seis-map{flex:1;width:100%}
.leaflet-control-zoom{border:1px solid var(--b1)!important}
.leaflet-control-zoom a{background:var(--panel)!important;color:var(--txt)!important;border-color:var(--b1)!important}
/* M-time plot */
#mt-wrap{flex:1;padding:6px;display:flex;flex-direction:column;gap:4px}
#mt-canvas{flex:1;width:100%;border-radius:4px;background:#0a0f15}
#mt-axis{height:18px;color:var(--dim);font-size:9px;display:flex;justify-content:space-between;padding:0 4px}
/* Station health */
#hlth-wrap{flex:1;overflow-y:auto;padding:4px}
.hgrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(110px,1fr));gap:3px}
.hcell{background:var(--bg);border-radius:3px;padding:4px 6px;border-left:3px solid var(--b1)}
.hcell.live{border-left-color:var(--grn)}.hcell.slow{border-left-color:var(--yel)}.hcell.dead{border-left-color:var(--red)}
.hcell .hn{font-size:10px;font-weight:600}.hcell .hs{font-size:9px;color:var(--mut);margin-top:1px}
/* Modal */
#mo{display:none;position:fixed;inset:0;background:rgba(0,0,0,.72);z-index:100;align-items:center;justify-content:center}
#mo.open{display:flex}
#md{background:var(--panel);border:1px solid var(--b1);border-radius:8px;width:820px;max-width:95vw;max-height:92vh;display:flex;flex-direction:column;box-shadow:0 20px 60px rgba(0,0,0,.7)}
#mhdr{padding:12px 18px;border-bottom:1px solid var(--b1);display:flex;justify-content:space-between;align-items:flex-start;flex-shrink:0}
#mt{font-size:15px;font-weight:700;color:#f0f6fc}
#ms{font-size:11px;color:var(--mut);margin-top:2px;font-family:'SF Mono',monospace}
#mc{background:none;border:none;color:var(--mut);font-size:20px;cursor:pointer;line-height:1}
#mc:hover{color:var(--txt)}
#mb{flex:1;overflow-y:auto;padding:14px 18px}
.mgrid{display:grid;grid-template-columns:1fr 1fr;gap:6px;margin-bottom:14px}
.ib{background:var(--bg);border-radius:5px;padding:8px 12px}
.ib .lbl{font-size:9px;color:var(--mut);text-transform:uppercase;letter-spacing:.8px;margin-bottom:3px}
.ib .val{font-size:14px;font-weight:600;color:#f0f6fc}.ib .un{font-size:10px;color:var(--mut)}
.mbs{display:flex;gap:6px;margin-bottom:14px}
.mb2{background:var(--bg);border-radius:5px;padding:8px 12px;flex:1;text-align:center}
.mb2 .mv{font-size:22px;font-weight:700}.mb2 .mty{font-size:10px;color:var(--mut)}
.ptbl-hdr{font-size:10px;color:var(--mut);text-transform:uppercase;letter-spacing:1px;margin-bottom:6px}
table.pt{width:100%;border-collapse:collapse;font-size:11px}
.pt th{padding:4px 8px;text-align:left;color:var(--dim);font-size:10px;border-bottom:1px solid var(--b1)}
.pt td{padding:3px 8px;border-bottom:1px solid var(--b2)}
.pp{color:var(--red);font-weight:600}.sp{color:var(--grn);font-weight:600}
.city-row{display:flex;gap:10px;margin-bottom:12px;flex-wrap:wrap}
.city-chip{background:var(--bg);border:1px solid var(--b1);border-radius:4px;padding:4px 10px;font-size:11px}
.city-chip .cn{color:#f0f6fc;font-weight:600}.city-chip .cd{color:var(--mut);font-size:10px}
.mactions{display:flex;gap:8px;margin-bottom:12px;align-items:center}
.mbtn{padding:5px 14px;border-radius:5px;font-size:11px;cursor:pointer;border:1px solid var(--b1);background:var(--b2);color:var(--txt);font-family:inherit}
.mbtn:hover{background:var(--b1)}
#bb-img{border-radius:6px;margin-bottom:14px;display:none}
/* Map */
#ev-map{width:100%;height:240px;border-radius:6px;margin-bottom:14px;border:1px solid var(--b1)}
.leaflet-container{background:#0a0f15!important;font-family:'SF Mono',monospace}
::-webkit-scrollbar{width:4px;height:4px}
::-webkit-scrollbar-track{background:transparent}
::-webkit-scrollbar-thumb{background:var(--b1);border-radius:2px}
/* scrttv-style trace view */
#scrttv-wrap{display:none;flex:1;overflow-y:auto;font-family:'SF Mono',monospace}
.scrttv-row{display:flex;align-items:stretch;border-bottom:1px solid var(--b2);height:52px}
.scrttv-lbl{width:110px;flex-shrink:0;display:flex;flex-direction:column;justify-content:center;
            padding:0 6px;border-right:1px solid var(--b2);background:var(--panel)}
.scrttv-lbl .sn{font-size:10px;font-weight:700}.scrttv-lbl .sc{font-size:9px;color:var(--dim)}
.scrttv-lbl .sa{font-size:8px;color:var(--mut);margin-top:1px}
canvas.scrttv-cv{flex:1;height:52px;background:#060a10;display:block;cursor:crosshair}
#scrttv-taxis{height:16px;background:var(--panel);border-top:1px solid var(--b2);
              display:flex;align-items:center;padding:0 4px;font-size:8px;color:var(--dim);
              position:sticky;bottom:0;z-index:10;flex-shrink:0}
</style></head><body>
<!-- Alert toast notification -->
<div id="alert-toast">
  <span class="at-close" onclick="dismissAlert()">✕</span>
  <div class="at-mag" id="at-mag">M?</div>
  <div class="at-loc" id="at-loc"></div>
  <div class="at-time" id="at-time"></div>
</div>
<div id="hdr">
  <h1>SeisComP Live Dashboard</h1>
  <span class="pill pg" id="blive">— / —</span>
  <span class="pill py" id="btrig" style="display:none">—</span>
  <span class="pill pr" id="balrm" style="display:none">—</span>
  <span id="lbanner"></span>
  <a href="/livemap" target="_blank"
    style="margin-left:auto;padding:3px 12px;border-radius:20px;font-size:10px;font-weight:600;cursor:pointer;border:1px solid #2d5a8e;background:#1f3a5f;color:#58a6ff;text-decoration:none;display:flex;align-items:center;gap:5px">
    &#9651; Live Map
  </a>
  <button id="voice-btn" onclick="toggleVoice()" title="Toggle scvoice audio alerts"
    style="padding:3px 12px;border-radius:20px;font-size:10px;font-weight:600;cursor:pointer;border:1px solid;font-family:inherit;transition:all .2s">
    voice …
  </button>
  <span id="utc" style="margin-left:8px">--:--:-- UTC</span>
</div>
<div id="main">
<div id="top-row">
  <div id="ep">
    <div class="ptitle">Event Catalog</div>
    <div id="el"><div style="padding:10px;color:var(--dim)">Loading…</div></div>
  </div>
  <div id="cp">
    <div id="wfs">
      <div id="wfh">
        <div class="ptitle" style="padding:0;border:none;white-space:nowrap">Live Waveforms</div>
        <div class="wfhdiv"></div>
        <button class="tbtn act" id="btn-top" onclick="setMode('top')">Top Active</button>
        <button class="tbtn"     id="btn-all" onclick="setMode('all')">All Stations</button>
        <div class="wfhdiv"></div>
        <span style="font-size:9px;color:var(--dim);white-space:nowrap">zoom:</span>
        <button class="tbtn" id="zoom-10"   onclick="setZoom(10)">10s</button>
        <button class="tbtn" id="zoom-30"   onclick="setZoom(30)">30s</button>
        <button class="tbtn act" id="zoom-60"  onclick="setZoom(60)">1min</button>
        <button class="tbtn" id="zoom-300"  onclick="setZoom(300)">5min</button>
        <button class="tbtn" id="zoom-600"  onclick="setZoom(600)">10min</button>
        <button class="tbtn" id="zoom-1800" onclick="setZoom(1800)">30min</button>
        <button class="tbtn" id="zoom-3600" onclick="setZoom(3600)">1hr</button>
        <button class="tbtn" id="zoom-7200" onclick="setZoom(7200)">2hr</button>
        <div class="wfhdiv"></div>
        <div class="wfhdiv"></div>
        <span style="font-size:9px;color:var(--dim);white-space:nowrap">filter:</span>
        <button class="tbtn act" id="filt-raw"  onclick="setFilter('raw')">Raw</button>
        <button class="tbtn"     id="filt-hp1"  onclick="setFilter('hp1')">HP&gt;1Hz</button>
        <button class="tbtn"     id="filt-bp1"  onclick="setFilter('bp1')">BP 1-10Hz</button>
        <button class="tbtn"     id="filt-bp2"  onclick="setFilter('bp2')">BP 2-20Hz</button>
        <button class="tbtn"     id="filt-lp"   onclick="setFilter('lp')">LP&lt;1Hz</button>
        <button class="tbtn" id="btn-scrttv" onclick="setMode('scrttv')">scrttv</button>
      </div>
      <div id="wfc"></div>
      <div id="wfag" style="display:none"></div>
      <div id="scrttv-wrap"></div>
    </div>
    <div id="ss">
      <div class="ptitle" style="padding:5px 10px;border-bottom:1px solid var(--b2)">STA/LTA Monitor</div>
      <div id="sst">Connecting…</div>
      <div id="sw">
        <table id="st"><thead><tr>
          <th>Net</th><th>Station</th><th>Chan</th>
          <th>STA/LTA</th><th style="width:100px">Level</th>
          <th>Status</th><th>Updated</th><th></th>
        </tr></thead><tbody id="stb"></tbody></table>
      </div>
    </div>
  </div>
  </div><!-- /top-row -->
  <!-- Bottom tab panel -->
  <div id="bp">
    <div id="btabs">
      <div class="btab" id="tab-map"    onclick="setTab('map')">Seismicity Map</div>
      <div class="btab" id="tab-mtime"  onclick="setTab('mtime')">M–Time</div>
      <div class="btab" id="tab-health" onclick="setTab('health')">Station Health</div>
      <div class="btab" id="tab-prelim" onclick="setTab('prelim')">Detections</div>
      <div id="bp-toggle" onclick="toggleDrawer()" title="Expand/collapse">▲</div>
    </div>
    <div id="bpanels">
      <div class="bpanel act" id="panel-map">
        <div id="seis-map"></div>
      </div>
      <div class="bpanel" id="panel-mtime">
        <div id="mt-wrap">
          <canvas id="mt-canvas"></canvas>
          <div id="mt-axis"></div>
        </div>
      </div>
      <div class="bpanel" id="panel-health">
        <div id="hlth-wrap"><div class="hgrid" id="hgrid"></div></div>
      </div>
      <div class="bpanel" id="panel-prelim">
        <div style="padding:5px 10px;font-size:10px;color:var(--mut);border-bottom:1px solid var(--b2);flex-shrink:0">
          Real-time STA/LTA multi-station detections — preliminary, before SeisComP confirmation
        </div>
        <div id="prelim-list" style="flex:1;overflow-y:auto;padding:4px"></div>
      </div>
    </div>
  </div>
</div>
<!-- Modal -->
<div id="mo" onclick="closeM(event)">
  <div id="md">
    <div id="mhdr">
      <div><div id="mt">—</div><div id="ms">—</div></div>
      <button id="mc" onclick="closeM()">&#x2715;</button>
    </div>
    <div id="mb"><div style="color:var(--dim)">Loading…</div></div>
  </div>
</div>
<script>
const socket=io({transports:['websocket','polling'],upgrade:true,reconnectionDelay:1000,reconnectionAttempts:Infinity});
const ALARM=4,TRIG=2.0;
let wfMode='top';
let wfFilter='raw';
let _lastWfData=[];  // cache last received waveform data for re-render on filter change

function setFilter(f){
  wfFilter=f;
  ['raw','hp1','bp1','bp2','lp'].forEach(k=>{
    document.getElementById('filt-'+k).classList.toggle('act',k===f);
  });
  // Re-draw current data with new filter
  if(_lastWfData.length) renderWaveforms(_lastWfData);
}

function applyFilter(data, filterType){
  if(!data||data.length<5||filterType==='raw') return data;
  // Simple IIR filter implementations (normalized 0-1 data, ~1sps downsampled)
  const n=data.length;
  const out=new Array(n);
  if(filterType==='hp1'){
    // High-pass > 1 Hz: DC-blocking IIR (first-order highpass, α=0.95)
    const alpha=0.95;
    let prev=data[0], prevOut=0;
    for(let i=0;i<n;i++){
      out[i]=alpha*(prevOut + data[i] - prev);
      prev=data[i]; prevOut=out[i];
    }
  } else if(filterType==='bp1'||filterType==='bp2'){
    // Bandpass: apply highpass then lowpass
    const hp_a = filterType==='bp1' ? 0.92 : 0.85;
    const lp_b = filterType==='bp1' ? 0.4  : 0.6;
    // High-pass pass first
    const tmp=new Array(n);
    let prev=data[0], prevOut=0;
    for(let i=0;i<n;i++){
      tmp[i]=hp_a*(prevOut + data[i] - prev);
      prev=data[i]; prevOut=tmp[i];
    }
    // Then low-pass
    out[0]=tmp[0];
    for(let i=1;i<n;i++){
      out[i]=out[i-1]+lp_b*(tmp[i]-out[i-1]);
    }
  } else if(filterType==='lp'){
    // Low-pass < 1 Hz: heavy smoothing IIR
    const b=0.05;
    out[0]=data[0];
    for(let i=1;i<n;i++) out[i]=out[i-1]+b*(data[i]-out[i-1]);
  } else {
    return data;
  }
  // Re-normalize
  const mx=Math.max(...out.map(Math.abs));
  if(mx>0) for(let i=0;i<n;i++) out[i]/=mx;
  return out;
}

function renderWaveforms(traces){
  const c=document.getElementById('wfc');
  if(!traces.length){c.innerHTML='<div style="padding:10px;color:var(--dim)">No active data yet…</div>';return;}
  const need=new Set(traces.map(t=>t.label));
  c.querySelectorAll('.trow').forEach(r=>{if(!need.has(r.dataset.label))r.remove();});
  traces.forEach(t=>{
    let row=c.querySelector(`.trow[data-label="${t.label}"]`);
    if(!row){
      row=document.createElement('div');row.className='trow';row.dataset.label=t.label;
      const lbl=document.createElement('div');lbl.className='tlbl';
      lbl.innerHTML=`<div class="ts" style="color:${t.color}">${t.label}</div><div class="tv"></div>`;
      const cv=document.createElement('canvas');cv.className='tc';cv.height=80;
      row.appendChild(lbl);row.appendChild(cv);cvs[t.label]=cv;c.appendChild(row);
    }
    const v=t.stalta;
    const ve=row.querySelector('.tv');
    ve.textContent=v.toFixed(2);
    ve.style.color=v>=ALARM?'var(--red)':v>=TRIG?'var(--yel)':'var(--dim)';
    if(wfZoom<=LONG_THRESHOLD&&cvs[t.label])
      drawW(cvs[t.label],applyFilter(t.data,wfFilter),t.color,t.t_end,t.secs||600);
  });
}
const cvs={};

// ── Bottom tab panel ───────────────────────────────────────────────────────────
let _seisMap=null,_seisMarkers=[];
let _mtEvents=[];
let _staltaRows=[];

let _drawerOpen=false,_activeTab='map';

function toggleDrawer(forceOpen){
  const bp=document.getElementById('bp');
  _drawerOpen = forceOpen===true ? true : !_drawerOpen;
  bp.classList.toggle('open',_drawerOpen);
  document.getElementById('bp-toggle').textContent=_drawerOpen?'▼':'▲';
  if(_drawerOpen){
    setTimeout(()=>{
      if(_activeTab==='map'&&_seisMap)_seisMap.invalidateSize();
      if(_activeTab==='mtime')drawMTime();
    },240);
  }
}

function setTab(name){
  // If clicking the already-active tab while open → collapse
  if(name===_activeTab && _drawerOpen){toggleDrawer();return;}
  _activeTab=name;
  ['map','mtime','health','prelim'].forEach(n=>{
    const el=document.getElementById('tab-'+n);
    if(el) el.classList.toggle('act',n===name);
  });
  document.querySelectorAll('.bpanel').forEach(p=>p.classList.remove('act'));
  document.getElementById('panel-'+name).classList.add('act');
  if(!_drawerOpen)toggleDrawer(true);
  else{
    if(name==='map')setTimeout(()=>{if(_seisMap)_seisMap.invalidateSize();},50);
    if(name==='mtime')drawMTime();
    if(name==='health')drawHealth();
  }
}

// ── Seismicity map ─────────────────────────────────────────────────────────────
function initSeisMap(){
  if(_seisMap)return;
  _seisMap=L.map('seis-map',{zoomControl:true,attributionControl:false}).setView([36,-118],6);
  L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png',
    {maxZoom:18,opacity:0.85}).addTo(_seisMap);
}
function magColor(m){
  if(m>=5)return'#f85149';
  if(m>=4)return'#ff9500';
  if(m>=3)return'#d29922';
  if(m>=2)return'#3fb950';
  return'#58a6ff';
}
function magRadius(m){return Math.max(4,Math.min(18,(m||0)*3.5));}
function updateSeisMap(evs){
  if(!_seisMap)return;
  _seisMarkers.forEach(m=>m.remove());_seisMarkers=[];
  const cutoff=Date.now()-48*3600*1000;
  evs.filter(ev=>{
    const d=new Date((ev.origin_time||'').replace(' ','T')+'Z');
    return d.getTime()>cutoff && ev.lat && ev.lon;
  }).forEach(ev=>{
    const m=parseFloat(ev.mag)||0;
    const circ=L.circleMarker([ev.lat,ev.lon],{
      radius:magRadius(m),
      fillColor:magColor(m),color:'#fff',
      weight:0.8,fillOpacity:0.8
    });
    const ot=(ev.origin_time||'').slice(0,19).replace('T',' ')+' UTC';
    circ.bindTooltip(`M${ev.mag||'?'}  ${ot}<br>${ev.lat}N ${Math.abs(ev.lon)}W  ${ev.depth}km`,{sticky:true});
    circ.on('click',()=>openEv(ev.evid));
    circ.addTo(_seisMap);
    _seisMarkers.push(circ);
  });
}

// ── M-Time plot ────────────────────────────────────────────────────────────────
function drawMTime(){
  const cv=document.getElementById('mt-canvas');
  const wrap=document.getElementById('mt-wrap');
  cv.width=wrap.clientWidth-12;
  cv.height=wrap.clientHeight-30;
  const w=cv.width,h=cv.height;
  if(!w||!h)return;
  const ctx=cv.getContext('2d');
  ctx.fillStyle='#0a0f15';ctx.fillRect(0,0,w,h);
  // Grid
  ctx.strokeStyle='#1c2333';ctx.lineWidth=0.5;
  [1,2,3,4,5].forEach(m=>{
    const y=h-(m/6)*h;
    ctx.beginPath();ctx.moveTo(0,y);ctx.lineTo(w,y);ctx.stroke();
    ctx.fillStyle='#6e7681';ctx.font='9px monospace';
    ctx.fillText('M'+m,2,y-2);
  });
  if(!_mtEvents.length)return;
  const now=Date.now();
  const span=48*3600*1000;
  const t0=now-span;
  _mtEvents.forEach(ev=>{
    const mag=parseFloat(ev.mag)||0;
    const t=new Date((ev.origin_time||'').replace(' ','T')+'Z').getTime();
    if(t<t0)return;
    const x=((t-t0)/span)*w;
    const y=h-(mag/6)*h;
    ctx.beginPath();
    ctx.arc(x,y,magRadius(mag)*0.8,0,Math.PI*2);
    ctx.fillStyle=magColor(mag);ctx.globalAlpha=0.85;ctx.fill();
    ctx.strokeStyle='#fff';ctx.lineWidth=0.5;ctx.globalAlpha=0.4;ctx.stroke();
    ctx.globalAlpha=1;
  });
  // Time axis
  const ax=document.getElementById('mt-axis');
  ax.innerHTML='';
  ['-48h','-36h','-24h','-12h','now'].forEach(l=>{
    const s=document.createElement('span');s.textContent=l;ax.appendChild(s);
  });
}

// ── Station health ─────────────────────────────────────────────────────────────
function drawHealth(){
  const grid=document.getElementById('hgrid');
  const now=Date.now();
  grid.innerHTML=_staltaRows.map(r=>{
    const updated=r.updated;
    let age=Infinity,cls='dead',hs='no data';
    if(updated){
      // parse "HH:MM:SS UTC" as today's UTC time
      const parts=updated.replace(' UTC','').split(':');
      const d=new Date();
      d.setUTCHours(+parts[0],+parts[1],+parts[2],0);
      age=(now-d.getTime())/1000;
      if(age>86400)age=now/1000%86400; // handle day rollover approximation
      cls=age<60?'live':age<300?'slow':'dead';
      hs=age<60?`${Math.round(age)}s ago`:`${Math.round(age/60)}m ago`;
    }
    const vc=r.stalta>=6?'var(--red)':r.stalta>=3?'var(--yel)':'var(--grn)';
    return `<div class="hcell ${cls}">
      <div class="hn" style="color:${r.color||'#aaa'}">${r.net}.${r.sta}</div>
      <div class="hs" style="color:${vc}">STA/LTA ${r.stalta.toFixed(1)}</div>
      <div class="hs">${hs}</div>
    </div>`;
  }).join('');
}

// Receive events and feed all panels
socket.on('events',msg=>{
  _mtEvents=msg.events||[];
  updateSeisMap(_mtEvents);
  if(document.getElementById('panel-mtime').classList.contains('act'))drawMTime();
});

// Feed health panel from stalta data
const _origStalta=socket.listeners?socket.listeners('stalta'):[];
socket.on('stalta',msg=>{_staltaRows=msg.rows||[];
  if(document.getElementById('panel-health').classList.contains('act'))drawHealth();
});

window.addEventListener('resize',()=>{
  if(document.getElementById('panel-mtime').classList.contains('act'))drawMTime();
  if(_seisMap)_seisMap.invalidateSize();
});

// Init map after first render
setTimeout(()=>{initSeisMap();},500);

// ── scvoice toggle ─────────────────────────────────────────────────────────────
function _setVoiceBtn(running){
  const b=document.getElementById('voice-btn');
  if(running){
    b.textContent='voice ON';
    b.style.background='#0d2318';b.style.color='#3fb950';b.style.borderColor='#1a4731';
  } else {
    b.textContent='voice OFF';
    b.style.background='#2d1f00';b.style.color='#8b949e';b.style.borderColor='#5a3e00';
  }
}
async function toggleVoice(){
  document.getElementById('voice-btn').textContent='…';
  const r=await fetch('/api/scvoice/toggle',{method:'POST'});
  const d=await r.json();
  _setVoiceBtn(d.running);
}
fetch('/api/scvoice').then(r=>r.json()).then(d=>_setVoiceBtn(d.running));

// ── Waveform zoom state ────────────────────────────────────────────────────────
let wfZoom=60;        // seconds to display (client-side zoom)
const wfStore={};     // label -> {data, t_end, secs, color}  — short (≤10min)
const wfLongStore={}; // label -> {data, t_end, secs, color}  — long  (≤2hr)
const LONG_THRESHOLD=600; // seconds above which we switch to the long buffer

function setMode(m){
  wfMode=m;
  document.getElementById('btn-top').classList.toggle('act',m==='top');
  document.getElementById('btn-all').classList.toggle('act',m==='all');
  document.getElementById('btn-scrttv').classList.toggle('act',m==='scrttv');
  document.getElementById('wfc').style.display=m==='top'?'block':'none';
  document.getElementById('wfag').style.display=m==='all'?'grid':'none';
  document.getElementById('scrttv-wrap').style.display=m==='scrttv'?'flex':'none';
  if(m==='scrttv'){
    document.getElementById('scrttv-wrap').style.flexDirection='column';
    drawScrttvAll();
  }
}

function setZoom(s){
  wfZoom=s;
  [10,30,60,300,600,1800,3600,7200].forEach(v=>{
    const b=document.getElementById('zoom-'+v);
    if(b)b.classList.toggle('act',v===s);
  });
  // Pick the right store: short buffer for ≤10min, long buffer for >10min
  const store=s>LONG_THRESHOLD?wfLongStore:wfStore;
  // Redraw all currently-visible waveforms
  Object.entries(store).forEach(([lbl,d])=>{
    const cv=cvs[lbl];
    if(cv)drawW(cv,applyFilter(d.data,wfFilter),d.color,d.t_end,d.secs);
  });
  // If switching to long range and long data not arrived yet, show placeholder
  if(s>LONG_THRESHOLD){
    Object.keys(wfStore).forEach(lbl=>{
      if(!wfLongStore[lbl]){
        const cv=cvs[lbl];
        if(!cv)return;
        const w=cv.offsetWidth||cv.width,h=cv.height;
        cv.width=w;
        const ctx=cv.getContext('2d');
        ctx.fillStyle='#0a0f15';ctx.fillRect(0,0,w,h);
        ctx.fillStyle='#3d4f66';ctx.font='9px monospace';ctx.textAlign='center';
        ctx.fillText('long-range data arriving (up to 30s)…',w/2,h/2);
        ctx.textAlign='left';
      }
    });
  }
}

// ── Core waveform draw — with UTC time axis ────────────────────────────────────
function drawW(canvas,data,color,t_end,totalSecs){
  const w=canvas.offsetWidth||canvas.width,h=canvas.height;
  if(!w||!h)return;
  canvas.width=w;
  const ctx=canvas.getContext('2d');

  const AXIS_H=14;
  const sigH=h-AXIS_H;

  // Background
  ctx.fillStyle='#06090f';ctx.fillRect(0,0,w,h);
  ctx.fillStyle='#0a0d14';ctx.fillRect(0,sigH,w,AXIS_H);

  // Clip zoom window
  const tSecs=totalSecs||300;
  const zSecs=Math.min(wfZoom,tSecs);
  const startFrac=1-zSecs/tSecs;
  const startIdx=data?Math.max(0,Math.floor(data.length*startFrac)):0;
  const displayData=data?data.slice(startIdx):[];

  // Amplitude grid lines (at ±0.5 and ±1.0 normalized)
  ctx.strokeStyle='#111827';ctx.lineWidth=0.5;
  [0.25,0.5,0.75].forEach(f=>{
    const yp=sigH/2 - f*(sigH/2 - 2);
    const yn=sigH/2 + f*(sigH/2 - 2);
    ctx.beginPath();ctx.moveTo(0,yp);ctx.lineTo(w,yp);ctx.stroke();
    ctx.beginPath();ctx.moveTo(0,yn);ctx.lineTo(w,yn);ctx.stroke();
  });
  // Centre line
  ctx.strokeStyle='#1c2a3a';ctx.lineWidth=0.8;
  ctx.beginPath();ctx.moveTo(0,sigH/2);ctx.lineTo(w,sigH/2);ctx.stroke();

  if(displayData&&displayData.length>=2){
    const mg=3;
    const pts=displayData;
    const n=pts.length;

    // Check for saturation (clipped signal)
    const maxAbs=Math.max(...pts.map(Math.abs));
    const saturated=maxAbs>=0.98;

    // Envelope fill (min/max per pixel column for visual width)
    if(n>w*1.5){
      // downsample: compute min/max per pixel
      const envTop=new Float32Array(w);
      const envBot=new Float32Array(w);
      envTop.fill(-2);envBot.fill(2);
      for(let i=0;i<n;i++){
        const px=Math.floor((i/(n-1))*(w-1));
        if(pts[i]>envTop[px])envTop[px]=pts[i];
        if(pts[i]<envBot[px])envBot[px]=pts[i];
      }
      ctx.fillStyle=color+'22';
      ctx.beginPath();
      for(let x=0;x<w;x++){
        const yT=sigH/2-envTop[x]*(sigH/2-mg);
        if(x===0)ctx.moveTo(x,yT);else ctx.lineTo(x,yT);
      }
      for(let x=w-1;x>=0;x--){
        const yB=sigH/2-envBot[x]*(sigH/2-mg);
        ctx.lineTo(x,yB);
      }
      ctx.closePath();ctx.fill();
    }

    // Waveform line
    ctx.beginPath();ctx.strokeStyle=saturated?'#ff6060':color;ctx.lineWidth=1.2;
    for(let i=0;i<n;i++){
      const x=(i/(n-1))*w;
      const y=sigH/2-pts[i]*(sigH/2-mg);
      i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
    }
    ctx.stroke();

    // Saturation label
    if(saturated){
      ctx.fillStyle='#ff606088';ctx.font='bold 8px monospace';
      ctx.fillText('SAT',4,11);
    }
  } else if(!displayData||!displayData.length){
    ctx.fillStyle='#3d4f66';ctx.font='9px monospace';
    ctx.textAlign='center';
    ctx.fillText('no data',w/2,sigH/2+3);
    ctx.textAlign='left';
  }

  // UTC time axis
  if(t_end){
    const tStart=t_end-zSecs;
    ctx.strokeStyle='#1c2333';ctx.lineWidth=0.5;
    ctx.beginPath();ctx.moveTo(0,sigH);ctx.lineTo(w,sigH);ctx.stroke();
    ctx.font='8px monospace';ctx.fillStyle='#3d4f66';
    const intervals=[2,5,10,15,30,60,120,300,600,1800,3600];
    const targetLabels=Math.max(3,Math.floor(w/80));
    let interval=intervals.find(iv=>zSecs/iv<=targetLabels)||3600;
    const firstTick=Math.ceil(tStart/interval)*interval;
    for(let t=firstTick;t<=t_end;t+=interval){
      const frac=(t-tStart)/zSecs;
      const x=frac*w;
      ctx.strokeStyle='#2a3a50';ctx.lineWidth=0.5;
      ctx.beginPath();ctx.moveTo(x,sigH-1);ctx.lineTo(x,sigH+3);ctx.stroke();
      const d=new Date(t*1000);
      const hh=String(d.getUTCHours()).padStart(2,'0');
      const mm=String(d.getUTCMinutes()).padStart(2,'0');
      const ss=String(d.getUTCSeconds()).padStart(2,'0');
      const lbl=zSecs>3600?`${hh}:${mm}`:`${hh}:${mm}:${ss}`;
      const tw=ctx.measureText(lbl).width;
      const lx=Math.min(Math.max(x-tw/2,0),w-tw);
      ctx.fillStyle='#4a6070';
      ctx.fillText(lbl,lx,sigH+11);
    }
    ctx.fillStyle='#2a3a50';ctx.font='7px monospace';
    ctx.fillText('UTC',w-20,sigH+11);
  }
}

// ── Waveform data receive ─────────────────────────────────────────────────────
socket.on('waveforms',msg=>{
  if(wfMode!=='top')return;
  _lastWfData=msg.traces||[];
  // Cache short buffer for zoom redraws
  _lastWfData.forEach(t=>{
    wfStore[t.label]={data:t.data,color:t.color,t_end:t.t_end,secs:t.secs||600};
  });
  renderWaveforms(_lastWfData);
});

socket.on('waveforms_long',msg=>{
  if(wfMode!=='top')return;
  msg.traces.forEach(t=>{
    wfLongStore[t.label]={data:t.data,color:t.color,t_end:t.t_end,secs:t.secs||7200};
    // If we're currently in a long zoom, update the canvas now
    if(wfZoom>LONG_THRESHOLD&&cvs[t.label])
      drawW(cvs[t.label],applyFilter(t.data,wfFilter),t.color,t.t_end,t.secs||7200);
  });
});

socket.on('waveforms_all',msg=>{
  if(wfMode!=='all')return;
  const g=document.getElementById('wfag');
  msg.traces.forEach(t=>{
    let box=g.querySelector(`[data-label="${t.label}"]`);
    if(!box){
      box=document.createElement('div');box.className='mt';box.dataset.label=t.label;
      const lbl=document.createElement('div');lbl.className='ml2';lbl.style.color=t.color;lbl.textContent=t.label;
      const cv=document.createElement('canvas');cv.height=36;
      box.appendChild(lbl);box.appendChild(cv);g.appendChild(box);
    }
    box.style.opacity=t.live?1:.3;
    if(t.data&&t.data.length)drawW(box.querySelector('canvas'),applyFilter(t.data,wfFilter),t.color,t.t_end,t.secs||60);
  });
});

socket.on('stalta',msg=>{
  document.getElementById('utc').textContent=msg.utc+' UTC';
  document.getElementById('blive').textContent=`${msg.active} / ${msg.total} live`;
  const bt=document.getElementById('btrig'),ba=document.getElementById('balrm');
  bt.style.display=msg.triggered>0?'':'none';ba.style.display=msg.alarmed>0?'':'none';
  if(msg.triggered>0)bt.textContent=`${msg.triggered} triggered`;
  if(msg.alarmed>0)ba.textContent=`${msg.alarmed} alarmed`;
  document.getElementById('sst').textContent=
    `${msg.active}/${msg.total} stations live  ·  ${msg.triggered} triggered  ·  ${msg.alarmed} alarmed`;
  const rows=msg.rows,tb=document.getElementById('stb');
  while(tb.children.length>rows.length)tb.removeChild(tb.lastChild);
  while(tb.children.length<rows.length)tb.appendChild(document.createElement('tr'));
  rows.forEach((r,i)=>{
    const tr=tb.children[i];
    const pct=Math.min((r.stalta/8)*100,100).toFixed(1);
    const bc=r.stalta>=ALARM?'br':r.stalta>=TRIG?'by':'bg';
    const vc=r.stalta>=ALARM?'var(--red)':r.stalta>=TRIG?'var(--yel)':'var(--grn)';
    const st=r.stalta>=ALARM?'ALARM':r.stalta>=TRIG?'TRIGGER':r.updated?'quiet':'waiting';
    const [net,sta,chan]=r.key.split('.');
    tr.innerHTML=`
      <td style="color:${r.color};font-weight:600">${r.net}</td>
      <td style="color:${r.color}">${r.sta}</td>
      <td style="color:var(--dim)">${r.chan}</td>
      <td style="color:${vc};font-weight:${r.stalta>=TRIG?700:400}">${r.stalta.toFixed(2)}</td>
      <td><div class="bw"><div class="bf ${bc}" style="width:${pct}%"></div></div></td>
      <td style="color:${vc}">${st}</td>
      <td style="color:var(--dim)">${r.updated?r.updated+' UTC':'—'}</td>
      <td><button class="obtn" onclick="window.open('/station/${net}/${sta}/${chan}','_blank')">waveform</button></td>`;
  });
});

let selEvid=null;
socket.on('events',msg=>{
  const evs=msg.events,el=document.getElementById('el');
  if(!evs.length){el.innerHTML='<div style="padding:10px;color:var(--dim)">No events</div>';return;}
  const lat=evs[0];
  if(lat.mag)document.getElementById('lbanner').textContent=`Latest  M${lat.mag}  ${(lat.origin_time||'').slice(11,19)} UTC`;
  el.innerHTML=evs.map(ev=>{
    const m=parseFloat(ev.mag)||0,mc=m>=4?'mh':m>=3?'mm':'ml',sel=ev.evid===selEvid?' sel':'';
    return `<div class="evrow${sel}" data-evid="${ev.evid}" onclick="openEv('${ev.evid}')">
      <div class="evtop"><span class="evmag ${mc}">M${ev.mag||'?'}</span><span class="evid">${ev.evid}</span></div>
      <div class="evloc">${ev.lat}N  ${Math.abs(ev.lon)}W  ·  ${ev.depth} km</div>
      <div class="evtime">${(ev.origin_time||'').slice(0,19)} UTC</div>
      <div class="evph">${ev.phases} phases</div>
    </div>`;
  }).join('');
});

let _alertTimer=null;
function dismissAlert(){
  document.getElementById('alert-toast').classList.remove('show');
  if(_alertTimer){clearTimeout(_alertTimer);_alertTimer=null;}
}
function showAlertToast(ev){
  const mag=parseFloat(ev.mag)||0;
  const mc=mag>=4?'mh':mag>=3?'mm':'ml';
  const ns=ev.lat>=0?'N':'S', ew=ev.lon<=0?'W':'E';
  document.getElementById('at-mag').className='at-mag '+mc;
  document.getElementById('at-mag').textContent='M'+(ev.mag||'?')+' Earthquake';
  document.getElementById('at-loc').textContent=
    Math.abs(ev.lat).toFixed(3)+'°'+ns+'  '+Math.abs(ev.lon).toFixed(3)+'°'+ew+
    '  ·  '+ev.depth+' km depth  ·  '+ev.phases+' phases';
  document.getElementById('at-time').textContent=(ev.origin_time||'').slice(0,19)+' UTC  |  '+ev.evid;
  const t=document.getElementById('alert-toast');
  t.classList.add('show');
  if(_alertTimer){clearTimeout(_alertTimer);_alertTimer=null;}
  // Toast stays until manually dismissed — click ✕ to close
}
socket.on('alert_event', ev => {
  showAlertToast(ev);
  // Also highlight this event in the list
  openEv(ev.evid);
});

let _evMap=null,_evMapEvid=null;
function _initMap(lat,lon,picks){
  const mapDiv=document.getElementById('ev-map');
  if(!mapDiv)return;
  if(_evMap){_evMap.remove();_evMap=null;}
  _evMap=L.map('ev-map',{zoomControl:true,attributionControl:false}).setView([lat,lon],7);
  L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png',
    {maxZoom:18,opacity:0.85}).addTo(_evMap);
  // Epicenter star
  const starIcon=L.divIcon({className:'',html:`<svg width="22" height="22" viewBox="-11 -11 22 22">
    <polygon points="0,-9 2.1,-2.9 8.6,-2.9 3.4,1.1 5.3,7.2 0,3.5 -5.3,7.2 -3.4,1.1 -8.6,-2.9 -2.1,-2.9"
      fill="#f85149" stroke="#fff" stroke-width="1"/>
    </svg>`,iconSize:[22,22],iconAnchor:[11,11]});
  L.marker([lat,lon],{icon:starIcon,zIndexOffset:1000})
    .bindTooltip(`Epicenter<br>${lat.toFixed(3)}°N ${Math.abs(lon).toFixed(3)}°W`,{sticky:true})
    .addTo(_evMap);
  // Station triangles
  const staSeen=new Set();
  picks.forEach(p=>{
    if(!p.sta_lat||!p.sta_lon)return;
    const k=`${p.net}.${p.sta}`;
    if(staSeen.has(k))return; staSeen.add(k);
    const isP=p.phase&&p.phase.startsWith('P');
    const col=isP?'#58a6ff':'#3fb950';
    const triIcon=L.divIcon({className:'',html:`<svg width="16" height="14" viewBox="0 0 16 14">
      <polygon points="8,1 15,13 1,13" fill="${col}" stroke="#fff" stroke-width="1" opacity="0.85"/>
      </svg>`,iconSize:[16,14],iconAnchor:[8,13]});
    L.marker([p.sta_lat,p.sta_lon],{icon:triIcon})
      .bindTooltip(`${k}<br>${p.phase||'?'}  ${p.pick_time?p.pick_time.slice(11,23)+' UTC':''}`,{sticky:true})
      .addTo(_evMap);
  });
  // Fit bounds to include all markers
  const pts=[[lat,lon],...picks.filter(p=>p.sta_lat).map(p=>[p.sta_lat,p.sta_lon])];
  if(pts.length>1)_evMap.fitBounds(L.latLngBounds(pts).pad(0.15));
  setTimeout(()=>_evMap.invalidateSize(),100);
}

async function openEv(evid){
  selEvid=evid;
  document.querySelectorAll('.evrow').forEach(r=>r.classList.toggle('sel',r.dataset.evid===evid));
  document.getElementById('mt').textContent=evid;
  document.getElementById('ms').textContent='Loading…';
  document.getElementById('mb').innerHTML='<div style="padding:20px;color:var(--dim)">Loading…</div>';
  document.getElementById('mo').classList.add('open');
  const res=await fetch(`/api/event/${evid}`);
  const ev=await res.json();
  if(ev.error){document.getElementById('mb').innerHTML=`<div style="color:var(--red)">${ev.error}</div>`;return;}
  const mags=ev.magnitudes||[],picks=ev.picks||[],cities=ev.cities||[];
  const lat=parseFloat(ev.lat||0),lon=parseFloat(ev.lon||0);
  const dep=parseFloat(ev.depth||0).toFixed(1);
  const rms=ev.rms?parseFloat(ev.rms).toFixed(2):'—';
  const gap=ev.az_gap?parseFloat(ev.az_gap).toFixed(0)+'°':'—';
  const dist=ev.min_dist?parseFloat(ev.min_dist).toFixed(1)+'°':'—';
  const mTag=mags.length?`M${mags[0].mag} ${mags[0].type}`:'—';
  document.getElementById('mt').textContent=`${evid}  —  ${mTag}`;
  document.getElementById('ms').textContent=(ev.origin_time||'').slice(0,19).replace('T',' ')+' UTC';
  const magBoxes=mags.map(m=>{
    const mc=m.mag>=4?'var(--red)':m.mag>=3?'var(--yel)':'var(--grn)';
    return `<div class="mb2"><div class="mv" style="color:${mc}">${m.mag}</div><div class="mty">${m.type} · ${m.stationCount||'?'} sta</div></div>`;
  }).join('');
  const cityChips=cities.map(c=>`<div class="city-chip"><span class="cn">${c.name}</span><span class="cd">, ${c.state} — ${c.km} km away</span></div>`).join('');
  const pickRows=picks.map(p=>{
    const pc=p.phase&&p.phase.startsWith('P')?'pp':'sp';
    const res2=p.residual!=null?parseFloat(p.residual).toFixed(2)+'s':'—';
    const d2=p.distance!=null?parseFloat(p.distance).toFixed(1)+'°':'—';
    const az=p.azimuth!=null?parseFloat(p.azimuth).toFixed(0)+'°':'—';
    const wt=p.weight!=null?parseFloat(p.weight).toFixed(2):'—';
    // Full UTC timestamp for pick time
    const t=p.pick_time?p.pick_time.slice(0,23).replace('T',' ')+' UTC':'—';
    const yellow=p.residual!=null&&Math.abs(p.residual)>1.5?'color:var(--yel)':'';
    return `<tr><td style="color:var(--txt)">${p.net}.${p.sta}</td><td style="color:var(--dim)">${p.cha}</td>
      <td class="${pc}">${p.phase||'—'}</td>
      <td style="color:var(--dim);font-size:10px;white-space:nowrap">${t}</td>
      <td>${d2}</td><td>${az}</td>
      <td style="${yellow}">${res2}</td><td style="color:var(--dim)">${wt}</td></tr>`;
  }).join('');
  document.getElementById('mb').innerHTML=`
    <div class="mactions">
      <button class="mbtn" onclick="window.open('/event/${evid}/waveforms','_blank')">View Waveforms</button>
      <a href="/api/beachball/${evid}.png" target="_blank" class="mbtn" style="text-decoration:none">Focal Mechanism</a>
      <img id="bb-img" src="/api/beachball/${evid}.png" height="100" style="border-radius:4px;display:none" onerror="this.style.display='none'" onload="this.style.display='block'">
    </div>
    ${cities.length?`<div class="city-row">${cityChips}</div>`:''}
    <div class="mbs">${magBoxes||'<div style="color:var(--dim)">No magnitude yet</div>'}</div>
    <div class="mgrid">
      <div class="ib"><div class="lbl">Origin Time (UTC)</div><div class="val" style="font-size:12px">${(ev.origin_time||'—').slice(0,19).replace('T',' ')}</div></div>
      <div class="ib"><div class="lbl">Location</div><div class="val" style="font-size:12px">${lat.toFixed(3)}°N  ${Math.abs(lon).toFixed(3)}°W</div></div>
      <div class="ib"><div class="lbl">Depth</div><div class="val">${dep}<span class="un"> km</span></div></div>
      <div class="ib"><div class="lbl">Phases Used</div><div class="val">${ev.phases||'—'}<span class="un"> arrivals</span></div></div>
      <div class="ib"><div class="lbl">RMS Residual</div><div class="val">${rms}<span class="un"> s</span></div></div>
      <div class="ib"><div class="lbl">Azimuthal Gap</div><div class="val">${gap}</div></div>
      <div class="ib"><div class="lbl">Min Station Dist</div><div class="val">${dist}</div></div>
      <div class="ib"><div class="lbl">Stations Used</div><div class="val">${ev.stations||'—'}</div></div>
    </div>
    <div id="ev-map"></div>
    <div class="ptbl-hdr" style="margin-top:14px">Peak STA/LTA by Station <span id="stalta-evt-note" style="font-weight:400;color:var(--dim)">loading…</span></div>
    <div id="stalta-evt-wrap" style="margin-bottom:14px">
      <div style="color:var(--dim);font-size:10px;padding:8px 0">Fetching station response data…</div>
    </div>
    <div class="ptbl-hdr">Phase Arrivals (${picks.length} picks)</div>
    <table class="pt"><thead><tr><th>Station</th><th>Chan</th><th>Phase</th><th>Pick Time (UTC)</th><th>Dist</th><th>Az</th><th>Residual</th><th>Weight</th></tr></thead>
    <tbody>${pickRows||'<tr><td colspan="8" style="color:var(--dim);padding:8px">No picks found</td></tr>'}</tbody></table>`;
  requestAnimationFrame(()=>_initMap(lat,lon,picks));

  // Async-load peak STA/LTA per station
  fetch(`/api/event_stalta/${evid}`).then(r=>r.json()).then(rows=>{
    const wrap=document.getElementById('stalta-evt-wrap');
    const note=document.getElementById('stalta-evt-note');
    if(!wrap)return;
    if(rows.error){note.textContent='(no history — event may be older than 4h)';wrap.innerHTML='';return;}
    if(!rows.length){note.textContent='(no data in window)';wrap.innerHTML='';return;}
    note.textContent=`(${rows.length} stations)`;
    wrap.innerHTML=`<table class="pt" style="font-size:11px">
      <thead><tr>
        <th>Station</th><th>Chan</th>
        <th>Peak STA/LTA</th><th style="width:120px">Level</th>
        <th>Status</th><th>Peak Time (UTC)</th><th>Delay</th>
      </tr></thead>
      <tbody>${rows.map(r=>{
        const vc=r.peak_stalta>=ALARM?'var(--red)':r.peak_stalta>=TRIG?'var(--yel)':'var(--grn)';
        const st=r.peak_stalta>=ALARM?'ALARM':r.peak_stalta>=TRIG?'TRIGGER':'quiet';
        const pct=Math.min((r.peak_stalta/8)*100,100).toFixed(1);
        const bc=r.peak_stalta>=ALARM?'br':r.peak_stalta>=TRIG?'by':'bg';
        const delay=r.delay_s>=0?`+${r.delay_s}s`:`${r.delay_s}s`;
        return `<tr>
          <td style="color:${r.color};font-weight:600">${r.net}.${r.sta}</td>
          <td style="color:var(--dim)">${r.chan}</td>
          <td style="color:${vc};font-weight:${r.peak_stalta>=TRIG?700:400}">${r.peak_stalta.toFixed(2)}</td>
          <td><div class="bw" style="width:120px"><div class="bf ${bc}" style="width:${pct}%"></div></div></td>
          <td style="color:${vc}">${st}</td>
          <td style="color:var(--dim)">${r.peak_time_utc} UTC</td>
          <td style="color:var(--dim)">${delay}</td>
        </tr>`;
      }).join('')}</tbody></table>`;
  }).catch(()=>{
    const note=document.getElementById('stalta-evt-note');
    if(note)note.textContent='(error loading)';
  });
}
function closeM(e){
  if(!e||e.target===document.getElementById('mo')){
    document.getElementById('mo').classList.remove('open');
    if(_evMap){_evMap.remove();_evMap=null;}
  }
}
document.addEventListener('keydown',e=>{if(e.key==='Escape')closeM();});

// ── scrttv-style helicorder view ─────────────────────────────────────────────
const scrtvCvs = {};
const SCRTTV_WIN = 600; // 10-minute window like scrttv default

function drawScrttvAll() {
  const wrap = document.getElementById('scrttv-wrap');
  if (!wrap || wfMode !== 'scrttv') return;
  const store = wfStore;
  const labels = Object.keys(store);
  if (!labels.length) { wrap.innerHTML = '<div style="padding:10px;color:var(--dim)">No data yet…</div>'; return; }

  // Keep existing rows, add missing ones
  labels.forEach(lbl => {
    if (!wrap.querySelector(`[data-scrttv="${lbl}"]`)) {
      const d = store[lbl];
      const row = document.createElement('div');
      row.className = 'scrttv-row'; row.dataset.scrttv = lbl;
      const [net, sta] = lbl.split('.');
      const staltaV = d.stalta ? d.stalta.toFixed(2) : '0.00';
      const staltaC = (d.stalta||0) >= ALARM ? 'var(--red)' : (d.stalta||0) >= TRIG ? 'var(--yel)' : 'var(--grn)';
      row.innerHTML = `
        <div class="scrttv-lbl">
          <div class="sn" style="color:${d.color||'#4fc3f7'}">${net}.${sta}</div>
          <div class="sc">BHZ</div>
          <div class="sa" style="color:${staltaC}" data-sl="${lbl}">STA/LTA ${staltaV}</div>
        </div>`;
      const cv = document.createElement('canvas');
      cv.className = 'scrttv-cv'; cv.height = 52;
      row.appendChild(cv);
      scrtvCvs[lbl] = cv;
      wrap.appendChild(row);
    }
  });

  // Draw each
  labels.forEach(lbl => {
    const d = store[lbl];
    if (!d || !d.data) return;
    const cv = scrtvCvs[lbl];
    if (!cv) return;
    drawScrttvTrace(cv, d.data, d.color, d.t_end, d.secs, d.stalta || 0);
    // Update STA/LTA label
    const sl = wrap.querySelector(`[data-sl="${lbl}"]`);
    if (sl) {
      const v = d.stalta || 0;
      sl.textContent = 'STA/LTA ' + v.toFixed(2);
      sl.style.color = v >= ALARM ? 'var(--red)' : v >= TRIG ? 'var(--yel)' : 'var(--grn)';
    }
  });
}

function drawScrttvTrace(canvas, data, color, t_end, totalSecs, stalta) {
  const w = canvas.offsetWidth || canvas.width || 600;
  const h = canvas.height;
  if (!w || !h || !data || data.length < 2) return;
  canvas.width = w;
  const ctx = canvas.getContext('2d');

  // Dark background
  ctx.fillStyle = '#060a10'; ctx.fillRect(0, 0, w, h);

  // Amplitude fill color: red if alarmed, orange if triggered, normal otherwise
  const fillCol = stalta >= ALARM ? '#f8514922' : stalta >= TRIG ? '#ff950022' : color + '18';
  const lineCol = stalta >= ALARM ? '#f85149' : stalta >= TRIG ? '#ff9500' : color;

  // Show last SCRTTV_WIN seconds from the buffer
  const zSecs = SCRTTV_WIN;
  const startFrac = 1 - zSecs / totalSecs;
  const startIdx = Math.max(0, Math.floor(data.length * startFrac));
  const pts = data.slice(startIdx);
  if (pts.length < 2) return;

  const mg = 3;
  const mid = h / 2;
  const amp = mid - mg;

  // Zero line
  ctx.strokeStyle = '#1c2a3a'; ctx.lineWidth = 0.5;
  ctx.beginPath(); ctx.moveTo(0, mid); ctx.lineTo(w, mid); ctx.stroke();

  // Amplitude fill
  ctx.fillStyle = fillCol;
  ctx.beginPath();
  for (let i = 0; i < pts.length; i++) {
    const x = (i / (pts.length - 1)) * w;
    const y = mid - pts[i] * amp;
    i === 0 ? ctx.moveTo(x, mid) : null;
    ctx.lineTo(x, y);
  }
  ctx.lineTo(w, mid); ctx.closePath(); ctx.fill();

  // Waveform line
  ctx.beginPath(); ctx.strokeStyle = lineCol; ctx.lineWidth = 1;
  for (let i = 0; i < pts.length; i++) {
    const x = (i / (pts.length - 1)) * w;
    const y = mid - pts[i] * amp;
    i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
  }
  ctx.stroke();

  // UTC time ticks (at right edge: show current time; at left: -10min)
  if (t_end) {
    const tStart = t_end - zSecs;
    ctx.font = '8px monospace'; ctx.fillStyle = '#3d4f66';
    const intervals = [60, 120, 300];
    const iv = intervals.find(v => zSecs / v <= 6) || 300;
    const first = Math.ceil(tStart / iv) * iv;
    for (let t = first; t <= t_end; t += iv) {
      const x = ((t - tStart) / zSecs) * w;
      ctx.strokeStyle = '#1c2a3a'; ctx.lineWidth = 0.5;
      ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, h); ctx.stroke();
      const d = new Date(t * 1000);
      const lbl = `${String(d.getUTCHours()).padStart(2,'0')}:${String(d.getUTCMinutes()).padStart(2,'0')}`;
      ctx.fillText(lbl, Math.min(x + 2, w - 28), 9);
    }
  }
}

// Hook waveforms event to also update scrttv
socket.on('waveforms', msg => {
  if (wfMode !== 'scrttv') return;
  msg.traces.forEach(t => {
    wfStore[t.label] = { data: t.data, color: t.color, t_end: t.t_end, secs: t.secs, stalta: t.stalta };
  });
  drawScrttvAll();
});

// ── Watchdog: reconnect if no stalta for 20 s ─────────────────────────────────
let _lastStaltaTs = Date.now();
socket.on('stalta', () => { _lastStaltaTs = Date.now(); });
socket.on('connect', () => {
  _lastStaltaTs = Date.now();
  const u = document.getElementById('utc');
  if (u) { u.style.color = 'var(--grn)'; }
});
socket.on('disconnect', () => {
  const u = document.getElementById('utc');
  if (u) { u.textContent = 'disconnected — reconnecting…'; u.style.color = 'var(--red)'; }
});
setInterval(() => {
  const stale = (Date.now() - _lastStaltaTs) / 1000;
  if (stale > 20 && socket.connected) {
    console.warn('[watchdog] stalta stale', stale.toFixed(0), 's — forcing reconnect');
    const u = document.getElementById('utc');
    if (u) { u.textContent = 'stale — reconnecting…'; u.style.color = 'var(--yel)'; }
    socket.disconnect();
    setTimeout(() => socket.connect(), 500);
  }
}, 5000);
</script></body></html>"""

# ── HTML: Station Waveform Page ────────────────────────────────────────────────
HTML_STATION = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>{{ net }}.{{ sta }}.{{ chan }} — Live</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
html,body{height:100%;background:#0d1117;color:#c9d1d9;font-family:'SF Mono',monospace;display:flex;flex-direction:column;overflow:hidden}
:root{--bg:#0d1117;--panel:#161b22;--b1:#30363d;--b2:#21262d;--txt:#c9d1d9;--mut:#8b949e;--dim:#6e7681;--grn:#3fb950;--yel:#d29922;--red:#f85149;--blu:#58a6ff}
#hdr{background:var(--panel);border-bottom:1px solid var(--b1);padding:0 16px;display:flex;align-items:center;gap:14px;flex-shrink:0;height:44px}
#hdr a{color:var(--mut);text-decoration:none;font-size:11px}
#hdr a:hover{color:var(--txt)}
#hdr h1{font-size:14px;font-weight:700;color:{{ color }};white-space:nowrap}
.stat{font-size:10px;color:var(--mut)}
#stalta-val{font-size:24px;font-weight:700;margin-left:auto;letter-spacing:-0.5px}
#zoom-bar{background:var(--b2);border-bottom:1px solid var(--b1);padding:4px 12px;display:flex;align-items:center;gap:6px;flex-shrink:0}
.zbtn{padding:2px 9px;border-radius:4px;font-size:10px;cursor:pointer;background:var(--bg);color:var(--mut);border:1px solid var(--b1);font-family:inherit;transition:all .15s}
.zbtn.act{background:#1f3a5f;color:var(--blu);border-color:#2d5a8e}
#cv-wrap{flex:1;padding:8px 8px 0 8px;min-height:0;display:flex;flex-direction:column}
canvas#wv{display:block;background:#0a0f15;border-radius:4px 4px 0 0;flex:1;width:100%;cursor:crosshair}
#info-bar{flex-shrink:0;padding:3px 10px;background:#0a0f15;border-radius:0 0 4px 4px;border-top:1px solid #1c2333;margin-bottom:8px;font-size:9px;color:var(--dim);display:flex;gap:16px}
::-webkit-scrollbar{width:4px}::-webkit-scrollbar-track{background:transparent}::-webkit-scrollbar-thumb{background:var(--b1);border-radius:2px}
</style></head><body>
<div id="hdr">
  <a href="/">← Dashboard</a>
  <h1>{{ net }}.{{ sta }}.{{ chan }}</h1>
  <div style="display:flex;flex-direction:column;gap:1px">
    <div class="stat" id="updated">connecting…</div>
    <div class="stat" id="buf-info"></div>
  </div>
  <div id="stalta-val" style="color:var(--grn)">—</div>
</div>
<div id="zoom-bar">
  <span style="font-size:9px;color:var(--dim)">zoom:</span>
  <button class="zbtn" id="z-10"   onclick="setZoom(10)">10s</button>
  <button class="zbtn" id="z-30"   onclick="setZoom(30)">30s</button>
  <button class="zbtn act" id="z-60" onclick="setZoom(60)">1min</button>
  <button class="zbtn" id="z-300"  onclick="setZoom(300)">5min</button>
  <button class="zbtn" id="z-600"  onclick="setZoom(600)">10min</button>
  <button class="zbtn" id="z-1800" onclick="setZoom(1800)">30min</button>
  <button class="zbtn" id="z-3600" onclick="setZoom(3600)">1hr</button>
  <button class="zbtn" id="z-7200" onclick="setZoom(7200)">2hr</button>
</div>
<div id="cv-wrap">
  <canvas id="wv"></canvas>
  <div id="info-bar">
    <span id="t-range">—</span>
    <span id="pts-info">—</span>
    <span style="margin-left:auto" id="stalta-lbl">STA/LTA: —</span>
  </div>
</div>
<script>
const KEY='{{ key }}',COLOR='{{ color }}';
const TRIG=3,ALARM=6,LONG_THRESH=600;
const canvas=document.getElementById('wv');
const ctx=canvas.getContext('2d');

let curZoom=60;
let shortData={values:null,t_end:null,secs:null};
let longData ={values:null,t_end:null,secs:null};
let lastStalta=0,lastUpdated='';

// ── Zoom ───────────────────────────────────────────────────────────────────────
function setZoom(s){
  curZoom=s;
  [10,30,60,300,600,1800,3600,7200].forEach(v=>{
    const b=document.getElementById('z-'+v);
    if(b)b.classList.toggle('act',v===s);
  });
  draw();
  // Trigger long refresh immediately if switching to long range
  if(s>LONG_THRESH&&!longData.values)refreshLong();
}

// ── Draw ───────────────────────────────────────────────────────────────────────
function resize(){
  const wrap=document.getElementById('cv-wrap');
  canvas.width =wrap.clientWidth  - 16;
  canvas.height=wrap.clientHeight - 44; // subtract info-bar+padding
  draw();
}

function draw(){
  const w=canvas.width,h=canvas.height;
  if(!w||!h)return;
  const AXIS_H=16;
  const sigH=h-AXIS_H;
  const src=curZoom>LONG_THRESH?longData:shortData;

  ctx.fillStyle='#0a0f15';ctx.fillRect(0,0,w,h);
  ctx.fillStyle='#080c12';ctx.fillRect(0,sigH,w,AXIS_H);

  // Horizontal grid lines (4 levels)
  ctx.strokeStyle='#141b25';ctx.lineWidth=0.5;
  [0.25,0.5,0.75].forEach(f=>{
    ctx.beginPath();ctx.moveTo(0,sigH*f);ctx.lineTo(w,sigH*f);ctx.stroke();
  });
  // Zero line
  ctx.strokeStyle='#1e2d40';ctx.lineWidth=1;
  ctx.beginPath();ctx.moveTo(0,sigH/2);ctx.lineTo(w,sigH/2);ctx.stroke();

  if(!src.values||!src.values.length){
    ctx.fillStyle='#3d4f66';ctx.font='11px monospace';ctx.textAlign='center';
    const msg=curZoom>LONG_THRESH?'Long-range data loading (up to 30s)…':'No data yet';
    ctx.fillText(msg,w/2,sigH/2+4);
    ctx.textAlign='left';
    return;
  }

  // Clip to zoom window
  const tSecs=src.secs||600;
  const zSecs=Math.min(curZoom,tSecs);
  const startFrac=1-zSecs/tSecs;
  const startIdx=Math.max(0,Math.floor(src.values.length*startFrac));
  const disp=src.values.slice(startIdx);

  // ── Waveform ────────────────────────────────────────────────────────────────
  if(disp.length>=2){
    ctx.beginPath();ctx.strokeStyle=COLOR;ctx.lineWidth=1.5;
    const mg=6;
    for(let i=0;i<disp.length;i++){
      const x=(i/(disp.length-1))*w;
      const y=sigH/2-disp[i]*(sigH/2-mg);
      i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
    }
    ctx.stroke();
  }

  // ── UTC time axis ─────────────────────────────────────────────────────────
  if(src.t_end){
    const tStart=src.t_end-zSecs;
    ctx.strokeStyle='#1a2535';ctx.lineWidth=0.5;
    ctx.beginPath();ctx.moveTo(0,sigH);ctx.lineTo(w,sigH);ctx.stroke();
    ctx.font='9px monospace';

    const intervals=[1,2,5,10,15,30,60,120,300,600,900,1800,3600];
    const target=6;
    let iv=intervals.find(i=>zSecs/i<=target)||3600;

    const firstTick=Math.ceil(tStart/iv)*iv;
    for(let t=firstTick;t<=src.t_end+0.5;t+=iv){
      const frac=(t-tStart)/zSecs;
      if(frac<0||frac>1)continue;
      const x=frac*w;
      ctx.strokeStyle='#243040';ctx.lineWidth=0.5;
      ctx.beginPath();ctx.moveTo(x,sigH-2);ctx.lineTo(x,sigH+4);ctx.stroke();
      // Vertical grid line through signal area
      ctx.strokeStyle='#0f1620';ctx.lineWidth=0.5;
      ctx.beginPath();ctx.moveTo(x,0);ctx.lineTo(x,sigH);ctx.stroke();
      const d=new Date(t*1000);
      const lbl=`${String(d.getUTCHours()).padStart(2,'0')}:${String(d.getUTCMinutes()).padStart(2,'0')}:${String(d.getUTCSeconds()).padStart(2,'0')}`;
      const tw=ctx.measureText(lbl).width;
      ctx.fillStyle='#3d5070';
      ctx.fillText(lbl,Math.min(Math.max(x-tw/2,0),w-tw),sigH+12);
    }
    // UTC label
    ctx.fillStyle='#243040';ctx.font='8px monospace';
    ctx.fillText('UTC',w-22,sigH+12);

    // Info bar
    const fmt=t=>{const d=new Date(t*1000);return d.toISOString().slice(11,19)+' UTC';};
    document.getElementById('t-range').textContent=`${fmt(tStart)} → ${fmt(src.t_end)}`;
  }
  document.getElementById('pts-info').textContent=`${disp.length} pts / ${zSecs}s`;
}

// ── STA/LTA display ───────────────────────────────────────────────────────────
function showStalta(v,updated){
  lastStalta=v;lastUpdated=updated||'';
  const el=document.getElementById('stalta-val');
  el.textContent=v.toFixed(2);
  el.style.color=v>=ALARM?'var(--red)':v>=TRIG?'var(--yel)':'var(--grn)';
  document.getElementById('stalta-lbl').textContent=`STA/LTA: ${v.toFixed(2)}`;
  if(updated)document.getElementById('updated').textContent=updated;
  const bi=document.getElementById('buf-info');
  bi.textContent=v>=ALARM?'⚠ ALARM':v>=TRIG?'▲ TRIGGERED':'quiet';
  bi.style.color=v>=ALARM?'var(--red)':v>=TRIG?'var(--yel)':'var(--grn)';
}

// ── Data fetch ────────────────────────────────────────────────────────────────
async function refreshShort(){
  try{
    const r=await fetch(`/api/live/${KEY}`);
    if(!r.ok)return;
    const d=await r.json();
    if(d.values&&d.values.length){
      shortData={values:d.values,t_end:d.t_end,secs:d.secs};
    }
    showStalta(d.stalta||0,d.updated);
    if(curZoom<=LONG_THRESH)draw();
  }catch(e){}
}

async function refreshLong(){
  try{
    const r=await fetch(`/api/live_long/${KEY}`);
    if(!r.ok)return;
    const d=await r.json();
    if(d.values&&d.values.length){
      longData={values:d.values,t_end:d.t_end,secs:d.secs};
    }
    showStalta(d.stalta||0,d.updated);
    if(curZoom>LONG_THRESH)draw();
  }catch(e){}
}

// ── Init ──────────────────────────────────────────────────────────────────────
window.addEventListener('load',()=>{
  resize();
  refreshShort();
  setInterval(refreshShort,2000);
  setInterval(refreshLong,30000);
});
window.addEventListener('resize',resize);
new ResizeObserver(()=>requestAnimationFrame(resize)).observe(document.getElementById('cv-wrap'));
</script></body></html>"""

# ── HTML: Event Waveforms Page ─────────────────────────────────────────────────
HTML_EVENT_WF = r"""<!DOCTYPE html>
<html><head><meta charset="UTF-8">
<title>{{ evid }} — Waveforms</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{background:#0d1117;color:#c9d1d9;font-family:'SF Mono',monospace;min-height:100vh;display:flex;flex-direction:column}
#hdr{background:#161b22;border-bottom:1px solid #30363d;padding:10px 16px}
#hdr h1{font-size:14px;font-weight:700;color:#f0f6fc}
#hdr .sub{font-size:11px;color:#8b949e;margin-top:3px}
#content{flex:1;padding:8px}
#status{padding:10px;color:#8b949e;font-size:11px}
.tr-wrap{display:flex;align-items:center;margin-bottom:4px;height:54px}
.tr-lbl{width:100px;flex-shrink:0;font-size:10px;color:#8b949e;text-align:right;padding-right:8px;line-height:1.5}
.tr-lbl .ls{font-weight:700;font-size:11px}.tr-lbl .ld{font-size:9px}
canvas.wc{flex:1;height:52px;background:#0a0f15;border-radius:2px}
.taxis{display:flex;margin-left:108px;color:#6e7681;font-size:9px;justify-content:space-between;padding:0 2px;margin-bottom:8px}
</style></head><body>
<div id="hdr">
  <h1 id="htitle">{{ evid }}</h1>
  <div class="sub" id="hsub">Loading waveforms…</div>
</div>
<div id="content">
  <div id="status">Fetching archive data…</div>
</div>
<script>
const EVID='{{ evid }}';
const SECS=60; // total window: 15s before pick + 45s after

function drawTrace(canvas,data,pickRel,color){
  const w=canvas.offsetWidth||canvas.width,h=canvas.height;
  canvas.width=w;
  const ctx=canvas.getContext('2d');
  ctx.fillStyle='#0a0f15';ctx.fillRect(0,0,w,h);
  ctx.strokeStyle='#1c2333';ctx.lineWidth=0.5;
  ctx.beginPath();ctx.moveTo(0,h/2);ctx.lineTo(w,h/2);ctx.stroke();
  if(!data||data.length<2)return;
  ctx.beginPath();ctx.strokeStyle=color;ctx.lineWidth=1.2;
  const mg=4;
  for(let i=0;i<data.length;i++){
    const x=(i/(data.length-1))*w,y=h/2-data[i]*(h/2-mg);
    i===0?ctx.moveTo(x,y):ctx.lineTo(x,y);
  }ctx.stroke();
  // Pick marker
  const px=(pickRel/SECS)*w;
  ctx.strokeStyle='#f85149';ctx.lineWidth=1.5;
  ctx.beginPath();ctx.moveTo(px,0);ctx.lineTo(px,h);ctx.stroke();
  ctx.fillStyle='#f85149';ctx.font='9px monospace';ctx.fillText('P',px+2,10);
}

async function load(){
  const res=await fetch(`/api/event/${EVID}`);
  const ev=await res.json();
  if(ev.origin_time){
    document.getElementById('htitle').textContent=`${EVID}  —  ${ev.magnitudes&&ev.magnitudes[0]?'M'+ev.magnitudes[0].mag+' '+ev.magnitudes[0].type:''}`;
    document.getElementById('hsub').textContent=`${ev.origin_time.slice(0,19)} UTC  ·  ${parseFloat(ev.lat).toFixed(3)}°N ${Math.abs(parseFloat(ev.lon)).toFixed(3)}°W  ·  ${parseFloat(ev.depth||0).toFixed(1)} km`;
  }
  document.getElementById('status').textContent='Loading archive waveforms…';
  const wres=await fetch(`/api/event/${EVID}/waveforms`);
  const traces=await wres.json();
  const content=document.getElementById('content');
  if(!traces.length){
    content.innerHTML='<div style="padding:20px;color:#8b949e">No archive data available yet. slarchive needs a few minutes to write waveform files.<br><br>Try again shortly or check /Users/OuOu/seiscomp/var/lib/archive/</div>';
    return;
  }
  content.innerHTML='';
  // Time axis
  const ax=document.createElement('div');ax.className='taxis';
  const labels=['-15s','-10s','-5s','P','+10s','+20s','+30s','+45s'];
  labels.forEach(l=>{const s=document.createElement('span');s.textContent=l;ax.appendChild(s);});
  content.appendChild(ax);

  traces.forEach(t=>{
    const wrap=document.createElement('div');wrap.className='tr-wrap';
    const lbl=document.createElement('div');lbl.className='tr-lbl';
    lbl.innerHTML=`<div class="ls" style="color:${t.color}">${t.label}</div><div class="ld">${t.chan}  ·  ${t.dist}°</div>`;
    const cv=document.createElement('canvas');cv.className='wc';cv.height=52;
    wrap.appendChild(lbl);wrap.appendChild(cv);content.appendChild(wrap);
    requestAnimationFrame(()=>drawTrace(cv,t.data,t.pick_rel,t.color));
  });
  document.getElementById('status').textContent='';
}
load();
window.addEventListener('resize',()=>{
  document.querySelectorAll('canvas.wc').forEach(cv=>{
    if(cv._drawn)drawTrace(cv,cv._data,cv._pickRel,cv._color);
  });
});
</script></body></html>"""

# ── HTML: Live Station Map ─────────────────────────────────────────────────────
HTML_LIVEMAP = r"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="UTF-8"><title>Station Live Map — SeisComP</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.socket.io/4.7.2/socket.io.min.js"></script>
<style>
*{box-sizing:border-box;margin:0;padding:0}
html,body{height:100%;background:#0d1117;color:#c9d1d9;font-family:'SF Mono',monospace;display:flex;flex-direction:column;overflow:hidden}
#hdr{background:#161b22;border-bottom:1px solid #30363d;padding:0 14px;display:flex;align-items:center;gap:10px;flex-shrink:0;height:36px}
#hdr h1{font-size:12px;font-weight:700;color:#f0f6fc;white-space:nowrap}
#hdr a{color:#58a6ff;text-decoration:none}#hdr a:hover{text-decoration:underline}
.pill{padding:2px 9px;border-radius:20px;font-size:10px;font-weight:600}
.pg{background:#0d2318;color:#3fb950;border:1px solid #1a4731}
.py{background:#2d1f00;color:#d29922;border:1px solid #5a3e00}
.pr{background:#2d0f0f;color:#f85149;border:1px solid #5a1a1a}
#utc-lbl{margin-left:auto;color:#6e7681;font-size:10px;white-space:nowrap}
#map{flex:1;position:relative}
.leaflet-container{background:#0a0e17!important;font-family:'SF Mono',monospace}
.leaflet-control-zoom{border:1px solid #30363d!important}
.leaflet-control-zoom a{background:#161b22!important;color:#c9d1d9!important;border-color:#30363d!important}
.leaflet-tooltip{background:#161b22;border:1px solid #30363d;color:#c9d1d9;font-family:'SF Mono',monospace;font-size:11px;padding:4px 8px;border-radius:4px;box-shadow:0 2px 8px rgba(0,0,0,.5)}
.leaflet-tooltip-top:before{border-top-color:#30363d}
/* Legend */
#legend{position:absolute;bottom:28px;right:12px;z-index:1000;background:rgba(22,27,34,.94);border:1px solid #30363d;border-radius:6px;padding:10px 14px;font-size:10px;backdrop-filter:blur(6px)}
#legend .lt{font-size:11px;font-weight:700;color:#f0f6fc;margin-bottom:8px}
.li{display:flex;align-items:center;gap:7px;margin-bottom:4px;color:#c9d1d9}
.lisym{width:14px;height:12px;display:flex;align-items:center;justify-content:center}
.ldiv{margin:6px 0;border-top:1px solid #30363d}
/* Status pill */
#sb{position:absolute;top:8px;left:50%;transform:translateX(-50%);z-index:1000;background:rgba(22,27,34,.92);border:1px solid #30363d;border-radius:20px;padding:3px 16px;font-size:10px;color:#8b949e;pointer-events:none;white-space:nowrap;transition:opacity .5s}
/* Stats dashboard */
#statpanel{position:absolute;top:10px;right:12px;z-index:1100;width:210px;background:rgba(13,17,23,.93);border:1px solid #30363d;border-radius:8px;padding:10px 12px;font-size:10px;font-family:'SF Mono',monospace;backdrop-filter:blur(8px);color:#c9d1d9;pointer-events:none}
#statpanel .sp-title{font-size:10px;font-weight:700;color:#f0f6fc;letter-spacing:.06em;text-transform:uppercase;margin-bottom:8px;display:flex;align-items:center;gap:5px}
#statpanel .sp-row{display:flex;justify-content:space-between;align-items:center;margin-bottom:5px}
#statpanel .sp-lbl{color:#6e7681}
#statpanel .sp-val{font-weight:700;color:#f0f6fc;text-align:right}
#statpanel .sp-div{border-top:1px solid #21262d;margin:7px 0}
#statpanel .sp-sub{font-size:9px;color:#6e7681;margin-bottom:4px;text-transform:uppercase;letter-spacing:.05em}
#statpanel .sp-ev{display:flex;justify-content:space-between;margin-bottom:3px;font-size:10px}
#statpanel .sp-ev-mag{font-weight:700}
#statpanel .sp-ev-info{color:#8b949e;text-align:right;font-size:9px}
.dot{width:6px;height:6px;border-radius:50%;display:inline-block;flex-shrink:0}
.dot-grn{background:#3fb950}.dot-yel{background:#d4b800}.dot-red{background:#f85149}.dot-dim{background:#444c56}
/* Prelim detection cards */
#prelim-stack{position:absolute;bottom:36px;left:12px;z-index:1100;display:flex;flex-direction:column-reverse;gap:8px;pointer-events:none;max-width:300px}
.prelim-card{background:rgba(13,17,23,.95);border-radius:8px;padding:10px 12px;font-size:11px;color:#c9d1d9;pointer-events:all;backdrop-filter:blur(8px);animation:card-in .3s ease;box-shadow:0 4px 20px rgba(0,0,0,.6);border-left:3px solid #555}
.prelim-card .pc-mag{font-size:22px;font-weight:700;line-height:1;margin-bottom:2px}
.prelim-card .pc-type{font-size:9px;font-weight:600;letter-spacing:.08em;text-transform:uppercase;opacity:.7;margin-bottom:6px}
.prelim-card .pc-city{font-size:12px;font-weight:600;color:#f0f6fc;margin-bottom:4px}
.prelim-card .pc-row{color:#8b949e;font-size:10px;margin-bottom:2px}
.prelim-card .pc-close{position:absolute;top:6px;right:8px;cursor:pointer;color:#8b949e;font-size:13px;line-height:1}
.prelim-card .pc-close:hover{color:#f0f6fc}
@keyframes card-in{from{opacity:0;transform:translateX(-18px)}to{opacity:1;transform:translateX(0)}}
@keyframes prelim-pulse{0%,100%{box-shadow:0 0 8px currentColor}50%{box-shadow:0 0 22px currentColor,0 0 40px currentColor}}
</style>
</head><body>
<div id="hdr">
  <h1><a href="/">&#8592; Dashboard</a> &nbsp;/&nbsp; Station Live Map</h1>
  <span class="pill pg" id="s-live">— / —</span>
  <span class="pill py" id="s-trig" style="display:none">0 triggered</span>
  <span class="pill pr" id="s-alrm" style="display:none">0 alarmed</span>
  <button id="pwave-btn" onclick="togglePwaves()"
    style="margin-left:6px;padding:2px 10px;border-radius:10px;font-size:10px;font-weight:600;cursor:pointer;
           background:#0d2318;color:#3fb950;border:1px solid #1a4731;font-family:inherit">
    P-waves ON</button>
  <button id="inactive-btn" onclick="toggleInactive()"
    style="margin-left:4px;padding:2px 10px;border-radius:10px;font-size:10px;font-weight:600;cursor:pointer;
           background:#161b22;color:#8b949e;border:1px solid #30363d;font-family:inherit">
    Show inactive</button>
  <button id="unmonitored-btn" onclick="toggleUnmonitored()"
    style="margin-left:4px;padding:2px 10px;border-radius:10px;font-size:10px;font-weight:600;cursor:pointer;
           background:#161b22;color:#8b949e;border:1px solid #30363d;font-family:inherit">
    Show unmonitored</button>
  <button id="epicenter-btn" onclick="toggleEpicenters()"
    style="margin-left:4px;padding:2px 10px;border-radius:10px;font-size:10px;font-weight:600;cursor:pointer;
           background:#0d2318;color:#3fb950;border:1px solid #1a4731;font-family:inherit">
    Epicenters ON</button>
  <span id="utc-lbl">Last updated: connecting…</span>
</div>
<div id="map">
  <div id="sb">Loading station inventory…</div>
  <div id="prelim-stack"></div>
  <div id="statpanel">
    <div class="sp-title"><span class="dot dot-grn" id="sp-dot"></span>SeisComP Monitor</div>
    <div class="sp-row"><span class="sp-lbl">Live stations</span><span class="sp-val" id="sp-live">—</span></div>
    <div class="sp-row"><span class="sp-lbl">Triggered</span><span class="sp-val" id="sp-trig" style="color:#d4b800">—</span></div>
    <div class="sp-row"><span class="sp-lbl">Alarm</span><span class="sp-val" id="sp-alrm" style="color:#f85149">—</span></div>
    <div class="sp-row"><span class="sp-lbl">Max STA/LTA</span><span class="sp-val" id="sp-max">—</span></div>
    <div class="sp-row"><span class="sp-lbl">UTC</span><span class="sp-val" id="sp-utc" style="color:#8b949e">—</span></div>
    <div class="sp-div"></div>
    <div class="sp-sub">Recent events (SeisComP)</div>
    <div id="sp-evlist"><span style="color:#6e7681">loading…</span></div>
    <div class="sp-div"></div>
    <div class="sp-sub">Latest detection</div>
    <div id="sp-prelim"><span style="color:#6e7681">none yet</span></div>
  </div>
  <div id="legend">
    <div class="lt">STA/LTA</div>
    <div class="li"><div class="lisym"><svg width="10" height="10"><rect width="10" height="10" rx="2" fill="#555"/></svg></div>Not monitored</div>
    <div class="li"><div class="lisym"><svg width="14" height="12"><polygon points="7,1 13,11 1,11" fill="#3a8fd4" stroke="#fff" stroke-width="0.7"/></svg></div>Quiet (&lt; 1.5)</div>
    <div class="li"><div class="lisym"><svg width="14" height="12"><polygon points="7,1 13,11 1,11" fill="#2ab5a0" stroke="#fff" stroke-width="0.7"/></svg></div>Elevated (1.5 – 3)</div>
    <div class="li"><div class="lisym"><svg width="14" height="12"><polygon points="7,1 13,11 1,11" fill="#d4b800" stroke="#fff" stroke-width="0.7"/></svg></div>Triggered (3 – 6)</div>
    <div class="li"><div class="lisym"><svg width="14" height="12"><polygon points="7,1 13,11 1,11" fill="#f85149" stroke="#fff" stroke-width="0.7"/></svg></div>Alarm (&ge; 6)</div>
    <div class="ldiv"></div>
    <div class="li"><div class="lisym"><svg width="10" height="10"><circle cx="5" cy="5" r="4" fill="#f85149" opacity="0.8"/></svg></div>Event (48h)</div>
    <div class="li"><div class="lisym"><svg width="14" height="14"><circle cx="7" cy="7" r="6" fill="none" stroke="#ff9500" stroke-width="2" stroke-dasharray="3,2"/></svg></div>PRELIM detection</div>
  </div>
</div>
<script>
const socket=io({transports:['websocket','polling'],upgrade:true,reconnectionDelay:1000,reconnectionAttempts:Infinity});
const ALARM=4,TRIG=2.0;

// ── Pacific Time helpers ───────────────────────────────────────────────────────
const _ptFmt = new Intl.DateTimeFormat('en-US', {
  timeZone:'America/Los_Angeles', year:'numeric', month:'2-digit', day:'2-digit',
  hour:'2-digit', minute:'2-digit', second:'2-digit', hour12:false
});
function utcToPT(utcStr) {
  // utcStr: "YYYY-MM-DD HH:MM:SS" (UTC)
  const d = new Date(utcStr.replace(' ','T')+'Z');
  if (isNaN(d)) return utcStr;
  const parts = _ptFmt.formatToParts(d);
  const p = {};
  parts.forEach(x => { p[x.type]=x.value; });
  // Determine PDT vs PST by checking if offset is -7 or -8
  const offsetH = -d.getTimezoneOffset ? 0 :
    Math.round((d - new Date(d.toLocaleString('en-US',{timeZone:'America/Los_Angeles'})))/3600000);
  const lbl = new Date().toLocaleString('en-US',{timeZone:'America/Los_Angeles',timeZoneName:'short'})
               .match(/(P[DS]T)/)?.[1] || 'PT';
  return `${p.year}-${p.month}-${p.day} ${p.hour}:${p.minute}:${p.second} ${lbl}`;
}
function utcHMStoPT(hms) {
  // hms: "HH:MM:SS" assumed today UTC
  const now = new Date();
  const [h,m,s] = hms.split(':').map(Number);
  const d = new Date(Date.UTC(now.getUTCFullYear(),now.getUTCMonth(),now.getUTCDate(),h,m,s));
  return utcToPT(d.toISOString().replace('T',' ').slice(0,19));
}

// ── Map ────────────────────────────────────────────────────────────────────────
const map=L.map('map',{zoomControl:true,attributionControl:false}).setView([37,-118],6);
L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png',
  {maxZoom:18,opacity:0.9}).addTo(map);
L.control.attribution({prefix:'<span style="color:#30363d">CartoDB | OSM</span>'}).addTo(map);

// ── Station markers ────────────────────────────────────────────────────────────
const _staMarkers={};   // "NET.STA" -> marker
let _totalStas=0,_monitoredStas=0;

function staltaColor(stalta,monitored){
  if(!monitored||stalta==null)return'#555555';
  if(stalta>=ALARM)return'#f85149';  // alarm   ≥4.0
  if(stalta>=TRIG) return'#d4b800';  // trigger ≥2.0
  if(stalta>=1.2)  return'#2ab5a0';  // elevated 1.2–2.0
  return'#3a8fd4';                    // normal
}

function makeTriIcon(col,size,zOff){
  const h=Math.round(size*0.866);
  return L.divIcon({
    className:'',
    html:`<svg width="${size}" height="${h}" viewBox="0 0 ${size} ${h}" style="display:block;filter:drop-shadow(0 1px 2px rgba(0,0,0,.5))">
      <polygon points="${size/2},1 ${size-1},${h-1} 1,${h-1}"
        fill="${col}" stroke="rgba(255,255,255,0.7)" stroke-width="0.8" opacity="0.95"/>
    </svg>`,
    iconSize:[size,h],
    iconAnchor:[size/2,h],
  });
}

// ── Inactive-station toggle ───────────────────────────────────────────────────
let _showInactive = true;
function toggleInactive() {
  _showInactive = !_showInactive;
  const btn = document.getElementById('inactive-btn');
  if (_showInactive) {
    btn.textContent = 'Show inactive';
    btn.style.background  = '#161b22'; btn.style.color = '#8b949e'; btn.style.borderColor = '#30363d';
  } else {
    btn.textContent = 'Hide inactive';
    btn.style.background  = '#0d1f2d'; btn.style.color = '#58a6ff'; btn.style.borderColor = '#1a4a6b';
  }
  Object.entries(_staMarkers).forEach(([key, entry]) => {
    if (entry.active === false) {
      if (_showInactive) entry.marker.addTo(map);
      else try { entry.marker.remove(); } catch(e){}
    }
  });
}

// ── Unmonitored-station toggle ────────────────────────────────────────────────
let _showUnmonitored = true;
function toggleUnmonitored() {
  _showUnmonitored = !_showUnmonitored;
  const btn = document.getElementById('unmonitored-btn');
  if (_showUnmonitored) {
    btn.textContent = 'Show unmonitored';
    btn.style.background  = '#161b22'; btn.style.color = '#8b949e'; btn.style.borderColor = '#30363d';
  } else {
    btn.textContent = 'Hide unmonitored';
    btn.style.background  = '#0d1f2d'; btn.style.color = '#58a6ff'; btn.style.borderColor = '#1a4a6b';
  }
  Object.entries(_staMarkers).forEach(([key, entry]) => {
    if (!entry.monitored) {
      if (_showUnmonitored) entry.marker.addTo(map);
      else try { entry.marker.remove(); } catch(e){}
    }
  });
}

// Fetch all station positions once
fetch('/api/stations').then(r=>r.json()).then(stas=>{
  _totalStas=stas.length;
  _monitoredStas=stas.filter(s=>s.monitored).length;
  stas.forEach(s=>{
    if(!s.lat||!s.lon)return;
    const key=`${s.net}.${s.sta}`;
    const isMonitored=s.monitored;
    const col=staltaColor(s.stalta,isMonitored);
    const sz=isMonitored?14:8;
    const icon=makeTriIcon(col,sz);
    const m=L.marker([s.lat,s.lon],{icon,zIndexOffset:isMonitored?200:0});
    const tip=isMonitored
      ?`<b>${key}</b><br>STA/LTA: ${s.stalta.toFixed(2)}<br>${s.lat.toFixed(3)}°N ${Math.abs(s.lon).toFixed(3)}°W`
      :`<b>${key}</b> (not monitored)<br>${s.lat.toFixed(3)}°N ${Math.abs(s.lon).toFixed(3)}°W`;
    m.bindTooltip(tip,{sticky:true});
    m.addTo(map);
    _staMarkers[key]={marker:m,monitored:isMonitored,active:null};
  });
  const sb=document.getElementById('sb');
  sb.textContent=`${_totalStas} stations loaded (${_monitoredStas} monitored)`;
  setTimeout(()=>{sb.style.opacity='0';setTimeout(()=>sb.style.display='none',600);},3000);
}).catch(()=>{
  document.getElementById('sb').textContent='Failed to load stations';
});

// ── Real-time STA/LTA from SocketIO ───────────────────────────────────────────
socket.on('stalta',msg=>{
  const rows=msg.rows||[];
  // Update header pills
  document.getElementById('utc-lbl').textContent='Last updated: '+utcToPT(msg.utc);
  document.getElementById('s-live').textContent=`${msg.active}/${msg.total} live`;
  const bt=document.getElementById('s-trig'),ba=document.getElementById('s-alrm');
  bt.style.display=msg.triggered>0?'':'none';
  ba.style.display=msg.alarmed>0?'':'none';
  if(msg.triggered>0)bt.textContent=`${msg.triggered} triggered`;
  if(msg.alarmed>0)  ba.textContent=`${msg.alarmed} alarmed`;

  // Update marker colors + activity state
  // Helper: "HH:MM:SS" UTC → "Xs ago" / "Xm Xs ago" / "Xh Xm ago"
  function agoStr(hms){
    if(!hms) return null;
    const now=new Date();
    const [h,m,s]=hms.split(':').map(Number);
    const ref=new Date(Date.UTC(now.getUTCFullYear(),now.getUTCMonth(),now.getUTCDate(),h,m,s));
    let diff=Math.round((now-ref)/1000);
    if(diff<0) diff+=86400; // crossed midnight
    if(diff<60) return `${diff}s ago`;
    if(diff<3600) return `${Math.floor(diff/60)}m ${diff%60}s ago`;
    return `${Math.floor(diff/3600)}h ${Math.floor((diff%3600)/60)}m ago`;
  }

  rows.forEach(r=>{
    const key=`${r.net}.${r.sta}`;
    const entry=_staMarkers[key];
    if(!entry)return;

    const wasActive = entry.active;
    entry.active = r.active !== false;  // default true if field absent

    // Handle hide/show based on active state change
    if(entry.active !== wasActive){
      if(!entry.active && !_showInactive){
        try { entry.marker.remove(); } catch(e){}
      } else if(entry.active && wasActive === false){
        try { entry.marker.addTo(map); } catch(e){}
      }
    }

    // Dim inactive markers (grey, small)
    if(!entry.active){
      entry.marker.setIcon(makeTriIcon('#444c56', 8));
      const updStr = r.updated ? utcHMStoPT(r.updated) : '—';
      const ago    = r.updated ? agoStr(r.updated) : null;
      entry.marker.setTooltipContent(
        `<b>${key}</b> — <span style="color:#6e7681">inactive</span><br>`+
        `Last updated: ${updStr}`+(ago?`<br><span style="color:#6e7681">${ago}</span>`:'')
      );
      entry.monitored = true;
      return;
    }

    const col=staltaColor(r.stalta,true);
    entry.marker.setIcon(makeTriIcon(col,14));
    const v=r.stalta?r.stalta.toFixed(2):'—';
    const st=r.stalta>=ALARM?'⚠ ALARM':r.stalta>=TRIG?'▲ TRIGGERED':'quiet';
    const updStr=r.updated?utcHMStoPT(r.updated):'—';
    const ago=r.updated?agoStr(r.updated):null;
    const agoLine=ago?`<br><span style="color:#8b949e">${ago}</span>`:'';
    if(r.lat&&r.lon){
      entry.marker.setTooltipContent(
        `<b>${key}</b> — <span style="color:${col}">${st}</span><br>`+
        `STA/LTA: <b>${v}</b><br>`+
        `${r.lat.toFixed(3)}°N ${Math.abs(r.lon).toFixed(3)}°W<br>`+
        `Last updated: ${updStr}${agoLine}`
      );
    } else {
      entry.marker.setTooltipContent(
        `<b>${key}</b> — <span style="color:${col}">${st}</span><br>`+
        `STA/LTA: <b>${v}</b><br>`+
        `Last updated: ${updStr}${agoLine}`
      );
    }
    entry.monitored=true;
  });
});

// ── GlobalQuake-style P-wave ring visualization ──────────────────────────────
const P_SPEED_KMpS = 6.0;          // km/s P-wave velocity
const triggers = {};                 // key -> {ts, lat, lon, circle, ring}
const KM_PER_DEG = 111.2;

// Epicenter (PRELIM) toggle
let _showEpicenters = true;
function toggleEpicenters() {
  _showEpicenters = !_showEpicenters;
  const btn = document.getElementById('epicenter-btn');
  if (_showEpicenters) {
    btn.textContent = 'Epicenters ON';
    btn.style.background  = '#0d2318';
    btn.style.color       = '#3fb950';
    btn.style.borderColor = '#1a4731';
  } else {
    btn.textContent = 'Epicenters OFF';
    btn.style.background  = '#2d1f00';
    btn.style.color       = '#8b949e';
    btn.style.borderColor = '#5a3e00';
  }
  Object.values(_prelims).forEach(d => {
    try {
      if (_showEpicenters) d.marker.addTo(map);
      else d.marker.remove();
    } catch(e){}
  });
}

// P-wave toggle
let _showPwaves = true;
function togglePwaves() {
  _showPwaves = !_showPwaves;
  const btn = document.getElementById('pwave-btn');
  if (_showPwaves) {
    btn.textContent = 'P-waves ON';
    btn.style.background   = '#0d2318';
    btn.style.color        = '#3fb950';
    btn.style.borderColor  = '#1a4731';
    // Re-add all existing rings to map
    Object.values(triggers).forEach(t => {
      if (t.ring) { try { t.ring.addTo(map); } catch(e){} }
    });
  } else {
    btn.textContent = 'P-waves OFF';
    btn.style.background   = '#2d1f00';
    btn.style.color        = '#8b949e';
    btn.style.borderColor  = '#5a3e00';
    // Remove all existing rings from map
    Object.values(triggers).forEach(t => {
      if (t.ring) { try { t.ring.remove(); } catch(e){} }
    });
  }
}

socket.on('trigger_new', msg => {
  (msg.triggers || []).forEach(t => {
    const k = t.net + '.' + t.sta;
    // Remove old ring if present
    if (triggers[k] && triggers[k].ring) {
      try { triggers[k].ring.remove(); } catch(e){}
    }
    // P-wave origin ring (expands from station location)
    const ring = L.circle([t.lat, t.lon], {
      radius: 1,
      color: '#ff9500',
      weight: 2,
      fill: false,
      opacity: 0.8,
      dashArray: '4,4',
    });
    if (_showPwaves) ring.addTo(map);
    triggers[k] = { ts: t.ts, lat: t.lat, lon: t.lon, ring, stalta: t.stalta };
    // Flash the station marker orange
    const entry = _staMarkers[k];
    if (entry) {
      entry.marker.setIcon(makeTriIcon('#ff9500', 14));
    }
  });
});

// Animate P-wave rings — smooth 60 fps via requestAnimationFrame
(function _pwaveLoop() {
  const now = Date.now() / 1000;
  Object.entries(triggers).forEach(([k, t]) => {
    const age = now - t.ts;
    if (age > 120) {
      try { t.ring.remove(); } catch(e){}
      delete triggers[k];
      return;
    }
    if (!_showPwaves) return;
    const radiusM = age * P_SPEED_KMpS * 1000;
    const opacity = Math.max(0, 1 - age / 120);
    t.ring.setRadius(radiusM);
    t.ring.setStyle({ opacity });
  });
  requestAnimationFrame(_pwaveLoop);
})();

// ── Event markers ──────────────────────────────────────────────────────────────
let _evMarkers=[];
function magColor(m){
  if(m>=5)return'#f85149';if(m>=4)return'#ff9500';
  if(m>=3)return'#d29922';if(m>=2)return'#3fb950';return'#58a6ff';
}
function magR(m){return Math.max(4,Math.min(18,(m||0)*3.5));}

socket.on('events',msg=>{
  _evMarkers.forEach(m=>m.remove());_evMarkers=[];
  const cutoff=Date.now()-48*3600*1000;
  (msg.events||[]).forEach(ev=>{
    if(!ev.lat||!ev.lon)return;
    const d=new Date((ev.origin_time||'').replace(' ','T')+'Z');
    if(d.getTime()<cutoff)return;
    const m=parseFloat(ev.mag)||0;
    const circ=L.circleMarker([ev.lat,ev.lon],{
      radius:magR(m),fillColor:magColor(m),color:'#fff',
      weight:1,fillOpacity:0.8,zIndexOffset:500,
    });
    const ot=(ev.origin_time||'').slice(0,19)+' UTC';
    circ.bindTooltip(
      `<b>M${ev.mag||'?'}</b> — ${ot}<br>`+
      `${parseFloat(ev.lat).toFixed(3)}°N ${Math.abs(parseFloat(ev.lon)).toFixed(3)}°W  ·  ${ev.depth} km`,
      {sticky:true}
    );
    circ.on('click',()=>window.opener?window.opener.openEv&&window.opener.openEv(ev.evid):null);
    circ.addTo(map);_evMarkers.push(circ);
  });
});

// ── PRELIM detection visualization ───────────────────────────────────────────
const _prelims = {};   // id → {marker, ts}
let _prelimId  = 0;

function _prelimColor(mag) {
  return mag>=5?'#f85149':mag>=4?'#ff6b35':mag>=3?'#d4b800':mag>=2?'#3fb950':'#58a6ff';
}
function _makePrelimIcon(mag, isTele) {
  const col  = _prelimColor(mag);
  const sz   = Math.max(26, Math.min(58, Math.round(mag * 9 + 10)));
  const fs   = Math.max(9,  Math.round(sz * 0.32));
  const label= `M${mag.toFixed(1)}`;
  // Large quake gets warning indicator, teleseism gets T prefix, otherwise blank
  const prefix = mag >= 4 ? '<span style="font-size:' + Math.max(7,fs-2) + 'px;opacity:.9">&#9888;</span>' :
                 isTele   ? '<span style="font-size:' + Math.max(7,fs-2) + 'px;opacity:.75">T</span>' : '';
  return L.divIcon({
    className:'',
    html:`<div style="
      width:${sz}px;height:${sz}px;border-radius:50%;
      background:${col};border:2px solid rgba(255,255,255,0.85);
      display:flex;flex-direction:column;align-items:center;justify-content:center;
      font-family:'SF Mono',monospace;font-weight:700;color:#fff;
      box-shadow:0 0 14px ${col},0 0 4px rgba(0,0,0,.8);
      animation:prelim-pulse 0.8s ease-in-out 4;
      line-height:1;gap:1px;cursor:pointer">
        ${prefix}
        <span style="font-size:${fs}px">${label}</span>
    </div>`,
    iconSize:[sz,sz], iconAnchor:[sz/2,sz/2],
  });
}

function _showPrelimCard(ev, id) {
  const mag    = ev.mag_est;
  const isTele = ev.teleseism||false;
  const col    = _prelimColor(mag);
  const lat    = ev.lat, lon = ev.lon;
  const ns     = lat>=0?'N':'S', ew = lon<=0?'W':'E';
  const loc    = `${Math.abs(lat).toFixed(3)}°${ns}  ${Math.abs(lon).toFixed(3)}°${ew}`;
  const typeStr= isTele ? 'Teleseism' : (mag>=4 ? '&#9888; Preliminary' : 'Preliminary');
  const mmiStr = ev.mmi_str || '?';
  const mmiDesc = ev.mmi_desc || '';
  const mmiCity = ev.mmi_city_str || '?';
  const rms    = ev.rms_sec != null ? ` · RMS ${ev.rms_sec}s` : '';

  const card = document.createElement('div');
  card.className = 'prelim-card';
  card.id = `pc-${id}`;
  card.style.borderLeftColor = col;
  card.style.position = 'relative';
  card.innerHTML = `
    <span class="pc-close" onclick="document.getElementById('pc-${id}')?.remove()">✕</span>
    <div class="pc-mag" style="color:${col}">${mag.toFixed(1)}</div>
    <div class="pc-type" style="color:${col}">${typeStr}</div>
    ${ev.city ? `<div class="pc-city">${ev.city}</div>` : ''}
    <div class="pc-row">${loc}</div>
    <div class="pc-row">${utcToPT(ev.ot_str||'')}</div>
    <div class="pc-row">${ev.n_stations} stations · max STA/LTA ${ev.max_stalta}${rms}</div>
    <div class="pc-row">Spread ${ev.spread_km||0} km · Depth ~${ev.depth_est||10} km</div>
    <div class="pc-row" style="color:#f0c040;font-weight:600">Max MMI ${mmiStr} — ${mmiDesc}${ev.city ? ' · ' + mmiCity + ' at ' + ev.city.split(',')[0] : ''}</div>
    <div class="pc-row" style="color:#6e7681;font-size:9px;margin-top:4px">${ev.stations||''}</div>
  `;
  // Pan-to button
  card.addEventListener('click', e => {
    if(e.target.classList.contains('pc-close')) return;
    map.flyTo([lat, lon], Math.min(map.getZoom(), isTele?5:8), {animate:true, duration:1.2});
  });
  const stack = document.getElementById('prelim-stack');
  stack.prepend(card);
  // Keep at most 4 cards
  while(stack.children.length > 4) stack.lastChild.remove();
  // Auto-dismiss after 60s
  setTimeout(() => { try { card.remove(); } catch(e){} }, 60000);
}

// Update stat panel "Latest detection" row
function _updateStatPanelPrelim(ev) {
  const mag  = ev.mag_est||0;
  const col  = mag>=5?'#f85149':mag>=4?'#ff6b35':mag>=3?'#d4b800':mag>=2?'#3fb950':'#58a6ff';
  const city = ev.city ? ev.city.slice(0,24) : `${(ev.lat||0).toFixed(2)}N ${Math.abs(ev.lon||0).toFixed(2)}W`;
  const isTele = ev.teleseism||false;
  const mmiStr = ev.mmi_str || '?';
  const mmiDesc = (ev.mmi_desc||'').slice(0,12);
  const prefix = isTele ? 'T ' : (mag>=4 ? '&#9888; ' : '');
  const el = document.getElementById('sp-prelim');
  if (el) el.innerHTML =
    `<div class="sp-ev">
      <span class="sp-ev-mag" style="color:${col}">${prefix}M${mag.toFixed(1)}</span>
      <span class="sp-ev-info">${city}<br>${utcToPT(ev.ot_str||'').replace(/\d{4}-\d{2}-\d{2} /,'')} · ${ev.n_stations} sta<br><span style="color:#f0c040">MMI ${mmiStr} — ${mmiDesc}</span></span>
    </div>`;
}

socket.on('preliminary', ev => {
  const id  = ++_prelimId;
  const lat = ev.lat, lon = ev.lon;
  const mag = ev.mag_est;
  const col = _prelimColor(mag);
  const isTele = ev.teleseism||false;
  const now = Date.now()/1000;

  // Epicenter marker
  const marker = L.marker([lat, lon], {
    icon: _makePrelimIcon(mag, isTele),
    zIndexOffset: 1200,
  });
  if (_showEpicenters) marker.addTo(map);

  const mmiStr  = ev.mmi_str  || '?';
  const mmiDesc = ev.mmi_desc || '';
  const typeLabel = isTele ? 'TELESEISM' : (mag>=4 ? '&#9888; PRELIMINARY' : 'PRELIMINARY');
  marker.bindTooltip(
    `<b>${typeLabel} — Est. M${mag.toFixed(1)}</b><br>`+
    `Near: ${ev.city||'—'}<br>`+
    `${lat.toFixed(3)}°N  ${Math.abs(lon).toFixed(3)}°W<br>`+
    `${ev.n_stations} stations · max STA/LTA ${ev.max_stalta}<br>`+
    `Spread: ${ev.spread_km||0} km<br>`+
    `Max MMI: <b>${mmiStr}</b> — ${mmiDesc}<br>`+
    `${utcToPT(ev.ot_str||'')}`,
    {sticky:true}
  );

  // Show info card + update stat panel
  _showPrelimCard(ev, id);
  _updateStatPanelPrelim(ev);

  // Pan map for significant events
  if (mag >= 3.0 || isTele) {
    map.flyTo([lat, lon], Math.min(map.getZoom(), isTele?5:7), {animate:true, duration:1.5});
  }

  _prelims[id] = {marker, ts: now};
});

// Replay historical preliminary detections sent on connect
socket.on('preliminary_history', msg => {
  const now = Date.now()/1000;
  (msg.events||[]).forEach(ev => {
    // Parse event timestamp — skip if too old (>10 min)
    let evTs = now;
    try {
      const d = new Date(ev.ot_str.replace(' ','T')+'Z');
      if (!isNaN(d)) evTs = d.getTime()/1000;
    } catch(e) {}
    const age = now - evTs;
    if (age > 600) return; // older than 10 min — skip

    const id  = ++_prelimId;
    const lat = ev.lat, lon = ev.lon;
    const mag = ev.mag_est;
    const col = _prelimColor(mag);
    const isTele = ev.teleseism||false;

    const marker = L.marker([lat, lon], {
      icon: _makePrelimIcon(mag, isTele),
      zIndexOffset: 1200,
      opacity: Math.max(0.3, 1 - age/600),
    });
    if (_showEpicenters) marker.addTo(map);

    const mmiStr  = ev.mmi_str  || '?';
    const mmiDesc = ev.mmi_desc || '';
    const typeLabel2 = isTele ? 'TELESEISM' : (mag>=4 ? '&#9888; PRELIMINARY' : 'PRELIMINARY');
    marker.bindTooltip(
      `<b>${typeLabel2} — Est. M${mag.toFixed(1)}</b><br>`+
      `Near: ${ev.city||'—'}<br>`+
      `${lat.toFixed(3)}°N  ${Math.abs(lon).toFixed(3)}°W<br>`+
      `${ev.n_stations} stations · max STA/LTA ${ev.max_stalta}<br>`+
      `Spread: ${ev.spread_km||0} km<br>`+
      `Max MMI: <b>${mmiStr}</b> — ${mmiDesc}<br>`+
      `${utcToPT(ev.ot_str||'')}`,
      {sticky:true}
    );

    _prelims[id] = {marker, ts: evTs};
  });
  // Update stat panel with most recent detection (first in array = newest)
  if ((msg.events||[]).length > 0) _updateStatPanelPrelim(msg.events[0]);
});

// Fade out and remove old PRELIM markers after 10 minutes
(function _prelimLoop() {
  const now    = Date.now()/1000;
  const MAX_AGE = 600;
  Object.entries(_prelims).forEach(([id, d]) => {
    const age = now - d.ts;
    if (age > MAX_AGE) {
      try { d.marker.remove(); } catch(e){}
      delete _prelims[id];
    }
  });
  requestAnimationFrame(_prelimLoop);
})();

// ── Stats dashboard panel ─────────────────────────────────────────────────────
let _spMaxSta='—', _spMaxVal=0;

socket.on('stalta', msg => {
  // Max STA/LTA
  _spMaxVal = 0; _spMaxSta = '—';
  (msg.rows||[]).forEach(r => {
    if((r.stalta||0) > _spMaxVal){ _spMaxVal=r.stalta; _spMaxSta=r.key||r.sta; }
  });
  const maxCol = _spMaxVal>=ALARM?'#f85149':_spMaxVal>=TRIG?'#d4b800':_spMaxVal>=1.2?'#2ab5a0':'#c9d1d9';
  const trig = msg.triggered||0, alrm = msg.alarmed||0;
  document.getElementById('sp-live').textContent = `${msg.active||0} / ${msg.total||0}`;
  document.getElementById('sp-trig').textContent = trig||'0';
  document.getElementById('sp-alrm').textContent = alrm||'0';
  document.getElementById('sp-max').innerHTML  =
    `<span style="color:${maxCol}">${_spMaxVal.toFixed(2)}</span> <span style="color:#6e7681;font-size:9px">${_spMaxSta}</span>`;
  document.getElementById('sp-utc').textContent = utcToPT(msg.utc||'').replace(/\d{4}-\d{2}-\d{2} /,'');
  // Status dot
  const dot = document.getElementById('sp-dot');
  dot.className = 'dot ' + (alrm>0?'dot-red':trig>0?'dot-yel':'dot-grn');
});

socket.on('events', msg => {
  const evs = (msg.events||[]).slice(0,4);
  const el  = document.getElementById('sp-evlist');
  if(!evs.length){ el.innerHTML='<span style="color:#6e7681">none in 48h</span>'; return; }
  el.innerHTML = evs.map(ev => {
    const m   = parseFloat(ev.mag)||0;
    const col = m>=5?'#f85149':m>=4?'#ff9500':m>=3?'#d29922':m>=2?'#3fb950':'#58a6ff';
    const ot  = utcToPT((ev.origin_time||'').slice(0,19)).replace(/\d{4}-\d{2}-\d{2} /,'');
    const loc = ev.place ? ev.place.slice(0,22) : `${parseFloat(ev.lat).toFixed(1)}N ${Math.abs(parseFloat(ev.lon)).toFixed(1)}W`;
    return `<div class="sp-ev">
      <span class="sp-ev-mag" style="color:${col}">M${m.toFixed(1)}</span>
      <span class="sp-ev-info">${loc}<br>${ot}</span>
    </div>`;
  }).join('');
});

// (stat panel preliminary update is handled by _updateStatPanelPrelim above)

// Connect
socket.on('connect',()=>{
  document.getElementById('utc-lbl').textContent='Last updated: connected';
  document.getElementById('utc-lbl').style.color='var(--grn)';
  _lastStaltaTs = Date.now();
});
socket.on('disconnect',()=>{
  document.getElementById('utc-lbl').textContent='Last updated: disconnected — reconnecting…';
  document.getElementById('utc-lbl').style.color='var(--red)';
});

// ── Watchdog: force reconnect if no stalta for 20 s ──────────────────────────
let _lastStaltaTs = Date.now();
socket.on('stalta', () => { _lastStaltaTs = Date.now(); });
setInterval(() => {
  const stale = (Date.now() - _lastStaltaTs) / 1000;
  if (stale > 20 && socket.connected) {
    console.warn('[watchdog] stalta stale for', stale.toFixed(0), 's — reconnecting');
    document.getElementById('utc-lbl').textContent = 'Last updated: stale — reconnecting…';
    document.getElementById('utc-lbl').style.color = 'var(--yel)';
    socket.disconnect();
    setTimeout(() => socket.connect(), 500);
  }
}, 5000);

// ── Preliminary real-time detections ─────────────────────────────────────────
function addPrelimRow(ev){
  const list = document.getElementById('prelim-list');
  if(!list) return;
  const mag = ev.mag_est;
  const col = mag>=3?'#f85149':mag>=2?'#ff9500':mag>=1.5?'#d29922':'#58a6ff';
  const row = document.createElement('div');
  row.style.cssText=`border-left:3px solid ${col};background:#0a0f15;padding:6px 10px;margin-bottom:3px;border-radius:2px;font-size:11px`;
  const ns = ev.lat>=0?'N':'S', ew = ev.lon<=0?'W':'E';
  const isTele2 = ev.teleseism||false;
  const typeTag = isTele2 ? 'TELESEISM' : (mag>=4 ? '&#9888; PRELIMINARY' : 'PRELIMINARY');
  const mmiStr2 = ev.mmi_str || '?';
  const mmiDesc2 = ev.mmi_desc || '';
  row.innerHTML = `
    <div style="display:flex;align-items:baseline;gap:8px;margin-bottom:3px">
      <span style="color:${col};font-weight:700;font-size:14px">Est. M${mag.toFixed(1)}</span>
      <span style="color:#8b949e;font-size:10px">${utcToPT(ev.ot_str||'')}</span>
      <span style="color:#58a6ff;font-size:10px;margin-left:auto">${typeTag}</span>
    </div>
    <div style="color:#c9d1d9">${Math.abs(ev.lat).toFixed(2)}°${ns}  ${Math.abs(ev.lon).toFixed(2)}°${ew}  ·  ${ev.city||'Unknown area'}</div>
    <div style="color:#f0c040;font-size:10px;margin-top:2px;font-weight:600">Max MMI ${mmiStr2} — ${mmiDesc2}</div>
    <div style="color:#6e7681;font-size:10px;margin-top:1px">${ev.n_stations} stations  ·  max STA/LTA ${ev.max_stalta}  ·  est. depth ~${ev.depth_est} km${ev.rms_sec != null ? '  ·  RMS ' + ev.rms_sec + 's' : ''}</div>
    <div style="color:#6e7681;font-size:9px;margin-top:1px">${ev.stations}</div>`;
  list.insertBefore(row, list.firstChild);
  // Flash the tab
  const tab = document.getElementById('tab-prelim');
  if(tab){tab.style.color='#f85149';setTimeout(()=>tab.style.color='',3000);}
}
socket.on('preliminary', ev => {
  addPrelimRow(ev);
  setTab('prelim');
  toggleDrawer(true);
});
// Load existing prelim events on startup
fetch('/api/preliminary').then(r=>r.json()).then(evs=>{
  evs.forEach(ev=>addPrelimRow(ev));
}).catch(()=>{});
</script></body></html>"""

# ── Startup ────────────────────────────────────────────────────────────────────
def _chain_plugin_cleanup():
    """
    Kill zombie chain_plugin processes every 5 minutes.
    SeisComP restarts chain_plugin on crash but never kills the old PID,
    so they accumulate and starve the real upstream SeedLink connection.
    Keep only the 2 newest PIDs (one per upstream: SCEDC + NCEDC).
    """
    import subprocess, signal
    while True:
        time.sleep(300)
        try:
            r = subprocess.run(["pgrep", "-f", "chain_plugin"],
                               capture_output=True, text=True)
            pids = sorted([int(p) for p in r.stdout.split() if p.strip()])
            if len(pids) > 4:
                zombies = pids[:-4]   # kill all but the 4 most recent
                for pid in zombies:
                    try: os.kill(pid, signal.SIGKILL)
                    except: pass
                print(f"[CLEANUP] killed {len(zombies)} zombie chain_plugin(s), "
                      f"{len(pids)-len(zombies)} kept", flush=True)
        except Exception as _ce:
            print(f"[CLEANUP] chain_plugin check error: {_ce}", flush=True)


if __name__=="__main__":
    _load_sta_coords()
    # Background inventory fetch — enriches station coords from IRIS (non-blocking)
    _inv_nets = list(FDSN_NETWORKS.keys()) + ["SB"]
    threading.Thread(target=_fetch_iris_inventory_for, args=(_inv_nets,), daemon=True).start()
    # Populate city labels (runs after coords loaded; re-run after inventory fetch completes)
    threading.Thread(target=_populate_cities, daemon=True).start()
    _load_prelim_log()   # restore saved detections before starting
    # Auto-kill zombie chain_plugin processes every 5 min to keep SeedLink healthy
    threading.Thread(target=_chain_plugin_cleanup, daemon=True).start()
    start_seedlink()
    threading.Thread(target=emit_loop,daemon=True).start()
    time.sleep(1)
    print(f"Dashboard → http://localhost:{PORT}")
    print(f"Live Map  → http://localhost:{PORT}/livemap")
    sio.run(app,host="0.0.0.0",port=PORT,allow_unsafe_werkzeug=True)
