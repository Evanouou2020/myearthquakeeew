#!/usr/bin/env python3
"""
seedlink_alert.py  —  Real-time seismic early-warning monitor  v4
=================================================================
  • Live OSM tile map (contextily) – zoom / pan – tiles auto-refresh on zoom
  • Expanding P & S wavefront rings once epicenter confirmed
  • Per-station suspected-epicenter rings (before solution)
  • Per-station sensitivity calibration (thresh_mult)
  • PGV (cm/s) + PGA (cm/s²) displayed on every waveform panel
  • Source ML  +  Max observed MMI  +  MMI @ San Ramon  shown separately
  • Single-station ML/depth (station dist as proxy; depth fixed at 10 km)
  • 3-D depth from Nelder-Mead when ≥2 stations detected
  • Predicted-P dashed line per station once epicenter is known
  • Time-since-last-update label under each station name
  • Test earthquake injection  →  type  t  + Enter
  • Full event log  →  seismic_alerts.log

Install:  pip install obspy matplotlib numpy scipy contextily pyproj
Run:      python seedlink_alert.py

Terminal commands
  s / status        – live station table
  t / test          – inject simulated M3.5 near SF Bay
  sm <usgs_id>      – fetch USGS ShakeMap + email  (e.g. sm nc75361631)
  q / quit          – exit
"""
import collections, io, math, os, platform, subprocess, sys, threading, time
import smtplib, email.mime.text, email.mime.multipart, email.mime.image
import numpy as np
from scipy.optimize import minimize
from scipy.signal  import butter, sosfilt, bilinear, lfilter

import matplotlib
matplotlib.use("MacOSX" if platform.system() == "Darwin" else "TkAgg")
import matplotlib.pyplot   as plt
import matplotlib.animation as animation

try:
    import contextily as ctx
    from pyproj import Transformer
    HAS_TILES = True
    _to_merc  = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    print("[MAP] contextily OK – live OSM tiles enabled")
except ImportError:
    HAS_TILES = False
    print("[MAP] pip install contextily pyproj  for OSM tiles (using static CA fallback)")

from obspy.clients.seedlink.easyseedlink import EasySeedLinkClient
from obspy.signal.trigger                import recursive_sta_lta

# ═══════════════════════════════════════════════════════════════════════════════
# PHYSICS
# ═══════════════════════════════════════════════════════════════════════════════
VP        = 6.0
VS        = 3.5
SP_FACTOR = 1/VS - 1/VP   # s/km ≈ 0.119

HOME_LAT, HOME_LON = 37.7644, -121.9540
HOME_LABEL         = "San Ramon, CA"

SEEDLINK_HOST = "rtserve.iris.washington.edu"
SEEDLINK_PORT = 18000

# ═══════════════════════════════════════════════════════════════════════════════
# STATIONS
# (net, sta, loc, cha, lat, lon, km_from_home, desc,
#  sensitivity_counts_per_m_s, thresh_mult)
#
#  thresh_mult < 1.0  → lower trigger bar  (quiet borehole / under-reports)
#  thresh_mult > 1.0  → raise trigger bar  (noisy coastal / over-triggers)
#
#  Calibration notes 2026-04:
#    B054 Berkeley Hills borehole – under-sensitive,     thresh_mult = 0.50
#    MCCM Point Reyes             – too sensitive,       thresh_mult = 1.90
#    B058 Watsonville             – slightly boosted,    thresh_mult = 0.80
#    CMB  Columbia Mine           – slightly boosted,    thresh_mult = 0.80
# ═══════════════════════════════════════════════════════════════════════════════
STATIONS = [
    # ── IRIS SeedLink real-time stations (confirmed streaming) ─────────────
    # Note: NC and most BK stations are not on the IRIS real-time feed.
    # Closer stations (RVRP, LLNL, WENL, BDM, BRIB) exist in the FDSN archive
    # but are only served by the NCEDC SeedLink server (ncedc.org:18000).
    ("PB","B054","","EHZ",  37.860199,-122.199501, 25,"Berkeley Hills",  7.81398e8, 0.50),
    ("BK","MCCM","00","HHZ",38.144779,-122.880180, 63,"Point Reyes",     2.51658e9, 1.90),
    ("PB","B057","","EHZ",  38.027302,-122.565498, 70,"Petaluma",        7.81398e8, 1.00),
    ("BK","SAO", "00","HHZ",36.764030,-121.447220, 85,"San Andreas Obs", 2.51558e9, 1.00),
    ("PB","B058","","EHZ",  36.799500,-121.580300,103,"Watsonville",     7.81398e8, 0.80),
    ("BK","CMB", "00","HHZ",38.034550,-120.386513,148,"Columbia Mine",   2.51995e9, 0.80),
]
N = len(STATIONS)

# Original thresh_mult per station — used to restore defaults
_DEFAULT_THRESH = {
    f"{net}.{sta}.{loc}.{cha}": tmult
    for net, sta, loc, cha, _la, _lo, _d, _desc, _sens, tmult in STATIONS
}

# ═══════════════════════════════════════════════════════════════════════════════
# DETECTION PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
STA_SEC         = 1.0
LTA_SEC         = 30.0
TRIGGER_ON      = 8.0      # raised from 6.5 — require a very clear signal onset
TRIGGER_OFF     = 2.5
P_STA_SEC       = 0.8    # raised from 0.5 s — longer window kills spike transients
P_LTA_SEC       = 10.0
P_THRESH        = 9.0    # raised from 6.0 — demands a prominent P onset, kills noise

# ── PhaseNet ML phase picker ──────────────────────────────────────────────────
PHASENET_ENABLED   = True   # set False to use STA/LTA+AIC only
PHASENET_PRETRAIN  = "stead"     # STEAD = Stanford EQ Dataset (N. California) — best match for BK/PB networks
PHASENET_P_THRESH  = 0.3    # minimum P probability to accept pick
PHASENET_S_THRESH  = 0.3    # minimum S probability to accept pick
PHASENET_MAX_SHIFT = 5.0    # reject PhaseNet P if it differs > this many s from STA/LTA

# ── Live WebSocket waveform server ────────────────────────────────────────────
WS_HOST = "localhost"
WS_PORT = 8765              # website connects to ws://localhost:8765

DISPLAY_SEC     = 300    # 5-minute rolling waveform buffer (was 3600 — 12× RAM reduction)
#   Detection only needs ~120 s (LTA=30 s + margin); 300 s gives a comfortable
#   display window.  Use the "Hide Waves" button to shrink further to 120 s.
_DISPLAY_SEC_FULL = DISPLAY_SEC   # buffer size when waveforms are visible
_DISPLAY_SEC_MIN  = 120           # buffer size when waveforms are hidden
_display_sec_active = [DISPLAY_SEC]  # mutable — Hide Waves shrinks to _DISPLAY_SEC_MIN
_WAVES_HIDDEN       = [False]        # True while waveform panels are collapsed
ALERT_COOLDOWN  = 30
EVENT_RESET     = 1800   # reset station after 30 min of quiet (was 5 min)
STALE_DATA_SEC  = 210    # data older than this (wall-clock) is backlog — skip detection
                         # 3.5 min: normal reconnect/buffer flush takes 30–90 s,
                         # so only flag truly stale data (missed events, not transient lag)
MIN_STA_CONFIRM = 2        # P-arrivals needed for "confirmed" event

P_FREQ_LO, P_FREQ_HI = 2.0, 8.0
S_FREQ_LO, S_FREQ_HI = 0.5, 2.5
PS_WIN_SEC  = 3.0
PS_S_THRESH = 0.65

# S-band dedicated STA/LTA (applied to S-filtered signal)
S_STA_SEC   = 1.5    # longer window suits lower S-wave frequencies
S_LTA_SEC   = 25.0
S_THRESH    = 3.5    # S-band STA/LTA trigger (lower than P — S weaker on Z)

# AIC refinement search window around STA/LTA trigger (seconds each side)
AIC_WIN_SEC = 4.0

_IS_MACOS    = platform.system() == "Darwin"
_SOUND       = next((s for s in ["/System/Library/Sounds/Sosumi.aiff",
                                  "/System/Library/Sounds/Glass.aiff"]
                     if os.path.exists(s)), "")
_sound_muted = False   # toggled by the Sound ON/Muted button in the banner

# ═══════════════════════════════════════════════════════════════════════════════
# LOGGING
# ═══════════════════════════════════════════════════════════════════════════════
LOG_FILE  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "seismic_alerts.log")
_log_lock = threading.Lock()

def _log(tag: str, msg: str):
    ts   = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
    line = f"[{ts}] [{tag:<12}] {msg}\n"
    with _log_lock:
        with open(LOG_FILE, "a", encoding="utf-8") as f:
            f.write(line)
    print(line, end="")

# ═══════════════════════════════════════════════════════════════════════════════
# EMAIL  CONFIGURATION  — fill in your details here
# ═══════════════════════════════════════════════════════════════════════════════
EMAIL_ENABLED   = True
EMAIL_FROM      = "evanouou2020@gmail.com"   # sender (Gmail recommended)
EMAIL_PASSWORD  = "ftjo qxsa yllu cjgr"      # Gmail App Password (16 chars)
EMAIL_TO        = "alt90228@gmail.com"   # recipient
EMAIL_SMTP_HOST = "smtp.gmail.com"
EMAIL_SMTP_PORT = 587

# Per-type cooldowns (seconds) — prevents duplicate emails for the same phase
EMAIL_COOLDOWN_P  = 30   # per station, per event
EMAIL_COOLDOWN_S  = 30   # per station, per event
EMAIL_COOLDOWN_ML = 60   # magnitude updates
EMAIL_COOLDOWN_EQ = 90   # confirmed-quake updates

# ═══════════════════════════════════════════════════════════════════════════════
# NTFY.SH  —  Free push notifications (no account needed)
#   1. Install the free ntfy app on iOS / Android
#   2. Subscribe to your unique topic name in the app
#   3. Set NTFY_ENABLED = True and fill in NTFY_TOPIC (make it unique!)
#
#   Priority levels: min | low | default | high | urgent
#   "urgent" + "high-priority" tag bypasses DND on most phones.
# ═══════════════════════════════════════════════════════════════════════════════
NTFY_ENABLED = True
NTFY_TOPIC   = "early_earthquake_warning_san_ramon"
NTFY_SERVER  = "https://ntfy.sh"

# ── Discord webhook ────────────────────────────────────────────────────────
# Paste your webhook URL here (Server Settings → Integrations → Webhooks).
# Set to None or "" to disable.
DISCORD_WEBHOOK_URLS = [
    "https://discord.com/api/webhooks/1499565828285141193/hoGzhuk2uwOLdZV4eZX-NhdUjkzFurAIee9AbMK6RkDMgRf-G8mq3JXz2uSokybUtAiu",
    "https://discord.com/api/webhooks/1508676731656077392/A_xCrUH11bb428dQu3ii6kNDS-JgaVd9msakiB7pI8pGHOarHA7PVcxu8qIP4WAmf6YN",
]

# Minimum ML to send a push notification.
# P-wave alerts fire before ML is known, so NTFY_MIN_ML_PWAVE controls whether
# the first P-wave detection pushes at all.  Set to None to always push.
# Confirmed-EQ and magnitude-update notifications obey NTFY_MIN_ML.
#   Examples:  0.0 = anything detected
#              1.5 = skip micro-quakes (saves your 250/day free quota)
#              2.5 = only felt events
NTFY_MIN_ML_PWAVE = 1.0   # push P-wave only if station ML estimate >= this
NTFY_MIN_ML       = 1.5   # push magnitude/confirmed alerts only if ML >= this (raised from 1.0)
NTFY_EEW_ML       = 2.5   # EARLY WARNING threshold — urgent alert with arrival countdown

# ── Alert radius filter ─────────────────────────────────────────────────────
ALERT_RADIUS_KM    = 300    # only send notifications if epicenter is within this distance
                            # from HOME. Set to None to always alert.

# ── Sound alarm ─────────────────────────────────────────────────────────────
SOUND_ALARM_ENABLED = True   # play audio alarm on HIGH/EXTREME detection (macOS afplay)
SOUND_ALARM_LEVEL   = "HIGH" # minimum level: "MEDIUM", "HIGH", or "EXTREME"

# ── NCEDC SeedLink (closer BK/NC stations) ──────────────────────────────────
NCEDC_ENABLED      = True    # connect to ncedc.org for additional nearby stations
NCEDC_HOST         = "ncedc.org"
NCEDC_PORT         = 18000
# NC/BK stations on NCEDC close to San Ramon (lat 37.78, lon -121.98):
NCEDC_STATIONS     = [
    ("BK","BRIB","","HHZ",  37.918900,-122.151800, 22,"Briones Hills",  2.51000e9, 1.00),
    ("BK","WENL","","HHZ",  37.622100,-121.757000, 26,"South Livermore",2.51000e9, 1.00),
    ("BK","BKS", "","HHZ",  37.876200,-122.235600, 27,"Berkeley",       2.51000e9, 1.20),
    ("BK","TESL","00","HHZ",37.614900,-121.590500, 39,"Tesla",          2.51000e9, 1.00),
    ("BK","JRSC","","HHZ",  37.403700,-122.238700, 48,"Jasper Ridge",   2.51000e9, 1.00),
]

# ── Aftershock probability ───────────────────────────────────────────────────
AFTERSHOCK_WINDOW_HOURS = 24   # forecast window for Omori-Utsu probability

# ── Alert cooldown (aftershock suppression) ──────────────────────────────────
AFTERSHOCK_COOLDOWN_SEC = 300  # suppress repeat EQ-confirmed alerts within this window
                               # for the same event cluster

# ── SQLite event catalog ─────────────────────────────────────────────────────
CATALOG_DB_PATH    = "seismic_catalog.db"   # path to SQLite database
CATALOG_ENABLED    = True

# ── Seismicity rate / b-value ────────────────────────────────────────────────
BVALUE_MIN_EVENTS  = 10   # minimum events needed to compute b-value

# ── Strong-signal alert (teleseism / large regional event) ───────────────────
# Fires when multiple stations simultaneously show an elevated general STA/LTA,
# even if the P-wave picker never triggers (teleseisms are long-period and
# don't excite the 2–8 Hz P-band filter the P-picker uses).
#
# How the level maps to raw STA/LTA ratio × trigger threshold:
#   LOW     =  1×– 3× TRIGGER_ON   →  barely above noise
#   MEDIUM  =  3×–10× TRIGGER_ON   →  clearly elevated
#   HIGH    = 10×–20× TRIGGER_ON   →  very strong shaking
#   EXTREME = ≥ 20×   TRIGGER_ON   →  exceptional
STRONG_SIG_ENABLED    = True     # set False to disable entirely
STRONG_SIG_MIN_LEVEL  = "MEDIUM" # minimum level to fire (MEDIUM = 3× trigger threshold)
STRONG_SIG_MIN_STA    = 2        # how many stations must be at that level simultaneously
STRONG_SIG_COOLDOWN   = 300      # seconds between repeated strong-signal alerts

# ── Discord false-detection guard ────────────────────────────────────────────
# P-wave Discord alerts fire only after this many stations have independently
# confirmed a P-arrival.  NTFY and email still fire on the first detection.
# Set to 1 to restore the original single-station Discord behaviour.
DISCORD_MIN_P_STA      = 2          # require 2 stations — eliminates single-station noise
DISCORD_PWAVE_MIN_LEVEL = "HIGH"    # require 10× threshold (was MEDIUM=3×)
DISCORD_MIN_P_ML        = 1.5       # don't Discord P-wave unless ML estimate ≥ this
DISCORD_EVERYONE_ML     = 4.0       # @everyone mention for earthquakes at or above this ML
DISCORD_SAFETY_ML       = 3.0       # attach Drop/Cover/Hold On card at or above this ML
DISCORD_SPAM_ML         = 5.0       # spam @everyone + safety card repeatedly above this ML
DISCORD_SPAM_COUNT      = 5         # how many rapid-fire @everyone safety pings to send
DISCORD_SPAM_DELAY      = 1.5       # seconds between each spam ping (Discord rate-limit safe)
# Path to the safety card image (relative to script directory)
import os as _os
_SAFETY_CARD_PATH = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                  "drop_cover_hold_on.png")

# ── STA/LTA @everyone alert ──────────────────────────────────────────────────
# Fires an @everyone Discord ping whenever ANY station's general STA/LTA ratio
# reaches TRIGGER_ON (8.0) — even before a P-wave pick or epicenter solve.
STA_LTA_EVERYONE_COOLDOWN = 120   # seconds between repeated STA/LTA @everyone pings
_sta_lta_everyone_last_t  = [0.0] # last send timestamp

# ── Thread state (all mutable via single-element lists) ──────────────────────
_email_lock          = threading.Lock()
_email_last_p        = {}     # station key → last P email timestamp
_email_last_s        = {}     # station key → last S email timestamp
_email_last_ml       = [0.0]  # last ML-update email timestamp
_email_last_eq       = [0.0]  # last confirmed-quake email timestamp
_strong_sig_last_t   = [0.0]  # last strong-signal alert timestamp (cooldown)
_event_thread_id     = [None] # Message-ID of first email in current event thread
_event_thread_subj   = [None] # Subject of first email (for display continuity)

def _eq_ensure_started():
    """
    Create a new EQ ID if none exists (atomic check-and-set).
    Returns the current EQ ID string.
    Safe to call from any thread.
    """
    with _eq_id_lock:
        if _current_eq_id[0] is None:
            _eq_id_counter[0] += 1
            _current_eq_id[0] = f"EQ-{_eq_id_counter[0]:04d}"
            _current_eq_upd[0] = 0
            _current_eq_usgs[0] = None
        return _current_eq_id[0]

def _eq_next_update():
    """
    Increment update counter.  Returns (update_num, eq_id_str).
    Call once per Discord milestone (ML update, confirmed, EEW, USGS, final).
    """
    with _eq_id_lock:
        _current_eq_upd[0] += 1
        return _current_eq_upd[0], (_current_eq_id[0] or "EQ-????")

def _eq_reset_id():
    """Clear active earthquake state on global reset."""
    with _eq_id_lock:
        _current_eq_id[0]   = None
        _current_eq_upd[0]  = 0
        _current_eq_usgs[0] = None

def _eq_discord_prefix(suffix: str = "") -> str:
    """
    Build the bracket prefix for Discord titles, e.g.:
      '[EQ-0042]'  or  '[EQ-0042 | Upd #3]'  or  '[EQ-0042 / USGS:nc12345 | Upd #3]'
    """
    with _eq_id_lock:
        eq   = _current_eq_id[0] or "EQ-????"
        usgs = _current_eq_usgs[0]
    usgs_s = f" / USGS:{usgs}" if usgs else ""
    inner  = eq + usgs_s
    if suffix:
        inner += f" | {suffix}"
    return f"[{inner}]"

def _new_event_thread():
    """Call when a new seismic event starts to open a fresh email thread."""
    import uuid
    with _email_lock:
        _event_thread_id[0]   = None   # will be set after first email sends
        _event_thread_subj[0] = None

def _reset_event_thread():
    """Call on global reset to close the thread.
    NOTE: does NOT clear the timeline here so _send_final_report can
    still read it.  Timeline is cleared explicitly by _animate after
    queuing the final report.
    """
    with _email_lock:
        _event_thread_id[0]   = None
        _event_thread_subj[0] = None
    _email_last_p.clear()
    _email_last_s.clear()
    _email_last_ml[0] = 0.0
    _email_last_eq[0] = 0.0
    # timeline intentionally NOT cleared here

# ── Event timeline (accumulated per-event, sent in final report) ──────────────
_event_timeline      = []
_event_timeline_lock = threading.Lock()

def _timeline_add(note: str):
    """Append a timestamped note to the current event timeline."""
    ts = time.strftime("%H:%M:%S UTC", time.gmtime())
    with _event_timeline_lock:
        _event_timeline.append((ts, note))

def _timeline_clear():
    with _event_timeline_lock:
        _event_timeline.clear()

def _timeline_snapshot():
    """Return a copy of the current timeline list."""
    with _event_timeline_lock:
        return list(_event_timeline)

# ── Subject icon helper (module-level so _send_final_report can use it) ───────
def _subj_icon(ml):
    """Return email-subject emoji prefix based on ML magnitude."""
    if ml is None or ml < 3.0: return ""
    if ml < 4.0:               return "[M3+] "
    if ml < 5.0:               return "[M4+] "
    return "[M5+] "

def _send_email(subject: str, body: str, img_bytes: bytes = None):
    """
    Send an alert email in a background daemon thread.

    • img_bytes — optional PNG screenshot attached as an inline image.
    All emails belonging to the same seismic event are threaded together:
    the first email sets Message-ID → stored as _event_thread_id[0].
    Every subsequent email sets In-Reply-To + References to that ID, which
    causes Gmail / Outlook / Apple Mail to group them as one conversation.
    """
    if not EMAIL_ENABLED:
        return
    if EMAIL_PASSWORD == "your-app-password-here":
        print("[EMAIL] Not configured — fill in EMAIL_FROM/PASSWORD/TO.")
        return

    import uuid
    msg_id  = f"<seismic-{int(time.time())}-{uuid.uuid4().hex[:8]}@quakemon>"
    with _email_lock:
        reply_to = _event_thread_id[0]

    def _worker(_mid=msg_id, _reply=reply_to):
        try:
            if img_bytes:
                # Outer mixed container (allows both related block + plain fallback)
                msg = email.mime.multipart.MIMEMultipart("mixed")
                msg["Subject"]    = subject
                msg["From"]       = EMAIL_FROM
                msg["To"]         = EMAIL_TO
                msg["Message-ID"] = _mid
                if _reply:
                    msg["In-Reply-To"] = _reply
                    msg["References"]  = _reply
                # Related block: HTML body + inline image
                related = email.mime.multipart.MIMEMultipart("related")
                html = (
                    f"<html><body style='background:#0d0d0d;color:#cccccc;"
                    f"font-family:monospace;font-size:13px'>"
                    f"<pre style='white-space:pre-wrap'>{body}</pre>"
                    f"<br><strong style='color:#88bbff'>Live waveforms at detection time:</strong><br>"
                    f"<img src='cid:waveform_cid' style='max-width:100%;"
                    f"border:1px solid #333;margin-top:6px'>"
                    f"</body></html>"
                )
                related.attach(email.mime.text.MIMEText(html, "html"))
                img_part = email.mime.image.MIMEImage(img_bytes)
                img_part.add_header("Content-ID", "<waveform_cid>")
                img_part.add_header("Content-Disposition", "inline",
                                    filename="waveforms.png")
                related.attach(img_part)
                msg.attach(related)
                msg.attach(email.mime.text.MIMEText(body, "plain"))
            else:
                msg = email.mime.multipart.MIMEMultipart("alternative")
                msg["Subject"]    = subject
                msg["From"]       = EMAIL_FROM
                msg["To"]         = EMAIL_TO
                msg["Message-ID"] = _mid
                if _reply:
                    msg["In-Reply-To"] = _reply
                    msg["References"]  = _reply
                msg.attach(email.mime.text.MIMEText(body, "plain"))

            with smtplib.SMTP(EMAIL_SMTP_HOST, EMAIL_SMTP_PORT, timeout=15) as srv:
                srv.ehlo(); srv.starttls()
                srv.login(EMAIL_FROM, EMAIL_PASSWORD)
                srv.sendmail(EMAIL_FROM, EMAIL_TO, msg.as_string())
            _log("EMAIL", f"sent{'  +screenshot' if img_bytes else ''}: {subject}")
            # Anchor the thread on the very first email of the event
            with _email_lock:
                if _event_thread_id[0] is None:
                    _event_thread_id[0]  = _mid
                    _event_thread_subj[0] = subject
        except Exception as exc:
            _log("EMAIL ERR", f"{exc}  |  {subject}")

    threading.Thread(target=_worker, daemon=True).start()


# ── Waveform screenshot buffer ────────────────────────────────────────────────
# Holds the most-recent PNG bytes captured by _capture_screenshot() (main thread).
# Background threads (P-wave emailer) read it safely via the GIL — one writer,
# many readers, bytes are immutable.
_screenshot_buf              = [None]   # [bytes | None]
_screenshot_email_requested  = [False]  # set True by ntfy "screenshot" cmd; cleared in _animate

# NTFY listener watchdog — track liveness so the watchdog can restart if stalled
_ntfy_listener_thread = [None]   # current listener Thread object
_ntfy_last_rx_t       = [0.0]    # epoch of last byte received from ntfy stream

# WebSocket live-waveform server state
_ws_connected = set()         # currently connected ServerConnection objects
_ws_loop      = [None]        # asyncio event loop running in background thread

# Aftershock cooldown state
_last_eq_confirmed_t  = [0.0]
_last_eq_confirmed_ml = [None]

# ── Earthquake identifier state ───────────────────────────────────────────────
# Each new seismic event gets a local ID like "EQ-0042" the moment the first
# P-wave is detected.  All Discord messages for that event share the same ID
# so it's always clear which updates belong together.
_eq_id_lock       = threading.Lock()
_eq_id_counter    = [0]     # global counter; never reset
_current_eq_id    = [None]  # e.g. "EQ-0042" — None when no active event
_current_eq_upd   = [0]     # update number within current event
_current_eq_usgs  = [None]  # USGS ComCat event ID once USGS cross-check matches

# PhaseNet ML picker state
_phasenet_model   = [None]    # loaded at startup in background thread
_phasenet_ready   = [False]   # True once model is loaded
_phasenet_pending = {}        # key → {p_time, p_prob, s_time, s_prob}; written by bg thread
_phasenet_pend_lk = threading.Lock()

def _capture_screenshot():
    """
    Render the current figure to PNG and store in _screenshot_buf.
    MUST be called from the main (matplotlib) thread only.
    Returns the PNG bytes, or None on failure.
    """
    try:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=90,
                    facecolor=fig.get_facecolor(), bbox_inches=None)
        buf.seek(0)
        data = buf.read()
        _screenshot_buf[0] = data
        return data
    except Exception as _e:
        _log("SCREENSHOT", f"capture failed: {_e}")
        return None

# ── Deferred screenshot emails ────────────────────────────────────────────────
# Emails that need an inline waveform screenshot are queued here with a
# fire_after timestamp.  _animate processes the queue each frame; when ready
# it captures a fresh screenshot on the main thread (matplotlib-safe) and
# dispatches the send in a background thread.  This guarantees the screenshot
# is taken AFTER the waveforms have had time to build up on screen, not at the
# moment the detection fires in a background on_data() thread.
SCREENSHOT_DELAY_SEC     = 5     # seconds after trigger before capturing

_pending_shot_emails      = []   # [(fire_after_t, subject, body)]
_pending_shot_emails_lock = threading.Lock()

def _queue_email_with_shot(subject, body, delay=SCREENSHOT_DELAY_SEC):
    """
    Schedule an email whose waveform screenshot will be captured ~delay
    seconds from now on the main thread, then sent in a background thread.
    Thread-safe — may be called from any thread.
    """
    fire_after = time.time() + delay
    with _pending_shot_emails_lock:
        _pending_shot_emails.append((fire_after, subject, body))
    _log("SCREENSHOT", f"deferred email queued (+{delay:.0f}s): {subject[:60]}")

def _send_ntfy(title: str, message: str, priority: str = "high",
               tags: str = ""):
    """
    Send a push notification via ntfy.sh in a daemon background thread.

    • No account required — install the free ntfy app and subscribe to NTFY_TOPIC.
    • priority="urgent" + tags="warning,rotating_light" bypasses iOS/Android DND.
    • Safe to call from any thread.
    """
    if not NTFY_ENABLED:
        return
    if NTFY_TOPIC == "seismic-alerts-CHANGEME":
        _log("NTFY", "Topic not configured — set NTFY_TOPIC in config and restart.")
        return

    def _worker():
        try:
            import urllib.request, json, ssl
            # Build an SSL context that works on macOS without the system cert store.
            # Try certifi first (pip install certifi); fall back to unverified context
            # (safe here — we're only sending outbound push notifications, not
            # receiving sensitive data).
            try:
                import certifi
                _ssl_ctx = ssl.create_default_context(cafile=certifi.where())
            except ImportError:
                _ssl_ctx = ssl._create_unverified_context()

            _pri_map = {"min": 1, "low": 2, "default": 3, "high": 4, "urgent": 5}
            payload = json.dumps({
                "topic":    NTFY_TOPIC,
                "title":    title,
                "message":  message,
                "priority": _pri_map.get(priority, 4),
                "tags":     [t.strip() for t in tags.split(",") if t.strip()],
            }).encode("utf-8")
            req = urllib.request.Request(
                NTFY_SERVER,
                data=payload,
                headers={"Content-Type": "application/json"},
                method="POST",
            )
            urllib.request.urlopen(req, timeout=10, context=_ssl_ctx)
            _log("NTFY", f"sent: {title}")
        except Exception as exc:
            _log("NTFY ERR", str(exc))

    threading.Thread(target=_worker, daemon=True).start()


def _send_discord(title: str, message: str, level: str = "", img_bytes: bytes = None,
                  bold_header: str = "", everyone: bool = False):
    """
    Send a Discord embed notification via webhook in a daemon background thread.

    level: "LOW" | "MEDIUM" | "HIGH" | "EXTREME" | "EEW" | ""
    img_bytes: optional PNG bytes — sent as a file attachment shown inside the embed.

    Embed colour:
      LOW     → grey   (#808080)
      MEDIUM  → yellow (#FFCC00)
      HIGH    → orange (#FF6600)
      EXTREME → red    (#FF0000)
      EEW     → bright red (#FF2222)
      (none)  → teal   (#1ABC9C)
    """
    if not DISCORD_WEBHOOK_URLS:
        return

    _colours = {
        "LOW":     0x808080,
        "MEDIUM":  0xFFCC00,
        "HIGH":    0xFF6600,
        "EXTREME": 0xFF0000,
        "EEW":     0xFF2222,
    }
    colour = _colours.get(level.upper(), 0x1ABC9C)

    # Capture snapshot of img_bytes for the closure (avoid late-binding issues)
    _img      = img_bytes
    _bhdr     = bold_header
    _everyone = everyone

    def _worker():
        import http.client, json as _json, ssl, urllib.parse
        _ssl_ctx = ssl.create_default_context()

        # Build embed + payload once — reused for every webhook URL
        if _bhdr:
            _hdr_md = f"**{_bhdr}**\n\n"
        else:
            _hdr_md = ""
        _max_code = 4090 - len(_hdr_md)
        _code_body = message[:_max_code] + "…" if len(message) > _max_code else message
        desc = _hdr_md + f"```\n{_code_body}\n```"

        embed = {
            "title":       title[:256],
            "description": desc[:4096],
            "color":       colour,
            "footer":      {"text": "This is a preliminary earthquake alert system "
                                    "compiled by Evan Li. For more information visit: "
                                    "https://myearthquake.dpdns.org/"},
        }
        if _img:
            embed["image"] = {"url": "attachment://epicenter_map.png"}

        # @everyone mention: put it in `content` (Discord only pings from that field)
        _mention_payload: dict = {}
        if _everyone:
            _mention_payload["content"]          = "@everyone"
            _mention_payload["allowed_mentions"] = {"parse": ["everyone"]}

        for _url in DISCORD_WEBHOOK_URLS:
            try:
                _parsed = urllib.parse.urlparse(_url)

                if _img:
                    # ── Multipart form-data: embed JSON + PNG file ─────────────
                    boundary    = "----QuakeAlertBoundary7f3a"
                    payload_json = _json.dumps({**_mention_payload, "embeds": [embed]})
                    parts = [
                        (f"--{boundary}\r\n"
                         f'Content-Disposition: form-data; name="payload_json"\r\n'
                         f"Content-Type: application/json\r\n\r\n"
                         f"{payload_json}\r\n").encode("utf-8"),
                        (f"--{boundary}\r\n"
                         f'Content-Disposition: form-data; name="files[0]";'
                         f' filename="epicenter_map.png"\r\n'
                         f"Content-Type: image/png\r\n\r\n").encode("utf-8")
                        + _img + b"\r\n",
                        f"--{boundary}--\r\n".encode("utf-8"),
                    ]
                    body = b"".join(parts)
                    conn = http.client.HTTPSConnection(_parsed.netloc, timeout=20,
                                                       context=_ssl_ctx)
                    conn.request("POST", _parsed.path, body=body, headers={
                        "Content-Type":   f"multipart/form-data; boundary={boundary}",
                        "Content-Length": str(len(body)),
                        "User-Agent":     "QuakeAlertBot/1.0",
                    })
                else:
                    # ── Plain JSON embed (no attachment) ───────────────────────
                    payload = _json.dumps({**_mention_payload, "embeds": [embed]}).encode("utf-8")
                    conn = http.client.HTTPSConnection(_parsed.netloc, timeout=10,
                                                       context=_ssl_ctx)
                    conn.request("POST", _parsed.path, body=payload, headers={
                        "Content-Type":   "application/json",
                        "Content-Length": str(len(payload)),
                        "User-Agent":     "QuakeAlertBot/1.0",
                    })

                resp = conn.getresponse()
                conn.close()
                if resp.status in (200, 204):
                    _log("DISCORD", f"sent → {_parsed.netloc}{_parsed.path[:40]}: "
                                    f"{title}" + (" (+map)" if _img else ""))
                else:
                    _log("DISCORD ERR", f"{_parsed.path[:40]} HTTP {resp.status}: "
                                        f"{resp.read(400)}")
            except Exception as exc:
                _log("DISCORD ERR", f"{_url[:60]}: {exc}")

    threading.Thread(target=_worker, daemon=True).start()


def _make_map_png(elat: float, elon: float, et0: float = None,
                  med_ml: float = None) -> bytes:
    """
    Generate an off-screen epicenter map PNG using the matplotlib Agg backend.

    Uses the Agg canvas directly so it never conflicts with the live MacOSX/Tk
    display.  Contextily basemap tiles are fetched if available; falls back to a
    clean dark background on network failure.

    Returns PNG bytes, or None on any error.
    """
    try:
        from matplotlib.figure       import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        import matplotlib.patches    as _mpatch
        import io as _io

        fig    = Figure(figsize=(8, 6), facecolor="#0d0d0d")
        canvas = FigureCanvasAgg(fig)
        ax     = fig.add_subplot(111)
        ax.set_facecolor("#111122")

        # ── Map extent: cover epicenter + all stations + home ──────────────────
        all_lats = [elat, HOME_LAT] + [s[4] for s in STATIONS + NCEDC_STATIONS]
        all_lons = [elon, HOME_LON] + [s[5] for s in STATIONS + NCEDC_STATIONS]
        pad_lat  = max(0.4, (max(all_lats) - min(all_lats)) * 0.20)
        pad_lon  = max(0.4, (max(all_lons) - min(all_lons)) * 0.20)
        xmin, xmax = min(all_lons) - pad_lon, max(all_lons) + pad_lon
        ymin, ymax = min(all_lats) - pad_lat, max(all_lats) + pad_lat
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        # ── Contextily basemap (CartoDB dark, no labels) ───────────────────────
        _basemap_ok = False
        if HAS_TILES:
            try:
                import contextily as _ctx
                _ctx.add_basemap(ax, crs="EPSG:4326", zoom=8,
                                 source=_ctx.providers.CartoDB.DarkMatterNoLabels,
                                 attribution=False)
                _basemap_ok = True
            except Exception as _bme:
                _log("MAP", f"basemap fetch failed: {_bme}  — using plain background")

        if not _basemap_ok:
            # Minimal grid lines as a fallback
            import numpy as _np_m
            for _gla in _np_m.arange(round(ymin, 0), round(ymax, 0) + 1, 0.5):
                ax.axhline(_gla, color="#222244", lw=0.4)
            for _glo in _np_m.arange(round(xmin, 0), round(xmax, 0) + 1, 0.5):
                ax.axvline(_glo, color="#222244", lw=0.4)

        # ── P-wave ring (elapsed travel time since origin) ─────────────────────
        if et0 is not None:
            _elapsed = time.time() - et0
            if _elapsed > 0:
                _p_deg = (VP * _elapsed) / 111.1   # km → degrees latitude
                _circ  = _mpatch.Circle(
                    (elon, elat), _p_deg,
                    color="#00ddff", fill=False, lw=1.5, ls="--",
                    alpha=0.7, zorder=4, transform=ax.transData)
                ax.add_patch(_circ)

        # ── Stations ───────────────────────────────────────────────────────────
        for _st in states.values():
            _sc = "#00ff88" if _st.p_time is not None else "#ffcc00"
            ax.plot(_st.lon, _st.lat, "^", color=_sc, ms=9, zorder=5,
                    markeredgecolor="#000000", markeredgewidth=0.5)
            ax.annotate(
                _st.description[:10], (_st.lon, _st.lat),
                xytext=(4, 4), textcoords="offset points",
                fontsize=6, color="white",
                bbox=dict(boxstyle="round,pad=0.1", fc="#000000", alpha=0.55,
                          ec="none"))

        # ── Home location ──────────────────────────────────────────────────────
        ax.plot(HOME_LON, HOME_LAT, "s", color="white", ms=9, zorder=6,
                markeredgecolor="#aaaaaa", markeredgewidth=0.5)
        ax.annotate(
            HOME_LABEL, (HOME_LON, HOME_LAT),
            xytext=(5, -13), textcoords="offset points",
            fontsize=7, color="white", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.15", fc="#000000", alpha=0.6, ec="none"))

        # ── Epicenter star ─────────────────────────────────────────────────────
        _mag_s = f"M{med_ml:+.1f}" if med_ml is not None else "M?"
        ax.plot(elon, elat, "*", color="#ff2222", ms=22, zorder=7,
                markeredgecolor="white", markeredgewidth=0.6)
        ax.annotate(
            f"Epicenter\n{_mag_s}", (elon, elat),
            xytext=(10, 8), textcoords="offset points",
            fontsize=8, color="#ff5555", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.2", fc="#000000", alpha=0.7, ec="none"))

        # ── Axes styling ───────────────────────────────────────────────────────
        ax.set_xlabel("Longitude", color="#888888", fontsize=8)
        ax.set_ylabel("Latitude",  color="#888888", fontsize=8)
        ax.tick_params(colors="#777777", labelsize=7)
        for _sp in ax.spines.values():
            _sp.set_edgecolor("#333333")

        _ts_s = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime(et0)) if et0 else ""
        ax.set_title(f"Seismic Alert — {_mag_s}  {_ts_s}",
                     color="white", fontsize=10, pad=6, loc="left")

        # Legend
        ax.plot([], [], "^", color="#00ff88", ms=7, label="Station (P detected)")
        ax.plot([], [], "^", color="#ffcc00", ms=7, label="Station")
        ax.plot([], [], "*", color="#ff2222", ms=11, label="Epicenter")
        ax.plot([], [], "s", color="white",   ms=7,  label=HOME_LABEL)
        if et0 is not None:
            ax.plot([], [], "--", color="#00ddff", lw=1.5, label="P-wave front")
        ax.legend(loc="lower right", fontsize=6.5, framealpha=0.8,
                  facecolor="#111111", edgecolor="#333333", labelcolor="white")

        buf = _io.BytesIO()
        canvas.print_figure(buf, format="png", dpi=130,
                            facecolor=fig.get_facecolor(), bbox_inches="tight")
        buf.seek(0)
        return buf.read()

    except Exception as _exc:
        _log("MAP", f"_make_map_png failed: {_exc}")
        return None


_alarm_last_t = [0.0]

def _discord_major_spam(med_ml: float, city: str, level: str):
    """
    For major earthquakes (ML ≥ DISCORD_SPAM_ML), rapidly send DISCORD_SPAM_COUNT
    @everyone pings with the Drop/Cover/Hold On safety card image to both servers.
    Runs in a daemon thread so it never blocks detection.
    """
    if med_ml < DISCORD_SPAM_ML:
        return

    def _spam_worker():
        try:
            _sc_bytes = open(_SAFETY_CARD_PATH, "rb").read()
        except OSError:
            _sc_bytes = None

        for i in range(1, DISCORD_SPAM_COUNT + 1):
            _spam_title = (f"🚨 MAJOR EARTHQUAKE M{med_ml:+.1f} — ALERT {i}/{DISCORD_SPAM_COUNT}"
                           f" | {city}")
            _spam_body  = (
                f"⚠️  MAJOR EARTHQUAKE DETECTED  ⚠️\n"
                f"Magnitude  : M{med_ml:+.1f}\n"
                f"Location   : {city}\n"
                f"Alert      : {i} of {DISCORD_SPAM_COUNT}\n"
                f"\n"
                f"  1. DROP   — get to the ground immediately\n"
                f"  2. COVER  — under a table or protect your head\n"
                f"  3. HOLD ON — until all shaking stops\n"
                f"\n"
                f"DO NOT run outside during shaking.\n"
                f"Expect strong aftershocks.\n"
                f"https://myearthquake.dpdns.org/")
            _send_discord(_spam_title, _spam_body,
                          level=level,
                          img_bytes=_sc_bytes,
                          bold_header="🚨 DROP!  COVER!  HOLD ON! 🚨",
                          everyone=True)
            if i < DISCORD_SPAM_COUNT:
                time.sleep(DISCORD_SPAM_DELAY)

    threading.Thread(target=_spam_worker, daemon=True).start()


def _play_alarm(level: str):
    """Play a system audio alarm on macOS/Linux for HIGH or EXTREME events."""
    if not SOUND_ALARM_ENABLED:
        return
    levels = ["LOW", "MEDIUM", "HIGH", "EXTREME"]
    if levels.index(level) < levels.index(SOUND_ALARM_LEVEL):
        return
    now = time.time()
    if now - _alarm_last_t[0] < 30:   # don't repeat within 30s
        return
    _alarm_last_t[0] = now
    def _play():
        try:
            import subprocess, platform
            if platform.system() == "Darwin":
                # Play system alert sound — repeat for EXTREME
                reps = 3 if level == "EXTREME" else 1
                for _ in range(reps):
                    subprocess.run(["afplay", "/System/Library/Sounds/Sosumi.aiff"],
                                   timeout=5, capture_output=True)
            elif platform.system() == "Linux":
                subprocess.run(["paplay", "/usr/share/sounds/freedesktop/stereo/alarm-clock-elapsed.oga"],
                               timeout=5, capture_output=True)
        except Exception as exc:
            _log("ALARM ERR", str(exc))
    threading.Thread(target=_play, daemon=True).start()


_mainshock_time         = [None]   # epoch of last confirmed mainshock
_mainshock_ml           = [None]   # magnitude of last mainshock
_as_forecast_sent       = [False]  # aftershock forecast sent for current event
_usgs_crosscheck_spawned = [False] # USGS crosscheck thread spawned for current event

def _omori_utsu_prob(ml_main, hours=AFTERSHOCK_WINDOW_HOURS):
    """
    Estimate probability of at least one M≥(ml_main-1) aftershock in the
    next `hours` hours using modified Omori-Utsu law (Reasenberg & Jones 1989).
    Returns probability 0.0–1.0.
    """
    if ml_main is None:
        return None
    # R&J parameters for California:
    a, b, p, c = -1.67, 0.91, 1.08, 0.05
    m_thresh = ml_main - 1.0    # forecast threshold
    # Expected number of aftershocks in [0, T] hours:
    #   N = 10^(a + b*(M_main - m_thresh)) * integral_0^T (t+c)^(-p) dt
    rate = 10 ** (a + b * (ml_main - m_thresh))
    if p == 1.0:
        integral = math.log((hours + c) / c)
    else:
        integral = ((hours + c) ** (1 - p) - c ** (1 - p)) / (1 - p)
    expected = rate * integral
    # Poisson probability of ≥1 event:
    prob = 1.0 - math.exp(-max(0.0, expected))
    return round(prob, 3)


_db_conn = [None]
_db_lock = threading.Lock()

def _init_catalog_db():
    """Create SQLite database and events table if not exists."""
    if not CATALOG_ENABLED:
        return
    import sqlite3
    try:
        conn = sqlite3.connect(CATALOG_DB_PATH, check_same_thread=False)
        conn.execute("""
            CREATE TABLE IF NOT EXISTS events (
                id          INTEGER PRIMARY KEY AUTOINCREMENT,
                detected_at REAL,
                origin_time REAL,
                lat         REAL,
                lon         REAL,
                depth_km    REAL,
                ml          REAL,
                ml_unc      REAL,
                n_stations  INTEGER,
                dist_home   REAL,
                location    TEXT,
                event_type  TEXT,
                rms_sec     REAL,
                confirmed   INTEGER DEFAULT 0
            )
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_origin ON events(origin_time)")
        conn.commit()
        _db_conn[0] = conn
        _log("CATALOG", f"Database ready: {CATALOG_DB_PATH}")
    except Exception as exc:
        _log("CATALOG ERR", str(exc))

def _catalog_log_event(origin_time, lat, lon, depth_km, ml, ml_unc,
                        n_stations, dist_home, location, event_type="EQ",
                        rms_sec=None, confirmed=False):
    """Insert or update an event in the catalog."""
    if not CATALOG_ENABLED or _db_conn[0] is None:
        return
    import sqlite3
    with _db_lock:
        try:
            _db_conn[0].execute("""
                INSERT INTO events
                (detected_at, origin_time, lat, lon, depth_km, ml, ml_unc,
                 n_stations, dist_home, location, event_type, rms_sec, confirmed)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
            """, (time.time(), origin_time, lat, lon, depth_km, ml, ml_unc,
                  n_stations, dist_home, location, event_type,
                  rms_sec, int(confirmed)))
            _db_conn[0].commit()
            _log("CATALOG", f"Logged {event_type} M{ml:+.1f} near {location}")
        except Exception as exc:
            _log("CATALOG ERR", str(exc))

def _catalog_recent(hours=168):
    """Return list of events from the last `hours` hours as dicts."""
    if not CATALOG_ENABLED or _db_conn[0] is None:
        return []
    import sqlite3
    cutoff = time.time() - hours * 3600
    with _db_lock:
        try:
            cur = _db_conn[0].execute(
                "SELECT * FROM events WHERE detected_at > ? ORDER BY origin_time DESC",
                (cutoff,))
            cols = [d[0] for d in cur.description]
            return [dict(zip(cols, row)) for row in cur.fetchall()]
        except Exception:
            return []

def _bvalue_and_rate(hours=168):
    """
    Compute Gutenberg-Richter b-value and event rate from recent catalog.
    Returns (b, rate_per_day, n_events) or (None, None, 0).
    """
    events = _catalog_recent(hours=hours)
    mls = [e["ml"] for e in events if e["ml"] is not None and e["ml"] >= 0.0]
    n = len(mls)
    rate = n / (hours / 24.0) if hours > 0 else 0
    if n < BVALUE_MIN_EVENTS:
        return None, round(rate, 2), n
    # Maximum likelihood b-value: b = log10(e) / (mean_M - Mc)
    # Use Mc = min observed magnitude
    mc   = min(mls)
    mean_m = sum(mls) / n
    if mean_m <= mc:
        return None, round(rate, 2), n
    import math as _m
    b = _m.log10(math.e) / (mean_m - mc)
    return round(b, 2), round(rate, 2), n


# ═══════════════════════════════════════════════════════════════════════════════
# USGS EARTHQUAKE API CROSS-CHECK
# Polls the USGS ComCat API ~60 s after a local detection to find any matching
# official event and posts official magnitude / location to Discord + NTFY.
# ═══════════════════════════════════════════════════════════════════════════════

_usgs_checked_events = set()   # set of (rounded_lat, rounded_lon, rounded_t0)

def _fetch_usgs_shakemap(usgs_id, max_attempts=4, retry_wait=90):
    """
    Download the USGS ShakeMap intensity image for a given event ID.
    Retries up to max_attempts times (ShakeMap may not exist immediately).
    Returns (image_bytes, filename) or (None, None) if unavailable.
    """
    import urllib.request as _ur, json as _js, ssl as _ssl
    _ctx = _ssl.create_default_context()
    _prefer = [
        "download/intensity.jpg",
        "download/intensity.png",
        "download/pga.jpg",
        "download/pgv.jpg",
    ]
    for attempt in range(max_attempts):
        try:
            _detail_url = (f"https://earthquake.usgs.gov/fdsnws/event/1/query"
                           f"?eventid={usgs_id}&format=geojson")
            _req = _ur.Request(_detail_url,
                               headers={"User-Agent": "QuakeAlertBot/1.0"})
            with _ur.urlopen(_req, timeout=20, context=_ctx) as _r:
                _ev = _js.loads(_r.read())
            _sms = _ev.get("properties", {}).get("products", {}).get("shakemap", [])
            if not _sms:
                _log("USGS", f"ShakeMap not published yet (attempt {attempt+1}/{max_attempts})")
            else:
                _contents = _sms[0].get("contents", {})
                for _key in _prefer:
                    if _key in _contents:
                        _img_url = _contents[_key]["url"]
                        _ir = _ur.Request(_img_url,
                                          headers={"User-Agent": "QuakeAlertBot/1.0"})
                        with _ur.urlopen(_ir, timeout=30, context=_ctx) as _img_r:
                            _bytes = _img_r.read()
                        _log("USGS", f"ShakeMap downloaded: {_key}"
                                     f"  ({len(_bytes)//1024} KB)")
                        return _bytes, _key
        except Exception as _fe:
            _log("USGS ERR", f"ShakeMap fetch attempt {attempt+1}: {_fe}")
        if attempt < max_attempts - 1:
            time.sleep(retry_wait)
    return None, None


def _usgs_crosscheck(est_lat, est_lon, est_t0, local_ml):
    """
    Called in a background thread after EQ confirmation.
    Waits 60 s then queries USGS ComCat for events near our estimate.
    Posts a comparison message to Discord + NTFY + Email if a match is found.
    Also fetches the ShakeMap intensity image and emails it.
    """
    time.sleep(60)   # give USGS time to process and publish
    try:
        import urllib.request, json as _json, ssl, urllib.parse
        _ssl_ctx = ssl.create_default_context()

        t_start = time.strftime("%Y-%m-%dT%H:%M:%S",
                                time.gmtime(est_t0 - 120))
        t_end   = time.strftime("%Y-%m-%dT%H:%M:%S",
                                time.gmtime(est_t0 + 300))
        params  = urllib.parse.urlencode({
            "format":    "geojson",
            "starttime": t_start,
            "endtime":   t_end,
            "latitude":  f"{est_lat:.3f}",
            "longitude": f"{est_lon:.3f}",
            "maxradius": 2.0,       # degrees (~220 km)
            "minmagnitude": 0.5,
            "orderby":  "time",
        })
        url = f"https://earthquake.usgs.gov/fdsnws/event/1/query?{params}"
        req = urllib.request.Request(url,
                                     headers={"User-Agent": "QuakeAlertBot/1.0"})
        with urllib.request.urlopen(req, timeout=20, context=_ssl_ctx) as r:
            data = _json.loads(r.read())

        feats = data.get("features", [])
        if not feats:
            _log("USGS", "No matching event found in catalog yet")
            return

        # Pick closest in time to our origin estimate
        best = min(feats,
                   key=lambda f: abs(f["properties"]["time"] / 1000.0 - est_t0))
        bp   = best["properties"]
        bc   = best["geometry"]["coordinates"]   # [lon, lat, depth]

        usgs_ml  = bp.get("mag")
        usgs_loc = bp.get("place", "unknown")
        usgs_t   = bp["time"] / 1000.0
        usgs_lat = bc[1]; usgs_lon = bc[0]; usgs_dep = bc[2]
        usgs_utc = time.strftime("%H:%M:%S UTC", time.gmtime(usgs_t))
        usgs_url = bp.get("url", "")
        mag_type = bp.get("magType", "M")

        # Deduplicate — don't re-notify for the same USGS event
        _evt_key = (round(usgs_lat, 1), round(usgs_lon, 1), round(usgs_t, -1))
        if _evt_key in _usgs_checked_events:
            return
        _usgs_checked_events.add(_evt_key)

        usgs_id = best.get("id")
        with _eq_id_lock:
            _current_eq_usgs[0] = usgs_id

        ml_diff  = (usgs_ml - local_ml) if (usgs_ml and local_ml) else None
        diff_str = (f"{ml_diff:+.1f}" if ml_diff is not None else "n/a")
        dist_off = haversine_km(est_lat, est_lon, usgs_lat, usgs_lon)

        title = (f"✅ USGS CONFIRMED | {mag_type}{usgs_ml:+.1f} | {usgs_loc}"
                 if usgs_ml else f"✅ USGS CONFIRMED | {usgs_loc}")
        body  = (
            f"---- USGS OFFICIAL ----\n"
            f"Magnitude  : {mag_type} {usgs_ml:+.1f}\n"
            f"Location   : {usgs_loc}\n"
            f"Origin time: {usgs_utc}\n"
            f"Depth      : {usgs_dep:.1f} km\n"
            f"Coords     : {usgs_lat:.3f}°N  {usgs_lon:.3f}°W\n"
            f"\n"
            f"---- VS LOCAL ESTIMATE ----\n"
            f"Local ML   : M{local_ml:+.1f}  (diff: {diff_str})\n"
            f"Loc offset : {dist_off:.1f} km\n"
            f"\n"
            f"Details    : {usgs_url}")

        _log("USGS", f"Match: {mag_type}{usgs_ml:+.1f}  {usgs_loc}  "
                     f"offset={dist_off:.1f}km  ml_diff={diff_str}")
        _send_ntfy(title, body, priority="default")
        _usgs_upd_n, _ = _eq_next_update()
        _usgs_dc_prefix = _eq_discord_prefix(f"Upd #{_usgs_upd_n} | USGS")
        _usgs_dc_title  = (f"{_usgs_dc_prefix} ✅ USGS Confirmed"
                           f" | {mag_type}{usgs_ml:+.1f} | {usgs_loc}"
                           if usgs_ml else
                           f"{_usgs_dc_prefix} ✅ USGS Confirmed | {usgs_loc}")
        _usgs_dc_bold   = (f"{mag_type}{usgs_ml:+.1f} | {usgs_loc}"
                           f" | Origin: {usgs_utc}"
                           if usgs_ml else
                           f"{usgs_loc} | Origin: {usgs_utc}")
        _ws_broadcast_event("USGS Confirmed", location=usgs_loc, ml=usgs_ml,
            detail=f"Local diff: {diff_str}  loc offset: {dist_off:.1f}km")

        # ── Fetch ShakeMap (retries up to 4× while USGS processes it) ────────
        _sm_bytes, _sm_key = _fetch_usgs_shakemap(usgs_id)
        _sm_note = (f"\nShakeMap: {'attached' if _sm_bytes else 'not yet available'}"
                    f" ({_sm_key or '—'})")

        # ── Email with ShakeMap attached ──────────────────────────────────────
        _email_subj = (f"USGS {mag_type}{usgs_ml:+.1f} — {usgs_loc}"
                       if usgs_ml else f"USGS Confirmed — {usgs_loc}")
        _send_email(_email_subj, body + _sm_note, img_bytes=_sm_bytes)
        _log("USGS", f"Email queued: ShakeMap={'yes' if _sm_bytes else 'no'}")

        # ── Discord with ShakeMap image ───────────────────────────────────────
        _send_discord(_usgs_dc_title, body + _sm_note,
                      bold_header=_usgs_dc_bold, img_bytes=_sm_bytes)

    except Exception as exc:
        _log("USGS ERR", str(exc))


# ═══════════════════════════════════════════════════════════════════════════════
# MYEARTHQUAKE QUAKE REPORT  —  iPhone 13-sized event summary image
# ═══════════════════════════════════════════════════════════════════════════════

def _generate_myearthquake_report(snap, timeline_copy):
    """
    Generate the MyEarthquake Quake Report — an iPhone 13-sized (390x844 pt)
    PNG earthquake summary image.  Returns raw PNG bytes, or None on error.

    Layout (bottom-up in matplotlib data coordinates 0,0 → 390,844):
        footer | cities | wave arrivals | MMI | details grid |
        origin | magnitude hero | app header | status bar
    """
    import io as _io
    import matplotlib.figure as _mfig
    import matplotlib.patches as _mpa
    from matplotlib.backends.backend_agg import FigureCanvasAgg as _Agg

    try:
        # ── Pull data ────────────────────────────────────────────────────────
        elat    = snap.get("elat");    elon    = snap.get("elon")
        edepth  = snap.get("edepth");  edist   = snap.get("edist_home")
        med_ml  = snap.get("med_ml");  en      = snap.get("en", 0)
        erms    = snap.get("erms");    eaz_gap = snap.get("eaz_gap")
        et0     = snap.get("et0")

        city_info = _nearest_city(elat, elon) if elat is not None else None
        city_lbl  = _city_label(elat, elon)   if elat is not None else None

        # Nearby cities sorted by distance (up to 6)
        nearby = []
        if elat is not None and elon is not None:
            for _cn, _clat, _clon in _CITIES:
                _d = haversine_km(elat, elon, _clat, _clon)
                if _d < 500:
                    nearby.append((_d, _cn))
            nearby.sort()
            nearby = nearby[:6]

        ml_val    = med_ml if med_ml is not None else 0.0
        ml_str_r  = f"M{med_ml:+.1f}" if med_ml is not None else "M ?"
        depth_str = f"{_eff_depth:.1f} km"  # use clamped depth (≥1 km)
        dist_str  = f"{edist:.0f} km"   if edist  is not None else "—"
        if elat is not None:
            lat_s = f"{abs(elat):.4f} {'N' if elat >= 0 else 'S'}"
            lon_s = f"{abs(elon):.4f} {'E' if elon >= 0 else 'W'}"
            loc_str = f"{lat_s}   {lon_s}"
        else:
            loc_str = "Not determined"

        time_utc   = (time.strftime("%Y-%m-%d  %H:%M:%S UTC", time.gmtime(et0))
                      if et0 else "Unknown")
        time_local = (time.strftime("%b %d, %Y  %I:%M:%S %p", time.localtime(et0))
                      if et0 else "Unknown")

        # Effective depth: clamp to 1 km minimum (solver can output 0 for very
        # shallow events; 0 km depth is unphysical for attenuation calculations)
        _eff_depth = max(edepth or 10.0, 1.0)

        # Max observed MMI — prefer ML-estimated PGV at station distance over the
        # raw P-onset peak counts, which systematically underestimate shaking by
        # missing the larger S-wave window.
        _obs_pgvs = []
        for _s in states.values():
            if _s.p_time is None:          # station didn't detect this event
                continue
            # Primary: ML estimate + actual station-to-hypocentre distance
            _dist_s = _s.event_dist_km or _s.dist_km or None
            if _s.ml_est is not None and _dist_s and _dist_s > 0:
                _pgv_s = pgv_at_dist(_s.ml_est, _dist_s)
                if _pgv_s:
                    _obs_pgvs.append(_pgv_s)
                    continue
            # Fallback: raw Wood-Anderson peak (less accurate)
            if _s.event_peak > 0 and _s.sensitivity > 0:
                _obs_pgvs.append(counts_to_pgv(_s.event_peak, _s.sensitivity))
        obs_mmi               = pgv_to_mmi(max(_obs_pgvs)) if _obs_pgvs else None
        mmi_obs_str, mmi_obs_col = mmi_label(obs_mmi)

        # MMI at home: use hypocentral distance (includes depth) for accuracy
        _hypo_home = (math.sqrt(edist**2 + _eff_depth**2) if edist else None)
        pgv_home   = (pgv_at_dist(ml_val, _hypo_home)
                      if (med_ml is not None and _hypo_home and _hypo_home > 0)
                      else None)
        mmi_home               = pgv_to_mmi(pgv_home) if pgv_home else None
        mmi_home_str, mmi_home_col = mmi_label(mmi_home)

        # Epicenter accuracy
        acc_label, _acc_col, pos_err = _epi_accuracy(en, erms, eaz_gap)
        acc_str = f"{acc_label}  ±{pos_err:.0f} km" if pos_err else acc_label

        # Station wave arrival rows
        sta_rows = []
        for _s in states.values():
            with _s.lock:
                _pt = _s.p_time; _st2 = _s.s_time
                _ml2 = _s.ml_est; _desc = _s.description
            _pt_s = time.strftime("%H:%M:%S", time.gmtime(_pt))  if _pt  else "—"
            _st_s = time.strftime("%H:%M:%S", time.gmtime(_st2)) if _st2 else "—"
            _sp_s = f"{_st2 - _pt:.1f}s"  if (_pt and _st2) else "—"
            _ml_s = f"M{_ml2:+.1f}"       if _ml2            else "—"
            sta_rows.append((_desc or _s.label, _pt_s, _st_s, _sp_s, _ml_s))

        # ── Canvas: iPhone 13 logical size 390x844 pt @ 200 dpi ─────────────
        _W, _H = 390, 844
        _PAD   = 12
        _R     = 7    # card corner radius

        fig_r  = _mfig.Figure(figsize=(3.90, 8.44), facecolor="#0b0b10")
        _cv    = _Agg(fig_r)
        _ax    = fig_r.add_axes([0, 0, 1, 1], facecolor="#0b0b10")
        _ax.set_xlim(0, _W); _ax.set_ylim(0, _H); _ax.axis("off")

        # ── Drawing helpers ──────────────────────────────────────────────────
        def _card(y_bot, height, x=_PAD, w=_W - 2 * _PAD, color="#181820"):
            """Rounded card rectangle (FancyBboxPatch compensated for pad)."""
            _ax.add_patch(_mpa.FancyBboxPatch(
                (x + _R, y_bot + _R), w - 2 * _R, height - 2 * _R,
                boxstyle=f"round,pad={_R}",
                facecolor=color, edgecolor="#26263a", linewidth=0.7, zorder=2))

        def _t(x, y, s, size=9, color="#e8e8f8", ha="left", va="bottom",
               weight="normal", zorder=4, alpha=1.0):
            _ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va,
                     fontweight=weight, alpha=alpha, zorder=zorder,
                     transform=_ax.transData)

        def _sec(x, y, label):
            _t(x, y, label, size=6.5, color="#4444bb", weight="bold")

        # ── Section geometry (bottom of each section, height) ────────────────
        _footer_y,   _footer_h   = 0,   20
        _cities_y,   _cities_h   = 24,  68
        _sta_y,      _sta_h      = 96,  154
        _mmi_y,      _mmi_h      = 254, 108
        _det_y,      _det_h      = 366, 130
        _loc_y,      _loc_h      = 500, 96
        _mag_y,      _mag_h      = 600, 146
        _hdr_y,      _hdr_h      = 750, 56
        _sb_y,       _sb_h       = 810, 34   # status bar → top = 844

        # ── STATUS BAR ───────────────────────────────────────────────────────
        _ax.add_patch(_mpa.Rectangle((0, _sb_y), _W, _sb_h,
                      facecolor="#070710", edgecolor="none", zorder=1))
        _t(_W / 2, _sb_y + _sb_h / 2, time_local,
           size=7.5, color="#666688", ha="center", va="center")

        # ── APP HEADER ───────────────────────────────────────────────────────
        _ax.add_patch(_mpa.Rectangle((0, _hdr_y), _W, _hdr_h,
                      facecolor="#0d0d18", edgecolor="none", zorder=1))
        # blue accent stripe at top of header
        _ax.add_patch(_mpa.Rectangle((0, _hdr_y + _hdr_h - 3), _W, 3,
                      facecolor="#2233cc", edgecolor="none", zorder=3))
        _t(_W / 2, _hdr_y + _hdr_h - 12, "MyEarthquake",
           size=16, color="#d0d0ff", ha="center", va="top", weight="bold")
        _t(_W / 2, _hdr_y + 10, "QUAKE REPORT",
           size=7.5, color="#3344bb", ha="center", va="bottom", weight="bold")

        # ── MAGNITUDE HERO ───────────────────────────────────────────────────
        if   ml_val >= 7.0: _mc = "#ff0033"
        elif ml_val >= 6.0: _mc = "#ff3300"
        elif ml_val >= 5.0: _mc = "#ff6600"
        elif ml_val >= 4.0: _mc = "#ffaa00"
        elif ml_val >= 3.0: _mc = "#eedd00"
        elif ml_val >= 2.0: _mc = "#88dd22"
        else:               _mc = "#44cc44"

        _card(_mag_y, _mag_h, color="#13131c")
        _t(_W / 2, _mag_y + _mag_h - 14,
           ml_str_r, size=58, color=_mc, ha="center", va="top", weight="bold")
        _loc_disp = city_lbl or (city_info[0] if city_info else "Location undetermined")
        _t(_W / 2, _mag_y + 40, _loc_disp,
           size=9.5, color="#ccccee", ha="center", va="bottom", weight="bold")
        _t(_W / 2, _mag_y + 22, "LOCAL MAGNITUDE  (Hutton & Boore 1987)",
           size=6, color="#44446a", ha="center", va="bottom")

        # ── ORIGIN CARD ──────────────────────────────────────────────────────
        _card(_loc_y, _loc_h, color="#121222")
        _sec(_PAD + 8, _loc_y + _loc_h - 10, "ORIGIN")
        _t(_PAD + 8,    _loc_y + _loc_h - 25, time_utc,  size=8.5, color="#e8e8f8", va="top")
        _t(_PAD + 8,    _loc_y + _loc_h - 42, loc_str,   size=7.5, color="#9999bb", va="top")
        _t(_W - _PAD - 8, _loc_y + _loc_h - 25,
           f"Depth  {depth_str}", size=8.5, color="#e8e8f8", ha="right", va="top")
        _t(_W - _PAD - 8, _loc_y + _loc_h - 42,
           f"{dist_str} from {HOME_LABEL}", size=7.5, color="#9999bb", ha="right", va="top")
        _ax.plot([_PAD + 8, _W - _PAD - 8],
                 [_loc_y + _loc_h - 52, _loc_y + _loc_h - 52],
                 color="#1e1e32", lw=0.6, zorder=4)
        _t(_PAD + 8, _loc_y + 10,
           f"Local:  {time.strftime('%b %d  %I:%M:%S %p', time.localtime(et0)) if et0 else '—'}",
           size=8, color="#7788bb", va="bottom")

        # ── DETAILS 2x2 GRID ─────────────────────────────────────────────────
        _card(_det_y, _det_h, color="#121222")
        _sec(_PAD + 8, _det_y + _det_h - 10, "SEISMIC DETAILS")
        _cw = (_W - 2 * _PAD) / 2
        _grid = [
            ("STATIONS USED",  str(en),
             _PAD + 8,          _det_y + _det_h - 30),
            ("RMS RESIDUAL",   f"{'%.2f s' % erms if erms else '—'}",
             _PAD + 8 + _cw,    _det_y + _det_h - 30),
            ("AZIMUTH GAP",    f"{'%.0f deg' % eaz_gap if eaz_gap else '—'}",
             _PAD + 8,          _det_y + 18),
            ("ACCURACY",       acc_str,
             _PAD + 8 + _cw,    _det_y + 18),
        ]
        for _gl, _gv, _gx, _gy in _grid:
            _t(_gx, _gy + 26, _gl, size=6.2, color="#44446a", va="bottom", weight="bold")
            _t(_gx, _gy,      _gv, size=9.5, color="#d8d8f0", va="bottom")
        _mx = _PAD + 8 + _cw - 4
        _my = (_det_y + 18 + _det_y + _det_h - 30) / 2 + 14
        _ax.plot([_mx, _mx], [_det_y + 8, _det_y + _det_h - 14],
                 color="#1e1e30", lw=0.7, zorder=4)
        _ax.plot([_PAD + 8, _W - _PAD - 8], [_my, _my],
                 color="#1e1e30", lw=0.7, zorder=4)

        # ── MMI INTENSITY CARD ───────────────────────────────────────────────
        _card(_mmi_y, _mmi_h, color="#121222")
        _sec(_PAD + 8, _mmi_y + _mmi_h - 10, "SHAKING INTENSITY (MMI)")

        _bx  = _PAD + 8
        _bw  = _W - 2 * _PAD - 16
        _bh  = 16           # bar height
        _by0 = _mmi_y + 50  # bar bottom (leaves room for 2 text rows below)
        _mmi_scale = [
            (1, "#333333"), (2, "#9999aa"), (3, "#99bbff"),
            (4, "#55ddff"), (5, "#eeee33"), (6, "#ffaa22"),
            (7, "#ff7700"), (8, "#ff3300"), (9, "#cc0000"),
        ]
        _sw = _bw / len(_mmi_scale)
        for _si, (_mv, _mc2) in enumerate(_mmi_scale):
            _ax.add_patch(_mpa.Rectangle(
                (_bx + _si * _sw, _by0), _sw, _bh,
                facecolor=_mc2, edgecolor="#0b0b10", linewidth=0.3, zorder=4))
            _t(_bx + _si * _sw + _sw / 2, _by0 + _bh + 2,
               str(_mv), size=5.5, color="#55557a", ha="center", va="bottom", zorder=5)

        # White downward-pointing triangle = Max Observed MMI (above bar)
        if obs_mmi is not None:
            _frac = min(max(obs_mmi - 1, 0), 8) / 8
            _ix = _bx + _frac * _bw + _sw / 2
            _ax.add_patch(_mpa.Polygon(
                [[_ix - 5, _by0 + _bh + 1],
                 [_ix + 5, _by0 + _bh + 1],
                 [_ix,     _by0 + _bh - 4]],
                closed=True, facecolor="#ffffff", edgecolor="none", zorder=6))

        # Cyan upward-pointing triangle = Home MMI (below bar)
        if mmi_home is not None:
            _frac_h = min(max(mmi_home - 1, 0), 8) / 8
            _ihx = _bx + _frac_h * _bw + _sw / 2
            _ax.add_patch(_mpa.Polygon(
                [[_ihx - 5, _by0 - 1],
                 [_ihx + 5, _by0 - 1],
                 [_ihx,     _by0 + 5]],
                closed=True, facecolor="#44aaff", edgecolor="none", zorder=6))

        # ── Two separate rows below the bar ────────────────────────────────────
        # Row 1 (upper): Max Observed   Row 2 (lower): Home MMI
        # Each row: coloured bullet + label + MMI string
        _row1_y = _mmi_y + 29
        _row2_y = _mmi_y + 12

        # White bullet = Max Observed (matches white triangle above bar)
        _ax.add_patch(_mpa.Circle((_PAD + 13, _row1_y + 3.5), 3,
                      facecolor="#cccccc", edgecolor="none", zorder=5))
        _t(_PAD + 20, _row1_y,
           "Max Observed:", size=7, color="#888899", va="bottom")
        _t(_PAD + 100, _row1_y,
           mmi_obs_str, size=7.5, color=mmi_obs_col, va="bottom", weight="bold")

        # Cyan bullet = Home (matches cyan triangle below bar)
        _ax.add_patch(_mpa.Circle((_PAD + 13, _row2_y + 3.5), 3,
                      facecolor="#44aaff", edgecolor="none", zorder=5))
        _home_short = HOME_LABEL.split(",")[0]   # "San Ramon" not "San Ramon, CA"
        _t(_PAD + 20, _row2_y,
           f"{_home_short}:", size=7, color="#888899", va="bottom")
        _t(_PAD + 100, _row2_y,
           mmi_home_str, size=7.5, color=mmi_home_col, va="bottom", weight="bold")

        # ── WAVE ARRIVALS TABLE ──────────────────────────────────────────────
        _card(_sta_y, _sta_h, color="#121222")
        _sec(_PAD + 8, _sta_y + _sta_h - 10, "WAVE ARRIVALS")

        _col_x   = [_PAD + 10, 108, 182, 258, 316]
        _col_lbl = ["STATION", "P-WAVE", "S-WAVE", "S-P", "ML"]
        _hdr_row_y = _sta_y + _sta_h - 26
        for _cx2, _cl in zip(_col_x, _col_lbl):
            _t(_cx2, _hdr_row_y, _cl, size=6, color="#33336a", va="bottom", weight="bold")
        _ax.plot([_PAD + 8, _W - _PAD - 8],
                 [_hdr_row_y - 1, _hdr_row_y - 1], color="#1e1e34", lw=0.6, zorder=4)

        _n_rows  = min(len(sta_rows), 6)
        _avail   = _sta_h - 34
        _rh      = _avail / max(_n_rows, 1)
        for _ri, (_sdesc, _pt_s2, _st_s2, _sp_s2, _ml_s2) in enumerate(sta_rows[:6]):
            _ry = _sta_y + _sta_h - 28 - (_ri + 0.5) * _rh
            if _ri % 2 == 0:
                _ax.add_patch(_mpa.Rectangle(
                    (_PAD + 4, _ry - _rh * 0.45), _W - 2 * _PAD - 8, _rh * 0.9,
                    facecolor="#15151e", edgecolor="none", zorder=2))
            _sn  = (_sdesc.split(",")[0] if "," in _sdesc else _sdesc)[:14]
            _pc  = "#4499ff" if _pt_s2 != "—" else "#333355"
            _sc  = "#ff9922" if _st_s2 != "—" else "#333355"
            _nc  = "#aaaacc" if _sp_s2 != "—" else "#333355"
            _mlc = "#88ee88" if _ml_s2 != "—" else "#333355"
            _t(_col_x[0], _ry, _sn,    size=7, color="#b0b0cc", va="center")
            _t(_col_x[1], _ry, _pt_s2, size=7, color=_pc,       va="center")
            _t(_col_x[2], _ry, _st_s2, size=7, color=_sc,       va="center")
            _t(_col_x[3], _ry, _sp_s2, size=7, color=_nc,       va="center")
            _t(_col_x[4], _ry, _ml_s2, size=7, color=_mlc,      va="center")

        # ── NEARBY CITIES ────────────────────────────────────────────────────
        _card(_cities_y, _cities_h, color="#121222")
        _sec(_PAD + 8, _cities_y + _cities_h - 10, "NEARBY CITIES")
        if nearby:
            _city_col_x = [_PAD + 10, _W // 2 + 4]
            for _ci, (_cd, _cn2) in enumerate(nearby[:4]):
                _cxc = _city_col_x[_ci % 2]
                _ryc = _cities_y + _cities_h - 26 - (_ci // 2) * 22
                _t(_cxc, _ryc,      _cn2,              size=7.5, color="#9999bb", va="bottom")
                _t(_cxc, _ryc - 12, f"  {_cd:.0f} km", size=6.5, color="#555577", va="bottom")
        else:
            _t(_PAD + 8, _cities_y + _cities_h / 2,
               "No nearby cities in database", size=8, color="#444455", va="center")

        # ── FOOTER ───────────────────────────────────────────────────────────
        _ax.add_patch(_mpa.Rectangle((0, _footer_y), _W, _footer_h,
                      facecolor="#070710", edgecolor="none", zorder=1))
        _t(_W / 2, _footer_h / 2,
           "MyEarthquake  |  Local Seismic Monitoring  |  Hutton & Boore ML",
           size=5.5, color="#2a2a44", ha="center", va="center")

        # ── Render at 2x iPhone scale → 780x1688 px ─────────────────────────
        _buf = _io.BytesIO()
        _cv.print_figure(_buf, format="png", dpi=200, facecolor="#0b0b10")
        _buf.seek(0)
        return _buf.read()

    except Exception as _exc:
        _log("MYEQ REPORT", f"Generation failed: {_exc}")
        import traceback as _tb; _tb.print_exc()
        return None


def _send_final_report(snap, timeline_copy):
    """
    Generate and send the post-event final report via ntfy AND email.

    ntfy: always sent (no email config required).
    Email: only sent when EMAIL_ENABLED and credentials are configured.

    Runs entirely in a daemon background thread.
    """
    # Guard: need at least a magnitude or epicenter — no snap = nothing to report
    if not snap:
        return

    def _worker():
        import io, uuid as _uuid
        import matplotlib.figure as mfig
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.patches import Polygon as _MplPoly

        elat    = snap.get("elat");    elon    = snap.get("elon")
        edepth  = snap.get("edepth");  edist   = snap.get("edist_home")
        med_ml  = snap.get("med_ml");  en      = snap.get("en", 0)
        erms    = snap.get("erms");    eaz_gap = snap.get("eaz_gap")
        et0     = snap.get("et0")

        city_info = _nearest_city(elat, elon)
        city_lbl_final = _city_label(elat, elon)   # "40 km E of Davis, CA"
        city_str  = (f"{city_info[0]} ({city_info[1]:.0f} km away)"
                     if city_info else "location unknown")

        ml_str    = f"M{med_ml:+.1f}" if med_ml is not None else "M?"
        depth_str = f"{edepth:.1f} km" if edepth is not None else "unknown"
        dist_str  = f"{edist:.0f} km"  if edist  is not None else "unknown"
        loc_str   = (f"{elat:.4f}°N  {elon:.4f}°W"
                     if elat is not None else "not determined")
        orig_str  = (time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime(et0))
                     if et0 else "unknown")

        acc_label, _, pos_err = _epi_accuracy(en, erms, eaz_gap)
        acc_str = (f"{acc_label} (±{pos_err:.0f} km)" if pos_err else acc_label)

        # ── Timeline text ──────────────────────────────────────────────────────
        tl_text = ("\n".join(f"  {ts}  {note}" for ts, note in timeline_copy)
                   if timeline_copy else "  (no events recorded)")

        # ── Station table ──────────────────────────────────────────────────────
        sta_lines = []
        for _s in states.values():
            with _s.lock:
                _pt2 = _s.p_time; _st2 = _s.s_time; _ml2 = _s.ml_est
            _pt_s = time.strftime("%H:%M:%S", time.gmtime(_pt2)) if _pt2 else "—"
            _st_s = time.strftime("%H:%M:%S", time.gmtime(_st2)) if _st2 else "—"
            _sp_s = f"{_st2-_pt2:.1f}s" if (_pt2 and _st2) else "—"
            _ml_s = f"ML{_ml2:+.1f}" if _ml2 else "—"
            sta_lines.append(
                f"  {_s.label:<22}  P={_pt_s}  S={_st_s}  S-P={_sp_s:<6}  {_ml_s}")
        sta_text = "\n".join(sta_lines)

        body = (
            f"FINAL EARTHQUAKE REPORT\n"
            f"{'='*58}\n"
            f"Magnitude     : {ml_str}\n"
            f"Epicenter     : {city_lbl_final or city_str}\n"
            f"Nearest city  : {city_str}\n"
            f"Coordinates   : {loc_str}\n"
            f"Depth         : {depth_str}\n"
            f"Distance      : {dist_str} from {HOME_LABEL}\n"
            f"Origin time   : {orig_str}\n"
            f"Accuracy      : {acc_str}\n"
            f"Stations used : {en}\n"
            f"RMS residual  : {'%.2f s' % erms if erms else '—'}\n"
            f"Azimuth gap   : {'%.0f°' % eaz_gap if eaz_gap else '—'}\n"
            f"\n"
            f"STATION DETECTIONS\n"
            f"{'-'*58}\n"
            f"{sta_text}\n"
            f"\n"
            f"EVENT TIMELINE\n"
            f"{'-'*58}\n"
            f"{tl_text}\n"
        )

        # ── ntfy push notification (always sent, no email config required) ────
        try:
            _ntfy_pri = "urgent" if (med_ml is not None and med_ml >= 4.0) else "high"
            _ntfy_loc = city_lbl_final or city_str
            _ntfy_body = (
                f"---- FINAL REPORT ----\n"
                f"Magnitude  : {ml_str}\n"
                f"Location   : {_ntfy_loc}\n"
                f"Depth      : {depth_str}\n"
                f"Dist home  : {dist_str}\n"
                f"Origin     : {orig_str}\n"
                f"Accuracy   : {acc_str}\n"
                f"Stations   : {en}\n"
                f"RMS        : {'%.2f s'%erms if erms else 'n/a'}\n"
                f"Az gap     : {'%.0f deg'%eaz_gap if eaz_gap else 'n/a'}\n"
                f"\n"
                f"---- STATIONS ----\n"
                f"{sta_text}\n"
                f"\n"
                f"---- TIMELINE ({len(timeline_copy)} events) ----\n"
                + "\n".join(f"{ts}  {note}" for ts, note in timeline_copy[-10:])
            )
            _send_ntfy(
                f"FINAL REPORT | {ml_str} | {_ntfy_loc}",
                _ntfy_body,
                priority=_ntfy_pri)
            _log("NTFY", f"Final report ntfy sent: {ml_str}")
        except Exception as _ne:
            _log("NTFY ERR", f"Final report ntfy failed: {_ne}")

        # ── Discord final report ───────────────────────────────────────────────
        try:
            _fr_eq_id  = snap.get("eq_id", "EQ-????")
            with _eq_id_lock:
                _fr_usgs = _current_eq_usgs[0]
            _fr_usgs_s  = f" / USGS:{_fr_usgs}" if _fr_usgs else ""
            _fr_upd_n, _fr_eq_raw = _eq_next_update()
            _fr_dc_title = (f"[{_fr_eq_raw}{_fr_usgs_s} | Final Report]"
                            f" {ml_str} | {city_lbl_final or city_str}")
            _fr_dc_bold  = (f"{ml_str} | {city_lbl_final or city_str}"
                            f" | Origin: {orig_str}")
            _fr_map = _make_map_png(elat, elon, et0=et0, med_ml=med_ml) if elat else None
            _send_discord(_fr_dc_title, _ntfy_body, bold_header=_fr_dc_bold,
                          img_bytes=_fr_map)
            _log("DISCORD", f"Final report sent to Discord: {ml_str}")
        except Exception as _fde:
            _log("DISCORD ERR", f"Final report Discord failed: {_fde}")

        # ── Generate map PNG (Agg — safe from background thread) ──────────────
        img_bytes = None
        try:
            fig_r    = mfig.Figure(figsize=(8, 6.5), facecolor="#0a0a0a")
            canvas_r = FigureCanvasAgg(fig_r)
            ax_r     = fig_r.add_subplot(111, facecolor="#08111e")
            fig_r.subplots_adjust(left=0.06, right=0.97, top=0.91, bottom=0.05)

            # California outline
            _CA_pts = np.array([
                (-124.25,42.00),(-120.00,42.00),(-120.00,39.00),(-119.32,38.50),
                (-119.00,37.50),(-118.20,36.50),(-116.50,35.75),(-114.63,35.00),
                (-114.62,32.73),(-117.13,32.53),(-118.19,33.73),(-119.05,34.04),
                (-120.66,34.58),(-121.33,35.79),(-121.90,36.96),(-122.17,37.20),
                (-122.51,37.73),(-122.53,38.01),(-122.98,38.10),(-123.38,38.56),
                (-123.70,38.85),(-123.97,39.84),(-124.24,40.30),(-124.25,42.00),
            ])
            ax_r.add_patch(_MplPoly(_CA_pts, closed=True,
                                    facecolor="#0d1f2d", edgecolor="#1e3a50",
                                    lw=1.5, zorder=1))
            # Determine map extent around epicenter + all stations
            all_lons = [lon_ for _,_,_,_,_,lon_,*_ in STATIONS] + [HOME_LON]
            all_lats = [lat_ for _,_,_,_,lat_,*_ in STATIONS] + [HOME_LAT]
            if elat is not None: all_lats.append(elat); all_lons.append(elon)
            pad_lo = 0.8; pad_hi = 0.8
            ax_r.set_xlim(min(all_lons) - pad_lo, max(all_lons) + pad_hi)
            ax_r.set_ylim(min(all_lats) - pad_lo, max(all_lats) + pad_hi)
            ax_r.set_aspect(1.0 / math.cos(math.radians(38.0)))
            ax_r.tick_params(colors="#666666", labelsize=7)
            ax_r.grid(True, color="#0e1a28", lw=0.5, zorder=0)
            ax_r.set_title(
                f"Epicenter Map — {ml_str}  |  {loc_str}  |  depth {depth_str}",
                fontsize=9, color="#cccccc", pad=5)

            # Station markers + S-P distance rings
            sta_cols2 = ["#1ab8e8","#f0a500","#2ecc71","#e056a0","#a78bfa","#fb923c"]
            for idx_, (net_, sta_, loc_, cha_, lat_, lon_, *_) in enumerate(STATIONS):
                col_ = sta_cols2[idx_ % len(sta_cols2)]
                ax_r.plot(lon_, lat_, "o", color=col_, ms=7, zorder=10,
                          markeredgecolor="#ffffff", markeredgewidth=0.5)
                ax_r.text(lon_, lat_, f"  {sta_}", color=col_,
                          fontsize=6.5, va="center", zorder=11)
                key_ = f"{net_}.{sta_}.{loc_}.{cha_}"
                st_  = states.get(key_)
                if st_ and st_.event_dist_km and st_.event_dist_km > 0:
                    theta_ = np.linspace(0, 2 * math.pi, 140)
                    dlat_  = st_.event_dist_km / 111.0
                    dlon_  = st_.event_dist_km / (111.0 * math.cos(math.radians(lat_)) + 1e-9)
                    ax_r.plot(lon_ + dlon_ * np.cos(theta_),
                              lat_ + dlat_ * np.sin(theta_),
                              color="#ffff44", lw=0.8, ls="--", alpha=0.55, zorder=8)

            # Home star
            ax_r.plot(HOME_LON, HOME_LAT, "*", color="#ffe033",
                      ms=14, zorder=14, markeredgecolor="#888800", markeredgewidth=0.5)
            ax_r.text(HOME_LON, HOME_LAT, "  San Ramon", color="#ffe033",
                      fontsize=6.5, va="center", zorder=15)

            # Epicenter cross + uncertainty ring + line to home
            if elat is not None and elon is not None:
                ax_r.plot(elon, elat, "+", color="#ff2222",
                          ms=22, mew=3.5, zorder=22)
                ax_r.text(elon, elat, f"  {ml_str}", color="#ff2222",
                          fontsize=9, fontweight="bold", va="bottom", zorder=23)
                ax_r.plot([elon, HOME_LON], [elat, HOME_LAT],
                          color="#ff2222", lw=0.7, ls="-", alpha=0.28, zorder=7)
                unc_km = max(15.0, (erms or 0.0) * VP)
                theta_u = np.linspace(0, 2 * math.pi, 140)
                dlat_u  = unc_km / 111.0
                dlon_u  = unc_km / (111.0 * math.cos(math.radians(elat)) + 1e-9)
                ax_r.plot(elon + dlon_u * np.cos(theta_u),
                          elat + dlat_u * np.sin(theta_u),
                          color="#ff3333", lw=1.5, ls="-", alpha=0.70, zorder=17)

            # Nearest city annotation
            if city_info:
                ax_r.text(0.02, 0.02,
                          f"Nearest: {city_info[0]} ({city_info[1]:.0f} km)",
                          transform=ax_r.transAxes, fontsize=7,
                          color="#aaaaaa", va="bottom")

            buf = io.BytesIO()
            canvas_r.print_figure(buf, format="png", dpi=130,
                                   facecolor="#0a0a0a")
            buf.seek(0)
            img_bytes = buf.read()
        except Exception as _me:
            _log("FINAL REPORT", f"Map generation error: {_me}")

        # ── Generate MyEarthquake Quake Report image ──────────────────────────
        myeq_bytes = None
        try:
            myeq_bytes = _generate_myearthquake_report(snap, timeline_copy)
            if myeq_bytes:
                _log("MYEQ REPORT", f"Generated MyEarthquake Quake Report "
                                    f"({len(myeq_bytes)//1024} KB)")
        except Exception as _mre:
            _log("MYEQ REPORT", f"Error: {_mre}")

        # ── Send email (only when email is configured) ────────────────────────
        if EMAIL_ENABLED and EMAIL_PASSWORD != "your-app-password-here":
            try:
                msg_id = (f"<seismic-final-{int(time.time())}"
                          f"-{_uuid.uuid4().hex[:8]}@quakemon>")
                with _email_lock:
                    reply_to = _event_thread_id[0]

                icon      = _subj_icon(med_ml)
                city_lbl  = _city_label(elat, elon) if elat is not None else None
                city_disp = city_lbl or (city_info[0] if city_info else "unknown")
                subj = f"{icon}FINAL REPORT — {ml_str}  {city_disp}"

                # Outer wrapper: mixed (allows attachments alongside related content)
                msg_outer = email.mime.multipart.MIMEMultipart("mixed")
                msg_outer["Subject"]    = subj
                msg_outer["From"]       = EMAIL_FROM
                msg_outer["To"]         = EMAIL_TO
                msg_outer["Message-ID"] = msg_id
                if reply_to:
                    msg_outer["In-Reply-To"] = reply_to
                    msg_outer["References"]  = reply_to

                # Inner related block: HTML body + inline image
                msg_related = email.mime.multipart.MIMEMultipart("related")

                # HTML body — epicenter map + MyEarthquake Quake Report inline
                html_body = (
                    f"<html><body style='background:#0d0d0d;color:#cccccc;"
                    f"font-family:monospace;font-size:13px'>"
                    f"<pre style='white-space:pre-wrap'>{body}</pre>"
                    + (f"<br><b style='color:#aaaaff'>Epicenter Map</b><br>"
                       f"<img src='cid:epicenter_map_cid' "
                       f"style='max-width:100%;border:1px solid #333'>"
                       if img_bytes else "")
                    + (f"<br><br><b style='color:#aaaaff'>MyEarthquake Quake Report</b><br>"
                       f"<img src='cid:myeq_report_cid' "
                       f"style='max-width:390px;border:1px solid #333'>"
                       if myeq_bytes else "")
                    + f"</body></html>"
                )
                msg_related.attach(email.mime.text.MIMEText(html_body, "html"))

                if img_bytes:
                    img_part = email.mime.image.MIMEImage(img_bytes)
                    img_part.add_header("Content-ID", "<epicenter_map_cid>")
                    img_part.add_header("Content-Disposition", "inline",
                                        filename="epicenter_map.png")
                    msg_related.attach(img_part)

                if myeq_bytes:
                    myeq_part = email.mime.image.MIMEImage(myeq_bytes)
                    myeq_part.add_header("Content-ID", "<myeq_report_cid>")
                    myeq_part.add_header("Content-Disposition", "inline",
                                         filename="myearthquake_quake_report.png")
                    msg_related.attach(myeq_part)

                msg_outer.attach(msg_related)

                # Also include plain-text fallback as an alternative
                msg_outer.attach(email.mime.text.MIMEText(body, "plain"))

                with smtplib.SMTP(EMAIL_SMTP_HOST, EMAIL_SMTP_PORT, timeout=25) as srv:
                    srv.ehlo(); srv.starttls()
                    srv.login(EMAIL_FROM, EMAIL_PASSWORD)
                    srv.sendmail(EMAIL_FROM, EMAIL_TO, msg_outer.as_string())
                _log("EMAIL",
                     f"Final report sent (map + MyEarthquake report): {ml_str}  {city_disp}")
            except Exception as exc:
                _log("EMAIL ERR", f"Final report send failed: {exc}")

    threading.Thread(target=_worker, daemon=True).start()


# ═══════════════════════════════════════════════════════════════════════════════
# COORDINATE HELPERS
# ═══════════════════════════════════════════════════════════════════════════════
def haversine_km(lat1, lon1, lat2, lon2):
    R    = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a    = (math.sin(dlat/2)**2
            + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
            * math.sin(dlon/2)**2)
    return 2 * R * math.asin(math.sqrt(min(a, 1.0)))

def _azimuth_deg(lat1, lon1, lat2, lon2):
    """Forward bearing from (lat1,lon1) to (lat2,lon2), degrees 0-360."""
    dlon = math.radians(lon2 - lon1)
    la1  = math.radians(lat1); la2 = math.radians(lat2)
    x    = math.sin(dlon) * math.cos(la2)
    y    = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(dlon)
    return (math.degrees(math.atan2(x, y)) + 360.0) % 360.0

def _az_gap_deg(epi_lat, epi_lon, sta_latlon_list):
    """
    Largest azimuthal gap (degrees) between any two consecutive stations
    as seen from the epicenter.  360 = all stations on one side (worst case).
    """
    if len(sta_latlon_list) < 2:
        return 360.0
    azs  = sorted(_azimuth_deg(epi_lat, epi_lon, sla, slo)
                  for sla, slo in sta_latlon_list)
    gaps = [azs[i+1] - azs[i] for i in range(len(azs) - 1)]
    gaps.append(360.0 - azs[-1] + azs[0])   # wrap-around gap
    return max(gaps)

def _epi_accuracy(n_sta, rms_sec, az_gap):
    """
    Return (label, hex_color, pos_err_km_or_None) describing epicenter quality.

    Scoring (0-6):
      n_sta ≥ 4             → +2
      n_sta == 3            → +1
      rms_sec < 0.30 s      → +2  | < 0.70 s → +1
      az_gap  < 150°        → +2  | < 220°   → +1

    Score 5-6 → HIGH  |  3-4 → MODERATE  |  1-2 → LOW  |  proxy → PROXY
    """
    if not n_sta:
        return "UNKNOWN", "#555555", None
    if n_sta == 1:
        return "PROXY", "#ff6644", None   # station used as geometric proxy

    pos_err = round(rms_sec * VP, 1) if rms_sec is not None else None

    score = 0
    if n_sta >= 4: score += 2
    elif n_sta >= 3: score += 1
    if rms_sec is not None:
        if   rms_sec < 0.30: score += 2
        elif rms_sec < 0.70: score += 1
    if az_gap is not None:
        if   az_gap < 150: score += 2
        elif az_gap < 220: score += 1

    if   score >= 5: return "HIGH",     "#44ff88", pos_err
    elif score >= 3: return "MODERATE", "#ffcc44", pos_err
    else:            return "LOW",      "#ff7744", pos_err

def _proj(lat, lon):
    """lat/lon → map plot coordinates (Web Mercator or pass-through)."""
    if HAS_TILES:
        x, y = _to_merc.transform(lon, lat)
        return float(x), float(y)
    return lon, lat

# ── Cities database (Bay Area · California · Pacific NW · Nevada) ─────────────
# Used for "nearest city" field in all alert emails.
_CITIES = [
    # Bay Area core
    ("San Francisco, CA",  37.7749, -122.4194),
    ("Oakland, CA",        37.8044, -122.2712),
    ("Berkeley, CA",       37.8716, -122.2727),
    ("San Jose, CA",       37.3382, -121.8863),
    ("San Ramon, CA",      37.7799, -121.9780),
    ("Fremont, CA",        37.5485, -121.9886),
    ("Hayward, CA",        37.6688, -122.0808),
    ("Concord, CA",        37.9780, -122.0311),
    ("Richmond, CA",       37.9358, -122.3478),
    ("Walnut Creek, CA",   37.9101, -122.0652),
    ("Livermore, CA",      37.6819, -121.7681),
    ("Pleasanton, CA",     37.6624, -121.8747),
    ("Antioch, CA",        38.0049, -121.8058),
    ("Pittsburg, CA",      38.0280, -121.8847),
    ("Daly City, CA",      37.6879, -122.4702),
    ("San Mateo, CA",      37.5630, -122.3255),
    ("Redwood City, CA",   37.4852, -122.2364),
    ("Palo Alto, CA",      37.4419, -122.1430),
    ("Mountain View, CA",  37.3861, -122.0839),
    ("Sunnyvale, CA",      37.3688, -122.0363),
    ("Santa Clara, CA",    37.3541, -121.9552),
    ("Milpitas, CA",       37.4323, -121.8996),
    # North Bay
    ("Santa Rosa, CA",     38.4404, -122.7141),
    ("Napa, CA",           38.2975, -122.2869),
    ("Petaluma, CA",       38.2324, -122.6366),
    ("Vallejo, CA",        38.1041, -122.2566),
    ("Fairfield, CA",      38.2494, -122.4000),
    ("Novato, CA",         38.1074, -122.5697),
    ("San Rafael, CA",     37.9735, -122.5311),
    ("Vacaville, CA",      38.3566, -121.9877),
    # South Bay / Peninsula
    ("Morgan Hill, CA",    37.1305, -121.6544),
    ("Gilroy, CA",         37.0058, -121.5683),
    ("Watsonville, CA",    36.9102, -121.7569),
    ("Santa Cruz, CA",     36.9741, -122.0308),
    ("Scotts Valley, CA",  37.0505, -122.0150),
    ("Salinas, CA",        36.6777, -121.6555),
    ("Monterey, CA",       36.6002, -121.8947),
    ("Seaside, CA",        36.6113, -121.8508),
    # Central Valley
    ("Stockton, CA",       37.9577, -121.2908),
    ("Modesto, CA",        37.6391, -120.9969),
    ("Turlock, CA",        37.4946, -120.8466),
    ("Sacramento, CA",     38.5816, -121.4944),
    ("Davis, CA",          38.5449, -121.7405),
    ("Lodi, CA",           38.1302, -121.2724),
    ("Manteca, CA",        37.7977, -121.2161),
    ("Tracy, CA",          37.7396, -121.4252),
    # Sierra Foothills
    ("Sonora, CA",         37.9835, -120.3824),
    ("Angels Camp, CA",    38.0682, -120.5399),
    ("Jackson, CA",        38.3485, -120.7730),
    ("Placerville, CA",    38.7296, -120.7985),
    ("Auburn, CA",         38.8966, -121.0769),
    # Southern CA
    ("Fresno, CA",         36.7378, -119.7871),
    ("Visalia, CA",        36.3302, -119.2921),
    ("Bakersfield, CA",    35.3733, -119.0187),
    ("Los Angeles, CA",    34.0522, -118.2437),
    ("Burbank, CA",        34.1808, -118.3090),
    # Pacific NW
    ("Portland, OR",       45.5051, -122.6750),
    ("Seattle, WA",        47.6062, -122.3321),
    ("Eugene, OR",         44.0521, -123.0868),
    # Nevada / Interior
    ("Reno, NV",           39.5296, -119.8138),
    ("Carson City, NV",    39.1638, -119.7674),
    ("Las Vegas, NV",      36.1699, -115.1398),
    # Coast
    ("Fort Bragg, CA",     39.4457, -123.8053),
    ("Eureka, CA",         40.8021, -124.1637),
    ("Crescent City, CA",  41.7558, -124.2026),
]

def _nearest_city(lat, lon, max_dist_km=600):
    """
    Return (city_name, dist_km) for the closest city in _CITIES within
    max_dist_km, or None if no city is within range.
    """
    if lat is None or lon is None:
        return None
    best_name, best_dist = None, float("inf")
    for name, clat, clon in _CITIES:
        d = haversine_km(lat, lon, clat, clon)
        if d < best_dist:
            best_dist, best_name = d, name
    if best_dist <= max_dist_km:
        return best_name, round(best_dist, 1)
    return None

def _city_label(lat, lon):
    """
    Return a human-readable epicenter label like '40 km E of Davis, CA'.
    Used in email subjects and headlines instead of raw distance-from-home.
    Returns None if no city found within 600 km.
    """
    city_info = _nearest_city(lat, lon)
    if city_info is None:
        return None
    name, dist = city_info
    # Look up the city's own coordinates
    try:
        city_lat, city_lon = next(
            (clat, clon) for cname, clat, clon in _CITIES if cname == name)
    except StopIteration:
        return f"{dist:.0f} km from {name}"
    # Bearing FROM city TO epicenter → cardinal direction
    az = _azimuth_deg(city_lat, city_lon, lat, lon)
    dirs = ["N", "NE", "E", "SE", "S", "SW", "W", "NW", "N"]
    dir_str = dirs[int((az + 22.5) / 45) % 8]
    return f"{dist:.0f} km {dir_str} of {name}"

def _ring_xy(lat, lon, radius_km, n=200):
    """(xs, ys) in map coordinates for a geodesic circle."""
    theta = np.linspace(0, 2 * math.pi, n)
    dlat  = radius_km / 111.0
    dlon  = radius_km / (111.0 * math.cos(math.radians(lat)) + 1e-9)
    lats  = lat + dlat * np.sin(theta)
    lons  = lon + dlon * np.cos(theta)
    if HAS_TILES:
        xs, ys = _to_merc.transform(lons, lats)
    else:
        xs, ys = lons, lats
    return np.asarray(xs), np.asarray(ys)

# ═══════════════════════════════════════════════════════════════════════════════
# SEISMIC MATH
# ═══════════════════════════════════════════════════════════════════════════════
# Maximum distance for local ML validity (Hutton & Boore 1987 calibrated to ~500 km)
ML_MAX_DIST_KM  = 500
# Maximum S-P time accepted as valid for local events (~1260 km at SP_FACTOR ≈ 0.119)
SP_MAX_SEC      = 55   # max S-P interval (~460km max event distance); rejects teleseism
                       # surface-wave tails being mis-picked as S-waves at distant stations
# Waveform window used for WA simulation — only first N seconds after P onset
ML_WAVE_WIN_SEC = 60

# Teleseism detection:
#   If ≥TELESEISM_MIN_STA stations trigger within TELESEISM_P_SPREAD seconds
#   of each other (nearly simultaneous → far planar wavefront) it's teleseismic.
TELESEISM_MIN_STA   = 3      # minimum stations that must trigger together
TELESEISM_P_SPREAD  = 25.0   # seconds — max inter-station P spread for teleseism
_teleseism_flag     = [False] # mutable flag shared across threads

def estimate_ml(peak_counts, sensitivity, dist_km, arr=None, sr=None,
                is_broadband=True, sta_lta_ratio=None, p_thr_eff=None,
                pre_p_arr=None):
    """
    Local magnitude — Hutton & Boore 1987 (NorCal).

    Valid range: ~10–500 km.  Returns None outside that range so that stale
    or teleseismic detections never produce absurd magnitudes.

    Parameters
    ----------
    arr, sr        : raw post-P waveform counts + sample rate for WA simulation.
    is_broadband   : True = HHZ flat sensor; False = EHZ ~1 Hz geophone.
                     Controls HP corner (0.10 vs 0.80 Hz) in WA pipeline.
    sta_lta_ratio  : P-band STA/LTA at detection time.  Acts as a quality gate:
                     if the ratio is < 1.5× the effective trigger threshold
                     AND no pre-P noise window is available, the SNR is too low
                     to produce a reliable amplitude and None is returned.
    p_thr_eff      : Effective P-picker threshold at this station
                     (= P_THRESH × thresh_mult).  Used with sta_lta_ratio.
    pre_p_arr      : Raw counts immediately before P onset (noise window).
                     Passed into _wa_amplitude_nm for quadrature noise
                     subtraction, removing the systematic over-estimation
                     that occurs when noise contaminates the amplitude peak.
    """
    if sensitivity <= 0 or dist_km < 1:
        return None
    # Hutton & Boore only valid for local distances — refuse teleseismic
    if dist_km > ML_MAX_DIST_KM:
        return None

    # ── STA/LTA quality gate ───────────────────────────────────────────────
    # A barely-triggered detection (ratio ≈ threshold) has low SNR; the
    # measured peak is significantly contaminated by noise.  Skip ML unless:
    #   (a) we have a pre-P noise window for correction, OR
    #   (b) the STA/LTA is clearly above threshold (ratio ≥ 1.5 × p_thr_eff)
    _has_noise_win = (pre_p_arr is not None and sr is not None
                      and len(pre_p_arr) >= max(int(sr * 1.5), 10))
    if (not _has_noise_win
            and sta_lta_ratio is not None
            and p_thr_eff     is not None
            and sta_lta_ratio < p_thr_eff * 1.5):
        return None   # SNR too low and no noise correction available

    # ── WA simulation (preferred path) ────────────────────────────────────
    # Limit to first ML_WAVE_WIN_SEC seconds after P onset so that very long
    # post-P segments (from stale events or 60-min buffers) don't inflate
    # the amplitude with coda or aftershock energy.
    A_nm = None
    if arr is not None and sr is not None:
        max_samp = int(ML_WAVE_WIN_SEC * sr)
        seg_wa   = arr[:max_samp] if len(arr) > max_samp else arr
        if len(seg_wa) >= int(sr * 1.5):
            A_nm = _wa_amplitude_nm(seg_wa, sensitivity, sr,
                                    is_broadband=is_broadband,
                                    pre_p_counts=pre_p_arr)

    # ── Fallback: velocity-peak approximation ─────────────────────────────
    if A_nm is None or A_nm <= 0:
        if peak_counts <= 0:
            return None
        # Approximate |H_WA_vel(f_dom)| = V₀·ω / |s²+2h·ω₀·s+ω₀²| at s=iω
        # Use station-type dominant frequency:
        #   HHZ broadband captures lower frequencies → 3 Hz typical
        #   EHZ geophone rolls off below 1 Hz → 5 Hz typical for local EQs
        f_dom = 3.0 if is_broadband else 5.0
        w     = 2.0 * math.pi * f_dom
        w0    = 2.0 * math.pi / 0.8      # WA natural freq
        h_wa  = 0.8
        denom = math.sqrt((w0 ** 2 - w ** 2) ** 2 + (2 * h_wa * w0 * w) ** 2)
        vel_nm_s = (peak_counts / sensitivity) * 1e9
        A_nm     = vel_nm_s * 2800.0 * w / denom

    A_mm = A_nm / 1e6
    if A_mm <= 0:
        return None
    try:
        # Hutton & Boore (1987) for Northern California, rearranged:
        #   ML = log10(A_mm) + 1.110*log10(R/100) + 0.00189*(R-100) + 3.0
        # Expanding the log10(R/100) term gives the constant +0.591:
        #   = log10(A_mm) + 1.110*log10(R) + 0.00189*R
        #     + (3.0 − 1.110×2 − 0.00189×100) = … + 0.591
        ml = (math.log10(A_mm)
              + 1.110 * math.log10(dist_km)
              + 0.00189 * dist_km + 0.591)
        # Sanity bounds: local network should never see ML < -3 or > 9
        if ml < -3.0 or ml > 9.0:
            return None
        return round(ml, 1)
    except ValueError:
        return None

def pgv_at_dist(ml, dist_km):
    """PGV (cm/s) at a given distance from ML source (NorCal attenuation)."""
    if ml is None or dist_km is None or dist_km < 1:
        return None
    try:
        return 10 ** (0.80*ml - 1.40*math.log10(dist_km) - 0.003*dist_km - 0.80)
    except ValueError:
        return None

def _sta_lta_level(ratio, threshold):
    """Classify a STA/LTA ratio relative to its effective trigger threshold.

    Returns one of four strings:
      "LOW"     — ratio is 1×–3× threshold  (barely above threshold)
      "MEDIUM"  — ratio is 3×–10× threshold (~triple threshold)
      "HIGH"    — ratio is 10×–20× threshold (very strong arrival)
      "EXTREME" — ratio is ≥ 20× threshold   (exceptional / major event)
    """
    if threshold and threshold > 0:
        mult = ratio / threshold
    else:
        mult = 0.0
    if mult >= 20.0:
        return "EXTREME"
    elif mult >= 10.0:
        return "HIGH"
    elif mult >= 3.0:
        return "MEDIUM"
    else:
        return "LOW"


def pgv_to_mmi(pgv_cm_s):
    """
    PGV (cm/s) → MMI using USGS ShakeMap threshold breakpoints with
    log-linear interpolation between them.  More accurate than Wald 1999
    at low PGV values (fixes "I Not felt" for ground motion that IS felt).

    Breakpoints (USGS ShakeMap Manual, Worden & Wald 2016, Table 5.1):
      < 0.1 cm/s → I   Not felt
      0.1 – 1.1  → II-III  Weak
      1.1 – 3.4  → IV  Light
      3.4 – 8.1  → V   Moderate
      8.1 – 16   → VI  Strong
      16  – 31   → VII Very Strong
      31  – 60   → VIII Severe
      60  – 116  → IX  Violent
      > 116       → X+
    """
    if pgv_cm_s is None or pgv_cm_s <= 0:
        return None
    if pgv_cm_s < 0.1:
        return 1.0
    # (pgv_upper_bound, mmi_at_upper_bound)
    _breaks = [
        (1.1,   3.0),
        (3.4,   4.0),
        (8.1,   5.0),
        (16.0,  6.0),
        (31.0,  7.0),
        (60.0,  8.0),
        (116.0, 9.0),
    ]
    pgv_lo, mmi_lo = 0.1, 2.0
    for pgv_hi, mmi_hi in _breaks:
        if pgv_cm_s <= pgv_hi:
            frac = math.log10(pgv_cm_s / pgv_lo) / math.log10(pgv_hi / pgv_lo)
            return round(min(12.0, mmi_lo + frac * (mmi_hi - mmi_lo)), 1)
        pgv_lo, mmi_lo = pgv_hi, mmi_hi
    return 12.0

def mmi_label(mmi):
    if mmi is None:  return "—",               "#444444"
    if mmi < 2:      return "I  Not felt",      "#666666"
    if mmi < 3:      return "II  Weak",         "#aaaaaa"
    if mmi < 4:      return "III  Felt",        "#99ccff"
    if mmi < 5:      return "IV  Light",        "#66ddff"
    if mmi < 6:      return "V  Moderate",      "#ffdd44"
    if mmi < 7:      return "VI  Strong",       "#ffaa22"
    if mmi < 8:      return "VII  V.Strong",    "#ff7700"
    if mmi < 9:      return "VIII  Severe",     "#ff4400"
    return               "IX+  Violent",        "#ff0000"

def counts_to_pgv(peak_counts, sensitivity):
    """Peak ground velocity in cm/s from raw counts."""
    if sensitivity <= 0: return 0.0
    return abs(peak_counts) / sensitivity * 100.0

def counts_to_pga(arr, sensitivity, sr):
    """Peak ground acceleration in cm/s² (differentiated velocity)."""
    if len(arr) < 2 or sensitivity <= 0: return 0.0
    vel = arr / sensitivity        # m/s
    acc = np.diff(vel) * sr        # m/s²
    return float(np.abs(acc).max()) * 100.0

# ── Advanced detection helpers ────────────────────────────────────────────────

def _aic_picker(x):
    """
    Vectorized AIC-CF onset picker — O(n) via prefix sums.
    Returns the sample index of the most likely P-wave onset.
    The AIC minimum separates pre-signal noise from the arriving wave.
    """
    n = len(x)
    if n < 8:
        return n // 2
    x  = x.astype(np.float64)
    s1 = np.cumsum(x)
    s2 = np.cumsum(x * x)
    k  = np.arange(1, n - 1, dtype=np.float64)
    # Left-window variance
    v1 = s2[:-2] / k       - (s1[:-2] / k) ** 2
    # Right-window variance
    k2 = (n - 1) - k
    rs1 = s1[-1] - s1[:-2]
    rs2 = s2[-1] - s2[:-2]
    v2  = rs2 / k2         - (rs1 / k2) ** 2
    aic = k * np.log(np.maximum(v1, 1e-30)) + k2 * np.log(np.maximum(v2, 1e-30))
    return int(np.argmin(aic)) + 1

def _wa_amplitude_nm(arr_counts, sensitivity, sr, is_broadband=True,
                     pre_p_counts=None):
    """
    Physically correct Wood-Anderson simulation on raw Z-velocity counts.
    Returns half peak-to-peak WA displacement in nanometres — the
    amplitude A used in the Hutton & Boore ML formula.

    WA seismometer transfer function (continuous-time, velocity → displacement):
        H(s) = V₀ · s / (s² + 2h·ω₀·s + ω₀²)
        T₀ = 0.8 s  →  ω₀ = 2π/0.8 = 7.854 rad/s
        h  = 0.8  (damping)
        V₀ = 2800 (static magnification)

    This is converted to a digital IIR filter via the bilinear transform,
    which replaces the old approach of:
       integrate → bandpass → ×2800
    The old approach applied a flat gain rather than the peaked WA response
    (WA peaks at ~1.25 Hz); this caused systematic amplitude errors.

    Station-type-aware HP pre-filter (applied to velocity before WA):
      HHZ broadband (flat ≥ 0.01 Hz): HP at 0.10 Hz — removes slow drift.
      EHZ geophone  (~1 Hz corner)  : HP at 0.80 Hz — must stay above the
        sensor corner; integrating counts below the geophone corner frequency
        amplifies noise that the scalar sensitivity doesn't account for.

    pre_p_counts : raw counts from the pre-P noise window.  When supplied,
        noise amplitude is estimated via the same pipeline and subtracted in
        quadrature:  A_corrected = sqrt(A_signal² − A_noise²).
        Eliminates the systematic over-estimation that occurs when the STA/LTA
        barely exceeds threshold (low-SNR detections).
    """
    min_samp = max(int(sr * 2), 10)
    if len(arr_counts) < min_samp or sensitivity <= 0 or sr <= 0:
        return None
    nyq = sr / 2.0

    # ── Station-type HP pre-filter ─────────────────────────────────────────
    # EHZ: 0.80 Hz keeps us above the geophone corner, preventing integration
    #      of sub-corner noise that the scalar sensitivity doesn't see.
    # HHZ: 0.10 Hz is sufficient for broadband sensors flat to ~0.01 Hz.
    # PB borehole EHZ sensors are flat well below 1 Hz — 0.50 Hz is safe
    # and avoids cutting into the S-wave coda that carries most ML energy.
    # (Old 0.80 Hz cut caused systematic underestimation at borehole stations.)
    hp_corner = 0.10 if is_broadband else 0.50
    sos_hp    = None
    if nyq > hp_corner * 1.5:
        sos_hp = butter(4, hp_corner / nyq, btype='high', output='sos')

    # ── WA IIR filter (continuous-time → digital via bilinear transform) ──
    # Numerator / denominator of H(s) in descending powers of s:
    #   b(s) = [V₀, 0]       →  V₀·s
    #   a(s) = [1, 2h·ω₀, ω₀²]
    T0, h_wa, V0 = 0.8, 0.8, 2800.0
    w0    = 2.0 * math.pi / T0            # 7.854 rad/s
    b_ct  = np.array([V0, 0.0])
    a_ct  = np.array([1.0, 2.0 * h_wa * w0, w0 ** 2])
    b_wa, a_wa = bilinear(b_ct, a_ct, sr)  # digital coefficients

    def _pipeline(cts):
        """HP → WA IIR → nm."""
        v = cts.astype(np.float64) / sensitivity
        if sos_hp is not None:
            v = sosfilt(sos_hp, v)
        return lfilter(b_wa, a_wa, v) * 1e9   # WA displacement, nm

    wa_sig = _pipeline(arr_counts)
    A_sig  = 0.5 * (wa_sig.max() - wa_sig.min())

    # ── Noise-floor correction (pre-P window) ──────────────────────────────
    if pre_p_counts is not None and len(pre_p_counts) >= max(int(sr * 1.5), 10):
        wa_noise = _pipeline(pre_p_counts)
        A_noise  = 0.5 * (wa_noise.max() - wa_noise.min())
        # Quadrature subtraction: A_true ≈ sqrt(A_meas² − A_noise²)
        A_sig = math.sqrt(max(0.0, A_sig ** 2 - A_noise ** 2))

    return float(A_sig) if A_sig > 0 else None

# ═══════════════════════════════════════════════════════════════════════════════
# RING BUFFER — fixed-capacity numpy circular buffer (~4× less RAM than deque)
# ═══════════════════════════════════════════════════════════════════════════════
class _RingBuf:
    """Fixed-capacity circular numpy buffer.
    Appending is O(m) array copy (no Python-object overhead).
    Reading returns a numpy view when possible (zero-copy) or a single
    np.concatenate when the data wraps around the end of the backing store.
    """
    __slots__ = ("_buf", "_ptr", "_n", "capacity")

    def __init__(self, capacity, dtype=np.float64):
        self._buf      = np.empty(capacity, dtype=dtype)
        self._ptr      = 0        # next-write index (mod capacity)
        self._n        = 0        # valid sample count (≤ capacity)
        self.capacity  = capacity

    def extend(self, arr):
        """Append a numpy array *arr* into the ring (oldest samples discarded)."""
        m = len(arr)
        if m == 0:
            return
        if m >= self.capacity:
            self._buf[:] = arr[-self.capacity:]
            self._ptr    = 0
            self._n      = self.capacity
            return
        end = self._ptr + m
        if end <= self.capacity:
            self._buf[self._ptr:end] = arr
        else:
            split = self.capacity - self._ptr
            self._buf[self._ptr:]  = arr[:split]
            self._buf[:m - split]  = arr[split:]
        self._ptr = end % self.capacity
        self._n   = min(self._n + m, self.capacity)

    def last(self):
        """Return the most-recently appended value, or None if empty."""
        if self._n == 0:
            return None
        return float(self._buf[(self._ptr - 1) % self.capacity])

    def to_array(self):
        """Return valid data in chronological order.
        Returns a *view* (zero-copy) when not yet wrapped or at ptr==0,
        otherwise a single np.concatenate copy."""
        if self._n == 0:
            return np.empty(0, dtype=self._buf.dtype)
        if self._n < self.capacity:
            return self._buf[:self._n]        # view
        if self._ptr == 0:
            return self._buf                  # view — full, not yet wrapped
        return np.concatenate((self._buf[self._ptr:], self._buf[:self._ptr]))

    def __len__(self):
        return self._n

    def __bool__(self):
        return self._n > 0


# ═══════════════════════════════════════════════════════════════════════════════
# PER-STATION STATE
# ═══════════════════════════════════════════════════════════════════════════════
class StationState:
    def __init__(self, net, sta, loc, cha, lat, lon, dist, desc, sens, tmult):
        self.net, self.sta, self.loc, self.cha = net, sta, loc, cha
        self.lat, self.lon   = lat, lon
        self.dist_km         = dist
        self.description     = desc
        self.sensitivity     = sens
        self.thresh_mult     = tmult
        self.alarm_enabled   = True    # toggled via Settings window
        self.label = f"{net}.{sta}.{loc}.{cha}" if loc else f"{net}.{sta}.{cha}"
        self.id    = f"{net}.{sta}.{loc}.{cha}"

        # Sensor type flag — affects WA simulation HP corner and bandpass:
        #   HHZ broadband (BK/STS-2/Trillium): flat to ~0.008 Hz → is_broadband=True
        #   EHZ borehole  (PB/Trident/L22):    corner ~1 Hz       → is_broadband=False
        self.is_broadband = cha.startswith("H")

        self.lock            = threading.Lock()
        self.sample_rate     = None
        # Waveform ring buffers — allocated lazily on first packet (SR unknown here)
        self.samples         = None   # _RingBuf once allocated
        self.times           = None   # _RingBuf once allocated
        self._buf_cap        = 0      # capacity of current ring buffers
        # CFT arrays — stored directly as numpy arrays (replaced each packet)
        self.cft             = None   # numpy float64 array or None
        self.p_cft           = None   # numpy float64 array or None

        self.alert           = False
        self.last_alert_t    = 0.0
        self.connected       = False
        self.last_data_time  = None    # wall-clock time of last on_data packet
        self.last_ratio      = 0.0
        self.last_p_ratio    = 0.0
        self.p_time          = None
        self.p_predicted     = None    # predicted P from confirmed epicenter
        self.s_time          = None
        self.s_predicted     = None
        self.p_cleared       = False
        self.event_dist_km   = None
        self.event_peak      = 0.0
        self.ml_est          = None
        self.pgv_cm_s        = 0.0     # rolling 5-s PGV
        self.pga_cm_s2       = 0.0     # rolling 5-s PGA
        self.last_ps_ratio   = None
        self.ps_hist         = collections.deque(maxlen=600)
        self.s_cft           = None   # numpy float64 array or None — S-band STA/LTA CF
        self.last_s_ratio    = 0.0
        self.ml_unc          = None   # std dev of ML across stations (set in _animate)
        self._sos_p          = None
        self._sos_s          = None
        self._sos_sr         = None

        self.gap_count       = 0      # number of data gaps detected
        self.clip_count      = 0      # number of clipped samples detected
        self.data_quality    = 100.0  # 0-100 health score
        self.bytes_received  = 0      # total samples received
        self.amplitude_suspect = False  # True if amplitude inconsistent with ML/distance
        self.freq_ratio      = None   # hi-band / lo-band power ratio (teleseism check)

        # Auto-calibration: collect background STA/LTA during first 90 s of
        # connection (only while no P-wave active) then auto-adjust thresh_mult
        self._calib_ratios   = collections.deque(maxlen=600)
        self._calib_done     = False   # True once calibration fires
        self._calib_start    = None    # wall-clock time of first packet

        # Sleep/gap recovery: True while processing backlog data (stale packets).
        # Set on the first stale packet; cleared when live data resumes.
        self._was_stale      = False

    def _ensure_bufs(self, sr):
        """Allocate (or resize) the waveform ring buffers for a given sample rate.
        Called once per on_data packet — no-op when capacity is unchanged.
        Reads _display_sec_active[0] so Hide Waves can shrink the buffer at runtime."""
        cap = int(_display_sec_active[0] * sr)
        if cap == self._buf_cap:
            return  # already the right size — nothing to do
        # First allocation or capacity changed (SR change or hide/show toggle):
        # Preserve existing data up to the new capacity.
        old_s = self.samples.to_array() if self.samples is not None else np.empty(0)
        old_t = self.times.to_array()   if self.times   is not None else np.empty(0)
        # Samples: float32 (seismometer counts are integers — float32 is plenty)
        # Times:   float64 (Unix epoch needs nanosecond-range precision)
        self.samples   = _RingBuf(cap, dtype=np.float32)
        self.times     = _RingBuf(cap, dtype=np.float64)
        self._buf_cap  = cap
        if len(old_s):
            self.samples.extend(old_s[-cap:].astype(np.float32))
            self.times.extend(old_t[-cap:])

states = {
    f"{net}.{sta}.{loc}.{cha}": StationState(
        net, sta, loc, cha, lat, lon, d, desc, sens, tmult)
    for net, sta, loc, cha, lat, lon, d, desc, sens, tmult in STATIONS
}

# ═══════════════════════════════════════════════════════════════════════════════
# EPICENTER ESTIMATOR  (3-D: lat, lon, depth_km, t_origin)
# ═══════════════════════════════════════════════════════════════════════════════
class EpicenterEstimator:
    def __init__(self):
        self.lock       = threading.Lock()
        self.arrivals   = {}
        self.lat        = None
        self.lon        = None
        self.depth_km   = None
        self.t_origin   = None
        self.dist_home  = None
        self.rms_sec    = None
        self.az_gap     = None   # largest azimuthal gap between stations (degrees)
        self.n_sta      = 0

    def add(self, key, slat, slon, tp):
        with self.lock:
            self.arrivals[key] = (slat, slon, tp)
            self._solve()

    @staticmethod
    def _tt(la, lo, dz, sl, so):
        h = haversine_km(la, lo, sl, so)
        return math.sqrt(h*h + max(dz, 0.1)**2) / VP

    def _solve(self):
        arr = list(self.arrivals.values())   # [(slat, slon, tp), ...]
        n   = len(arr)
        if n < 1:
            return

        # ── Collect S-wave arrivals from station states ────────────────────
        # S arrivals tightly constrain depth and provide extra equations.
        s_arr = []   # [(slat, slon, ts), ...]
        for key, (slat, slon, _) in self.arrivals.items():
            st = states.get(key)
            if st and st.s_time is not None:
                s_arr.append((slat, slon, st.s_time))

        # ── Initial guess: use the earliest-arriving station as location seed ──
        # The first P arrival comes from the station nearest the epicenter,
        # so it's a much better seed than the plain centroid of all stations.
        tp_min    = min(a[2] for a in arr)
        arr_by_t  = sorted(arr, key=lambda x: x[2])
        la0, lo0  = arr_by_t[0][0], arr_by_t[0][1]
        # Estimate origin time: earliest P minus short travel-time guess (3 s).
        # For near-source events (< 30 km) this is within ±3 s of the truth.
        t00       = tp_min - 3.0

        if n == 1:
            # Single station: use station as proxy epicenter, fixed 10 km depth
            slat, slon, tp0 = arr[0]
            key0 = list(self.arrivals.keys())[0]
            st0  = states.get(key0)
            dist = (st0.event_dist_km if st0 and st0.event_dist_km else
                    st0.dist_km       if st0 else 50.0)
            self.lat       = slat
            self.lon       = slon
            self.depth_km  = 10.0
            self.t_origin  = tp0 - dist / VP
            self.dist_home = haversine_km(slat, slon, HOME_LAT, HOME_LON)
            self.rms_sec   = None
            self.az_gap    = 360.0   # proxy: no geometric constraint
            self.n_sta     = 1
            return

        # ── Objective function: P + S residuals (weighted) ─────────────────
        def _ts(la, lo, dz, sl, so):
            h = haversine_km(la, lo, sl, so)
            return math.sqrt(h*h + max(dz, 0.1)**2) / VS

        def obj(p):
            la, lo, dz, t0 = p
            # P residuals
            res  = [(t0 + self._tt(la, lo, dz, sl, so) - tp)**2
                    for sl, so, tp in arr]
            # S residuals (same weight — S arrivals constrain depth well)
            res += [(t0 + _ts(la, lo, dz, sl, so) - ts)**2
                    for sl, so, ts in s_arr]
            return sum(res)

        n_obs = n + len(s_arr)
        # Multi-start Nelder-Mead: try several t0 offsets and take the best.
        # This prevents the solver from getting stuck in a local minimum when
        # the initial origin-time guess is wrong (e.g., deep events, slow P).
        _best_r = None
        for _t0_off in (3.0, 8.0, 15.0, 1.0, 25.0):
            _r = minimize(obj, [la0, lo0, 10.0, tp_min - _t0_off],
                          method="Nelder-Mead",
                          options={"xatol": 1e-4, "fatol": 1e-8, "maxiter": 30000})
            if _best_r is None or _r.fun < _best_r.fun:
                _best_r = _r
        r = _best_r
        # Accept any solution; use RMS as quality indicator rather than hard cutoff.
        # A large RMS just means POOR accuracy — still better than a proxy.
        if r.fun < 500.0 * n_obs:
            la, lo, dz, t0 = r.x
            self.lat       = la
            self.lon       = lo
            self.depth_km  = max(1.0, round(dz, 1))  # minimum 1 km — 0 km is unphysical
            self.t_origin  = t0
            self.dist_home = haversine_km(la, lo, HOME_LAT, HOME_LON)
            self.rms_sec   = math.sqrt(r.fun / n_obs)
            self.az_gap    = _az_gap_deg(la, lo, [(sl, so) for sl, so, _ in arr])
            self.n_sta     = n
            # ── Update every triggered station's event_dist_km to the actual
            #    hypocentral distance and recompute ML with the correct distance.
            #    This eliminates the systematic error from using S-P proxy distance
            #    or the fixed station-to-home distance.
            # ── P-wave arrival-order consistency ─────────────────────────────
            # Stations closer to the solved epicentre must have gotten P first.
            # Count inversions (farther station arriving before closer station).
            # If >50% of pairs are inverted the solution is suspect.
            _order_pairs = 0; _order_ok = 0
            _arr_sorted = sorted(
                [(haversine_km(la, lo, sl, so), tp)
                 for sl, so, tp in arr],
                key=lambda x: x[0])   # sort by distance
            for _pi in range(len(_arr_sorted)):
                for _pj in range(_pi+1, len(_arr_sorted)):
                    _di, _ti = _arr_sorted[_pi]
                    _dj, _tj = _arr_sorted[_pj]
                    _order_pairs += 1
                    if _ti <= _tj:   # closer station arrived first
                        _order_ok += 1
            if _order_pairs > 0:
                _order_frac = _order_ok / _order_pairs
                if _order_frac < 0.5:
                    _log("EPICENTER",
                         f"WARNING: P-arrival order inconsistent with solution "
                         f"({_order_ok}/{_order_pairs} pairs OK = {_order_frac:.0%})"
                         f" — location may be unreliable")
            dz_clip = max(dz, 0.1)
            for _key, (_slat, _slon, _tp) in self.arrivals.items():
                _st = states.get(_key)
                if _st is None or _st.p_time is None:
                    continue
                _h  = haversine_km(la, lo, _slat, _slon)
                _hypo = math.sqrt(_h * _h + dz_clip * dz_clip)
                _st.event_dist_km = _hypo
                # Recompute ML with the corrected distance if we have the waveform
                if (_st.event_peak > 0 and _st.sensitivity > 0):
                    _ml_new = estimate_ml(
                        _st.event_peak, _st.sensitivity, _hypo,
                        arr=None, sr=None,
                        is_broadband=_st.is_broadband,
                        sta_lta_ratio=_st.last_p_ratio,
                        p_thr_eff=P_THRESH * _st.thresh_mult)
                    if _ml_new is not None:
                        _st.ml_est = _ml_new
                    # Amplitude-distance consistency: expected PGV at this distance
                    # vs. observed amplitude.  If station shows < 8% of expected,
                    # it is likely a false pick or badly calibrated sensitivity.
                    if _ml_new is not None and _hypo > 0:
                        _pgv_exp = pgv_at_dist(_ml_new, _hypo)   # cm/s
                        _pgv_obs = counts_to_pgv(_st.event_peak, _st.sensitivity)
                        if (_pgv_exp and _pgv_obs is not None and _pgv_exp > 0
                                and _pgv_obs < 0.08 * _pgv_exp):
                            _st.amplitude_suspect = True
                            _log("AMP-CHECK",
                                 f"{_key}  pgv_obs={_pgv_obs:.4f} < 8% of"
                                 f" pgv_exp={_pgv_exp:.4f} at {_hypo:.0f}km"
                                 f" — flagged as amplitude-suspect")
                        else:
                            _st.amplitude_suspect = False

            # ── Post-solve S-pick validation ──────────────────────────────────
            # Reject any S-pick whose implied S-P distance is inconsistent
            # with the solved station-to-hypocentre distance.  False S-picks
            # (e.g. coda noise picked as S-wave on a 130 km station after a
            # 10 km event) are the #1 cause of blown epicentre solutions.
            _bad_s_keys = []
            for _key, (_slat, _slon, _tp) in self.arrivals.items():
                _st = states.get(_key)
                if _st is None or _st.s_time is None:
                    continue
                _sp_obs  = _st.s_time - _st.p_time          # observed S-P (s)
                _sp_dist = _sp_obs / SP_FACTOR               # S-P implied distance (km)
                _epi_dist = haversine_km(la, lo, _slat, _slon)
                # Expected S-P at the solved epicentre-to-station distance
                _sp_exp  = _epi_dist * SP_FACTOR
                # Reject if implied distance differs by more than 50 % from solved distance
                if _sp_exp > 0 and abs(_sp_dist - _epi_dist) / max(_epi_dist, 1.0) > 0.50:
                    _bad_s_keys.append(_key)
                    _log("S-REJECT",
                         f"{_key}  S-P={_sp_obs:.1f}s → {_sp_dist:.0f}km"
                         f"  but epi-dist={_epi_dist:.0f}km  (>50% off) — dropping S-pick")
            if _bad_s_keys:
                for _bk in _bad_s_keys:
                    _st = states.get(_bk)
                    if _st:
                        _st.s_time       = None
                        _st.event_dist_km = None   # will be re-estimated next packet
                # Re-collect s_arr without bad picks and re-solve once
                s_arr = []
                for _key2, (_slat2, _slon2, _) in self.arrivals.items():
                    _st2 = states.get(_key2)
                    if _st2 and _st2.s_time is not None:
                        s_arr.append((_slat2, _slon2, _st2.s_time))
                # Quick re-solve with cleaned S-picks
                _best_r2 = None
                for _t0_off2 in (3.0, 8.0, 15.0):
                    _r2 = minimize(obj, [la, lo, max(dz, 1.0), t0],
                                   method="Nelder-Mead",
                                   options={"xatol": 1e-4, "fatol": 1e-8,
                                            "maxiter": 20000})
                    if _best_r2 is None or _r2.fun < _best_r2.fun:
                        _best_r2 = _r2
                if _best_r2 is not None and _best_r2.fun < 500.0 * max(n + len(s_arr), 1):
                    la, lo, dz, t0 = _best_r2.x
                    self.lat      = la
                    self.lon      = lo
                    self.depth_km = max(1.0, round(dz, 1))
                    self.t_origin = t0
                    self.dist_home = haversine_km(la, lo, HOME_LAT, HOME_LON)
                    self.rms_sec   = math.sqrt(_best_r2.fun / max(n + len(s_arr), 1))
                    self.az_gap    = _az_gap_deg(la, lo, [(sl, so) for sl, so, _ in arr])
                    _log("EPICENTER", f"Re-solved after S-pick cleanup: "
                                      f"{la:.4f}°N {lo:.4f}°W  RMS={self.rms_sec:.2f}s")

            s_note         = f" + {len(s_arr)}S" if s_arr else ""
            if n >= 3:
                _log("EPICENTER",
                     f"{la:.4f}°N  {lo:.4f}°W  depth={self.depth_km:.1f}km"
                     f"  dist_home={self.dist_home:.0f}km"
                     f"  RMS={self.rms_sec:.2f}s  az_gap={self.az_gap:.0f}°"
                     f"  ({n}P{s_note})")
                _city_epi = _nearest_city(la, lo)
                _city_epi_s = f"  near {_city_epi[0]}" if _city_epi else ""
                _timeline_add(
                    f"Epicenter updated: {la:.3f}°N {lo:.3f}°W"
                    f"  depth={self.depth_km:.1f}km"
                    f"  RMS={self.rms_sec:.2f}s"
                    f"  ({n}P{s_note}){_city_epi_s}")
        else:
            # Solver failed to converge — keep existing location (proxy or previous
            # solution) but update n_sta so the display shows the real arrival count.
            self.n_sta = n
            _log("EPICENTER",
                 f"solver did not converge with {n} arrivals (fun={r.fun:.3f}>"
                 f"{500.0*n_obs:.1f}) — keeping previous position")

    def remove_arrival(self, key):
        """Remove one station from the solution; re-solve with remaining stations."""
        with self.lock:
            if key not in self.arrivals:
                return
            del self.arrivals[key]
            if self.arrivals:
                self._solve()
            else:
                self.lat = self.lon = self.depth_km = self.t_origin = None
                self.dist_home = self.rms_sec = self.az_gap = None
                self.n_sta = 0

    def reset(self):
        with self.lock:
            self.arrivals.clear()
            self.lat = self.lon = self.depth_km = self.t_origin = None
            self.dist_home = self.rms_sec = self.az_gap = None
            self.n_sta = 0

epicenter = EpicenterEstimator()

# ═══════════════════════════════════════════════════════════════════════════════
# AUDIO / VOICE
# ═══════════════════════════════════════════════════════════════════════════════
def _play_sound():
    if not _sound_muted and _IS_MACOS and _SOUND:
        for _ in range(4):
            subprocess.Popen(["afplay", _SOUND]); time.sleep(0.18)

def _speak(txt):
    if not _sound_muted and _IS_MACOS:
        subprocess.Popen(["say", "-r", "200", txt])

# ── Manual alarm (terminal 'alarm' command) ───────────────────────────────────
_manual_alarm_active = [False]

def _trigger_manual_alarm():
    """Start a continuous alarm loop (runs until _stop_manual_alarm is called)."""
    if _manual_alarm_active[0]:
        print("  [ALARM] Already active — type 'alarm off' to stop")
        return
    _manual_alarm_active[0] = True
    def _loop():
        _speak("Manual alarm activated.")
        while _manual_alarm_active[0]:
            if _IS_MACOS and _SOUND:
                for _ in range(6):
                    if not _manual_alarm_active[0]:
                        break
                    subprocess.Popen(["afplay", _SOUND])
                    time.sleep(0.22)
            time.sleep(0.4)
    threading.Thread(target=_loop, daemon=True).start()
    print("  [ALARM] ACTIVE — type  alarm off  to stop")

def _stop_manual_alarm():
    """Stop the continuous alarm loop."""
    if not _manual_alarm_active[0]:
        print("  Alarm is not active")
        return
    _manual_alarm_active[0] = False
    print("  [ALARM] Stopped")

def _fire_alert(st, ratio):
    if not st.alarm_enabled:
        return
    _play_sound()
    _speak(f"Earthquake warning! Elevated seismic activity at station {st.sta}. "
           f"Ratio {ratio:.1f}.")

# ═══════════════════════════════════════════════════════════════════════════════
# SEEDLINK CLIENT
# ═══════════════════════════════════════════════════════════════════════════════
class QuakeClient(EasySeedLinkClient):
    def on_data(self, trace):
        key  = trace.id
        st   = states.get(key)
        if st is None:
            return
        sr    = trace.stats.sampling_rate
        data  = trace.data.astype(np.float64)
        t0    = trace.stats.starttime.timestamp
        times = t0 + np.arange(len(data)) / sr

        with st.lock:
            st.sample_rate    = sr
            st.connected      = True
            st.last_data_time = time.time()
            st._ensure_bufs(sr)          # allocate / resize ring buffers if needed
            mx = st._buf_cap             # convenient alias used later

            # ── Drop samples that overlap with the tail of the buffer ─────────
            # SeedLink regularly re-sends the last 1-2 seconds of the previous
            # miniSEED record.  Non-monotonic timestamps produce the zig-zag /
            # crossing-line artifact on the waveform plots.
            if st.times:
                last_t = st.times.last()
                valid  = times > (last_t + 1e-4)   # strict forward-only; 0.1 ms tolerance
                if valid.any():
                    data  = data[valid]
                    times = times[valid]
                else:
                    # Entire packet already covered — skip it
                    return

            # Ring buffers handle capacity trimming automatically — no popleft loop
            st.samples.extend(data.astype(np.float32))  # float32 — halves sample RAM
            st.times.extend(times)                       # float64 — timestamps need precision
            arr   = st.samples.to_array()   # view (zero-copy) most of the time
            t_arr = st.times.to_array()     # view (zero-copy) most of the time
            now   = time.time()

            # ── Clip detection (ADC saturation) ──────────────────────────────
            _max_val = np.abs(data).max()
            if _max_val > 0.95 * 2**23:   # within 5% of 24-bit ADC limit
                st.clip_count += 1
            st.bytes_received += len(data)
            # Update health score (degrades with gaps/clips, recovers over time)
            st.data_quality = max(0.0, min(100.0,
                100.0 - st.gap_count * 2.0 - st.clip_count * 5.0))

            # ── Broadcast new samples to WebSocket clients (website live feed) ─
            if _ws_connected:
                _ws_broadcast_data({
                    "key":    key,
                    "t_end":  float(times[-1]),
                    "sr":     float(sr),
                    "stalta": round(float(st.last_ratio), 2),
                    "v":      data.astype(int).tolist(),
                })

            # ── Sleep / gap recovery guard ────────────────────────────────────
            # When the computer wakes from sleep (or reconnects after a network
            # drop), the SeedLink server sends a burst of buffered packets whose
            # seismic timestamps are many minutes in the past.  Without this
            # guard the STA/LTA spikes on the sudden data burst and produces
            # false P detections (the "B054 P at -1333s" scenario).
            #
            # Wall-clock age of the newest sample in this packet:
            _data_age  = now - t_arr[-1]
            _now_stale = _data_age > STALE_DATA_SEC

            if _now_stale:
                if not st._was_stale:
                    # First stale packet after a gap — clear any lingering
                    # detections exactly once and restart STA/LTA calibration.
                    st._was_stale = True
                    st.gap_count += 1
                    st.data_quality = max(0.0, min(100.0,
                        100.0 - st.gap_count * 2.0 - st.clip_count * 5.0))
                    _log("GAP",
                         f"{st.label} backlog data {_data_age:.0f}s old "
                         f"— clearing detection state (sleep/reconnect recovery)")
                    st.alert      = False
                    st.p_time     = None;  st.s_time      = None
                    st.p_predicted = None; st.s_predicted  = None
                    st.event_peak = 0.0;   st.ml_est       = None
                    st.p_cleared  = False; st.event_dist_km = None
                    st.amplitude_suspect = False
                    st._calib_done = False
                    st._calib_ratios.clear()
                    st._calib_start = None
                    epicenter.remove_arrival(key)
                # Ring buffer already updated — skip all detection this packet
                return

            # Data is live; mark recovery complete
            if st._was_stale:
                _log("GAP", f"{st.label} caught up to live data — detection resumed")
                st._was_stale = False

            # ── Apply pending PhaseNet pick refinement ────────────────────────
            # PhaseNet background thread stores results here; we apply them on
            # the next on_data() call so inference never blocks the SeedLink thread.
            if _phasenet_ready[0]:
                with _phasenet_pend_lk:
                    _pn = _phasenet_pending.pop(key, None)
                if _pn is not None:
                    _pn_p = _pn.get("p_time")
                    # Refine P: accept if within PHASENET_MAX_SHIFT of STA/LTA pick
                    if (_pn_p is not None
                            and st.p_time is not None
                            and abs(_pn_p - st.p_time) <= PHASENET_MAX_SHIFT):
                        _old_pt   = st.p_time
                        st.p_time = _pn_p
                        _log("PHASENET",
                             f"{st.label}  P refined"
                             f"  {_old_pt - now:+.3f}s → {st.p_time - now:+.3f}s"
                             f"  (prob={_pn['p_prob']:.2f})")
                        epicenter.add(key, st.lat, st.lon, st.p_time)
                    # Apply S pick if we don't have one yet
                    _pn_s = _pn.get("s_time")
                    if (_pn_s is not None
                            and st.s_time is None
                            and st.p_time is not None
                            and _pn_s > st.p_time + 1.0):
                        _sp_pn = _pn_s - st.p_time
                        if 1.0 < _sp_pn <= SP_MAX_SEC:
                            st.s_time        = _pn_s
                            st.event_dist_km = _sp_pn / SP_FACTOR
                            _log("PHASENET",
                                 f"{st.label}  S from ML"
                                 f"  S-P={_sp_pn:.1f}s"
                                 f"  dist={st.event_dist_km:.0f}km"
                                 f"  (prob={_pn['s_prob']:.2f})")

            thr_on  = TRIGGER_ON  * st.thresh_mult
            thr_off = TRIGGER_OFF * st.thresh_mult
            p_thr   = P_THRESH    * st.thresh_mult

            # ── Adaptive noise floor tracking ─────────────────────────────────
            # Track rolling RMS of the last 60s of data outside any event window.
            # Adjust thresh_mult if noise floor drifts significantly from baseline.
            if st.p_time is None and st.alert is False:
                _noise_rms = float(np.sqrt(np.mean(arr[-int(sr*10):]**2))) if len(arr) >= int(sr*10) else None
                if _noise_rms is not None and _noise_rms > 0:
                    if not hasattr(st, '_noise_baseline'):
                        st._noise_baseline = _noise_rms
                        st._noise_samples  = 1
                    else:
                        # Exponential moving average
                        alpha = 0.01
                        st._noise_baseline = (1-alpha)*st._noise_baseline + alpha*_noise_rms
                        st._noise_samples  = min(st._noise_samples + 1, 10000)
                        # If noise is 3× baseline, raise threshold temporarily
                        if st._noise_samples > 100:
                            noise_ratio = _noise_rms / max(st._noise_baseline, 1.0)
                            if noise_ratio > 3.0:
                                st.thresh_mult = min(st.thresh_mult * 1.05,
                                                    _DEFAULT_THRESH.get(key, 1.0) * 5.0)
                            elif noise_ratio < 1.5 and st.thresh_mult > _DEFAULT_THRESH.get(key, 1.0):
                                st.thresh_mult = max(st.thresh_mult * 0.99,
                                                    _DEFAULT_THRESH.get(key, 1.0))

            # ── Alert STA/LTA ──────────────────────────────────────────────
            ns, nl = int(STA_SEC*sr), int(LTA_SEC*sr)
            if len(arr) >= ns + nl:
                cft    = recursive_sta_lta(arr, ns, nl)
                st.cft = cft          # numpy array — replaced each packet, no deque needed
                ratio  = float(cft[-1]); st.last_ratio = ratio

                # ── Background auto-calibration (fires once per station) ──────
                if not st._calib_done:
                    if st._calib_start is None:
                        st._calib_start = now
                    # Collect only while quiet and within calibration window
                    if (st.p_time is None and not st.alert
                            and (now - st._calib_start) < 90.0):
                        st._calib_ratios.append(ratio)
                    elif len(st._calib_ratios) >= 30:
                        st._calib_done = True
                        p95 = float(np.percentile(list(st._calib_ratios), 95))
                        # Borehole EHZ stations are quieter → tighter target
                        # Surface HHZ broadband stations have more cultural noise
                        target = TRIGGER_OFF * (0.45 if not st.is_broadband else 0.60)
                        if p95 > 0.05:
                            correction = max(0.40, min(2.50, target / p95))
                            new_mult   = round(
                                max(0.20, min(4.0, st.thresh_mult * correction)), 2)
                            if abs(new_mult - st.thresh_mult) > 0.08:
                                old_mult       = st.thresh_mult
                                st.thresh_mult = new_mult
                                _log("AUTOCAL",
                                     f"{st.label}  thresh_mult "
                                     f"{old_mult:.2f}→{new_mult:.2f}"
                                     f"  P95_bg={p95:.2f}"
                                     f"  target={target:.2f}")

                if not st.alert and ratio >= thr_on:
                    st.alert = True
                    _log("ALERT ON",
                         f"{st.label}  STA/LTA={ratio:.2f}"
                         f"  [{_sta_lta_level(ratio, thr_on)}]"
                         f"  dist={st.dist_km}km  ({st.description})")
                    if now - st.last_alert_t > ALERT_COOLDOWN:
                        st.last_alert_t = now
                        threading.Thread(target=_fire_alert, args=(st, ratio),
                                         daemon=True).start()
                elif st.alert and ratio < thr_off:
                    st.alert = False
                    _log("ALERT OFF", f"{st.label}  STA/LTA={ratio:.2f}")

            # ── Rolling PGV / PGA (5-second window) ───────────────────────
            n5 = max(int(5.0 * sr), 10)
            if len(arr) >= n5:
                chunk        = arr[-n5:]
                st.pgv_cm_s  = counts_to_pgv(float(np.abs(chunk).max()),
                                              st.sensitivity)
                st.pga_cm_s2 = counts_to_pga(chunk, st.sensitivity, sr)

            # ── ML from P onset — uses WA simulation when waveform available ─
            # Skip ML if P detection is stale (older than EVENT_RESET) to
            # prevent pairing with a completely different event's S-wave.
            # Also cap the amplitude measurement window at 60 s after P onset
            # so late surface waves / coda / noise don't inflate the estimate.
            _p_age = (now - st.p_time) if st.p_time is not None else None
            if st.p_time is not None and (_p_age is None or _p_age <= EVENT_RESET):
                mask = t_arr >= st.p_time
                if mask.sum() > 0:
                    seg = arr[mask]
                    # Only measure peak within first 60 s of P onset
                    _ml_win = int(ML_WAVE_WIN_SEC * sr)
                    seg_ml  = seg[:_ml_win] if len(seg) > _ml_win else seg
                    pk  = float(np.abs(seg_ml).max())
                    if pk > st.event_peak:
                        st.event_peak = pk
                        dist = st.event_dist_km or st.dist_km
                        # Pre-P noise window: up to 10 s before P onset
                        # Used by WA simulation for quadrature noise correction
                        _pre_mask = ((t_arr >= st.p_time - 10.0)
                                     & (t_arr < st.p_time))
                        _pre_p = (arr[_pre_mask]
                                  if _pre_mask.sum() >= max(int(sr * 2), 10)
                                  else None)
                        st.ml_est = estimate_ml(
                            pk, st.sensitivity, dist,
                            arr=seg_ml, sr=sr,
                            is_broadband=st.is_broadband,
                            sta_lta_ratio=st.last_p_ratio,
                            p_thr_eff=p_thr,
                            pre_p_arr=_pre_p)

            if st.event_peak > 0 and st.p_time is None:
                st.event_peak = 0.0

            # ── P-picker: STA/LTA trigger + AIC onset refinement ──────────
            nps, npl = int(P_STA_SEC*sr), int(P_LTA_SEC*sr)
            if len(arr) >= nps + npl:
                pcft   = recursive_sta_lta(arr, nps, npl)
                st.p_cft = pcft       # numpy array — replaced each packet, no deque needed
                pr     = float(pcft[-1]); st.last_p_ratio = pr

                if st.p_time is None and pr >= p_thr:
                    # ── Sustained-signal guard ─────────────────────────────────
                    # Real P-waves remain elevated in the 2–8 Hz band for ≥2 s.
                    # Impulse noise (vehicles, taps, blasts) spikes and drops
                    # within one STA window.  Reject if P-band RMS in the last
                    # 2 s is not clearly above the preceding quiet background.
                    _P_SUSTAIN_SEC  = 1.5   # 1.5 s sustained signal — balance speed vs noise
                    _P_SUSTAIN_MULT = 3.5   # must be 3.5× above background
                    _ns_p    = max(int(_P_SUSTAIN_SEC * sr), 5)
                    _p_accept = True
                    if st._sos_p is not None and len(arr) >= _ns_p * 4:
                        _tail_p  = sosfilt(st._sos_p, arr[-_ns_p:])
                        _sig_rms = float(np.sqrt(np.mean(_tail_p**2)))
                        _bg_lo   = max(0, len(arr) - int(60*sr))
                        _bg_hi   = max(_bg_lo + _ns_p, len(arr) - int(10*sr))
                        if _bg_hi > _bg_lo + _ns_p:
                            _bg_p   = sosfilt(st._sos_p, arr[_bg_lo:_bg_hi])
                            _bg_rms = float(np.sqrt(np.mean(_bg_p**2))) + 1e-12
                            if _sig_rms < _P_SUSTAIN_MULT * _bg_rms:
                                _p_accept = False
                                _log("P-REJECT",
                                     f"{st.label}  pr={pr:.2f} (need {p_thr:.2f})"
                                     f"  P-band {_sig_rms:.1f}/{_bg_rms:.1f}"
                                     f"={_sig_rms/_bg_rms:.1f}× < {_P_SUSTAIN_MULT}×"
                                     f" — transient noise, skipping")

                    # ── P-wave time-coherence gate ─────────────────────────────
                    # If other stations already have P-picks and this station's
                    # trigger is >40 s *later*, it is almost certainly a separate
                    # noise event — not the same earthquake.  Reject it so we
                    # don't mix stale noise picks into the ongoing event solution.
                    if _p_accept:
                        _other_pts = [_os.p_time for _os in states.values()
                                      if _os.p_time is not None and _os is not st]
                        if _other_pts:
                            _earliest_p = min(_other_pts)
                            _this_t     = float(t_arr[-1])  # approximate trigger time
                            if (_this_t - _earliest_p) > 40.0:
                                _p_accept = False
                                _log("P-REJECT",
                                     f"{st.label}  pr={pr:.2f}  trigger"
                                     f" {_this_t - _earliest_p:.1f}s after earliest P"
                                     f" — time-incoherent with existing picks, skipping")

                    # Step 1: coarse onset from STA/LTA crossing (always computed)
                    above    = pcft >= p_thr
                    cross    = np.where(above & ~np.roll(above, 1))[0]
                    trig_idx = max(0, (cross[-1] if len(cross)
                                       else len(pcft) - 1) - nps // 2)
                    # Step 2: AIC refinement — search ±AIC_WIN_SEC around trigger
                    if st._sos_p is not None:
                        aic_win = int(AIC_WIN_SEC * sr)
                        lo = max(0, trig_idx - aic_win)
                        hi = min(len(arr), trig_idx + aic_win)
                        if hi - lo >= 8:
                            seg_p   = sosfilt(st._sos_p, arr[lo:hi])
                            aic_off = _aic_picker(seg_p)
                            aic_idx = lo + aic_off
                            if abs(aic_idx - trig_idx) <= aic_win:
                                trig_idx = aic_idx

                    # Only commit the P-wave pick if it passed the quality check
                    if not _p_accept:
                        pass   # transient — discard, wait for next packet
                    else:
                        idx       = trig_idx
                        st.p_time = float(t_arr[min(idx, len(t_arr) - 1)])
                        _log("P-WAVE",
                             f"{st.label}  rel={st.p_time-now:+.1f}s"
                             f"  picker={pr:.2f}  dist={st.dist_km}km")
                        _timeline_add(
                            f"P-wave @ {st.description} ({st.label})"
                            f"  ratio={pr:.2f}")
                        epicenter.add(key, st.lat, st.lon, st.p_time)
                        # ── Frequency content check (teleseism vs local) ───────
                        # Teleseisms dominate at 0.1–1 Hz; local events at 2–10 Hz.
                        # Compute ratio of power in these bands.  If P-band power
                        # is < 10% of low-frequency power, flag as likely teleseism.
                        if len(arr) >= int(sr * 20):
                            _fc_arr = arr[-int(sr * 20):].astype(float)
                            from scipy.signal import butter as _bt, sosfilt as _sf, welch as _wl
                            try:
                                nyq_fc = sr / 2.0
                                _sos_lo = _bt(4, [0.05/nyq_fc, 1.0/nyq_fc],
                                              btype='band', output='sos')
                                _sos_hi = _bt(4, [2.0/nyq_fc, min(8.0, nyq_fc*0.95)/nyq_fc],
                                              btype='band', output='sos')
                                _pwr_lo = float(np.mean(_sf(_sos_lo, _fc_arr)**2)) + 1e-20
                                _pwr_hi = float(np.mean(_sf(_sos_hi, _fc_arr)**2)) + 1e-20
                                _freq_ratio = _pwr_hi / _pwr_lo
                                st.freq_ratio = _freq_ratio   # store for display
                                if _freq_ratio < 0.05:   # teleseism: nearly all energy at low freq
                                    _teleseism_flag[0] = True
                                    _log("TELESEISM",
                                         f"{st.label}  P-band/LF-band={_freq_ratio:.3f}"
                                         f" < 0.05 — likely teleseism, suppressing alerts")
                            except Exception:
                                pass
                        # ── PhaseNet refinement ────────────────────────────────
                        if _phasenet_ready[0]:
                            _snap_len = min(len(arr), int(60.0 * sr))
                            _arr_snap = arr[-_snap_len:].copy()
                            threading.Thread(
                                target=_run_phasenet_pick,
                                args=(key, _arr_snap, sr, float(t_arr[-1])),
                                daemon=True).start()
                        # ── P-wave email ───────────────────────────────────────
                        with _email_lock:
                            _p_ok = (now - _email_last_p.get(key, 0.0)) >= EMAIL_COOLDOWN_P
                            if _p_ok:
                                _email_last_p[key] = now
                        if _p_ok:
                            _pt_utc  = time.strftime("%Y-%m-%d %H:%M:%S UTC",
                                                      time.gmtime(st.p_time))
                            _pgv_str = f"{st.pgv_cm_s:.4f} cm/s" if st.pgv_cm_s else "—"
                            _ncity_p = _nearest_city(st.lat, st.lon, max_dist_km=200)
                            _city_p  = (f"near {_ncity_p[0]} ({_ncity_p[1]:.0f} km)"
                                        if _ncity_p else "—")
                            _city_p_lbl = _city_label(st.lat, st.lon) or _city_p
                            _queue_email_with_shot(
                                f"P-Wave Detected — {st.description}  ({_city_p_lbl})",
                                f"P-WAVE DETECTED\n"
                                f"{'='*52}\n"
                                f"Station     : {st.label}  ({st.description})\n"
                                f"Location    : {_city_p_lbl}  (station proxy)\n"
                                f"P-arrival   : {_pt_utc}\n"
                                f"STA/LTA     : {pr:.2f}  (threshold {p_thr:.2f})"
                                f"  [{_sta_lta_level(pr, p_thr)}]\n"
                                f"Dist to {HOME_LABEL}: {st.dist_km} km\n"
                                f"PGV (5 s)   : {_pgv_str}\n"
                                f"\n"
                                f"Early detection — epicenter not yet confirmed.\n")
                            _pgv_val  = st.pgv_cm_s  or 0.0
                            _pga_val  = st.pga_cm_s2 or 0.0
                            _peak_cts = int(st.event_peak) if st.event_peak else 0
                            _ml_now   = (f"M{st.ml_est:+.2f}" if st.ml_est is not None
                                         else "not yet estimated")
                            _sens_str = f"{st.sensitivity:.3e} cts/(m/s)"
                            _cha_type = "broadband HH" if st.is_broadband else "borehole EH"
                            _ntfy_ml_ok = (NTFY_MIN_ML_PWAVE is None or
                                           st.ml_est is None or
                                           st.ml_est >= NTFY_MIN_ML_PWAVE)
                            _p_level = _sta_lta_level(pr, p_thr)
                            if _ntfy_ml_ok:
                                _pwave_title = (f"[{_p_level}] P-WAVE DETECTED"
                                                f" | {st.description} | {_city_p_lbl}")
                                _pwave_body  = (
                                    f"---- STATION ----\n"
                                    f"ID         : {st.label}\n"
                                    f"Type       : {_cha_type}\n"
                                    f"Location   : {_city_p_lbl}\n"
                                    f"Dist home  : {st.dist_km} km\n"
                                    f"\n"
                                    f"---- DETECTION ----\n"
                                    f"P-arrival  : {_pt_utc}\n"
                                    f"Rel to now : {st.p_time - now:+.1f} s\n"
                                    f"STA/LTA    : {pr:.3f}"
                                    f"  [{_sta_lta_level(pr, p_thr)}]\n"
                                    f"Threshold  : {p_thr:.2f}"
                                    f"  (base {P_THRESH} × mult {st.thresh_mult:.2f})\n"
                                    f"Level      : {_sta_lta_level(pr, p_thr)}"
                                    f"  ({pr/p_thr:.1f}× threshold)\n"
                                    f"\n"
                                    f"---- WAVEFORM ----\n"
                                    f"PGV        : {_pgv_val:.5f} cm/s\n"
                                    f"PGA        : {_pga_val:.5f} cm/s2\n"
                                    f"Peak amp   : {_peak_cts:,} counts\n"
                                    f"Sensitivity: {_sens_str}\n"
                                    f"ML estimate: {_ml_now}  (single station)\n"
                                    f"\n"
                                    f"STATUS: Early detection — epicenter not yet confirmed")
                                _send_ntfy(_pwave_title, _pwave_body, priority="high")
                                # Discord: require 2 stations + HIGH level + ML floor
                                _n_p_now   = sum(
                                    1 for _s in states.values() if _s.p_time is not None)
                                _levels_ord = ["LOW", "MEDIUM", "HIGH", "EXTREME"]
                                _p_sig = (_levels_ord.index(_p_level)
                                          >= _levels_ord.index(DISCORD_PWAVE_MIN_LEVEL))
                                _ml_sig = (DISCORD_MIN_P_ML is None or
                                           st.ml_est is None or
                                           st.ml_est >= DISCORD_MIN_P_ML)
                                if _n_p_now >= DISCORD_MIN_P_STA and _p_sig and _ml_sig:
                                    _eq_ensure_started()
                                    _p_prefix = _eq_discord_prefix()
                                    _pd_title = (f"{_p_prefix} P-Wave Detected"
                                                 f" | {st.description} | {_city_p_lbl}")
                                    _pd_bold  = (f"STA/LTA: {pr:.1f}×"
                                                 f" ({pr/p_thr:.1f}× threshold)"
                                                 f" | {_pt_utc}")
                                    _send_discord(_pd_title, _pwave_body,
                                                  level=_p_level,
                                                  bold_header=_pd_bold)
                                _play_alarm(_p_level)
                                _ws_broadcast_event("P-Wave", level=_p_level,
                                    location=_city_p_lbl,
                                    ml=st.ml_est,
                                    detail=f"STA/LTA={pr:.2f}  dist={st.dist_km}km")

                if st.p_time is not None and not st.p_cleared and pr < p_thr * 0.55:
                    st.p_cleared = True

                # S via STA/LTA re-trigger
                if (st.p_time is not None and st.s_time is None
                        and st.p_cleared and pr >= p_thr
                        and (now - st.p_time) >= 2.0):
                    st.s_time = float(t_arr[-1])
                    sp = st.s_time - st.p_time
                    _sp_min_dist = max(1.2, st.dist_km * SP_FACTOR * 0.25)
                    if sp < _sp_min_dist:
                        _log("S-REJECT",
                             f"{st.label}  S-P={sp:.1f}s < min {_sp_min_dist:.1f}s"
                             f" for {st.dist_km}km station — impossible geometry")
                        st.s_time = None
                    elif sp > SP_MAX_SEC:
                        _log("S-WAVE REJECT",
                             f"{st.label}  method=re-trig  S-P={sp:.0f}s > {SP_MAX_SEC}s"
                             "  (stale P?) — ignoring")
                        st.s_time = None
                    else:
                        st.event_dist_km = sp / SP_FACTOR
                        _log("S-WAVE",
                             f"{st.label}  method=re-trig  S-P={sp:.1f}s"
                             f"  dist≈{st.event_dist_km:.0f}km")
                        _timeline_add(
                            f"S-wave @ {st.description}  S-P={sp:.1f}s"
                            f"  dist≈{st.event_dist_km:.0f}km")
                        # ── S-wave email ───────────────────────────────────────
                        with _email_lock:
                            _s_ok = (now - _email_last_s.get(key, 0.0)) >= EMAIL_COOLDOWN_S
                            if _s_ok:
                                _email_last_s[key] = now
                        if _s_ok:
                            _st_utc  = time.strftime("%H:%M:%S UTC", time.gmtime(st.s_time))
                            _pt_utc  = time.strftime("%H:%M:%S UTC", time.gmtime(st.p_time))
                            _ncity_s = _nearest_city(st.lat, st.lon, max_dist_km=200)
                            _city_s  = (f"near {_ncity_s[0]} ({_ncity_s[1]:.0f} km)"
                                        if _ncity_s else "—")
                            _send_email(
                                f"S-Wave Detected — {st.description}  (S-P={sp:.1f}s  ~{st.event_dist_km:.0f}km)",
                                f"S-WAVE DETECTED\n"
                                f"{'='*52}\n"
                                f"Station     : {st.label}  ({st.description})\n"
                                f"Location    : {_city_s}  (station proxy)\n"
                                f"P-arrival   : {_pt_utc}\n"
                                f"S-arrival   : {_st_utc}\n"
                                f"S-P interval: {sp:.2f} s\n"
                                f"Est. distance: {st.event_dist_km:.0f} km from station\n"
                                f"Dist to {HOME_LABEL}: {st.dist_km} km\n")

                # Predicted P (before arrival) and predicted S (after P)
                with epicenter.lock:
                    if epicenter.lat is not None and epicenter.t_origin is not None:
                        h   = haversine_km(epicenter.lat, epicenter.lon,
                                           st.lat, st.lon)
                        dz  = epicenter.depth_km or 10.0
                        d3d = math.sqrt(h*h + dz*dz)
                        if st.p_time is None:
                            st.p_predicted = epicenter.t_origin + d3d / VP
                        if st.p_time is not None and st.s_time is None:
                            st.s_predicted = epicenter.t_origin + d3d / VS

                # Reset after quiet period
                if (st.p_time is not None and not st.alert
                        and (now - st.p_time) > EVENT_RESET):
                    _log("RESET",
                         f"{st.label}  quiet for {EVENT_RESET}s"
                         + (f"  final ML={st.ml_est:+.1f}" if st.ml_est else ""))
                    st.p_time = st.s_time = st.s_predicted = st.p_predicted = None
                    st.event_dist_km = None
                    st.p_cleared = False; st.event_peak = 0.0; st.ml_est = None
                    epicenter.remove_arrival(key)   # drop from epicenter solution

            # ── Build / cache bandpass filters ─────────────────────────────
            if sr != st._sos_sr and sr > 10:
                nyq        = sr / 2.0
                p_hi       = min(P_FREQ_HI, nyq * 0.90)
                s_hi       = min(S_FREQ_HI, nyq * 0.90)
                st._sos_p  = butter(4, [P_FREQ_LO/nyq, p_hi/nyq],
                                    btype='band', output='sos')
                st._sos_s  = butter(4, [S_FREQ_LO/nyq, s_hi/nyq],
                                    btype='band', output='sos')
                st._sos_sr = sr

            # ── Method 1: P/S frequency-ratio S detector ──────────────────
            if st._sos_p is not None and len(arr) >= int((PS_WIN_SEC + 6.0) * sr):
                wu    = max(int(6 * sr), int(PS_WIN_SEC * sr))
                chunk = arr[-wu - int(PS_WIN_SEC * sr):]
                fp    = sosfilt(st._sos_p, chunk)[-int(PS_WIN_SEC * sr):]
                fs    = sosfilt(st._sos_s, chunk)[-int(PS_WIN_SEC * sr):]
                p_rms = float(np.sqrt(np.mean(fp**2))) + 1e-12
                s_rms = float(np.sqrt(np.mean(fs**2))) + 1e-12
                ps_r  = p_rms / s_rms
                st.last_ps_ratio = ps_r
                st.ps_hist.append(ps_r)

                if (st.p_time is not None and st.s_time is None
                        and (now - st.p_time) >= 2.0
                        and ps_r < PS_S_THRESH
                        and len(st.ps_hist) >= 3
                        and st.ps_hist[-2] > st.ps_hist[-1]):
                    st.s_time = float(t_arr[-1])
                    sp = st.s_time - st.p_time
                    _sp_min_dist = max(1.2, st.dist_km * SP_FACTOR * 0.25)
                    if sp < _sp_min_dist:
                        _log("S-REJECT",
                             f"{st.label}  S-P={sp:.1f}s < min {_sp_min_dist:.1f}s"
                             f" for {st.dist_km}km station — impossible geometry")
                        st.s_time = None
                    elif sp > SP_MAX_SEC:
                        _log("S-WAVE REJECT",
                             f"{st.label}  method=freq-ratio  S-P={sp:.0f}s > {SP_MAX_SEC}s"
                             "  (stale P?) — ignoring")
                        st.s_time = None
                    else:
                        st.event_dist_km = sp / SP_FACTOR
                        _log("S-WAVE",
                             f"{st.label}  method=freq-ratio  S-P={sp:.1f}s"
                             f"  dist≈{st.event_dist_km:.0f}km  P/S={ps_r:.2f}")
                        _timeline_add(
                            f"S-wave @ {st.description}  method=freq-ratio"
                            f"  S-P={sp:.1f}s  dist≈{st.event_dist_km:.0f}km")

            # ── Method 2: dedicated S-band STA/LTA + AIC onset ────────────
            nss, nsl = int(S_STA_SEC * sr), int(S_LTA_SEC * sr)
            if (st._sos_s is not None
                    and st.p_time is not None and st.s_time is None
                    and (now - st.p_time) >= 2.0
                    and len(arr) >= nss + nsl):
                s_filt = sosfilt(st._sos_s, arr)
                s_cft  = recursive_sta_lta(s_filt, nss, nsl)
                st.s_cft = s_cft      # numpy array — replaced each packet, no deque needed
                sr_val = float(s_cft[-1]); st.last_s_ratio = sr_val
                if sr_val >= S_THRESH * st.thresh_mult:
                    # AIC refinement on S-filtered signal
                    above_s  = s_cft >= S_THRESH * st.thresh_mult
                    cross_s  = np.where(above_s & ~np.roll(above_s, 1))[0]
                    trig_s   = (cross_s[-1] if len(cross_s) else len(s_cft) - 1)
                    aic_win  = int(AIC_WIN_SEC * sr)
                    lo_s = max(0, trig_s - aic_win)
                    hi_s = min(len(s_filt), trig_s + aic_win)
                    if hi_s - lo_s >= 8:
                        aic_off_s = _aic_picker(s_filt[lo_s:hi_s])
                        aic_s_idx = lo_s + aic_off_s
                        if abs(aic_s_idx - trig_s) <= aic_win:
                            trig_s = aic_s_idx
                    s_t = float(t_arr[min(trig_s, len(t_arr) - 1)])
                    sp_candidate = s_t - st.p_time
                    _sp_min_dist = max(1.2, st.dist_km * SP_FACTOR * 0.25)
                    if sp_candidate < _sp_min_dist:
                        _log("S-REJECT",
                             f"{st.label}  S-P={sp_candidate:.1f}s < min {_sp_min_dist:.1f}s"
                             f" for {st.dist_km}km station — impossible geometry")
                        st.s_time = None
                    elif s_t > st.p_time + 1.0 and sp_candidate <= SP_MAX_SEC:
                        st.s_time        = s_t
                        sp               = sp_candidate
                        st.event_dist_km = sp / SP_FACTOR
                    elif sp_candidate > SP_MAX_SEC:
                        _log("S-WAVE REJECT",
                             f"{st.label}  method=S-STA/LTA+AIC  S-P={sp_candidate:.0f}s"
                             f" > {SP_MAX_SEC}s (stale P?) — ignoring")
                        _log("S-WAVE",
                             f"{st.label}  method=S-STA/LTA+AIC  S-P={sp:.1f}s"
                             f"  dist≈{st.event_dist_km:.0f}km  S-ratio={sr_val:.2f}")
                        _timeline_add(
                            f"S-wave @ {st.description}  method=S-STA/LTA+AIC"
                            f"  S-P={sp:.1f}s  dist≈{st.event_dist_km:.0f}km")


def _run_seedlink():
    url = f"{SEEDLINK_HOST}:{SEEDLINK_PORT}"
    try:
        c = QuakeClient(url)
        for net, sta, loc, cha, *_ in STATIONS:
            c.select_stream(net, sta, f"{loc}{cha}" if loc else cha)
            print(f"[SeedLink] ✓  {net}.{sta}.{loc or '--'}.{cha}")
        print(f"[SeedLink] → {url}\n")
        c.run()
    except Exception as e:
        print(f"[SeedLink] Error: {e}")


def _run_ncedc_seedlink():
    """Connect to NCEDC SeedLink server for additional nearby BK/NC stations."""
    if not NCEDC_ENABLED or not NCEDC_STATIONS:
        return
    import time as _time

    # Register NCEDC stations in the states dict using same StationState class
    for net, sta, loc, cha, lat, lon, dist, desc, sens, tmult in NCEDC_STATIONS:
        key = f"{net}.{sta}.{loc}.{cha}"
        if key not in states:
            _st = StationState(net, sta, loc, cha, lat, lon, dist, desc, sens, tmult)
            states[key] = _st
            _log("NCEDC", f"Registered station {key}")

    # Exponential back-off: 30s → 60s → 120s → … → 1800s cap
    # Only log the first failure and then once every 10 retries so the
    # log isn't flooded when ncedc.org is temporarily unreachable.
    _retry_delay   = 30
    _retry_count   = 0
    _max_delay     = 1800   # 30 minutes maximum between attempts

    while True:
        try:
            class _NCEDCClient(QuakeClient):
                pass  # reuse QuakeClient.on_data handler

            client = _NCEDCClient(f"{NCEDC_HOST}:{NCEDC_PORT}",
                                  autoconnect=False)
            client.connect()
            for net, sta, loc, cha, *_ in NCEDC_STATIONS:
                try:
                    client.select_stream(net, sta, loc + cha if loc else cha)
                except Exception:
                    try:
                        client.select_stream(net, sta, cha)
                    except Exception as e2:
                        _log("NCEDC WARN", f"{net}.{sta}: {e2}")
            _log("NCEDC", f"Connected to {NCEDC_HOST}:{NCEDC_PORT}  "
                           f"({len(NCEDC_STATIONS)} stations)")
            # Reset back-off on successful connection
            _retry_delay = 30
            _retry_count = 0
            client.run()
        except Exception as exc:
            _retry_count += 1
            # Log first failure, then only every 10th attempt
            if _retry_count == 1 or _retry_count % 10 == 0:
                _log("NCEDC", f"Unavailable (attempt {_retry_count}) — "
                               f"retrying every {_retry_delay}s "
                               f"(will back off to {_max_delay}s)")
            _time.sleep(_retry_delay)
            _retry_delay = min(_retry_delay * 2, _max_delay)


# ═══════════════════════════════════════════════════════════════════════════════
# NTFY REMOTE COMMAND LISTENER
# ═══════════════════════════════════════════════════════════════════════════════
# Command topic = alert topic + "_cmd"
# e.g.  early_earthquake_warning_san_ramon_cmd
#
# Supported commands (case-insensitive, sent as the ntfy message body):
#   reconnect  — drop and re-spawn the SeedLink thread
#   reset      — full global reset (clears P/S detections + epicenter)
#   status     — replies with a status ntfy message
#
# To send from your phone: open ntfy app → tap the (+) icon → enter the _cmd
# topic → type "reconnect" → publish.  Or from any terminal:
#   curl -d "reconnect" https://ntfy.sh/early_earthquake_warning_san_ramon_cmd
# ─────────────────────────────────────────────────────────────────────────────
NTFY_CMD_TOPIC = NTFY_TOPIC + "_cmd"   # listen on this topic for commands

# ═══════════════════════════════════════════════════════════════════════════════
# PHASENET ML PHASE PICKER
# ═══════════════════════════════════════════════════════════════════════════════

def _do_broadcast(msg):
    """Called on the asyncio event loop — broadcast msg to all WS clients."""
    if _ws_connected:
        import websockets
        websockets.broadcast(_ws_connected, msg)


def _ws_broadcast_data(data_dict):
    """Thread-safe: schedule a broadcast from any thread."""
    loop = _ws_loop[0]
    if loop is None or not _ws_connected:
        return
    import json
    try:
        msg = json.dumps(data_dict, separators=(',', ':'))
        loop.call_soon_threadsafe(_do_broadcast, msg)
    except Exception:
        pass


def _ws_broadcast_event(etype, level="", location="", ml=None, detail=""):
    """Broadcast a seismic event to the web dashboard event log."""
    _ws_broadcast_data({
        "type":     "event",
        "time":     time.time(),
        "etype":    etype,
        "level":    level,
        "location": location,
        "ml":       ml,
        "detail":   detail,
    })


def _run_ws_server():
    """Run a WebSocket server in its own asyncio event loop / daemon thread.

    Each connected client receives every SeedLink packet as a JSON message:
      { "key":  "PB.B054..EHZ",
        "t_end": 1234567890.5,   // epoch of last sample
        "sr":    100.0,
        "v":     [1234, 1235, …] // raw counts, newest packet only }

    The browser maintains its own ring buffer and draws from it.
    """
    try:
        import asyncio
        import websockets

        async def _handler(ws):
            _ws_connected.add(ws)
            _log("WS", f"Client connected  total={len(_ws_connected)}")
            try:
                await ws.wait_closed()
            finally:
                _ws_connected.discard(ws)
                _log("WS", f"Client disconnected  total={len(_ws_connected)}")

        async def _serve():
            async with websockets.serve(_handler, WS_HOST, WS_PORT):
                _log("WS", f"Live waveform server ready — ws://{WS_HOST}:{WS_PORT}")
                await asyncio.Future()   # run forever

        loop = asyncio.new_event_loop()
        _ws_loop[0] = loop
        asyncio.set_event_loop(loop)
        loop.run_until_complete(_serve())

    except ImportError:
        _log("WS", "websockets not installed — pip install websockets")
    except Exception as exc:
        _log("WS", f"Server error: {exc}")


def _run_http_server():
    """Serve index.html on http://localhost:8080.

    Accessing the site via http:// (not https://) means the browser allows
    ws://localhost:8765 without any mixed-content or certificate issues.
    Open http://localhost:8080 in your browser for the live waveform.
    """
    import http.server, os
    _dir = os.path.dirname(os.path.abspath(__file__))

    class _Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *a, **kw):
            super().__init__(*a, directory=_dir, **kw)
        def log_message(self, fmt, *args):
            pass   # silence per-request access log

    srv = http.server.HTTPServer(("localhost", 8080), _Handler)
    _log("HTTP", "Website served at  http://localhost:8080  (open this in Chrome)")
    print("[HTTP] Open your browser to:  http://localhost:8080")
    srv.serve_forever()


def _load_phasenet():
    """Load PhaseNet from seisbench in a background thread at startup.

    Sets _phasenet_ready[0] = True when the model is ready.
    Falls back silently if seisbench / torch are not installed.
    """
    if not PHASENET_ENABLED:
        return
    try:
        _log("PHASENET", f"Loading model '{PHASENET_PRETRAIN}' (first run downloads ~50 MB)...")
        import seisbench.models as sbm   # noqa
        model = sbm.PhaseNet.from_pretrained(PHASENET_PRETRAIN)
        model.eval()
        _phasenet_model[0] = model
        _phasenet_ready[0] = True
        _log("PHASENET", "Model ready — ML-assisted phase picking active")
    except ImportError:
        _log("PHASENET", "seisbench not installed — run: pip install seisbench")
    except Exception as exc:
        _log("PHASENET", f"Model load failed: {exc} — falling back to STA/LTA+AIC")


def _run_phasenet_pick(key, arr_snap, sr, t_end):
    """Run PhaseNet on a waveform snapshot in a background thread.

    Stores refined P/S picks in _phasenet_pending[key] so on_data()
    can apply them on the next packet without blocking the SeedLink thread.

    Parameters
    ----------
    key      : station key (e.g. 'BK.CMB.00.HHZ')
    arr_snap : np.ndarray of raw counts (most recent ≤60 s)
    sr       : original sample rate (Hz)
    t_end    : wall-clock epoch of the last sample in arr_snap
    """
    try:
        import torch
        from fractions import Fraction
        from scipy.signal import resample_poly

        model = _phasenet_model[0]
        if model is None:
            return

        TARGET_SR = 100.0
        WIN_LEN   = 3001          # 30.00 s at 100 Hz

        # ── Resample to 100 Hz ────────────────────────────────────────────
        if abs(sr - TARGET_SR) > 0.5:
            frac = Fraction(int(TARGET_SR), int(sr)).limit_denominator(20)
            data = resample_poly(arr_snap.astype(np.float64),
                                 frac.numerator, frac.denominator)
            sr_r = float(sr * frac.numerator / frac.denominator)
        else:
            data = arr_snap.astype(np.float64)
            sr_r = float(sr)

        if len(data) < WIN_LEN:
            return

        # Most recent 30 s window, amplitude-normalised
        window = data[-WIN_LEN:].astype(np.float32)
        std = float(window.std())
        if std < 1e-10:
            return
        window /= std

        # PhaseNet input: (batch=1, channels=3, WIN_LEN)
        # component_order = ZNE → Z at index 0; N and E filled with zeros
        x = np.zeros((1, 3, WIN_LEN), dtype=np.float32)
        x[0, 0] = window

        with torch.no_grad():
            out = model(torch.from_numpy(x))   # (1, 3, WIN_LEN)

        p_prob = out[0, 1].cpu().numpy()   # P probability trace
        s_prob = out[0, 2].cpu().numpy()   # S probability trace

        # Epoch of the first sample in the window
        t_start = t_end - (WIN_LEN - 1) / sr_r

        # P — global peak
        pi     = int(p_prob.argmax())
        pp_val = float(p_prob[pi])
        p_time = (t_start + pi / sr_r) if pp_val >= PHASENET_P_THRESH else None

        # S — search only after P + 1 s to avoid picking P coda as S
        s0     = (pi + int(sr_r * 1.0)) if p_time else (WIN_LEN // 2)
        sp_arr = s_prob.copy()
        sp_arr[:s0] = 0.0
        si     = int(sp_arr.argmax())
        sp_val = float(s_prob[si])
        s_time = (t_start + si / sr_r) if sp_val >= PHASENET_S_THRESH else None

        result = {"p_time": p_time, "p_prob": pp_val,
                  "s_time": s_time, "s_prob": sp_val}

        with _phasenet_pend_lk:
            _phasenet_pending[key] = result

        _log("PHASENET",
             f"{key}"
             f"  P={('%+.2f s' % (p_time - time.time())) if p_time else 'none':>10}"
             f"  Pprob={pp_val:.2f}"
             f"  S={('%+.2f s' % (s_time - time.time())) if s_time else 'none':>10}"
             f"  Sprob={sp_val:.2f}")

    except Exception as exc:
        _log("PHASENET", f"Inference error ({key}): {exc}")


def _run_ntfy_listener():
    """
    Stream the ntfy command topic for remote control messages.

    Uses ntfy's persistent streaming API (no poll=1) so messages are
    delivered immediately when published — no polling delay.
    The 'since' parameter is set to the program start time so commands
    sent before this session are never replayed.

    Runs in its own daemon thread; auto-reconnects on any network error.
    Only active when NTFY_ENABLED = True.

    Supported commands (send as the ntfy message body, case-insensitive):
      reconnect — re-spawn the SeedLink thread
      reset     — full global reset (clears detections + epicenter)
      status    — reply with live station/ML status notification
    """
    if not NTFY_ENABLED:
        return

    import urllib.request, json, ssl

    def _ssl_ctx():
        try:
            import certifi
            return ssl.create_default_context(cafile=certifi.where())
        except ImportError:
            return ssl._create_unverified_context()

    # Only act on commands sent after this program started
    _start_ts = str(int(time.time()))
    # After receiving a message, advance the cursor so reconnects don't replay
    _since = [_start_ts]

    _log("NTFY CMD", f"Streaming commands on topic: {NTFY_CMD_TOPIC}")
    print(f"[NTFY CMD] Remote control active — topic: {NTFY_CMD_TOPIC}")
    print(f"[NTFY CMD] Commands: reconnect | reset | status")

    def _handle(cmd):
        """Execute a validated command string on the main flag variables."""
        _log("NTFY CMD", f"Executing command: '{cmd}'")
        print(f"[NTFY CMD] Command received: '{cmd}'")

        if cmd == "reconnect":
            global _reconnect_requested
            _reconnect_requested = True
            _send_ntfy(
                "REMOTE COMMAND ACCEPTED: reconnect",
                "Reconnect command received.\n"
                "Re-spawning SeedLink thread now.\n\n"
                f"Stations: {N}  |  Server: {SEEDLINK_HOST}:{SEEDLINK_PORT}",
                priority="default")

        elif cmd == "reset":
            _animate._force_reset_requested = True
            _send_ntfy(
                "REMOTE COMMAND ACCEPTED: reset",
                "Global reset command received.\n"
                "Clearing all P/S detections and epicenter on next frame.",
                priority="default")

        elif cmd.startswith("sensitivity"):
            _parts = cmd.split()
            _sta_list = "\n".join(
                f"  {_s.sta:<6}  current={_s.thresh_mult:.2f}"
                f"  default={_DEFAULT_THRESH.get(_k, '?')}"
                f"  ({_s.description})"
                for _k, _s in states.items())

            if len(_parts) == 1:
                # "sensitivity" — reset all to defaults
                _changed = []
                for _key, _st in states.items():
                    _def = _DEFAULT_THRESH.get(_key)
                    if _def is not None and _st.thresh_mult != _def:
                        _changed.append(f"  {_st.description}: {_st.thresh_mult:.2f}→{_def:.2f}")
                        _st.thresh_mult = _def
                _log("NTFY CMD", "Sensitivity reset to defaults"
                                 + (f" ({len(_changed)} changed)" if _changed else " (already at defaults)"))
                _body = ("All sensitivities reset to defaults.\n\n"
                         + ("\n".join(_changed) if _changed else "  (all already at defaults)"))
                _send_ntfy("REMOTE: sensitivity reset to defaults", _body, priority="default")

            elif len(_parts) == 2 and _parts[1] == "list":
                # "sensitivity list" — show current values
                _send_ntfy("SENSITIVITY VALUES",
                           f"Station sensitivity (thresh_mult):\n{_sta_list}\n\n"
                           f"Lower = more sensitive  |  Range 0.10–5.00\n"
                           f"To set: sensitivity B054 0.7\n"
                           f"To reset all: sensitivity",
                           priority="default")

            elif len(_parts) == 3:
                # "sensitivity B054 0.7" — set one station
                _sta_name = _parts[1].upper()
                try:
                    _new_val = round(float(_parts[2]), 2)
                    if not (0.10 <= _new_val <= 5.00):
                        raise ValueError("out of range 0.10–5.00")
                    _matched = [(k, s) for k, s in states.items()
                                if s.sta.upper() == _sta_name]
                    if not _matched:
                        _avail = ", ".join(s.sta for s in states.values())
                        _send_ntfy("SENSITIVITY ERROR",
                                   f"Station '{_sta_name}' not found.\n"
                                   f"Available: {_avail}\n\n"
                                   f"Tip: send 'sensitivity list' to see all",
                                   priority="default")
                    else:
                        _k0, _s0 = _matched[0]
                        _old = _s0.thresh_mult
                        _s0.thresh_mult = _new_val
                        _log("NTFY CMD", f"sensitivity {_k0}: {_old:.2f}→{_new_val:.2f}")
                        _def0 = _DEFAULT_THRESH.get(_k0, "?")
                        _send_ntfy(
                            f"REMOTE: {_s0.description} sensitivity → {_new_val:.2f}",
                            f"Station   : {_s0.description} ({_s0.label})\n"
                            f"Old value : {_old:.2f}\n"
                            f"New value : {_new_val:.2f}\n"
                            f"Default   : {_def0}\n\n"
                            f"Lower = more sensitive  |  Higher = less sensitive\n"
                            f"Send 'sensitivity' (no args) to reset all to defaults",
                            priority="default")
                except ValueError as _ve:
                    _send_ntfy("SENSITIVITY ERROR",
                               f"Invalid value '{_parts[2]}': {_ve}\n"
                               f"Must be a number between 0.10 and 5.00.\n\n"
                               f"Example: sensitivity B054 0.7",
                               priority="default")

            else:
                # Unrecognised — show help
                _avail = ", ".join(s.sta for s in states.values())
                _send_ntfy("SENSITIVITY HELP",
                           f"Usage:\n"
                           f"  sensitivity              — reset all to defaults\n"
                           f"  sensitivity list         — show current values\n"
                           f"  sensitivity B054 0.7     — set B054 to 0.7\n\n"
                           f"Range: 0.10 (most sensitive) – 5.00 (least sensitive)\n"
                           f"Stations: {_avail}",
                           priority="default")

        elif cmd == "status":
            _n_live  = sum(1 for s in states.values() if s.connected)
            _n_p     = sum(1 for s in states.values() if s.p_time)
            _n_alert = sum(1 for s in states.values() if s.alert)
            _ml_vals_st = [s.ml_est for s in states.values()
                           if s.ml_est is not None]
            _ml_st   = (f"M{sorted(_ml_vals_st)[len(_ml_vals_st)//2]:+.1f}"
                        if _ml_vals_st else "none")
            _uptime  = (time.time() - _animate._last_auto_reset) / 60
            _sta_lines = "\n".join(
                f"  {s.label}: {'LIVE' if s.connected else 'DEAD'}"
                f"  P={'yes' if s.p_time else 'no'}"
                f"  ML={f'{s.ml_est:+.1f}' if s.ml_est else '—'}"
                f"  sens={s.thresh_mult:.2f}"
                for s in states.values()
            )
            _send_ntfy(
                f"REMOTE STATUS | {_n_live}/{N} live | ML={_ml_st}",
                f"Stations live : {_n_live}/{N}\n"
                f"In alert      : {_n_alert}\n"
                f"P detections  : {_n_p}\n"
                f"Current ML    : {_ml_st}\n"
                f"Since reset   : {_uptime:.0f} min\n"
                f"\n{_sta_lines}\n\n"
                f"Commands: reconnect | reset | status | sensitivity | alarm | alarm off | screenshot",
                priority="default")

        elif cmd in ("alarm", "alarm on"):
            _trigger_manual_alarm()
            _send_ntfy(
                "REMOTE COMMAND: alarm activated",
                "Continuous alarm started on the monitor.\n"
                "Send  alarm off  to stop it.",
                priority="urgent")

        elif cmd in ("alarm off", "alarm stop"):
            _stop_manual_alarm()
            _send_ntfy(
                "REMOTE COMMAND: alarm stopped",
                "Alarm has been stopped.",
                priority="default")

        elif cmd == "screenshot":
            _screenshot_email_requested[0] = True
            _send_ntfy(
                "REMOTE COMMAND: screenshot requested",
                "Capturing live waveform screenshot now.\n"
                "It will be emailed to you within a few seconds.\n\n"
                f"Destination: {EMAIL_TO if EMAIL_ENABLED else '(email not configured)'}",
                priority="default")

        else:
            _log("NTFY CMD", f"Unknown command ignored: '{cmd}'")
            print(f"[NTFY CMD] Unknown command: '{cmd}'")

    _log("NTFY CMD", "Streaming commands (real-time, no poll delay)")
    print("[NTFY CMD] Remote control active — real-time streaming")

    # ntfy sends a keepalive event every ~55 s.
    # We use a 70 s socket timeout so a stalled/dead connection is detected
    # within one keepalive interval.  The watchdog (started below) checks
    # _ntfy_last_rx_t and force-restarts if no bytes arrive within 90 s.
    _STREAM_TIMEOUT = 70
    _backoff = 5   # seconds to wait before reconnect; doubles on repeated failures

    while True:
        try:
            url = (f"{NTFY_SERVER}/{NTFY_CMD_TOPIC}/json"
                   f"?since={_since[0]}")
            req = urllib.request.Request(
                url, headers={"Accept": "application/x-ndjson"})
            with urllib.request.urlopen(req, timeout=_STREAM_TIMEOUT,
                                        context=_ssl_ctx()) as resp:
                _backoff = 5   # reset backoff on successful connection
                for raw_line in resp:
                    # Heartbeat — updated on every byte received (keepalives too)
                    _ntfy_last_rx_t[0] = time.time()

                    line = raw_line.strip()
                    if not line:
                        continue
                    try:
                        msg = json.loads(line)
                    except Exception:
                        continue

                    # Advance cursor so reconnects never replay old messages
                    if msg.get("time"):
                        _since[0] = str(int(msg["time"]) + 1)

                    if msg.get("event") != "message":
                        continue   # skip keepalives / open events

                    cmd = (msg.get("message") or "").strip().lower()
                    _handle(cmd)

        except Exception as exc:
            err = str(exc).lower()
            if "timed out" in err or "timeout" in err:
                # Normal — no keepalive in 70 s (network blip / sleep/wake)
                _log("NTFY CMD", f"Stream timeout — reconnecting in {_backoff}s")
            else:
                _log("NTFY CMD", f"Stream error: {exc} — reconnecting in {_backoff}s")
            time.sleep(_backoff)
            _backoff = min(_backoff * 2, 120)   # cap at 2 minutes


def _start_ntfy_listener():
    """Spawn a fresh ntfy listener daemon thread and record it."""
    t = threading.Thread(target=_run_ntfy_listener, daemon=True)
    t.start()
    _ntfy_listener_thread[0] = t
    return t


def _ntfy_watchdog():
    """
    Watchdog for the NTFY command listener.

    Wakes every 30 s and checks two conditions:
      1. Is the listener thread still alive?
      2. Has it received any data (keepalive or message) in the last 90 s?

    ntfy sends keepalives every ~55 s, so a 90 s silence means the
    connection is genuinely dead (sleep/wake, network drop, etc.).
    Either failure spawns a fresh listener thread.
    """
    if not NTFY_ENABLED:
        return
    _WATCHDOG_INTERVAL = 30   # check every 30 s
    _STALE_THRESHOLD   = 90   # no data for this long → restart

    while True:
        time.sleep(_WATCHDOG_INTERVAL)
        thr  = _ntfy_listener_thread[0]
        age  = time.time() - _ntfy_last_rx_t[0]
        dead = thr is None or not thr.is_alive()
        stale = _ntfy_last_rx_t[0] > 0 and age > _STALE_THRESHOLD

        if dead or stale:
            reason = "thread died" if dead else f"no data for {age:.0f}s"
            _log("NTFY CMD", f"Watchdog: {reason} — restarting listener")
            _start_ntfy_listener()


# ═══════════════════════════════════════════════════════════════════════════════
# TEST EARTHQUAKE INJECTION
# ═══════════════════════════════════════════════════════════════════════════════
def _inject_test_quake(epi_lat=37.85, epi_lon=-122.24, depth_km=8.0, ml=3.5):
    """Inject synthetic P/S arrivals to verify the detection pipeline."""
    t_now    = time.time()
    t_origin = t_now - 5.0     # event happened 5 s ago

    _log("TEST",
         f"SIMULATED M{ml} at {epi_lat:.3f}°N {epi_lon:.3f}°W "
         f"depth={depth_km}km")
    print(f"\n{'='*64}")
    print(f"  TEST EARTHQUAKE  M{ml}  {epi_lat:.3f}°N  {epi_lon:.3f}°W")
    print(f"  depth={depth_km}km   simulated t_origin = now − 5 s")
    print(f"{'='*64}\n")

    epicenter.reset()
    for key, st in states.items():
        horiz    = haversine_km(epi_lat, epi_lon, st.lat, st.lon)
        dist_3d  = math.sqrt(horiz**2 + depth_km**2)
        tp       = t_origin + dist_3d / VP
        ts       = t_origin + dist_3d / VS
        pgv_sim  = pgv_at_dist(ml, dist_3d) or 0.001
        peak_cts = pgv_sim / 100.0 * st.sensitivity   # cm/s → m/s → counts

        with st.lock:
            st.alert        = True
            st.last_alert_t = t_now
            if st.p_time is None and tp <= t_now:
                st.p_time = tp
            if st.s_time is None and ts <= t_now:
                st.s_time        = ts
                st.event_dist_km = (ts - tp) / SP_FACTOR
            st.event_peak = max(st.event_peak, peak_cts)
            st.ml_est     = estimate_ml(peak_cts, st.sensitivity, dist_3d)
            epicenter.add(key, st.lat, st.lon, tp)

    threading.Thread(
        target=lambda: _speak(
            f"Test earthquake. Simulated magnitude {ml} "
            f"near the San Francisco Bay Area."
        ), daemon=True).start()

# ═══════════════════════════════════════════════════════════════════════════════
# TERMINAL
# ═══════════════════════════════════════════════════════════════════════════════
def _print_status():
    ts = time.strftime("%Y-%m-%d  %H:%M:%S UTC", time.gmtime())
    print(f"\n{'─'*104}\n  {ts}")
    print(f"  {'Station':<22}  {'Alert':>7}  {'P-pick':>7}  "
          f"{'P/S':>6}  {'ML':>5}  {'PGV cm/s':>10}  {'PGA cm/s²':>10}  State")
    print(f"  {'─'*22}  {'─'*7}  {'─'*7}  "
          f"{'─'*6}  {'─'*5}  {'─'*10}  {'─'*10}  {'─'*12}")
    for st in states.values():
        with st.lock:
            r, pr  = st.last_ratio, st.last_p_ratio
            alrt   = st.alert;   pt  = st.p_time
            ps_r   = st.last_ps_ratio; ml = st.ml_est
            pgv, pga = st.pgv_cm_s, st.pga_cm_s2
        state  = "!! ALERT !!" if alrt else ("P-detected" if pt else "normal")
        ps_str = f"{ps_r:.2f}" if ps_r is not None else "  — "
        ml_str = f"{ml:+.1f}"  if ml  is not None else "  —"
        print(f"  {st.label:<22}  {r:>7.3f}  {pr:>7.3f}  "
              f"{ps_str:>6}  {ml_str:>5}  {pgv:>10.5f}  {pga:>10.3f}  {state}")
    with epicenter.lock:
        if epicenter.lat is not None:
            dstr = f"  depth={epicenter.depth_km:.1f}km" if epicenter.depth_km else ""
            rstr = f"  RMS={epicenter.rms_sec:.2f}s"     if epicenter.rms_sec  else ""
            print(f"\n  Epicenter: {epicenter.lat:.4f}°N  {epicenter.lon:.4f}°W"
                  f"{dstr}  ≈{epicenter.dist_home:.0f}km{rstr}"
                  f"  ({epicenter.n_sta} station{'s' if epicenter.n_sta>1 else ''})")
    print(f"{'─'*104}\n")

def _cmd():
    for line in sys.stdin:
        parts = line.strip().split()
        if not parts or parts[0] in ("s", "status"):
            _print_status()
        elif parts[0] in ("t", "test"):
            # t [lat] [lon] [depth_km] [ml]
            # e.g.  t 37.85 -122.24 8.0 3.5
            defaults = dict(epi_lat=37.85, epi_lon=-122.24, depth_km=8.0, ml=3.5)
            keys     = ["epi_lat", "epi_lon", "depth_km", "ml"]
            kw       = {}
            try:
                for k, v in zip(keys, parts[1:]):
                    kw[k] = float(v)
            except ValueError as e:
                print(f"  Bad value: {e}")
                print("  Usage: t [lat] [lon] [depth_km] [ml]")
                print(f"  Example: t 37.85 -122.24 8.0 4.2")
                continue
            params = {**defaults, **kw}
            print(f"  Injecting test quake: "
                  f"lat={params['epi_lat']}  lon={params['epi_lon']}  "
                  f"depth={params['depth_km']}km  M{params['ml']}")
            _inject_test_quake(**params)
        elif parts[0] in ("q", "quit", "exit"):
            os._exit(0)
        elif parts[0] in ("sm", "shakemap"):
            # sm <usgs_event_id>  — fetch ShakeMap + email for any USGS event
            # e.g.  sm nc75361631
            if len(parts) < 2:
                print("  Usage: sm <usgs_event_id>   e.g. sm nc75361631")
            else:
                _sm_id = parts[1].strip()
                print(f"  Fetching ShakeMap for {_sm_id} ...")
                def _sm_manual():
                    try:
                        import urllib.request as _ur, json as _js, ssl as _ssl
                        _ctx  = _ssl.create_default_context()
                        _durl = (f"https://earthquake.usgs.gov/fdsnws/event/1/query"
                                 f"?eventid={_sm_id}&format=geojson")
                        _req  = _ur.Request(_durl,
                                            headers={"User-Agent": "QuakeAlertBot/1.0"})
                        with _ur.urlopen(_req, timeout=20, context=_ctx) as _r:
                            _ev   = _js.loads(_r.read())
                        _bp   = _ev.get("properties", {})
                        _mag  = _bp.get("mag")
                        _loc  = _bp.get("place", "unknown")
                        _t0   = _bp.get("time", 0) / 1000.0
                        _utc  = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime(_t0))
                        _url  = _bp.get("url", "")
                        _magt = _bp.get("magType", "M")
                        _body = (
                            f"---- USGS EVENT ----\n"
                            f"ID         : {_sm_id}\n"
                            f"Magnitude  : {_magt} {_mag:+.1f}\n"
                            f"Location   : {_loc}\n"
                            f"Origin     : {_utc}\n"
                            f"Details    : {_url}\n"
                        )
                        _sm_bytes, _sm_key = _fetch_usgs_shakemap(_sm_id, max_attempts=2)
                        _subj = f"ShakeMap: {_magt}{_mag:+.1f} — {_loc}"
                        _send_email(_subj, _body, img_bytes=_sm_bytes)
                        _send_discord(f"🗺 ShakeMap | {_magt}{_mag:+.1f} | {_loc}",
                                      _body,
                                      bold_header=f"{_magt}{_mag:+.1f} | {_loc} | {_utc}",
                                      img_bytes=_sm_bytes)
                        print(f"  Sent: ShakeMap={'yes' if _sm_bytes else 'none'} ({_sm_key})")
                    except Exception as _sme:
                        print(f"  Error: {_sme}")
                threading.Thread(target=_sm_manual, daemon=True).start()
        elif parts[0] in ("c", "config", "settings"):
            global _settings_requested
            _settings_requested = True   # handled on main thread in _animate
        elif parts[0] == "alarm":
            # alarm        → start continuous alarm
            # alarm off / alarm stop / alarm 0 → stop alarm
            if len(parts) > 1 and parts[1].lower() in ("off", "stop", "0"):
                _stop_manual_alarm()
            else:
                _trigger_manual_alarm()
        elif cmd in ("te", "test-everyone", "testeveryone"):
            # Send a test Discord @everyone ping to verify the mention works
            _te_ml = 5.0
            print(f"  Sending test @everyone Discord message (simulated M{_te_ml:+.1f}) ...")
            _send_discord(
                title="🔔 TEST — @everyone mention check",
                message=(f"This is a test of the @everyone Discord mention.\n"
                         f"In a real alert this fires for ML ≥ {DISCORD_EVERYONE_ML:+.1f}.\n"
                         f"If you see a ping notification, it is working correctly."),
                level="HIGH",
                bold_header=f"TEST @everyone | M{_te_ml:+.1f} | simulated event",
                everyone=True)
            print("  Test message sent — check Discord for the @everyone ping.")

        else:
            print("  Commands: [s]tatus  [t]est [lat] [lon] [depth_km] [ml]"
                  "  [sm] shakemap <id>  [te]st-everyone  [c]onfig  alarm  [q]uit")
            print("  Example:  t 37.85 -122.24 8.0 4.2")
            print("  Example:  sm nc75361631")
            print("  Example:  te   (test @everyone Discord ping)")

# ═══════════════════════════════════════════════════════════════════════════════
# SETTINGS OVERLAY  (in-figure — avoids macOS second-window SIGSEGV crash)
# ═══════════════════════════════════════════════════════════════════════════════
from matplotlib.widgets import CheckButtons, Slider, Button as MplButton

_settings_overlay_axes = []   # every axes added by the overlay
_settings_refs         = []   # widget refs (prevent GC)
_settings_requested    = False

def _close_settings(_event=None):
    global _settings_overlay_axes
    for ax in _settings_overlay_axes:
        try:
            ax.remove()
        except Exception:
            pass
    _settings_overlay_axes.clear()
    _settings_refs.clear()
    try:
        fig.canvas.draw_idle()
    except Exception:
        pass

def _reset_sensitivity(_event=None):
    """Restore every station's thresh_mult to its original value from STATIONS."""
    changed = []
    for key, st in states.items():
        default = _DEFAULT_THRESH.get(key)
        if default is not None:
            if st.thresh_mult != default:
                _log("SETTINGS", f"{key}  thresh_mult reset {st.thresh_mult:.2f}→{default:.2f}")
                changed.append(f"{st.description}: {st.thresh_mult:.2f}→{default:.2f}")
                st.thresh_mult = default
    if not changed:
        _log("SETTINGS", "Reset sensitivity: already at defaults")
    else:
        _log("SETTINGS", f"Sensitivity reset to defaults ({len(changed)} stations changed)")
    # Close the overlay so the sliders re-open at the correct default positions
    _close_settings()

def _open_settings(_event=None):
    global _settings_overlay_axes, _settings_refs
    # Toggle: close if already visible
    if _settings_overlay_axes:
        _close_settings()
        return

    print("[SETTINGS] Opening overlay…")
    n = len(STATIONS)

    try:
        # Overlay covers the waveform-panel area of the main figure
        OX, OY, OW, OH = 0.032, 0.022, 0.572, 0.910

        # ── Background panel ──────────────────────────────────────────────────
        ax_bg = fig.add_axes([OX, OY, OW, OH])
        ax_bg.set_facecolor("#0b160b")
        ax_bg.set_zorder(30)
        for sp in ax_bg.spines.values():
            sp.set_edgecolor("#2a4a2a")
            sp.set_linewidth(1.8)
        ax_bg.set_xticks([]); ax_bg.set_yticks([])
        ax_bg.set_xlim(0, 1);  ax_bg.set_ylim(0, 1)
        _settings_overlay_axes.append(ax_bg)

        # Title
        ax_bg.text(0.5, 0.984, "SETTINGS",
                   transform=ax_bg.transAxes, fontsize=12, fontweight="bold",
                   ha="center", va="top", color="#88cc88")

        # Column headers
        hdr_y = 0.945
        ax_bg.text(0.04,  hdr_y, "Station",
                   transform=ax_bg.transAxes, fontsize=8,
                   color="#556655", va="top", style="italic")
        ax_bg.text(0.515, hdr_y, "Alert",
                   transform=ax_bg.transAxes, fontsize=8,
                   color="#556655", va="top", style="italic", ha="center")
        ax_bg.text(0.62,  hdr_y, "Sensitivity  (↓ = more sensitive)",
                   transform=ax_bg.transAxes, fontsize=8,
                   color="#556655", va="top", style="italic")
        # Header divider — use data coords (xlim/ylim = 0..1 matches axes fraction)
        ax_bg.plot([0, 1], [0.924, 0.924], color="#1e3a1e", lw=1.0, clip_on=False)

        # ── Per-station rows ──────────────────────────────────────────────────
        ROW_ZONE_H = 0.76          # fraction of ax_bg data height used for rows
        row_fh     = ROW_ZONE_H / n

        for i, (net, sta, loc, cha, *_) in enumerate(STATIONS):
            key = f"{net}.{sta}.{loc}.{cha}"
            st  = states[key]

            row_top_f = 0.918 - i * row_fh
            row_ctr_f = row_top_f - row_fh * 0.42

            # Row separator (data coords, same as axes fraction since xlim/ylim=0..1)
            if i > 0:
                y_sep = row_top_f + 0.004
                ax_bg.plot([0, 1], [y_sep, y_sep], color="#162816",
                           lw=0.8, clip_on=False)

            # Station name + description
            lbl_short = f"{net}.{sta}" if not loc else f"{net}.{sta}.{loc}"
            ax_bg.text(0.04, row_top_f - 0.008, lbl_short,
                       transform=ax_bg.transAxes, fontsize=9, fontweight="bold",
                       color="#dddddd", va="top")
            ax_bg.text(0.04, row_top_f - row_fh * 0.42, st.description,
                       transform=ax_bg.transAxes, fontsize=6.5,
                       color="#556655", va="center")

            # ── Alarm checkbox ────────────────────────────────────────────────
            chk_fig_x = OX + 0.445 * OW
            chk_fig_y = OY + (row_ctr_f - row_fh * 0.28) * OH
            chk_fig_w = 0.12  * OW
            chk_fig_h = row_fh * 0.55 * OH

            ax_chk = fig.add_axes([chk_fig_x, chk_fig_y, chk_fig_w, chk_fig_h])
            ax_chk.set_facecolor("#101a10")
            ax_chk.set_zorder(31)
            for sp in ax_chk.spines.values(): sp.set_edgecolor("#1e3a1e")
            chk = CheckButtons(ax_chk, ["  On"], [st.alarm_enabled])
            chk.labels[0].set_color("#cccccc")
            chk.labels[0].set_fontsize(8.5)

            def _on_alarm(_lbl, k=key):
                states[k].alarm_enabled = not states[k].alarm_enabled
                print(f"[SETTINGS] {k} alarm={'ON' if states[k].alarm_enabled else 'OFF'}")
                _log("SETTINGS", f"{k}  alarm={'ON' if states[k].alarm_enabled else 'OFF'}")

            chk.on_clicked(_on_alarm)
            _settings_refs.append(chk)
            _settings_overlay_axes.append(ax_chk)

            # ── Sensitivity slider ────────────────────────────────────────────
            sl_fig_x = OX + 0.615 * OW
            sl_fig_y = OY + (row_ctr_f - row_fh * 0.22) * OH
            sl_fig_w = 0.360 * OW
            sl_fig_h = row_fh * 0.42 * OH

            ax_sl = fig.add_axes([sl_fig_x, sl_fig_y, sl_fig_w, sl_fig_h])
            ax_sl.set_facecolor("#0b160b")
            ax_sl.set_zorder(31)
            for sp in ax_sl.spines.values(): sp.set_edgecolor("#1e3a1e")
            sl = Slider(ax_sl, "", 0.10, 3.0,
                        valinit=round(st.thresh_mult, 2),
                        color="#1ab8e8", initcolor="none")
            sl.valtext.set_color("#dddddd")
            sl.valtext.set_fontsize(8)

            def _on_mult(val, k=key):
                states[k].thresh_mult = round(val, 2)
                print(f"[SETTINGS] {k} thresh_mult={states[k].thresh_mult:.2f}")
                _log("SETTINGS", f"{k}  thresh_mult={states[k].thresh_mult:.2f}")

            sl.on_changed(_on_mult)
            _settings_refs.append(sl)
            _settings_overlay_axes.append(ax_sl)

        # ── Footer — Reset Defaults + Close (side by side) ───────────────────
        ax_bg.text(0.5, 0.038, "Changes apply instantly  —  no restart needed",
                   transform=ax_bg.transAxes, fontsize=6.5,
                   ha="center", color="#3a4a3a", style="italic")

        # ↺ Reset to Defaults button (left)
        rd_fig_x = OX + 0.04 * OW
        rd_fig_y = OY + 0.016 * OH
        rd_fig_w = 0.43 * OW
        rd_fig_h = 0.062 * OH
        ax_rd = fig.add_axes([rd_fig_x, rd_fig_y, rd_fig_w, rd_fig_h])
        ax_rd.set_zorder(31)
        btn_rd = MplButton(ax_rd, "Reset to Defaults",
                           color="#111a18", hovercolor="#1a2a30")
        btn_rd.label.set_color("#88aacc")
        btn_rd.label.set_fontsize(8.5)
        btn_rd.on_clicked(_reset_sensitivity)
        _settings_refs.append(btn_rd)
        _settings_overlay_axes.append(ax_rd)

        # ✕ Close button (right)
        cl_fig_x = OX + 0.53 * OW
        cl_fig_y = OY + 0.016 * OH
        cl_fig_w = 0.43 * OW
        cl_fig_h = 0.062 * OH
        ax_cl = fig.add_axes([cl_fig_x, cl_fig_y, cl_fig_w, cl_fig_h])
        ax_cl.set_zorder(31)
        btn_cl = MplButton(ax_cl, "✕  Close Settings",
                           color="#111a11", hovercolor="#1e2e1e")
        btn_cl.label.set_color("#778877")
        btn_cl.label.set_fontsize(8.5)
        btn_cl.on_clicked(_close_settings)
        _settings_refs.append(btn_cl)
        _settings_overlay_axes.append(ax_cl)

        fig.canvas.draw_idle()
        print("[SETTINGS] Overlay drawn")

    except Exception as e:
        print(f"[SETTINGS] Error building overlay: {e}")
        _close_settings()   # clean up any partially-added axes

# ═══════════════════════════════════════════════════════════════════════════════
# COLOUR PALETTE
# ═══════════════════════════════════════════════════════════════════════════════
C = dict(
    bg       ="#0c0c0c", bg_alert="#250404", bg_p="#0f1808",
    wave     ="#1ab8e8", cft="#f0a500",      pkl="#22e8cc",
    on       ="#ff2222", off="#2ecc71",
    p        ="#5599ff", s_det="#ff7744",    s_pred="#7dff9a",
    home     ="#ffe033", epi="#ff2222",
    dim      ="#444444", label="#888888",    bright="#dddddd",
    ring_p   ="#5599ff", ring_s="#ff7744",
    ring_unc ="#ff3333", ring_sus="#ffff44",
)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE  LAYOUT
# ═══════════════════════════════════════════════════════════════════════════════
FIG_W, FIG_H = 20.0, 14.5
BANNER_H     = 0.065
BOTTOM, TOP  = 0.025, 1.0 - BANNER_H - 0.010
LX, LW       = 0.035, 0.560
RX, RW       = 0.620, 0.365

unit_h = (TOP - BOTTOM) / N
wave_h = unit_h * 0.59
cft_h  = unit_h * 0.27
gap    = 0.005

# Speed up all line rendering — simplify short line segments during draw
plt.rcParams.update({
    "path.simplify":           True,
    "path.simplify_threshold": 1.0,
    "agg.path.chunksize":      0,    # disable chunking (faster for rasterized lines)
})

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="#080808")
try:
    fig.canvas.manager.set_window_title("Seismic Early-Warning Monitor v4")
except Exception:
    pass

# ── Banner ────────────────────────────────────────────────────────────────────
ax_ban = fig.add_axes([0.0, 1.0 - BANNER_H, 1.0, BANNER_H])
ax_ban.set_facecolor("#0a2a0a"); ax_ban.axis("off")
_ban_main = ax_ban.text(0.5, 0.63,
    "MONITORING  —  NO SEISMIC EVENT DETECTED",
    ha="center", va="center", fontsize=13, fontweight="bold",
    color=C["off"], transform=ax_ban.transAxes)
_ban_sub = ax_ban.text(0.5, 0.13, "",
    ha="center", va="center", fontsize=8.5,
    color=C["label"], transform=ax_ban.transAxes)

# ── Zoom state (shared across animation frames) ───────────────────────────────
_wave_window   = [DISPLAY_SEC]   # waveform time-window (seconds shown)
_wave_amp_zoom = [1.0]           # waveform amplitude zoom multiplier (>1 = zoomed in)

# Helper: button axes factory
def _ban_btn(x, w, label, tooltip=""):
    _bh = BANNER_H - 0.005
    _ax = fig.add_axes([x, 1.0 - BANNER_H + 0.0025, w, _bh])
    _bt = MplButton(_ax, label, color="#131f13", hovercolor="#1e311e")
    _bt.label.set_color("#99bb99"); _bt.label.set_fontsize(7.5)
    return _bt

# Time-window zoom  W+ (halve window)  /  W− (double window)
_wt_in_btn  = _ban_btn(0.686, 0.034, "W+")
_wt_out_btn = _ban_btn(0.722, 0.034, "W−")
_wt_in_btn.label.set_fontsize(8.5)
_wt_out_btn.label.set_fontsize(8.5)

# Predefined zoom steps (seconds) — clicking W+/W- snaps to nearest step
_ZOOM_STEPS = [10, 30, 60, 120, 300, 600, 1200, 1800, 3600]

def _wt_zoom_in(_e=None):
    cur = _wave_window[0]
    # Move to the next smaller step
    smaller = [s for s in _ZOOM_STEPS if s < cur]
    _wave_window[0] = smaller[-1] if smaller else _ZOOM_STEPS[0]

def _wt_zoom_out(_e=None):
    cur = _wave_window[0]
    # Move to the next larger step (capped at active buffer size)
    larger = [s for s in _ZOOM_STEPS if s > cur and s <= _display_sec_active[0]]
    _wave_window[0] = larger[0] if larger else _display_sec_active[0]

_wt_in_btn.on_clicked(_wt_zoom_in)
_wt_out_btn.on_clicked(_wt_zoom_out)

# Amplitude zoom  A+ (zoom in)  /  A− (zoom out)
_wa_in_btn  = _ban_btn(0.758, 0.034, "A+")
_wa_out_btn = _ban_btn(0.794, 0.034, "A−")
_wa_in_btn.label.set_fontsize(8.5)
_wa_out_btn.label.set_fontsize(8.5)

def _wa_zoom_in(_e=None):
    _wave_amp_zoom[0] = min(32.0, _wave_amp_zoom[0] * 2.0)
def _wa_zoom_out(_e=None):
    _wave_amp_zoom[0] = max(0.125, _wave_amp_zoom[0] / 2.0)
_wa_in_btn.on_clicked(_wa_zoom_in)
_wa_out_btn.on_clicked(_wa_zoom_out)

# Reset Sensitivity button — restores all thresh_mult to STATIONS defaults
_rss_btn = _ban_btn(0.635, 0.044, "Rst Sens")
_rss_btn.label.set_fontsize(7.5)
_rss_btn.on_clicked(_reset_sensitivity)

# Reconnect button — re-spawns the SeedLink thread without restarting the program
_reconnect_requested = False
_rcn_btn = _ban_btn(0.830, 0.042, "Reconnect")
_rcn_btn.label.set_fontsize(7.5)

def _request_reconnect(_e=None):
    global _reconnect_requested
    _reconnect_requested = True
_rcn_btn.on_clicked(_request_reconnect)

# Sound toggle button
ax_snd_btn = fig.add_axes([0.876, 1.0 - BANNER_H + 0.002, 0.060, BANNER_H - 0.004])
_snd_btn   = MplButton(ax_snd_btn, "Sound ON", color="#1a2a1a", hovercolor="#2a3a2a")
_snd_btn.label.set_color("#aaaaaa")
_snd_btn.label.set_fontsize(7.5)

def _toggle_sound(_event=None):
    global _sound_muted
    _sound_muted = not _sound_muted
    if _sound_muted:
        _snd_btn.label.set_text("Sound OFF")
        _snd_btn.label.set_color("#666666")
    else:
        _snd_btn.label.set_text("Sound ON")
        _snd_btn.label.set_color("#aaaaaa")
    _log("SOUND", "muted" if _sound_muted else "unmuted")
    fig.canvas.draw_idle()

_snd_btn.on_clicked(_toggle_sound)

# Settings button
ax_set_btn = fig.add_axes([0.940, 1.0 - BANNER_H + 0.002, 0.054, BANNER_H - 0.004])
_set_btn   = MplButton(ax_set_btn, "Settings", color="#1a2a1a", hovercolor="#2a3a2a")
_set_btn.label.set_color("#aaaaaa")
_set_btn.label.set_fontsize(7.5)
_set_btn.on_clicked(_open_settings)

# Map zoom +/- label text  (will be added after ax_map is created below)
_map_zoom_label = None

# ── Waveform + CFT axes ───────────────────────────────────────────────────────
_wax, _cax                    = [], []
_wln, _cln, _pkl_ln           = [], [], []
_pln, _sln, _spln, _ppln      = [], [], [], []
_atx, _rtx, _sptx             = [], [], []
_pgv_txt, _age_txt, _psr_txt  = [], [], []
_yscale_txt, _utc_txt         = [], []   # count-scale label + UTC timestamp
# P–S window shade (blue) and S-to-now shade (orange) — updated each frame
_p_shd, _s_shd                = [], []
# "P" and "S" text labels that float at the top of each detection line
_p_lbl, _s_lbl                = [], []

for i, (net, sta, loc, cha, lat, lon, dist, desc, _sens, _tm) in enumerate(STATIONS):
    row_idx = N - 1 - i
    cb      = BOTTOM + row_idx * unit_h + gap * 0.5
    wb      = cb + cft_h + gap

    aw = fig.add_axes([LX, wb, LW, wave_h])
    ac = fig.add_axes([LX, cb, LW, cft_h])
    for ax in (aw, ac):
        ax.set_facecolor(C["bg"])
        ax.tick_params(colors=C["label"], labelsize=6)
        for sp in ax.spines.values(): sp.set_edgecolor("#1a1a1a")

    lbl = f"{net}.{sta}.{loc}.{cha}" if loc else f"{net}.{sta}.{cha}"
    wl, = aw.plot([], [], color=C["wave"], lw=0.45,
                  antialiased=False, rasterized=True)
    # Fixed xlim — keeps tick positions stable; updated dynamically to _wave_window
    aw.set_xlim(-DISPLAY_SEC, 0); aw.xaxis.set_visible(False)
    aw.set_ylabel("")
    # Y-axis: show compact counts scale (3 ticks max, scientific notation off)
    aw.yaxis.set_major_locator(plt.MaxNLocator(3, symmetric=True))
    aw.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, _: (f"{v/1e6:.1f}M" if abs(v)>=5e5 else
                      f"{v/1e3:.0f}k"  if abs(v)>=500   else f"{v:.0f}")))
    aw.tick_params(axis="y", labelsize=4.5, colors="#555555", length=2, pad=1)
    aw.yaxis.get_offset_text().set_visible(False)
    aw.text(0.004, 0.97, f"{lbl}  ·  {desc}  ({dist} km)",
            transform=aw.transAxes, fontsize=7.5, fontweight="bold",
            color=C["bright"], va="top")

    # Age indicator (top-left, line 2)
    age_t = aw.text(0.004, 0.64, "",
                    transform=aw.transAxes, fontsize=5, color="#88aa88",
                    va="top", fontfamily="monospace")
    # UTC timestamp (bottom-left)
    utc_t = aw.text(0.004, 0.03, "",
                    transform=aw.transAxes, fontsize=4.5, color="#446644",
                    va="bottom", fontfamily="monospace")
    rt    = aw.text(0.998, 0.97, "—",
                    transform=aw.transAxes, fontsize=7, color=C["cft"],
                    va="top", ha="right", fontfamily="monospace")
    pgv_t = aw.text(0.998, 0.64, "",
                    transform=aw.transAxes, fontsize=5.5, color="#888888",
                    va="top", ha="right", fontfamily="monospace")
    # Amplitude zoom scale indicator (top-right, line 3)
    ys_t  = aw.text(0.998, 0.30, "",
                    transform=aw.transAxes, fontsize=4.5, color="#555555",
                    va="top", ha="right", fontfamily="monospace")
    at    = aw.text(0, 0, "", transform=aw.transAxes, visible=False)
    spt   = aw.text(0.004, 0.06, "",
                    transform=aw.transAxes, fontsize=6, color=C["s_det"],
                    va="bottom", fontfamily="monospace")

    # ── P/S background shades (axvspan with blended transform) ──────────────
    # P-to-S window: translucent blue; S-to-now: translucent orange
    _ps = aw.axvspan(0, 0, alpha=0.13, color=C["p"],     visible=False, zorder=0)
    _ss = aw.axvspan(0, 0, alpha=0.10, color=C["s_det"], visible=False, zorder=0)
    _p_shd.append(_ps); _s_shd.append(_ss)

    # ── Detection lines (made thicker and more solid for easier reading) ──────
    pl   = aw.axvline(-9999, color=C["p"],      lw=2.0, ls="-",
                       alpha=0.95, visible=False, zorder=6)
    sl   = aw.axvline(-9999, color=C["s_det"],  lw=2.0, ls="-",
                       alpha=0.95, visible=False, zorder=6)
    spl  = aw.axvline(-9999, color=C["s_pred"], lw=1.0, ls=":",
                       alpha=0.75, visible=False)
    ppl  = aw.axvline(-9999, color=C["p"],      lw=0.9, ls=":",
                       alpha=0.45, visible=False)   # predicted P

    # "P" and "S" floating labels at top of each detection line
    _pt = aw.text(0, 0.97, "P", transform=aw.get_xaxis_transform(),
                  fontsize=7, fontweight="bold", color=C["p"],
                  ha="center", va="top", visible=False, zorder=7,
                  bbox=dict(boxstyle="round,pad=0.15", fc="#0a1a2a", ec=C["p"],
                            lw=0.8, alpha=0.85))
    _st = aw.text(0, 0.97, "S", transform=aw.get_xaxis_transform(),
                  fontsize=7, fontweight="bold", color=C["s_det"],
                  ha="center", va="top", visible=False, zorder=7,
                  bbox=dict(boxstyle="round,pad=0.15", fc="#1a0a00", ec=C["s_det"],
                            lw=0.8, alpha=0.85))
    _p_lbl.append(_pt); _s_lbl.append(_st)

    cl,   = ac.plot([], [], color=C["cft"],   lw=0.9,
                    antialiased=False, rasterized=True, zorder=3)
    pkl2, = ac.plot([], [], color=C["pkl"],   lw=0.8, ls="--",
                    antialiased=False, rasterized=True, alpha=0.85, zorder=2)
    ac.axhline(TRIGGER_ON,  color=C["on"],  lw=0.8, ls="--", zorder=1)
    ac.axhline(TRIGGER_OFF, color=C["off"], lw=0.8, ls="--", zorder=1)
    ac.axhline(P_THRESH,    color=C["p"],   lw=0.6, ls=":",  zorder=1)
    # Fixed xlim on CFT — keeps tick spacing stable
    ac.set_xlim(-DISPLAY_SEC, 0)
    ac.set_ylim(0, TRIGGER_ON * 2.4)
    ac.set_ylabel("")
    ac.yaxis.set_ticklabels([])
    if i == N - 1:
        # Pick a "nice" tick step that rounds to the nearest clean time interval
        _raw_step = max(10, DISPLAY_SEC // 8)
        _tick_step = _raw_step
        for _ns in (10, 15, 20, 30, 60, 120, 180, 300, 600, 900, 1200, 1800, 3600):
            if _ns >= _raw_step:
                _tick_step = _ns
                break
        ac.set_xlabel("time ago", fontsize=6, color=C["label"])
        ac.tick_params(axis="x", labelsize=5.5, length=3)
        ac.xaxis.set_major_locator(plt.MultipleLocator(_tick_step))
        # Formatter: show minutes for ≥60 s values, seconds otherwise; 0 → "now"
        ac.xaxis.set_major_formatter(plt.FuncFormatter(
            lambda v, _: ("now"         if v == 0     else
                          f"{int(-v//60)}m" if -v >= 60 else
                          f"{int(-v)}s")))
    else:
        ac.xaxis.set_ticklabels([])
        ac.tick_params(axis="x", length=0)   # hide tick marks on non-bottom panels

    # Threshold labels only on the first (top) panel to avoid repetition
    if i == 0:
        yn = TRIGGER_ON * 2.4
        for val, col, lab2 in [(TRIGGER_ON,  C["on"],  f"ON={TRIGGER_ON}"),
                                (TRIGGER_OFF, C["off"], f"OFF={TRIGGER_OFF}"),
                                (P_THRESH,   C["p"],   f"P={P_THRESH}")]:
            ac.text(0.002, val/yn + 0.01, lab2, transform=ac.transAxes,
                    fontsize=4.5, color=col, va="bottom")
    psr_t = ac.text(0.998, 0.04, "",
                    transform=ac.transAxes, fontsize=5,
                    ha="right", va="bottom", color="#888888",
                    fontfamily="monospace")

    _wax.append(aw);     _cax.append(ac)
    _wln.append(wl);     _cln.append(cl);    _pkl_ln.append(pkl2)
    _pln.append(pl);     _sln.append(sl);    _spln.append(spl); _ppln.append(ppl)
    _atx.append(at);     _rtx.append(rt);    _sptx.append(spt)
    _pgv_txt.append(pgv_t); _age_txt.append(age_t); _psr_txt.append(psr_t)
    _yscale_txt.append(ys_t); _utc_txt.append(utc_t)

# ── Right panel ───────────────────────────────────────────────────────────────
total_rh = TOP - BOTTOM
tbl_h    = total_rh * 0.25
mag_h    = total_rh * 0.14
map_h    = total_rh * 0.57
rsp      = 0.007

map_bot = BOTTOM
mag_bot = map_bot + map_h + rsp
tbl_bot = mag_bot + mag_h + rsp

# Station table
ax_tbl = fig.add_axes([RX, tbl_bot, RW, TOP - tbl_bot])
ax_tbl.set_facecolor("#090909"); ax_tbl.axis("off")
ax_tbl.set_title("Station Status", fontsize=8, color=C["bright"],
                 pad=4, fontweight="bold")

TCOL = [0.02, 0.30, 0.50, 0.68, 0.84]
for ci, (hdr, col) in enumerate(zip(
        ["Station", "Alert", "P-arriv", "Dist", "P/S"],
        [C["bright"], C["cft"], C["p"], C["s_det"], "#aa88ff"])):
    ax_tbl.text(TCOL[ci], 0.97, hdr,
                transform=ax_tbl.transAxes, fontsize=6.5,
                color=col, va="top", fontweight="bold", fontfamily="monospace")
ax_tbl.plot([0, 1], [0.92, 0.92],
            color="#1e1e1e", lw=0.8, transform=ax_tbl.transAxes)

_row_h = 0.82 / (N + 1)
_tl, _tr, _tp, _td, _tps = [], [], [], [], []
for j in range(N):
    y  = 0.88 - j * _row_h
    kw = dict(transform=ax_tbl.transAxes, fontsize=6.8,
              va="top", fontfamily="monospace")
    _tl.append(ax_tbl.text(TCOL[0], y, "", color=C["label"],  **kw))
    _tr.append(ax_tbl.text(TCOL[1], y, "", color=C["label"],  **kw))
    _tp.append(ax_tbl.text(TCOL[2], y, "", color=C["p"],      **kw))
    _td.append(ax_tbl.text(TCOL[3], y, "", color=C["s_det"],  **kw))
    _tps.append(ax_tbl.text(TCOL[4], y, "", color="#aa88ff",  **kw))

_epi_txt = ax_tbl.text(0.5, 0.13, "Awaiting P-arrivals …",
    transform=ax_tbl.transAxes, fontsize=7, ha="center", va="bottom",
    color=C["dim"], fontfamily="monospace")

_acc_txt = ax_tbl.text(0.5, 0.01, "",
    transform=ax_tbl.transAxes, fontsize=7.5, fontweight="bold",
    ha="center", va="bottom", color=C["dim"], fontfamily="monospace")

# Magnitude / MMI panel  —  three columns: Source ML | Max MMI | MMI @ SanRamon
ax_mag = fig.add_axes([RX, mag_bot, RW, mag_h])
ax_mag.set_facecolor("#080f08"); ax_mag.axis("off")
ax_mag.set_title("Magnitude  &  Intensity",
                 fontsize=7.5, color=C["bright"], pad=3, fontweight="bold")

for xc, lbl2 in [(0.17, "Source ML"),
                  (0.50, "Max observed MMI"),
                  (0.83, "MMI @ San Ramon")]:
    ax_mag.text(xc, 0.12, lbl2,
                transform=ax_mag.transAxes, fontsize=5,
                ha="center", va="bottom", color=C["label"], style="italic")
for xd in [0.335, 0.665]:
    ax_mag.plot([xd, xd], [0.06, 0.94],
                color="#333333", lw=0.8, transform=ax_mag.transAxes)

_ml_src_txt    = ax_mag.text(0.17, 0.55, "M —",
    transform=ax_mag.transAxes, fontsize=16, fontweight="bold",
    color=C["dim"], ha="center", va="center")
_mmi_quake_txt = ax_mag.text(0.50, 0.55, "—",
    transform=ax_mag.transAxes, fontsize=10, fontweight="bold",
    color=C["dim"], ha="center", va="center")
_mmi_local_txt = ax_mag.text(0.83, 0.55, "—",
    transform=ax_mag.transAxes, fontsize=10, fontweight="bold",
    color=C["dim"], ha="center", va="center")
_mag_note = ax_mag.text(0.5, 0.01, "± approx; requires P-wave detection",
    transform=ax_mag.transAxes, fontsize=5,
    ha="center", va="bottom", color=C["dim"], style="italic")

_countdown_txt = ax_mag.text(0.5, 0.97, "",
    transform=ax_mag.transAxes, fontsize=9, fontweight="bold",
    ha="center", va="top", color=C["home"], fontfamily="monospace")

# ── Map panel ─────────────────────────────────────────────────────────────────
ax_map = fig.add_axes([RX, map_bot, RW, map_h])
ax_map.set_facecolor("#08111e")
for sp in ax_map.spines.values(): sp.set_edgecolor("#1e2e3e")
# Hide Web-Mercator tick labels (metres would show as "×10^7")
ax_map.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax_map.set_title(
    "Epicenter Map  (OSM – use toolbar to zoom & pan)"
    if HAS_TILES else "Epicenter Map",
    fontsize=7, color=C["bright"], pad=3)

# Initial map extent: Bay Area + all stations tightly framed
if HAS_TILES:
    _mx0, _my0 = _to_merc.transform(-124.2, 36.2)   # SW corner
    _mx1, _my1 = _to_merc.transform(-119.8, 39.8)   # NE corner
else:
    _mx0, _my0 = -124.2, 36.2
    _mx1, _my1 = -119.8, 39.8

ax_map.set_xlim(_mx0, _mx1)
ax_map.set_ylim(_my0, _my1)
# Equal aspect so Web Mercator x/y (both metres) scale identically on zoom
ax_map.set_aspect('equal', adjustable='datalim')
# Critical: prevent ax_map.plot() calls in the animation loop from
# auto-scaling and silently resetting the user's zoom level
ax_map.set_autoscale_on(False)

# ── Map zoom buttons (+/−) ─────────────────────────────────────────────────────
# Placed in the top-right corner of the map panel area
_mzoom_frac = 0.028    # button height as fraction of figure
_mzoom_bw   = 0.030    # button width
_mzoom_y    = map_bot + map_h - _mzoom_frac - 0.002
_mzin_ax  = fig.add_axes([RX + RW - 2*_mzoom_bw - 0.004, _mzoom_y, _mzoom_bw, _mzoom_frac])
_mzout_ax = fig.add_axes([RX + RW - _mzoom_bw,            _mzoom_y, _mzoom_bw, _mzoom_frac])
_mzin_btn  = MplButton(_mzin_ax,  "+", color="#081408", hovercolor="#0f220f")
_mzout_btn = MplButton(_mzout_ax, "−", color="#081408", hovercolor="#0f220f")
for _b in (_mzin_btn, _mzout_btn):
    _b.label.set_color("#55aa55"); _b.label.set_fontsize(13); _b.label.set_fontweight("bold")

def _map_zoom_in(_e=None):
    global _tile_pending, _last_tile_t
    xl = ax_map.get_xlim(); yl = ax_map.get_ylim()
    cx = (xl[0]+xl[1])/2;  cy = (yl[0]+yl[1])/2
    dx = (xl[1]-xl[0])/4;  dy = (yl[1]-yl[0])/4
    ax_map.set_xlim(cx-dx, cx+dx); ax_map.set_ylim(cy-dy, cy+dy)
    _tile_pending = True; _last_tile_t = 0.0   # trigger immediate tile refresh

def _map_zoom_out(_e=None):
    global _tile_pending, _last_tile_t
    xl = ax_map.get_xlim(); yl = ax_map.get_ylim()
    cx = (xl[0]+xl[1])/2;  cy = (yl[0]+yl[1])/2
    dx = (xl[1]-xl[0]);     dy = (yl[1]-yl[0])
    ax_map.set_xlim(cx-dx, cx+dx); ax_map.set_ylim(cy-dy, cy+dy)
    _tile_pending = True; _last_tile_t = 0.0

_mzin_btn.on_clicked(_map_zoom_in)
_mzout_btn.on_clicked(_map_zoom_out)

# Static CA polygon fallback when no tiles
if not HAS_TILES:
    _CA = np.array([
        (-124.25,42.00),(-120.00,42.00),(-120.00,39.00),(-119.32,38.50),
        (-119.00,37.50),(-118.20,36.50),(-116.50,35.75),(-114.63,35.00),
        (-114.62,32.73),(-117.13,32.53),(-118.19,33.73),(-119.05,34.04),
        (-120.66,34.58),(-121.33,35.79),(-121.90,36.96),(-122.17,37.20),
        (-122.51,37.73),(-122.53,38.01),(-122.98,38.10),(-123.38,38.56),
        (-123.70,38.85),(-123.97,39.84),(-124.24,40.30),(-124.25,42.00),
    ])
    from matplotlib.patches import Polygon as _MplPoly
    ax_map.add_patch(_MplPoly(_CA, closed=True,
                               facecolor="#0d1f2d", edgecolor="#1e3a50",
                               lw=1.2, zorder=1))
    ax_map.set_aspect(1.0 / math.cos(math.radians(38.0)))
    ax_map.grid(True, color="#0e1a28", lw=0.5, zorder=0)

# Station markers + home star (static)
_STA_COLS = ["#1ab8e8","#f0a500","#2ecc71","#e056a0","#a78bfa","#fb923c"]
_map_static_artists = []
_sta_map_xy = {}   # key → (px, py, lat, lon)

for idx, (net, sta, loc, cha, lat, lon, dist, desc, _s, _t) in enumerate(STATIONS):
    col       = _STA_COLS[idx]
    px, py    = _proj(lat, lon)
    key_      = f"{net}.{sta}.{loc}.{cha}"
    _sta_map_xy[key_] = (px, py, lat, lon)
    mk, = ax_map.plot(px, py, "o", color=col, ms=7, zorder=12,
                      markeredgecolor="#ffffff", markeredgewidth=0.5)
    lb  = ax_map.text(px, py, f"  {sta}", fontsize=5.5, color=col,
                      va="center", zorder=13)
    _map_static_artists.extend([mk, lb])

_hx, _hy = _proj(HOME_LAT, HOME_LON)
hm, = ax_map.plot(_hx, _hy, "*", color=C["home"], ms=12, zorder=14,
                  markeredgecolor="#888800", markeredgewidth=0.4)
hl  = ax_map.text(_hx, _hy, "  San Ramon", fontsize=5.5,
                  color=C["home"], va="center", zorder=15)
_map_static_artists.extend([hm, hl])

# Dynamic artists (epicenter marker, rings) — rebuilt every frame
_epi_mk,   = ax_map.plot([], [], "+", color=C["epi"], ms=20,
                          mew=3.0, zorder=22)
_epi_lbl   = ax_map.text(0, 0, "", fontsize=6.5, color=C["epi"],
                          va="top", zorder=23, visible=False,
                          fontweight="bold")
_epi_home, = ax_map.plot([], [], color=C["epi"], lw=0.8, ls="-",
                          alpha=0.30, zorder=6)
_dyn_art   = []    # rebuilt every frame

# Tile management (zoom/pan refresh)
# Tiles are fetched in a background thread to avoid blocking the GUI.
_tile_pending   = False
_tile_ready     = False    # set True by bg thread when new tiles are fetched
_tile_fetching  = False    # guard against overlapping fetches
_last_tile_t    = 0.0
_tile_xl_req    = [None]   # xlim at the time the fetch was requested
_tile_yl_req    = [None]   # ylim at the time the fetch was requested
TILE_DEBOUNCE   = 2.0      # seconds of stillness before refetching

def _on_map_xlim_change(_ax):
    global _tile_pending
    _tile_pending = True

def _fetch_tiles_bg():
    """Background thread: fetch OSM tiles for the current map extent."""
    global _tile_fetching, _tile_ready
    if not HAS_TILES or _tile_fetching:
        return
    _tile_fetching = True
    try:
        xl = tuple(_tile_xl_req[0]) if _tile_xl_req[0] else None
        yl = tuple(_tile_yl_req[0]) if _tile_yl_req[0] else None
        if xl is None or yl is None:
            return
        # Remove old tile images and add new ones on main thread via flag
        # (contextily is not thread-safe for the actual Axes manipulation,
        # so we only fetch here and let _animate apply them)
        # ── actually contextily downloads tiles independently of Axes ──
        # We set _tile_ready and let _animate call ctx.add_basemap() itself,
        # but on the debounced schedule.  The real benefit of this function
        # is enforcing that only one fetch happens at a time.
        _tile_ready = True   # signal _animate to call ctx.add_basemap
    except Exception as _te:
        print(f"[MAP] Tile-fetch thread error: {_te}")
    finally:
        _tile_fetching = False

if HAS_TILES:
    try:
        # Load initial tiles for the framed NorCal extent (zoom 9 ≈ 150 m/px)
        ctx.add_basemap(ax_map, crs="EPSG:3857",
                        source=ctx.providers.CartoDB.DarkMatter,
                        zoom=9, reset_extent=False, attribution=False)
        for a in _map_static_artists: a.set_zorder(20)
        print("[MAP] Initial OSM tiles loaded ✓")
    except Exception as _te:
        print(f"[MAP] Initial tile fetch failed (offline?): {_te}")
    # Connect AFTER initial load so the first tile-add doesn't set _tile_pending
    ax_map.callbacks.connect("xlim_changed", _on_map_xlim_change)
    _last_tile_t = time.time()   # suppress immediate re-trigger

# Waveform y-axis smooth-scaling state
_ylw = [[1.0, -1.0] for _ in range(N)]
_ylc = [TRIGGER_ON * 2.4] * N

# Maximum display points per waveform panel (downsampling target).
# At 100 sps × 3 600 s = 360 000 raw points → downsample to 3 000.
# 3 000 pts gives ~2 s/pixel at full zoom-out; zooming in raises effective res.
_WDISP_MAX = 3000

# ═══════════════════════════════════════════════════════════════════════════════
# ANIMATION
# ═══════════════════════════════════════════════════════════════════════════════
def _vline(vl, t_ep, now):
    if t_ep is not None:
        xv = t_ep - now
        if -DISPLAY_SEC <= xv <= 0:
            vl.set_xdata([xv, xv]); vl.set_visible(True); return
    vl.set_visible(False)

def _event_quality_score(en, erms, eaz_gap, ml_unc, sp_coherent=True):
    """
    Compute a 0–100 event quality score from four orthogonal indicators.
    Higher = more confident location and magnitude estimate.
    """
    score = 0
    # Station count (up to 30 pts)
    if   en >= 4: score += 30
    elif en >= 3: score += 22
    elif en >= 2: score += 12
    # RMS residual (up to 30 pts)
    if erms is not None:
        if   erms < 0.15: score += 30
        elif erms < 0.40: score += 20
        elif erms < 1.00: score += 10
    else:
        score += 8   # proxy solution, partial credit
    # Azimuthal gap (up to 25 pts)
    if eaz_gap is not None:
        if   eaz_gap < 90:  score += 25
        elif eaz_gap < 180: score += 15
        elif eaz_gap < 270: score += 6
    # ML uncertainty (up to 15 pts)
    if ml_unc is not None:
        if   ml_unc < 0.25: score += 15
        elif ml_unc < 0.40: score += 10
        elif ml_unc < 0.60: score += 5
    else:
        score += 5   # single station
    # S-P coherence penalty
    if not sp_coherent:
        score = max(0, score - 20)
    return min(100, score)


def _animate(_frame):
    global _tile_pending, _last_tile_t, _settings_requested, _reconnect_requested
    now       = time.time()
    any_alert = False
    any_p     = False
    n_p_sta   = 0

    # ── Confirmed-event timeout: force final report + reset if event persists
    #    too long without a natural end (e.g. surface-wave coda keeps STA/LTA
    #    elevated so per-station resets never fire).
    #    EVENT_RESET (30 min) quiet timer + 5 min grace = 35 min hard cap.
    _EVENT_HARD_CAP = EVENT_RESET + 300   # 35 minutes
    if _animate._was_confirmed:
        _snap = _animate._last_confirmed_snap
        _et0  = _snap.get("et0") or _snap.get("_snap_time", now)
        if not hasattr(_animate, "_confirmed_since"):
            _animate._confirmed_since = now
        if (now - _animate._confirmed_since) >= _EVENT_HARD_CAP:
            _log("EVENT TIMEOUT",
                 f"Confirmed event exceeded {_EVENT_HARD_CAP}s — forcing final report")
            _global_reset(reason=f"event timeout {_EVENT_HARD_CAP}s")
            # After reset epicenter is gone → confirmed=False → report fires below
    else:
        _animate._confirmed_since = now   # reset the clock when no event

    # ── Auto global reset every AUTO_RESET_SEC ────────────────────────────
    if now - _animate._last_auto_reset >= AUTO_RESET_SEC:
        _animate._last_auto_reset = now
        _global_reset(reason=f"timer {AUTO_RESET_SEC}s")

    # ── Open settings window on main thread (safe for macOS) ─────────────
    if _settings_requested:
        _settings_requested = False
        _open_settings()

    # ── Remote global reset (triggered by ntfy "reset" command) ─────────
    if getattr(_animate, "_force_reset_requested", False):
        _animate._force_reset_requested = False
        _global_reset(reason="remote ntfy command")

    # ── SeedLink reconnect (triggered by Reconnect button OR ntfy "reconnect") ──
    if _reconnect_requested:
        _reconnect_requested = False
        _animate._last_auto_rc_time   = now   # suppress auto-reconnect overlap
        _animate._last_seedlink_spawn = now
        _log("RECONNECT", "User-requested SeedLink reconnect — resetting station state")
        for st in states.values():
            with st.lock:
                st.connected      = False
                st.last_data_time = None
                # Reset calibration so it re-runs cleanly on fresh data
                st._calib_done    = False
                if hasattr(st, "_calib_ratios"):  st._calib_ratios.clear()
                if hasattr(st, "_calib_start"):   st._calib_start = None
        threading.Thread(target=_run_seedlink, daemon=True).start()
        _log("RECONNECT", "New SeedLink thread spawned")

    # ── Per-station dead detection — retry every 60 s until back online ──────
    # The all-dead check below only fires when EVERY station is silent.
    # This block handles the common case of one station going offline while
    # others remain alive: mark it dead immediately and reconnect every 60 s
    # until it recovers.  last_data_time is NOT reset so the check keeps
    # firing on the 60 s cadence regardless of intermediate reconnects.
    _STATION_DEAD_SEC = 60
    _startup_elapsed  = now - getattr(_animate, "_last_seedlink_spawn", now)
    if _startup_elapsed > 30:   # grace period after spawn
        _dead_sta = [(k, st) for k, st in states.items()
                     if (st.last_data_time is not None
                         and (now - st.last_data_time) > _STATION_DEAD_SEC)]
        if _dead_sta:
            for _dk, _dst in _dead_sta:
                if _dst.connected:
                    _log("DEAD",
                         f"{_dk} silent >{_STATION_DEAD_SEC}s — marking offline")
                    _dst.connected = False
            _rc_gap = now - getattr(_animate, "_last_auto_rc_time", 0)
            if _rc_gap >= _STATION_DEAD_SEC:
                _dead_labels = [st.label for _, st in _dead_sta]
                _log("RECONNECT",
                     f"Dead station(s) {_dead_labels} — retrying SeedLink")
                _animate._last_auto_rc_time   = now
                _animate._last_seedlink_spawn  = now
                # Do NOT reset last_data_time here: keeping the old timestamp
                # ensures the check re-fires 60 s later if the station is
                # still dead after the reconnect.
                threading.Thread(target=_run_seedlink, daemon=True).start()

    # ── All-stations-dead reconnect: sleep/network-drop recovery ─────────────
    # Grace period: don't check until 30s after last spawn so the initial
    # connection has time to deliver its first packets.
    _silence_threshold = 90   # seconds since last packet before declaring dead
    _startup_grace     = 30   # seconds to wait before first liveness check
    _retry_cooldown    = 120  # minimum seconds between reconnect attempts
    _time_since_spawn  = now - getattr(_animate, "_last_seedlink_spawn", now)
    if _time_since_spawn > _startup_grace:
        _any_live = any(
            st.last_data_time is not None
            and (now - st.last_data_time) < _silence_threshold
            for st in states.values()
        )
        if not _any_live:
            _since_last_rc = now - getattr(_animate, "_last_auto_rc_time", 0)
            if _since_last_rc >= _retry_cooldown:
                _animate._last_auto_rc_time = now
                _animate._last_seedlink_spawn = now
                _animate._auto_reconnect_pending = True
                _log("RECONNECT",
                     f"No live stations for >{_silence_threshold}s "
                     f"(sleep/network drop?) — auto-reconnecting")
                for st in states.values():
                    with st.lock:
                        st.connected      = False
                        st.last_data_time = None
                        st._calib_done    = False
                        if hasattr(st, "_calib_ratios"): st._calib_ratios.clear()
                        if hasattr(st, "_calib_start"):  st._calib_start = None
                threading.Thread(target=_run_seedlink, daemon=True).start()
        else:
            _animate._auto_reconnect_pending = False

    # ── Refresh OSM tiles after zoom/pan (debounced, non-blocking) ───────
    if HAS_TILES and _tile_pending and (now - _last_tile_t) > TILE_DEBOUNCE:
        _tile_pending    = False
        _last_tile_t     = now
        _tile_xl_req[0]  = ax_map.get_xlim()
        _tile_yl_req[0]  = ax_map.get_ylim()
        # Kick off a background thread so tile download doesn't freeze the GUI.
        # The thread sets _tile_ready; we apply the result on the NEXT frame.
        if not _tile_fetching:
            threading.Thread(target=_fetch_tiles_bg, daemon=True).start()

    # ── P-wave corroboration watchdog ────────────────────────────────────────
    # If only a single station has a P-pick and no second station has confirmed
    # within 35 s, the pick is almost certainly noise — auto-clear it so it
    # cannot accumulate into a false event or keep the alert state hot.
    # We do NOT clear during a confirmed event (≥2 stations) to avoid
    # disrupting a legitimately confirmed earthquake.
    _P_WATCHDOG_SEC = 35.0
    _p_sta_with_pick = [(k, st) for k, st in states.items()
                        if st.p_time is not None and not st.alert]
    _n_p_picks = len(_p_sta_with_pick)
    if _n_p_picks == 1:
        _solo_key, _solo_st = _p_sta_with_pick[0]
        _solo_age = now - _solo_st.p_time
        if _solo_age > _P_WATCHDOG_SEC:
            _log("P-WATCHDOG",
                 f"{_solo_st.label}  P-pick unconfirmed after {_solo_age:.0f}s"
                 f" — no second station — auto-clearing as likely noise")
            with _solo_st.lock:
                _solo_st.p_time    = None
                _solo_st.s_time    = None
                _solo_st.p_cleared = False
                _solo_st.alert     = False

    # ── STA/LTA @everyone alert ──────────────────────────────────────────────
    # Fires an @everyone Discord ping the moment ANY station's general STA/LTA
    # ratio reaches TRIGGER_ON (≥8.0), before P-wave pick or epicenter solve.
    # Also fires when EQ is confirmed by ≥2 stations (the confirmed block below
    # already does this via everyone=True, so this covers the pre-confirm case).
    if (now - _sta_lta_everyone_last_t[0]) >= STA_LTA_EVERYONE_COOLDOWN:
        _sle_triggered = [(k, st) for k, st in states.items()
                          if st.connected and st.last_ratio is not None
                          and st.last_ratio >= TRIGGER_ON * st.thresh_mult]
        if _sle_triggered:
            _sta_lta_everyone_last_t[0] = now
            # Pick station with highest ratio
            _sle_top = max(_sle_triggered,
                           key=lambda x: x[1].last_ratio / (TRIGGER_ON * x[1].thresh_mult))
            _sle_key, _sle_st = _sle_top
            _sle_ratio = _sle_st.last_ratio
            _sle_thr   = TRIGGER_ON * _sle_st.thresh_mult
            _sle_level = _sta_lta_level(_sle_ratio, _sle_thr)
            _sle_lines = "\n".join(
                f"  {k:<28}  STA/LTA={st.last_ratio:.2f}"
                f"  ({st.last_ratio/(TRIGGER_ON*st.thresh_mult):.1f}× thr)  [{_sta_lta_level(st.last_ratio, TRIGGER_ON*st.thresh_mult)}]"
                for k, st in _sle_triggered)
            _sle_body = (
                f"⚠️  SEISMIC SIGNAL DETECTED  ⚠️\n"
                f"Station(s) triggered : {len(_sle_triggered)}\n"
                f"Highest ratio        : {_sle_key}  {_sle_ratio:.2f}"
                f"  ({_sle_ratio/_sle_thr:.1f}× threshold)\n"
                f"Level                : {_sle_level}\n"
                f"\n"
                f"---- TRIGGERED STATIONS ----\n"
                f"{_sle_lines}\n"
                f"\n"
                f"Epicenter not yet confirmed.\n"
                f"Stay alert — this may indicate a nearby earthquake.\n"
                f"https://myearthquake.dpdns.org/")
            _send_discord(
                f"[{_sle_level}] ⚡ Seismic Trigger | {_sle_st.description}",
                _sle_body,
                level=_sle_level,
                bold_header=f"STA/LTA {_sle_ratio:.1f} ({_sle_ratio/_sle_thr:.1f}× thr) | {_sle_st.description}",
                everyone=True)
            _log("STA/LTA @EVERYONE",
                 f"{_sle_level}  {_sle_key}  ratio={_sle_ratio:.2f}"
                 f"  ({_sle_ratio/_sle_thr:.1f}× thr)  @everyone sent")

    # ── Strong-signal alert (catches teleseisms + large regionals) ───────────
    # Fires when STRONG_SIG_MIN_STA stations simultaneously hit
    # STRONG_SIG_MIN_LEVEL without requiring a P-wave detection.
    # Teleseisms are long-period and never excite the 2–8 Hz P-band filter,
    # so they never reach P_THRESH.  This parallel path catches them.
    if STRONG_SIG_ENABLED and (now - _strong_sig_last_t[0]) >= STRONG_SIG_COOLDOWN:
        _lvl_ord    = ["LOW", "MEDIUM", "HIGH", "EXTREME"]
        _min_idx    = _lvl_ord.index(STRONG_SIG_MIN_LEVEL)
        _ss_sta     = []   # stations currently at or above the threshold level
        for _ssk, _sst in states.items():
            if not _sst.connected or _sst.last_ratio is None:
                continue
            _ss_thr = TRIGGER_ON * _sst.thresh_mult
            _ss_lvl = _sta_lta_level(_sst.last_ratio, _ss_thr)
            if _lvl_ord.index(_ss_lvl) >= _min_idx:
                _ss_sta.append((_ssk, _sst, _ss_lvl, _sst.last_ratio, _ss_thr))
        if len(_ss_sta) >= STRONG_SIG_MIN_STA:
            # Find the highest level across all triggered stations
            _ss_peak_lvl = max(_ss_sta, key=lambda x: _lvl_ord.index(x[2]))[2]
            _strong_sig_last_t[0] = now
            _ss_lines = "\n".join(
                f"  {_ssk:<28}  STA/LTA={_sr:.2f}  ({_sr/_st:.1f}× thr)  [{_sl}]"
                for _ssk, _, _sl, _sr, _st in _ss_sta)
            _ss_title = (f"[{_ss_peak_lvl}] STRONG SIGNAL"
                         f" — {len(_ss_sta)} stations")
            _ss_body  = (
                f"---- STRONG SIGNAL DETECTED ----\n"
                f"Level        : {_ss_peak_lvl}\n"
                f"Stations hit : {len(_ss_sta)} of {len(states)}\n"
                f"Min threshold: {STRONG_SIG_MIN_LEVEL}  ({STRONG_SIG_MIN_STA}+ stations)\n"
                f"\n"
                f"NOTE: This alert fires on raw STA/LTA energy without requiring\n"
                f"a P-wave pick.  Likely causes: teleseism (distant M5+), regional\n"
                f"M4+ event, or strong local shaking with slow onset.\n"
                f"\n"
                f"---- STATION LEVELS ----\n"
                f"{_ss_lines}\n"
                f"\n"
                f"P-wave detection did NOT trigger (normal for teleseisms —\n"
                f"long-period energy is outside the 2–8 Hz P-band filter).")
            _log("STRONG-SIGNAL",
                 f"{_ss_peak_lvl}  {len(_ss_sta)} stations  — "
                 f"raw STA/LTA alert (no P-pick needed)")
            _timeline_add(f"Strong signal: {_ss_peak_lvl}  {len(_ss_sta)} sta")
            _send_ntfy(_ss_title, _ss_body, priority="high")
            _send_discord(_ss_title, _ss_body, level=_ss_peak_lvl)
            threading.Thread(
                target=lambda b=_ss_body, t=_ss_title: _send_email(t, b),
                daemon=True).start()
            _play_alarm(_ss_peak_lvl)

    global _tile_ready
    if HAS_TILES and _tile_ready and not _tile_fetching:
        _tile_ready = False
        try:
            _xl = _tile_xl_req[0] or ax_map.get_xlim()
            _yl = _tile_yl_req[0] or ax_map.get_ylim()
            for img in list(ax_map.images):
                img.remove()
            ctx.add_basemap(ax_map, crs="EPSG:3857",
                            source=ctx.providers.CartoDB.DarkMatter,
                            reset_extent=False, attribution=False)
            ax_map.set_xlim(_xl)
            ax_map.set_ylim(_yl)
            ax_map.set_aspect('equal', adjustable='datalim')
            for a in _map_static_artists:
                a.set_zorder(20)
        except Exception:
            pass

    # ── Waveform panels ───────────────────────────────────────────────────
    for i, (net, sta, loc, cha, _la, _lo, _dist, _desc, _s, _t) in enumerate(STATIONS):
        key = f"{net}.{sta}.{loc}.{cha}"
        st  = states[key]
        with st.lock:
            if st.samples is None or len(st.samples) < 2:
                _rtx[i].set_text("no data"); continue
            t_arr  = st.times.to_array()    # zero-copy view when possible
            samp   = st.samples.to_array()  # zero-copy view when possible
            cft    = st.cft                 # already numpy array or None
            pcft   = st.p_cft              # already numpy array or None
            alert  = st.alert;        ratio  = st.last_ratio
            pr     = st.last_p_ratio; ps_r   = st.last_ps_ratio
            pt     = st.p_time;       pp     = st.p_predicted
            st2    = st.s_time;       sp2    = st.s_predicted
            edist  = st.event_dist_km; conn  = st.connected
            ml     = st.ml_est
            pgv    = st.pgv_cm_s;     pga    = st.pga_cm_s2
            ldt    = st.last_data_time

        rel = t_arr - now

        # ── Apply time window — clip to current zoom window ───────────────
        win = _wave_window[0]
        mask_win = rel >= -win
        if mask_win.sum() >= 2:
            samp_w = samp[mask_win]
            rel_w  = rel[mask_win]
        else:
            samp_w, rel_w = samp, rel

        # ── Downsample for display (keeps rendering fast at any window size) ──
        n_pts = len(samp_w)
        if n_pts > _WDISP_MAX:
            step   = n_pts // _WDISP_MAX
            samp_d = samp_w[::step]
            rel_d  = rel_w[::step]
        else:
            samp_d, rel_d = samp_w, rel_w

        # Waveform — fixed xlim keeps tick positions stable
        _wln[i].set_data(rel_d, samp_d)
        _wax[i].set_xlim(-win, 0)

        # Amplitude limits with zoom multiplier applied
        lo  = float(np.percentile(samp_d, 1)); hi = float(np.percentile(samp_d, 99))
        pad = max((hi - lo) * 0.18, 1.0); lo -= pad; hi += pad
        pl0, ph0 = _ylw[i]
        _ylw[i] = [min(lo, pl0*.995+lo*.005), max(hi, ph0*.995+hi*.005)]
        # Apply amplitude zoom: narrow the y-range around the centre
        if _wave_amp_zoom[0] != 1.0:
            cy  = (_ylw[i][0] + _ylw[i][1]) / 2
            half = (_ylw[i][1] - _ylw[i][0]) / 2 / _wave_amp_zoom[0]
            _wax[i].set_ylim(cy - half, cy + half)
        else:
            _wax[i].set_ylim(*_ylw[i])

        # Y-scale indicator (peak half-range in counts)
        cur_ylim = _wax[i].get_ylim()
        half_rng = (cur_ylim[1] - cur_ylim[0]) / 2
        if   half_rng >= 5e5: _yscale_txt[i].set_text(f"±{half_rng/1e6:.1f}M cts")
        elif half_rng >= 500: _yscale_txt[i].set_text(f"±{half_rng/1e3:.0f}k cts")
        else:                 _yscale_txt[i].set_text(f"±{half_rng:.0f} cts")

        # ── UTC timestamp (right edge = now) ─────────────────────────────
        _utc_txt[i].set_text(time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime(now)))

        # ── Last-update age ───────────────────────────────────────────────
        if ldt:
            age_s = now - ldt
            if age_s < 2:      age_str = "● live"
            elif age_s < 60:   age_str = f"upd {age_s:.0f}s ago"
            else:              age_str = f"upd {age_s/60:.1f}min ago"
            _age_txt[i].set_text(age_str)
        else:
            _age_txt[i].set_text("no data")

        # Ratio + ML label
        ml_str = f"  ML{ml:+.1f}" if ml is not None else ""
        _rtx[i].set_text(f"A:{ratio:.2f}  P:{pr:.2f}{ml_str}")
        _rtx[i].set_color(C["on"] if alert else (C["p"] if pt else C["cft"]))

        # PGV / PGA label
        _pgv_txt[i].set_text(f"PGV {pgv:.4f} cm/s   PGA {pga:.3f} cm/s²")

        _atx[i].set_visible(alert)
        _vline(_pln[i],  pt,  now)
        _vline(_sln[i],  st2, now)
        _vline(_spln[i], sp2, now)

        # ── P/S background shades + floating text labels ──────────────────────
        # axvspan returns a Rectangle in matplotlib 3.5+; update via set_x/set_width
        # (x coords are "seconds relative to now" — negative = past, 0 = now)
        _win_lo = -_wave_window[0]
        if pt is not None and (pt - now) >= _win_lo:
            xp = pt - now                         # e.g. -30.0 s
            _p_lbl[i].set_x(xp); _p_lbl[i].set_visible(True)
            if st2 is not None and (st2 - now) >= _win_lo:
                xs = st2 - now                    # e.g. -22.0 s
                # Blue shade: P → S  (width = xs - xp > 0)
                _p_shd[i].set_x(xp);  _p_shd[i].set_width(xs - xp)
                _p_shd[i].set_visible(True)
                # Orange shade: S → now  (width = 0 - xs = -xs > 0)
                _s_shd[i].set_x(xs);  _s_shd[i].set_width(-xs)
                _s_shd[i].set_visible(True)
                _s_lbl[i].set_x(xs);  _s_lbl[i].set_visible(True)
            else:
                # Only P — blue shade from P arrival to now
                _p_shd[i].set_x(xp);  _p_shd[i].set_width(-xp)
                _p_shd[i].set_visible(True)
                _s_shd[i].set_visible(False)
                _s_lbl[i].set_visible(False)
        else:
            _p_shd[i].set_visible(False)
            _s_shd[i].set_visible(False)
            _p_lbl[i].set_visible(False)
            _s_lbl[i].set_visible(False)

        # Predicted-P line (dotted, shown up to 8 s before expected arrival)
        if pp is not None and pt is None:
            xv = pp - now
            if -DISPLAY_SEC <= xv <= 8:
                _ppln[i].set_xdata([xv, xv]); _ppln[i].set_visible(True)
            else:
                _ppln[i].set_visible(False)
        else:
            _ppln[i].set_visible(False)

        # P / S arrival time annotation
        if pt:
            pt_str = time.strftime("%H:%M:%S", time.gmtime(pt))
            if st2:
                st2_str = time.strftime("%H:%M:%S", time.gmtime(st2))
                sp_dt   = st2 - pt
                _sptx[i].set_text(
                    f"P {pt_str}  S {st2_str}  S−P={sp_dt:.1f}s  ≈{sp_dt/SP_FACTOR:.0f}km")
                _sptx[i].set_color(C["s_det"])
            elif sp2:
                sp2_str = time.strftime("%H:%M:%S", time.gmtime(sp2))
                _sptx[i].set_text(
                    f"P {pt_str}  pred-S {sp2_str}  (in {max(0.0, sp2-now):.0f}s)")
                _sptx[i].set_color(C["s_pred"])
            else:
                _sptx[i].set_text(f"P {pt_str}")
                _sptx[i].set_color(C["p"])
        else:
            _sptx[i].set_text("")

        # CFT plots (clipped to window + downsampled for performance)
        if cft is not None and len(cft) == len(t_arr):
            cft_w = cft[mask_win] if mask_win.sum() >= 2 else cft
            rel_cw = rel_w if mask_win.sum() >= 2 else rel
            if len(cft_w) > _WDISP_MAX:
                _step_c = len(cft_w) // _WDISP_MAX
                _cln[i].set_data(rel_cw[::_step_c], cft_w[::_step_c])
            else:
                _cln[i].set_data(rel_cw, cft_w)
            _cax[i].set_xlim(-win, 0)
            pk = max(float(cft.max()), TRIGGER_ON * 2.4)
            _ylc[i] = max(pk, _ylc[i]*.998 + TRIGGER_ON*2.4*.002)
            _cax[i].set_ylim(0, _ylc[i] * 1.05)
        if pcft is not None and len(pcft) == len(t_arr):
            pcft_w = pcft[mask_win] if mask_win.sum() >= 2 else pcft
            if len(pcft_w) > _WDISP_MAX:
                _step_p2 = len(pcft_w) // _WDISP_MAX
                _pkl_ln[i].set_data(rel_w[::_step_p2], pcft_w[::_step_p2])
            else:
                _pkl_ln[i].set_data(rel_w, pcft_w)

        # P/S ratio label
        psc = "#888888"
        if ps_r is not None:
            if   ps_r >= 2.0:        psc, pst = C["pkl"],      f"P/S={ps_r:.1f} P↑"
            elif ps_r <= PS_S_THRESH: psc, pst = C["s_det"],   f"P/S={ps_r:.1f} S↑"
            else:                     psc, pst = "#888888",    f"P/S={ps_r:.1f}"
            _psr_txt[i].set_text(pst); _psr_txt[i].set_color(psc)
        else:
            _psr_txt[i].set_text("")

        bg = C["bg_alert"] if alert else (C["bg_p"] if pt else C["bg"])
        _wax[i].set_facecolor(bg); _cax[i].set_facecolor(bg)
        if alert: any_alert = True
        if pt:    any_p = True; n_p_sta += 1

        # Table row
        if   alert: sta_tag, row_col = f"{sta}!", C["on"]
        elif pt:    sta_tag, row_col = f"{sta}p", C["p"]
        elif conn:  sta_tag, row_col = sta,        C["label"]
        else:       sta_tag, row_col = sta,        C["dim"]
        _tl[i].set_text(sta_tag);              _tl[i].set_color(row_col)
        _tr[i].set_text(f"{ratio:.2f}");       _tr[i].set_color(row_col)
        _tp[i].set_text(f"{pt-now:+.0f}s" if pt else "—")
        _td[i].set_text(f"{edist:.0f}km"   if edist else "—")
        if ps_r is not None:
            _tps[i].set_text(f"{ps_r:.2f}"); _tps[i].set_color(psc)
        else:
            _tps[i].set_text("—")

    # ── Deferred screenshot email flush ──────────────────────────────────────
    # Emails queued via _queue_email_with_shot() fire here once their delay
    # has elapsed.  Screenshot is captured now (main thread) so matplotlib is
    # safe, then the actual SMTP send runs in a daemon thread.
    _now_shot = time.time()
    with _pending_shot_emails_lock:
        _still_pend, _ready_now = [], []
        for _fat, _subj, _bdy in _pending_shot_emails:
            (_ready_now if _now_shot >= _fat else _still_pend).append(
                (_fat, _subj, _bdy))
        _pending_shot_emails[:] = _still_pend
    for _, _subj, _bdy in _ready_now:
        _shot = _capture_screenshot()
        threading.Thread(
            target=_send_email,
            args=(_subj, _bdy),
            kwargs={"img_bytes": _shot},
            daemon=True).start()
        _log("SCREENSHOT", f"deferred email sent with fresh screenshot: {_subj[:60]}")

    # ── On-demand screenshot email (ntfy "screenshot" command) ───────────────
    if _screenshot_email_requested[0]:
        _screenshot_email_requested[0] = False
        _od_shot  = _capture_screenshot()
        _od_ts    = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())
        _od_subj  = f"Live Waveforms — {time.strftime('%H:%M:%S UTC', time.gmtime())}"
        _od_body  = (
            f"ON-DEMAND WAVEFORM SCREENSHOT\n"
            f"{'='*46}\n"
            f"Captured : {_od_ts}\n"
            f"Requested via ntfy remote command.\n"
            f"\n"
            f"STATION STATUS\n"
            f"{'-'*46}\n"
            + "\n".join(
                f"  {s.label}: {'LIVE' if s.connected else 'DEAD'}"
                f"  P={'yes' if s.p_time else 'no'}"
                f"  ML={f'M{s.ml_est:+.1f}' if s.ml_est else '—'}"
                for s in states.values())
        )
        threading.Thread(
            target=_send_email,
            args=(_od_subj, _od_body),
            kwargs={"img_bytes": _od_shot},
            daemon=True).start()
        _send_ntfy(
            "SCREENSHOT CAPTURED",
            f"Live waveform screenshot emailed.\n"
            f"Captured: {_od_ts}\n"
            f"To: {EMAIL_TO if EMAIL_ENABLED else '(email not configured)'}",
            priority="default")
        _log("SCREENSHOT", f"On-demand screenshot emailed at {_od_ts}")

    # ── Epicenter snapshot ─────────────────────────────────────────────────
    with epicenter.lock:
        have_epi   = epicenter.lat is not None
        elat       = epicenter.lat;  elon  = epicenter.lon
        edepth     = epicenter.depth_km
        edist_home = epicenter.dist_home
        erms       = epicenter.rms_sec; en = epicenter.n_sta
        eaz_gap    = epicenter.az_gap
        et0        = epicenter.t_origin

    # ── Inter-station S-P consistency check ──────────────────────────────────
    # If multiple stations have S-picks, their implied distances must agree.
    # A spread > 3× the median implies at least one false S-pick; reduce
    # confidence to prevent confirming on a bad location.
    _sp_dists = []
    for _sst in states.values():
        if _sst.p_time is not None and _sst.s_time is not None:
            _sp = _sst.s_time - _sst.p_time
            if 1.0 < _sp <= SP_MAX_SEC:
                _sp_dists.append(_sp / SP_FACTOR)
    _sp_coherent = True
    if len(_sp_dists) >= 2:
        _sp_med = float(np.median(_sp_dists))
        _sp_max = max(_sp_dists); _sp_min = min(_sp_dists)
        if _sp_med > 0 and (_sp_max / _sp_med > 3.5 or _sp_min / _sp_med < 0.15):
            _sp_coherent = False
            _log("S-COHERENCE",
                 f"S-P distances incoherent: {[f'{d:.0f}km' for d in _sp_dists]}"
                 f"  spread={_sp_max/_sp_med:.1f}× — suppressing confirmation")
    confirmed = have_epi and en >= MIN_STA_CONFIRM and _sp_coherent

    # ── Earthquake event log (fires once when confirmed event ends) ────────
    if confirmed:
        # Keep snapshot fresh every frame while confirmed so we always have
        # the most recent data even if a global reset fires next frame.
        _ml_vals = [st.ml_est for st in states.values() if st.ml_est is not None]
        _med_ml  = sorted(_ml_vals)[len(_ml_vals)//2] if _ml_vals else None
        # Save timeline snapshot HERE (while confirmed) so it survives global reset
        _animate._last_confirmed_snap = dict(
            elat=elat, elon=elon, edepth=edepth, edist_home=edist_home,
            en=en, med_ml=_med_ml, erms=erms, eaz_gap=eaz_gap, et0=et0,
            eq_id=_eq_discord_prefix(),
            _timeline=_timeline_snapshot())   # ← saved inside snap

    if _animate._was_confirmed and not confirmed:
        snap      = _animate._last_confirmed_snap
        _med_ml   = snap.get("med_ml")
        _elat     = snap.get("elat"); _elon = snap.get("elon")
        _edepth   = snap.get("edepth"); _edist = snap.get("edist_home")
        _en       = snap.get("en", 0)
        _loc_str  = (f"{_elat:.4f}°N  {_elon:.4f}°W  depth≈{_edepth:.1f}km"
                     if _elat is not None else "unknown")
        _dist_str = f"  {_edist:.0f}km from San Ramon" if _edist else ""
        _ml_str   = f"ML={_med_ml:+.1f}" if _med_ml is not None else "ML=?"
        _pgv_end  = pgv_at_dist(_med_ml, _edist) if (_med_ml and _edist) else None
        _mmi_str  = mmi_label(pgv_to_mmi(_pgv_end))[0] if _pgv_end else "?"
        # Max observed MMI across stations
        _obs_pgvs = [counts_to_pgv(st.event_peak, st.sensitivity)
                     for st in states.values()
                     if st.event_peak > 0]
        _max_mmi_str = mmi_label(pgv_to_mmi(max(_obs_pgvs)))[0] if _obs_pgvs else "?"
        _city_end = _nearest_city(_elat, _elon)
        _city_end_str = f"  near {_city_end[0]}" if _city_end else ""
        _log("EARTHQUAKE",
             f"{_ml_str}  loc={_loc_str}{_dist_str}{_city_end_str}"
             f"  MaxMMI={_max_mmi_str}  MMI@SanRamon={_mmi_str}"
             f"  ({_en} stations)")
        # Use the timeline saved inside the snap (survives global reset)
        _saved_tl = snap.get("_timeline") or []
        _saved_tl.append((time.strftime("%H:%M:%S UTC", time.gmtime()),
                          f"Event ended — final {_ml_str}  {_loc_str}"))
        # Send final report (map + full timeline) in background thread
        _send_final_report(snap, _saved_tl)
        _timeline_clear()   # clear AFTER report is queued (not during global reset)
    _animate._was_confirmed = confirmed

    # ── Teleseism detection ────────────────────────────────────────────────
    # A teleseism shows nearly-simultaneous P arrivals at all stations
    # (planar wavefront from very far away) but produces absurdly large S-P
    # distances.  Detect it and suppress local-ML to avoid fake magnitudes.
    _p_times_active = [st.p_time for st in states.values()
                       if st.p_time is not None]
    if len(_p_times_active) >= TELESEISM_MIN_STA:
        _p_spread = max(_p_times_active) - min(_p_times_active)
        # All stations triggered within TELESEISM_P_SPREAD seconds AND
        # no valid local epicenter (teleseismic → epicenter far away)
        _local_epi = (have_epi and edist_home is not None
                      and edist_home < ML_MAX_DIST_KM)
        if _p_spread <= TELESEISM_P_SPREAD and not _local_epi:
            if not _teleseism_flag[0]:
                _teleseism_flag[0] = True
                _timeline_add(
                    f"TELESEISM detected: {len(_p_times_active)} sta triggered "
                    f"within {_p_spread:.1f}s (planar wavefront — distant event)")
                _log("TELESEISM",
                     f"{len(_p_times_active)} sta  P-spread={_p_spread:.1f}s"
                     f"  (ML calculation suppressed — not a local event)")
                _tele_title = "TELESEISM DETECTED — distant earthquake"
                _tele_body  = (
                    f"TELESEISM\n\n"
                    f"All {len(_p_times_active)} stations triggered within "
                    f"{_p_spread:.1f}s of each other.\n"
                    f"This is a DISTANT (teleseismic) earthquake, not local.\n"
                    f"Local ML formula does not apply.\n\n"
                    f"No shaking expected at {HOME_LABEL}.")
                _send_ntfy(_tele_title, _tele_body, priority="default")
                _send_discord(_tele_title, _tele_body)
        else:
            _teleseism_flag[0] = False
    else:
        _teleseism_flag[0] = False

    # ── Magnitude / MMI ────────────────────────────────────────────────────
    # Suppress ML display for teleseisms — formula is invalid at those distances
    ml_vals = ([] if _teleseism_flag[0] else
               [st.ml_est for st in states.values()
                if st.ml_est is not None and st.p_time is not None])
    med_ml  = None
    ml_unc  = None   # std-dev uncertainty across stations
    if ml_vals:
        _ml_arr = np.array(sorted(ml_vals), dtype=float)
        if len(_ml_arr) >= 4:
            # Trimmed mean: drop the highest and lowest reading
            _ml_trim = _ml_arr[1:-1]
        elif len(_ml_arr) == 3:
            # With 3 stations, drop whichever single value is furthest from median
            _med3 = _ml_arr[1]
            _diffs = np.abs(_ml_arr - _med3)
            _ml_trim = np.delete(_ml_arr, np.argmax(_diffs))
        else:
            _ml_trim = _ml_arr
        med_ml = round(float(np.mean(_ml_trim)), 1)
        if len(ml_vals) >= 2:
            ml_unc = round(float(np.std(_ml_arr)), 2)

        # Source ML with uncertainty
        unc_str = f" ±{ml_unc:.1f}" if ml_unc is not None else ""
        _ml_src_txt.set_text(f"M {med_ml:+.1f}{unc_str}")
        _ml_src_txt.set_color(C["on"])

        # Max observed MMI — best available estimate per detecting station:
        #   1st choice: ML estimate + actual hypocentral distance (physics-based)
        #   2nd choice: rolling 5-s PGV from the waveform stream
        #   3rd choice: peak count / sensitivity (P-wave onset only — least accurate)
        obs_pgvs = []
        for st in states.values():
            if st.p_time is None:
                continue
            _dist_ml = st.event_dist_km or st.dist_km or None
            if st.ml_est is not None and _dist_ml and _dist_ml > 0:
                _pgv_ml = pgv_at_dist(st.ml_est, _dist_ml)
                if _pgv_ml:
                    obs_pgvs.append(_pgv_ml); continue
            candidates = []
            if st.pgv_cm_s > 0:
                candidates.append(st.pgv_cm_s)
            if st.event_peak > 0:
                candidates.append(counts_to_pgv(st.event_peak, st.sensitivity))
            if candidates:
                obs_pgvs.append(max(candidates))
        max_pgv     = max(obs_pgvs) if obs_pgvs else None
        obs_mmi     = pgv_to_mmi(max_pgv)
        qstr, qcol  = mmi_label(obs_mmi)
        _mmi_quake_txt.set_text(qstr); _mmi_quake_txt.set_color(qcol)

        # MMI at San Ramon via attenuation — use hypocentral distance (3-D)
        # to correctly account for depth, especially for shallow nearby events.
        _eff_dep = max(edepth or 10.0, 1.0) if edepth is not None else 10.0
        _hypo_sr = (math.sqrt(edist_home**2 + _eff_dep**2)
                    if (have_epi and edist_home) else None)
        pgv_sr   = pgv_at_dist(med_ml, _hypo_sr) if _hypo_sr else None
        mmi_sr   = pgv_to_mmi(pgv_sr)
        lstr, lcol = mmi_label(mmi_sr)
        _mmi_local_txt.set_text(lstr); _mmi_local_txt.set_color(lcol)

        dstr    = f"  depth≈{edepth:.0f}km" if edepth is not None else ""
        pgv_str = f"{pgv_sr:.3f}cm/s" if pgv_sr else "need epi"
        _mag_note.set_text(
            f"({len(ml_vals)} sta{unc_str}){dstr}  PGV@SanRamon≈{pgv_str}")
        _mag_note.set_color(C["label"])
        ax_mag.set_facecolor("#1a0606" if obs_mmi and obs_mmi >= 5 else "#080f08")

        prev = _animate._last_ml
        ml_changed = (prev is None or abs(med_ml - prev) >= 0.1)
        if ml_changed:
            unc_log = f"  unc=±{ml_unc:.2f}" if ml_unc is not None else ""
            city_ml = _nearest_city(elat, elon) if have_epi else None
            city_ml_s = f"  near {city_ml[0]}" if city_ml else ""
            _log("MAGNITUDE",
                 f"ML={med_ml:+.1f}{unc_log}  MMI_obs={qstr}  MMI_SanRamon={lstr}"
                 + (f"  depth={edepth:.1f}km" if edepth else "")
                 + city_ml_s)
            _timeline_add(
                f"Magnitude updated: M{med_ml:+.1f}"
                f"{'  ±%.1f'%ml_unc if ml_unc else ''}  ({len(ml_vals)} sta)"
                + city_ml_s)
            _animate._last_ml = med_ml

        # ── Build shared station table (used by both ML and EQ emails) ──
        def _station_table():
            lines = []
            for _s in states.values():
                with _s.lock:
                    _pt = _s.p_time; _st2 = _s.s_time; _ml2 = _s.ml_est
                _pt_s = time.strftime("%H:%M:%S", time.gmtime(_pt)) if _pt else "—"
                _st_s = time.strftime("%H:%M:%S", time.gmtime(_st2)) if _st2 else "—"
                _sp_s = f"{_st2-_pt:.1f}s" if (_pt and _st2) else "—"
                _ml_s = f"ML{_ml2:+.1f}" if _ml2 else "—   "
                lines.append(
                    f"  {_s.label:<22} P={_pt_s}  S={_st_s}  S-P={_sp_s:<6}  {_ml_s}")
            return "\n".join(lines)

        # ── ML-update email (fires when magnitude changes by ≥0.1) ──────
        with _email_lock:
            _ml_ok = ml_changed and (now - _email_last_ml[0]) >= EMAIL_COOLDOWN_ML
            if _ml_ok:
                _email_last_ml[0] = now
        if _ml_ok:
            _dep_s   = f"{edepth:.1f} km"            if edepth is not None else "not yet estimated"
            _pgv_s   = f"{pgv_sr:.4f} cm/s"          if pgv_sr             else "need epicenter"
            _loc_s   = f"{elat:.4f}°N  {elon:.4f}°W" if elat is not None   else "not yet estimated"
            _dist_s  = f"{edist_home:.0f} km"         if edist_home         else "—"
            _orig_s  = (time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime(et0))
                        if et0 else "not yet estimated")
            _acc_lbl2, _, _pos_err2 = _epi_accuracy(en, erms, eaz_gap)
            _acc_s   = (f"{_acc_lbl2}  (±{_pos_err2:.0f} km)" if _pos_err2
                        else _acc_lbl2 if en else "—")
            _conf_s  = f"confirmed ({en} sta)" if confirmed else f"unconfirmed ({en} sta)"
            _city_ml = _nearest_city(elat, elon) if have_epi else None
            _city_ml_s = (f"{_city_ml[0]} ({_city_ml[1]:.0f} km away)"
                          if _city_ml else "not yet determined")
            _ml_city_lbl = _city_label(elat, elon) if have_epi else None
            _send_email(
                f"{_subj_icon(med_ml)}Magnitude Update — M{med_ml:+.1f}"
                + (f"  {_ml_city_lbl}" if _ml_city_lbl else ""),
                f"MAGNITUDE UPDATE\n"
                f"{'='*54}\n"
                f"Source ML          : M {med_ml:+.1f}"
                f"{'  ±'+str(ml_unc) if ml_unc else ''}\n"
                f"ML station spread  : {'±%.2f (N=%d)'%(ml_unc,len(ml_vals)) if ml_unc else 'single station'}\n"
                f"Max observed MMI   : {qstr}\n"
                f"MMI @ {HOME_LABEL:<13} : {lstr}\n"
                f"PGV @ {HOME_LABEL:<13} : {_pgv_s}\n"
                f"Stations w/ P      : {len(ml_vals)}\n"
                f"\n"
                f"LOCATION\n"
                f"{'-'*54}\n"
                f"Nearest city : {_city_ml_s}\n"
                f"Coordinates  : {_loc_s}\n"
                f"Depth        : {_dep_s}\n"
                f"Distance     : {_dist_s} from {HOME_LABEL}\n"
                f"Origin time  : {_orig_s}\n"
                f"Epicenter    : {_conf_s}\n"
                f"Accuracy     : {_acc_s}\n"
                f"RMS residual : {'%.2f s'%erms if erms else '—'}\n"
                f"Azimuth gap  : {'%.0f°'%eaz_gap if eaz_gap else '—'}\n"
                f"\n"
                f"STATION DETECTIONS\n"
                f"{'-'*54}\n"
                f"{_station_table()}\n")
            _ntfy_ml_loc = _ml_city_lbl or _city_ml_s
            _ntfy_ml_pri = "urgent" if med_ml >= 4.0 else "high"
            _ml_per_sta  = "\n".join(
                f"  {s.label:<22} M{s.ml_est:+.2f}  peak {int(s.event_peak):>10,} cts"
                f"  PGV {s.pgv_cm_s:.5f} cm/s"
                for s in states.values() if s.ml_est is not None)
            if NTFY_MIN_ML is None or med_ml >= NTFY_MIN_ML:
                _ml_title = f"MAGNITUDE UPDATE | M{med_ml:+.1f} | {_ntfy_ml_loc}"
                _ml_body  = (
                    f"---- MAGNITUDE ----\n"
                    f"Source ML  : M {med_ml:+.1f}"
                    f"{(' +/-' + str(ml_unc)) if ml_unc else ''}\n"
                    f"Stations   : {len(ml_vals)} P-arrivals\n"
                    f"Status     : {_conf_s}\n"
                    f"\n"
                    f"Per-station ML / waveform peaks:\n"
                    f"{_ml_per_sta if _ml_per_sta else '  none'}\n"
                    f"\n"
                    f"---- INTENSITY ----\n"
                    f"MMI max obs: {qstr}\n"
                    f"MMI local  : {lstr}  ({HOME_LABEL})\n"
                    f"PGV local  : {_pgv_s}  ({HOME_LABEL})\n"
                    f"\n"
                    f"---- LOCATION ----\n"
                    f"Nearest    : {_city_ml_s}\n"
                    f"Coordinates: {_loc_s}\n"
                    f"Depth      : {_dep_s}\n"
                    f"Dist home  : {_dist_s}\n"
                    f"Origin time: {_orig_s}\n"
                    f"\n"
                    f"---- QUALITY ----\n"
                    f"Accuracy   : {_acc_s}\n"
                    f"RMS resid  : {'%.3f s' % erms if erms else 'n/a'}\n"
                    f"Azimuth gap: {'%.1f deg' % eaz_gap if eaz_gap else 'n/a'}")
                _send_ntfy(_ml_title, _ml_body, priority=_ntfy_ml_pri)
                # Attach map if epicenter is known
                _ml_map = (_make_map_png(elat, elon, et0=et0, med_ml=med_ml)
                           if have_epi else None)
                _ml_upd_n, _ml_eq_id = _eq_next_update()
                _ml_dc_prefix = _eq_discord_prefix(f"Upd #{_ml_upd_n}")
                _ml_dc_title  = f"{_ml_dc_prefix} Magnitude Update | {_ntfy_ml_loc}"
                _ml_dc_bold   = f"M{med_ml:+.1f} | Origin: {_orig_s}"
                _send_discord(_ml_dc_title, _ml_body, img_bytes=_ml_map,
                              bold_header=_ml_dc_bold)

        # ── Confirmed-earthquake alert (epicenter known, ≥2 stations) ───
        # Require ML ≥ 1.5 to send any confirmed-event notification.
        # Only re-send when something actually changed:
        #   • first time this event becomes confirmed (station count was 0)
        #   • station count increased (more arrivals added)
        #   • epicenter moved > 5 km (solver revised the location)
        # ML changes are already handled by the magnitude-update block above —
        # no need to repeat the full confirmed message just because ML ticked.
        _CONFIRMED_ML_MIN = 1.5
        if confirmed and (med_ml is None or med_ml >= _CONFIRMED_ML_MIN):
            _eq_en_changed  = (en != _animate._last_eq_sent_en)
            _eq_epi_moved   = (
                _animate._last_eq_sent_elat is None or
                elat is None or
                haversine_km(elat, elon,
                             _animate._last_eq_sent_elat,
                             _animate._last_eq_sent_elon) > 5.0
            )
            _eq_something_changed = _eq_en_changed or _eq_epi_moved
            # Also honour the email cooldown so we never send faster than 30 s
            # even if multiple things update in quick succession.
            with _email_lock:
                _eq_ok = ((_eq_something_changed)
                          and (now - _email_last_eq[0]) >= 30.0)
                if _eq_ok:
                    _email_last_eq[0] = now
            if _eq_ok:
                _utc_o = (time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime(et0))
                          if et0 else "unknown")
                _acc_lbl, _, _pos_err = _epi_accuracy(en, erms, eaz_gap)
                _city_eq = _nearest_city(elat, elon)
                _city_eq_s = (f"{_city_eq[0]} ({_city_eq[1]:.0f} km away)"
                              if _city_eq else "unknown")
                _timeline_add(
                    f"Earthquake confirmed: M{med_ml:+.1f}  {en} stations"
                    f"  near {_city_eq[0] if _city_eq else '?'}")
                _eq_city_lbl = _city_label(elat, elon) or f"{edist_home:.0f} km from {HOME_LABEL}"
                # Defer screenshot so waveforms + S-waves have time to develop.
                _queue_email_with_shot(
                    f"{_subj_icon(med_ml)}Earthquake Confirmed — M{med_ml:+.1f}  {_eq_city_lbl}",
                    f"EARTHQUAKE CONFIRMED\n"
                    f"{'='*54}\n"
                    f"Origin time  : {_utc_o}\n"
                    f"Epicenter    : {_eq_city_lbl}\n"
                    f"Nearest city : {_city_eq_s}\n"
                    f"Coordinates  : {elat:.4f}°N  {elon:.4f}°W\n"
                    f"Depth        : {'%.1f km'%edepth if edepth is not None else '—'}\n"
                    f"Distance     : {'%.0f km'%edist_home if edist_home else '—'} from {HOME_LABEL}\n"
                    f"\n"
                    f"MAGNITUDE & INTENSITY\n"
                    f"{'-'*54}\n"
                    f"Source ML          : M {med_ml:+.1f}"
                    f"{'  ±'+str(ml_unc) if ml_unc else ''}\n"
                    f"Max observed MMI   : {qstr}\n"
                    f"MMI @ {HOME_LABEL:<13} : {lstr}\n"
                    f"PGV @ {HOME_LABEL:<13} : {'%.4f cm/s'%pgv_sr if pgv_sr else '—'}\n"
                    f"\n"
                    f"LOCATION QUALITY\n"
                    f"{'-'*54}\n"
                    f"Accuracy     : {_acc_lbl}  "
                    f"({'±%.0f km'%_pos_err if _pos_err else '—'})\n"
                    f"Stations     : {en}\n"
                    f"RMS residual : {'%.2f s'%erms if erms else '—'}\n"
                    f"Azimuth gap  : {'%.0f°'%eaz_gap if eaz_gap else '—'}\n"
                    f"\n"
                    f"STATION DETECTIONS\n"
                    f"{'-'*54}\n"
                    f"{_station_table()}\n")
                if med_ml >= 5.0:
                    _ntfy_eq_pri  = "urgent"
                    _ntfy_eq_tags = "warning,rotating_light,sos"
                elif med_ml >= 3.0:
                    _ntfy_eq_pri  = "high"
                    _ntfy_eq_tags = "warning,rotating_light"
                else:
                    _ntfy_eq_pri  = "default"
                    _ntfy_eq_tags = "seismograph"
                _sta_detail   = "\n".join(
                    f"  {s.label:<22}"
                    f"  P={time.strftime('%H:%M:%S UTC', time.gmtime(s.p_time)) if s.p_time else 'no P':>15}"
                    f"  ML={('M%+.2f' % s.ml_est) if s.ml_est is not None else '  n/a':>7}"
                    f"  STA/LTA={s.last_p_ratio:.1f}"
                    f"  [{_sta_lta_level(s.last_p_ratio, P_THRESH * s.thresh_mult)}]"
                    f"  peak={int(s.event_peak):>10,} cts"
                    f"  PGV={s.pgv_cm_s:.5f} cm/s"
                    f"  PGA={s.pga_cm_s2:.5f} cm/s2"
                    for s in states.values()
                    if s.p_time is not None or s.event_peak)
                _eq_level = max(
                    (_sta_lta_level(s.last_p_ratio, P_THRESH * s.thresh_mult)
                     for s in states.values()
                     if s.p_time is not None or s.event_peak),
                    key=lambda l: ["LOW", "MEDIUM", "HIGH", "EXTREME"].index(l),
                    default="LOW")
                # Alert radius filter
                _radius_ok = True
                if ALERT_RADIUS_KM is not None and edist_home is not None:
                    if edist_home > ALERT_RADIUS_KM:
                        _log("ALERT FILTER",
                             f"Epicenter {edist_home:.0f}km > radius {ALERT_RADIUS_KM}km — suppressed")
                        _radius_ok = False

                # Aftershock cooldown
                _now_eq = time.time()
                _is_aftershock = (
                    _last_eq_confirmed_ml[0] is not None
                    and (_now_eq - _last_eq_confirmed_t[0]) < AFTERSHOCK_COOLDOWN_SEC
                    and med_ml is not None
                    and med_ml < _last_eq_confirmed_ml[0] - 0.5
                )

                if NTFY_MIN_ML is None or med_ml >= NTFY_MIN_ML:
                    _eq_title = f"[{_eq_level}] EARTHQUAKE CONFIRMED | M{med_ml:+.1f} | {_eq_city_lbl}"
                    _eq_body  = (
                        f"---- MAGNITUDE ----\n"
                        f"Source ML  : M {med_ml:+.1f}"
                        f"{(' +/-' + str(ml_unc)) if ml_unc else ''}\n"
                        f"Stations   : {en} confirmed\n"
                        f"\n"
                        f"---- LOCATION ----\n"
                        f"Epicenter  : {_eq_city_lbl}\n"
                        f"Nearest    : {_city_eq_s}\n"
                        f"Coordinates: {elat:.4f} N  {elon:.4f} W\n"
                        f"Depth      : {'%.2f km' % edepth if edepth is not None else 'n/a'}\n"
                        f"Dist home  : {'%.1f km' % edist_home if edist_home else 'n/a'}\n"
                        f"Origin time: {_utc_o}\n"
                        f"\n"
                        f"---- INTENSITY ----\n"
                        f"MMI max obs: {qstr}\n"
                        f"MMI local  : {lstr}  ({HOME_LABEL})\n"
                        f"PGV local  : {'%.5f cm/s' % pgv_sr if pgv_sr else 'n/a'}\n"
                        f"\n"
                        f"---- QUALITY ----\n"
                        f"Quality    : {_event_quality_score(en, erms, eaz_gap, ml_unc)}/100\n"
                        f"Accuracy   : {_acc_lbl}  ({'+-%.0f km' % _pos_err if _pos_err else 'n/a'})\n"
                        f"RMS resid  : {'%.3f s' % erms if erms else 'n/a'}\n"
                        f"Azimuth gap: {'%.1f deg' % eaz_gap if eaz_gap else 'n/a'}\n"
                        f"\n"
                        f"---- DETECTION STRENGTH ----\n"
                        + "\n".join(
                            f"  {s.label:<22}  STA/LTA={s.last_p_ratio:.1f}"
                            f"  [{_sta_lta_level(s.last_p_ratio, P_THRESH * s.thresh_mult)}]"
                            f"  ({s.last_p_ratio / (P_THRESH * s.thresh_mult):.1f}× thr)"
                            for s in states.values()
                            if s.p_time is not None or s.event_peak) +
                        f"\n\n"
                        f"---- PER-STATION DETAIL ----\n"
                        f"{_sta_detail if _sta_detail else '  none'}")
                    if _radius_ok and not _is_aftershock:
                        _last_eq_confirmed_t[0]  = _now_eq
                        _last_eq_confirmed_ml[0] = med_ml
                        # Record what we sent so next frame can diff against it
                        _animate._last_eq_sent_en   = en
                        _animate._last_eq_sent_elat = elat
                        _animate._last_eq_sent_elon = elon
                        _send_ntfy(_eq_title, _eq_body, priority=_ntfy_eq_pri)
                        # Generate epicenter map and attach to Discord alert
                        _eq_map = _make_map_png(elat, elon, et0=et0, med_ml=med_ml)
                        _eq_upd_n, _eq_uid = _eq_next_update()
                        _eq_dc_prefix = _eq_discord_prefix(f"Upd #{_eq_upd_n}")
                        _eq_dc_title  = (f"{_eq_dc_prefix} ✅ Earthquake Confirmed"
                                         f" | {_eq_city_lbl}")
                        _eq_dc_bold   = (f"M{med_ml:+.1f}"
                                         f" | {_eq_city_lbl}"
                                         f" | Origin: {_utc_o}")
                        # Gate Discord on ML uncertainty — if spread is too high,
                        # append a caution note rather than blocking the alert
                        if ml_unc is not None and ml_unc > 0.5:
                            _eq_dc_bold += f"  ⚠️ HIGH UNCERTAINTY (±{ml_unc:.1f})"
                            _eq_body += (f"\n\n⚠️ WARNING: ML uncertainty = ±{ml_unc:.1f}"
                                         f" (spread across stations is large)."
                                         f" Magnitude estimate may be unreliable.")
                        _send_discord(_eq_dc_title, _eq_body, level=_eq_level,
                                      img_bytes=_eq_map, bold_header=_eq_dc_bold,
                                      everyone=True)  # @everyone on all confirmed EQs

                        # ── Safety card (M ≥ 3.0) ─────────────────────────────
                        # Send a second Discord message with the Drop/Cover/Hold On
                        # image and actionable instructions for felt earthquakes.
                        if (med_ml is not None and med_ml >= DISCORD_SAFETY_ML
                                and _os.path.exists(_SAFETY_CARD_PATH)):
                            try:
                                _sc_bytes = open(_SAFETY_CARD_PATH, "rb").read()
                            except OSError:
                                _sc_bytes = None
                            if _sc_bytes:
                                _sc_mmi   = lstr if lstr else "Unknown"
                                _sc_body  = (
                                    f"---- WHAT TO DO NOW ----\n"
                                    f"A M{med_ml:+.1f} earthquake was just detected"
                                    f" near {_eq_city_lbl}.\n"
                                    f"Expected shaking at San Ramon: {_sc_mmi}\n"
                                    f"\n"
                                    f"  1. DROP  — get down on hands and knees\n"
                                    f"             (protects you from being knocked over)\n"
                                    f"\n"
                                    f"  2. COVER — get under a sturdy desk or table,\n"
                                    f"             or cover your head/neck with your arms\n"
                                    f"             Stay away from windows and heavy objects\n"
                                    f"\n"
                                    f"  3. HOLD ON — hold on until shaking stops\n"
                                    f"               Do NOT run outside during shaking\n"
                                    f"\n"
                                    f"---- AFTER SHAKING STOPS ----\n"
                                    f"  • Check for injuries and hazards (gas leaks, fire)\n"
                                    f"  • Expect aftershocks — be prepared to DROP again\n"
                                    f"  • Stay away from damaged buildings and power lines\n"
                                    f"  • Text, don't call — keeps phone lines clear for emergency services\n"
                                    f"  • Tune to local emergency broadcasts for updates\n"
                                    f"\n"
                                    f"More info: https://myearthquake.dpdns.org/")
                                _send_discord(
                                    f"⚠️ Safety Instructions | M{med_ml:+.1f} | {_eq_city_lbl}",
                                    _sc_body,
                                    level=_eq_level,
                                    img_bytes=_sc_bytes,
                                    bold_header="DROP!  COVER!  HOLD ON!")

                        _play_alarm(_eq_level)
                        # Major EQ spam — rapid-fire @everyone pings
                        if med_ml is not None:
                            _discord_major_spam(med_ml, _eq_city_lbl, _eq_level)
                    _ws_broadcast_event("EQ Confirmed", level=_eq_level,
                        location=_eq_city_lbl, ml=med_ml,
                        detail=f"{en} stations  depth={edepth:.0f}km" if edepth else f"{en} stations")
                    # Spawn USGS cross-check ONCE per event (60 s later)
                    if (elat is not None and elon is not None and et0 is not None
                            and not _usgs_crosscheck_spawned[0]):
                        _usgs_crosscheck_spawned[0] = True
                        threading.Thread(
                            target=_usgs_crosscheck,
                            args=(elat, elon, et0, med_ml),
                            daemon=True).start()
                    # Catalog logging + aftershock forecast (catalog updates each time;
                    # forecast sent ONCE per event to prevent repeat spam)
                    if elat is not None and elon is not None:
                        _catalog_log_event(
                            origin_time=et0, lat=elat, lon=elon, depth_km=edepth,
                            ml=med_ml, ml_unc=ml_unc, n_stations=en,
                            dist_home=edist_home, location=_eq_city_lbl,
                            event_type="EQ", rms_sec=erms, confirmed=True)
                        _mainshock_time[0] = et0
                        _mainshock_ml[0]   = med_ml
                        # Aftershock probability — sent once per event only
                        if not _as_forecast_sent[0]:
                            _as_prob = _omori_utsu_prob(med_ml)
                            if _as_prob is not None:
                                _as_forecast_sent[0] = True
                                _log("AFTERSHOCK",
                                     f"P(M≥{med_ml-1:.1f} in {AFTERSHOCK_WINDOW_HOURS}h)"
                                     f" = {_as_prob*100:.0f}%")
                                _send_discord(
                                    f"Aftershock Forecast | M{med_ml:+.1f} | {_eq_city_lbl}",
                                    f"Omori-Utsu forecast ({AFTERSHOCK_WINDOW_HOURS}h window)\n"
                                    f"P(M≥{med_ml-1:.1f} aftershock in next {AFTERSHOCK_WINDOW_HOURS}h)"
                                    f" = {_as_prob*100:.0f}%\n"
                                    f"(Reasenberg & Jones 1989  —  sent once per event)")
    else:
        _ml_src_txt.set_text("M —");     _ml_src_txt.set_color(C["label"])
        _mmi_quake_txt.set_text("—");    _mmi_quake_txt.set_color(C["label"])
        _mmi_local_txt.set_text("—");    _mmi_local_txt.set_color(C["label"])
        _mag_note.set_text("± approx; requires P-wave detection")
        _mag_note.set_color(C["label"])
        ax_mag.set_facecolor("#080f08")
        _animate._last_ml = None

    # ── P/S-wave countdown to San Ramon ───────────────────────────────────
    if have_epi and et0 is not None:
        h3d      = math.sqrt((edist_home or 0.0)**2 + (edepth or 10.0)**2)
        t_p_home = et0 + h3d / VP
        t_s_home = et0 + h3d / VS
        p_secs   = t_p_home - now
        s_secs   = t_s_home - now
        if p_secs > 0:
            p_str = f"P in {p_secs:.1f}s"
        else:
            p_str = f"P {abs(p_secs):.0f}s ago"
        if s_secs > 0:
            s_str = f"S in {s_secs:.1f}s"
        else:
            s_str = f"S {abs(s_secs):.0f}s ago"
        _countdown_txt.set_text(f"!! San Ramon:  {p_str}   {s_str}")

        # ── EARLY EARTHQUAKE WARNING (fires once per event, M >= NTFY_EEW_ML) ─
        if (not _animate._eew_sent
                and med_ml is not None
                and med_ml >= NTFY_EEW_ML
                and s_secs > -30):        # still relevant (S not yet arrived >30s ago)
            _animate._eew_sent = True
            _eew_city  = _city_label(elat, elon) or f"{edist_home:.0f} km from {HOME_LABEL}"
            _eew_mmi   = pgv_to_mmi(pgv_at_dist(med_ml, edist_home) if edist_home else None)
            _eew_mmistr, _ = mmi_label(_eew_mmi)
            _p_utc     = time.strftime("%H:%M:%S UTC", time.gmtime(t_p_home))
            _s_utc     = time.strftime("%H:%M:%S UTC", time.gmtime(t_s_home))
            _eew_pri   = "urgent" if med_ml >= 4.0 else "high"
            _log("EEW", f"M{med_ml:+.1f}  {_eew_city}  P_home={p_secs:.0f}s  S_home={s_secs:.0f}s")
            _eew_title = f"*** EARLY EARTHQUAKE WARNING *** M{med_ml:+.1f} | {_eew_city}"
            _eew_body  = (
                f"**** EARLY EARTHQUAKE WARNING ****\n"
                f"\n"
                f"MAGNITUDE : M {med_ml:+.1f}"
                f"{(' +/-%.1f' % ml_unc) if ml_unc else ''}\n"
                f"LOCATION  : {_eew_city}\n"
                f"DEPTH     : {'%.1f km' % edepth if edepth is not None else 'n/a'}\n"
                f"DIST HOME : {'%.1f km' % edist_home if edist_home else 'n/a'}\n"
                f"\n"
                f"ARRIVAL TIMES AT {HOME_LABEL.upper()}\n"
                f"  P-wave : {_p_utc}  ({('in %.0fs' % p_secs) if p_secs > 0 else ('%.0fs ago' % abs(p_secs))})\n"
                f"  S-wave : {_s_utc}  ({('in %.0fs' % s_secs) if s_secs > 0 else ('%.0fs ago' % abs(s_secs))})\n"
                f"\n"
                f"PREDICTED SHAKING AT {HOME_LABEL.upper()}\n"
                f"  MMI    : {_eew_mmistr}\n"
                f"  PGV    : {'%.4f cm/s' % pgv_at_dist(med_ml, edist_home) if edist_home else 'n/a'}\n"
                f"\n"
                f"Stations: {en}  RMS: {'%.2fs' % erms if erms else 'n/a'}"
                f"  Az gap: {'%.0fdeg' % eaz_gap if eaz_gap else 'n/a'}")
            _eew_radius_ok = (ALERT_RADIUS_KM is None or edist_home is None
                              or edist_home <= ALERT_RADIUS_KM)
            if _eew_radius_ok:
                _send_ntfy(_eew_title, _eew_body, priority=_eew_pri)
                _eew_map = _make_map_png(elat, elon, et0=et0, med_ml=med_ml)
                _eew_upd_n, _ = _eq_next_update()
                _eew_dc_prefix = _eq_discord_prefix(f"Upd #{_eew_upd_n}")
                _eew_orig_s = (time.strftime("%H:%M:%S UTC", time.gmtime(et0))
                               if et0 else "unknown")
                _eew_dc_title = (f"{_eew_dc_prefix} ⚠️ EARLY WARNING"
                                 f" | M{med_ml:+.1f} | {_eew_city}")
                _eew_dc_bold  = (f"M{med_ml:+.1f} | {_eew_city}"
                                 f" | P in {p_secs:.0f}s | Origin: {_eew_orig_s}")
                _send_discord(_eew_dc_title, _eew_body, level="EEW",
                              img_bytes=_eew_map, bold_header=_eew_dc_bold,
                              everyone=True)  # EEW always @everyone
            _ws_broadcast_event("EEW", level="EXTREME",
                location=_eew_city, ml=med_ml,
                detail=f"S-wave in {s_secs:.0f}s" if s_secs > 0 else "S-wave arrived")
        if s_secs > 0:
            col = "#ff4444" if p_secs < 10 else C["home"]
        else:
            col = C["label"]
        _countdown_txt.set_color(col)
    else:
        _countdown_txt.set_text("")

    # ── Epicenter table text ───────────────────────────────────────────────
    if have_epi:
        dstr = f"  depth≈{edepth:.0f}km" if edepth is not None else ""
        rstr = f"  RMS={erms:.2f}s"       if erms is not None   else ""
        conf_str = "confirmed" if confirmed else "suspected"
        _epi_txt.set_text(
            f"Epicenter: {elat:.3f}°N  {elon:.3f}°W{dstr}\n"
            f"≈{edist_home:.0f}km from San Ramon{rstr}\n"
            f"({conf_str}, {en} station{'s' if en > 1 else ''})")
        _epi_txt.set_color(C["s_pred"] if confirmed else "#ffff88")

        acc_label, acc_color, pos_err = _epi_accuracy(en, erms, eaz_gap)
        pos_str = f"  ±{pos_err:.0f}km" if pos_err is not None else ""
        gap_str = f"  az.gap={eaz_gap:.0f}°" if eaz_gap is not None else ""
        _acc_txt.set_text(f"Loc. accuracy: {acc_label}{pos_str}{gap_str}")
        _acc_txt.set_color(acc_color)
    else:
        _epi_txt.set_text("Epicenter: awaiting P-arrivals …")
        _epi_txt.set_color(C["dim"])
        _acc_txt.set_text("")

    # ── Banner ─────────────────────────────────────────────────────────────
    utc    = time.strftime("%Y-%m-%d  %H:%M:%S UTC", time.gmtime(now))
    conn_n = sum(1 for s in states.values() if s.connected)
    if any_alert:
        if _teleseism_flag[0]:
            ax_ban.set_facecolor("#061828")
            _ban_main.set_text("TELESEISM — DISTANT EARTHQUAKE (no local shaking)")
            _ban_main.set_color("#44aaff")
            _p_tms = [st.p_time for st in states.values() if st.p_time is not None]
            _spd   = (max(_p_tms) - min(_p_tms)) if len(_p_tms) >= 2 else 0
            _ban_sub.set_text(
                f"{len(_p_tms)} stations triggered  |  P-spread {_spd:.1f}s"
                f"  |  No local ML  |  {utc}")
            _ban_sub.set_color("#88ccff")
        else:
            ax_ban.set_facecolor("#280606")
            _ban_main.set_text("***   EARTHQUAKE DETECTED   ***")
            _ban_main.set_color(C["on"])
            depth_tag = (f"  depth≈{edepth:.0f}km"
                         if (have_epi and edepth is not None) else "")
            if confirmed and edist_home is not None:
                rms_tag = f"  RMS {erms:.2f}s" if erms else ""
                acc_label, _ac, _pe = _epi_accuracy(en, erms, eaz_gap)
                sub = (f"Epicenter ≈{elat:.3f}°N {elon:.3f}°W{depth_tag}"
                       f"  — {edist_home:.0f}km from San Ramon"
                       f"  ({en} sta{rms_tag}  acc={acc_label})")
            else:
                sub = "Locating … need ≥2 confirmed P-arrivals"
            _q_score = _event_quality_score(en, erms, eaz_gap, ml_unc,
                                             sp_coherent=_sp_coherent)
            _q_color = ("#44ff44" if _q_score >= 70
                        else "#ffcc00" if _q_score >= 40
                        else "#ff6644")
            _ban_sub.set_text(f"{sub}    {utc}    Quality: {_q_score}%")
            _ban_sub.set_color("#ffcc88")
    elif any_p:
        ax_ban.set_facecolor("#1a1a06")
        _ban_main.set_text("P-WAVE DETECTED  —  Analyzing …")
        _ban_main.set_color("#ffff44")
        _ban_sub.set_text(
            f"{n_p_sta}/{N} P-arrivals  |  {conn_n}/{N} live    {utc}")
        _ban_sub.set_color(C["label"])
    else:
        ax_ban.set_facecolor("#0a260a")
        _ban_main.set_text("MONITORING  —  NO SEISMIC EVENT DETECTED")
        _ban_main.set_color(C["off"])
        _ban_sub.set_text(f"{conn_n}/{N} stations live    {utc}")
        _ban_sub.set_color(C["label"])

    # ── Map rings + epicenter ──────────────────────────────────────────────
    for a in _dyn_art:
        try: a.remove()
        except Exception: pass
    _dyn_art.clear()

    # Per-station P+S confirmed rings — drawn always, regardless of epicenter
    # Only when both P and S are confirmed so event_dist_km is reliable
    for st in states.values():
        with st.lock:
            ed2  = st.event_dist_km
            sla, slo = st.lat, st.lon
            has_ps   = st.p_time is not None and st.s_time is not None
        if has_ps and ed2 and ed2 > 0:
            xs, ys = _ring_xy(sla, slo, ed2)
            r, = ax_map.plot(xs, ys, color=C["ring_sus"], lw=1.4,
                             ls="--", alpha=0.80, zorder=15)
            _dyn_art.append(r)

    if have_epi:
        ex, ey = _proj(elat, elon)
        _epi_mk.set_data([ex], [ey]); _epi_mk.set_visible(True)
        _epi_lbl.set_position((ex, ey))
        _epi_lbl.set_text(
            f"  M{med_ml:+.1f}" if med_ml is not None else "  EPI")
        _epi_lbl.set_visible(True)
        _epi_home.set_data([ex, _hx], [ey, _hy])

        # Uncertainty ring (solid red — positional error envelope)
        unc = max(15.0, (erms or 0.0) * VP)
        xs, ys = _ring_xy(elat, elon, unc)
        r, = ax_map.plot(xs, ys, color=C["ring_unc"], lw=2.2,
                         ls="-", alpha=0.85, zorder=17)
        _dyn_art.append(r)

        # Expanding P-wavefront and S-wavefront rings
        if et0 is not None:
            p_dist = VP * (now - et0)
            s_dist = VS * (now - et0)
            for wdist, wcol, wlw, wls, walpha in [
                    (p_dist, C["ring_p"], 1.6, "--", 0.50),
                    (s_dist, C["ring_s"], 1.6, "--", 0.50)]:
                if 1 < wdist < 1800:
                    xs, ys = _ring_xy(elat, elon, wdist)
                    r, = ax_map.plot(xs, ys, color=wcol, lw=wlw,
                                     ls=wls, alpha=walpha, zorder=16)
                    _dyn_art.append(r)

    else:
        _epi_mk.set_visible(False)
        _epi_lbl.set_visible(False)
        _epi_home.set_data([], [])


def _stats_broadcast_loop():
    while True:
        time.sleep(60)
        b, rate, n = _bvalue_and_rate(hours=168)
        events = _catalog_recent(hours=168)
        mls = [e["ml"] for e in events if e["ml"] is not None]
        max_ml = max(mls) if mls else None
        as_prob = _omori_utsu_prob(_mainshock_ml[0]) if _mainshock_ml[0] else None
        _ws_broadcast_data({
            "type":   "stats",
            "n7d":    n,
            "rate":   rate,
            "bval":   b,
            "maxml":  max_ml,
            "asprob": as_prob,
        })


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
_animate._force_reset_requested   = False        # remote ntfy "reset" command
_animate._last_ml                 = None         # ML change tracking
_animate._was_confirmed           = False        # event→no-event transition
_animate._last_auto_reset         = time.time()
_animate._last_confirmed_snap     = {}
_animate._eew_sent                = False        # early warning fired for current event
_animate._auto_reconnect_pending  = False        # sleep-recovery reconnect in progress
_animate._last_seedlink_spawn     = time.time()  # when last SeedLink thread was started
_animate._last_auto_rc_time       = 0.0          # last auto-reconnect timestamp
_animate._confirmed_since         = time.time()  # when current confirmed event started
# ── EQ-confirmed change tracking — only re-send when something actually changed
_animate._last_eq_sent_en         = 0            # station count at last confirmed send
_animate._last_eq_sent_elat       = None         # epicenter at last confirmed send
_animate._last_eq_sent_elon       = None

AUTO_RESET_SEC = 3600  # forced global reset every 60 min (matches display window)

def _global_reset(reason="auto"):
    """Clear all P/S detections, epicenter estimate, and email thread state."""
    for st in states.values():
        with st.lock:
            st.p_time = st.s_time = st.s_predicted = st.p_predicted = None
            st.event_dist_km = None
            st.p_cleared = False
            st.event_peak = 0.0
            st.ml_est = None
            st.alert = False
    epicenter.reset()
    _reset_event_thread()
    _eq_reset_id()
    _animate._eew_sent           = False
    _teleseism_flag[0]           = False
    _as_forecast_sent[0]         = False
    _usgs_crosscheck_spawned[0]  = False
    _animate._last_eq_sent_en    = 0
    _animate._last_eq_sent_elat  = None
    _animate._last_eq_sent_elon  = None
    _sta_lta_everyone_last_t[0]  = 0.0
    _log("GLOBAL RESET", f"reason={reason}  all P/S times, epicenter, and email thread cleared")

if __name__ == "__main__":
    _log("SESSION",
         f"Monitor started  —  {N} stations"
         f"  server={SEEDLINK_HOST}:{SEEDLINK_PORT}"
         f"  tiles={'OSM' if HAS_TILES else 'static'}"
         f"  log={LOG_FILE}")

    print("═" * 72)
    print("  Seismic Monitor v4  ·  OSM Map  ·  P/S Rings  ·  ML/MMI/PGV/PGA")
    print(f"  Server    : {SEEDLINK_HOST}:{SEEDLINK_PORT}")
    print(f"  Thresholds: Alert≥{TRIGGER_ON}  P≥{P_THRESH}  (×thresh_mult per station)")
    print(f"  Stations  : {N}")
    _map_str = 'OSM tiles (contextily) — scroll/drag to zoom+pan' if HAS_TILES else 'static CA polygon (install contextily)'
    print(f"  Map       : {_map_str}")
    print(f"  Log       : {LOG_FILE}")
    print(f"\n  Commands: [s]tatus   [t]est earthquake   [q]uit")
    print("═" * 72 + "\n")

    # Prevent macOS from sleeping while the monitor runs
    if platform.system() == "Darwin":
        try:
            subprocess.Popen(["caffeinate", "-i", "-w", str(os.getpid())])
            _log("SYSTEM", "caffeinate active — Mac will not sleep while monitor runs")
            print("[SYSTEM] caffeinate active — Mac will not sleep while monitor runs")
        except Exception as _caf_exc:
            _log("SYSTEM", f"caffeinate failed: {_caf_exc}")

    _init_catalog_db()
    threading.Thread(target=_run_seedlink,        daemon=True).start()
    threading.Thread(target=_run_ncedc_seedlink,  daemon=True).start()
    threading.Thread(target=_stats_broadcast_loop, daemon=True).start()
    threading.Thread(target=_cmd,                 daemon=True).start()
    _start_ntfy_listener()
    threading.Thread(target=_ntfy_watchdog,       daemon=True).start()
    threading.Thread(target=_load_phasenet,       daemon=True).start()
    threading.Thread(target=_run_ws_server,       daemon=True).start()
    threading.Thread(target=_run_http_server,     daemon=True).start()
    _send_discord(
        "🟢 Quake Monitor Online",
        f"Seismic monitor started.\n"
        f"Stations : {N}\n"
        f"Server   : {SEEDLINK_HOST}:{SEEDLINK_PORT}\n"
        f"Discord alerts are active.")
    ani = animation.FuncAnimation(fig, _animate, interval=200,
                                  blit=False, cache_frame_data=False)
    plt.show()