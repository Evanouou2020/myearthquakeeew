#!/usr/bin/env python3
"""
SeismoPhone Computer Server
────────────────────────────
Run this on your computer, then in the phone app tap ⇆ and enter:
    ws://<your-computer-ip>:8765

Dashboard opens automatically at http://localhost:8080

Requires: pip install websockets
"""

import asyncio
import json
import socket
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer

WS_PORT   = 8765
HTTP_PORT = 8080

# ── Shared state ──────────────────────────────────────────────────────────
viewers: set = set()
latest:  dict = {}
phone_connected = False

# ── Embedded dashboard HTML ───────────────────────────────────────────────
DASHBOARD = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SeismoPhone — Computer Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/qrcodejs@1.0.0/qrcode.min.js"></script>
<style>
  *{box-sizing:border-box;margin:0;padding:0}
  :root{
    --bg:#080812;--surface:#10101e;--border:#1e1e3a;
    --accent:#44aaff;--green:#2ecc71;--red:#e74c3c;--amber:#f39c12;
    --text:#c8d0e8;--dim:#5a607a;
  }
  html,body{background:var(--bg);color:var(--text);font-family:'SF Mono','Fira Code','Courier New',monospace;height:100vh;overflow:hidden}
  header{display:flex;align-items:center;justify-content:space-between;padding:10px 20px;background:var(--surface);border-bottom:1px solid var(--border)}
  .logo{font-size:1.1rem;font-weight:700;color:var(--accent);letter-spacing:.04em}
  .logo small{color:var(--dim);font-weight:400;font-size:.68rem;margin-left:8px}
  .hdr-right{display:flex;align-items:center;gap:12px;font-size:.72rem;color:var(--dim)}
  #phone-dot{width:9px;height:9px;border-radius:50%;background:var(--dim);transition:background .3s,box-shadow .3s}
  #phone-dot.on{background:var(--green);box-shadow:0 0 8px var(--green);animation:blink 1.5s infinite}
  @keyframes blink{0%,100%{opacity:1}50%{opacity:.4}}

  .layout{display:grid;grid-template-columns:1fr 300px;grid-template-rows:1fr auto;height:calc(100vh - 48px);gap:0}
  .chart-col{display:flex;flex-direction:column;border-right:1px solid var(--border)}
  .chart-wrap{flex:1;min-height:0;padding:10px 8px 0;background:var(--surface)}
  .legend{display:flex;gap:14px;padding:4px 14px;font-size:.65rem;color:var(--dim)}
  .legend span{display:flex;align-items:center;gap:5px}
  .dot{width:8px;height:8px;border-radius:50%}

  .metrics-row{display:grid;grid-template-columns:1fr 1fr 1fr 1fr;border-top:1px solid var(--border)}
  .metric{display:flex;flex-direction:column;align-items:center;justify-content:center;padding:10px 4px 8px;border-right:1px solid var(--border);gap:2px}
  .metric:last-child{border-right:none}
  .metric-label{font-size:.58rem;color:var(--dim);text-transform:uppercase;letter-spacing:.09em}
  .metric-value{font-size:1.5rem;font-weight:700;color:var(--accent);transition:color .25s}
  .metric-sub{font-size:.65rem;color:var(--dim)}
  .metric-unit{font-size:.52rem;color:#2a2a4a}
  .metric-value.warn{color:var(--amber)}
  .metric-value.high{color:var(--red)}

  .side-col{display:flex;flex-direction:column;overflow:hidden}
  .side-section{padding:14px 16px;border-bottom:1px solid var(--border)}
  .side-title{font-size:.65rem;text-transform:uppercase;letter-spacing:.09em;color:var(--dim);margin-bottom:10px}

  /* QR code */
  #qrcode{display:flex;justify-content:center;margin:8px 0}
  #qrcode canvas,#qrcode img{border-radius:6px}
  .qr-url{font-size:.62rem;color:var(--dim);text-align:center;word-break:break-all;margin-top:4px}

  /* Axis bars */
  .axis-row{display:flex;align-items:center;gap:8px;font-size:.65rem;margin-bottom:5px}
  .axis-name{width:14px;color:var(--dim)}
  .axis-track{flex:1;height:6px;background:var(--border);border-radius:3px;overflow:hidden;position:relative}
  .axis-bar{position:absolute;top:0;bottom:0;border-radius:3px;transition:width .07s,left .07s}
  .axis-val{width:52px;text-align:right;font-size:.62rem;color:var(--dim)}

  /* Scale */
  .scale-row{display:flex;gap:2px;padding:0 0 2px}
  .scale-seg{flex:1;height:5px;border-radius:2px}
  .scale-lbl{display:flex;justify-content:space-between;font-size:.52rem;color:var(--dim)}

  /* Event log */
  .evt-log{flex:1;overflow-y:auto;padding:0}
  .evt-item{padding:5px 14px;border-bottom:1px solid var(--border);font-size:.62rem;animation:fi .3s}
  @keyframes fi{from{opacity:0;background:#1a1a00}to{opacity:1;background:transparent}}
  .evt-ts{color:var(--dim)}
  .evt-pga{color:var(--amber);margin:0 6px}
  .evt-si{font-weight:700}

  /* Download button */
  #dl-btn{display:block;width:calc(100% - 28px);margin:10px 14px;padding:8px;background:var(--surface);border:1px solid var(--border);color:var(--text);font-family:inherit;font-size:.68rem;border-radius:6px;cursor:pointer;text-align:center;text-transform:uppercase;letter-spacing:.06em}
  #dl-btn:hover{border-color:var(--accent);color:var(--accent)}
</style>
</head>
<body>
<header>
  <div class="logo">SeismoPhone <small>computer dashboard</small></div>
  <div class="hdr-right">
    <span id="phone-label">Phone: disconnected</span>
    <div id="phone-dot"></div>
  </div>
</header>

<div class="layout">
  <!-- Left: chart + metrics -->
  <div class="chart-col">
    <div class="chart-wrap"><canvas id="chart"></canvas></div>
    <div class="legend">
      <span><span class="dot" style="background:#44aaff"></span>X</span>
      <span><span class="dot" style="background:#2ecc71"></span>Y</span>
      <span><span class="dot" style="background:#e74c3c"></span>Z</span>
      <span><span class="dot" style="background:#f39c12"></span>Vector</span>
    </div>
    <div class="metrics-row">
      <div class="metric">
        <span class="metric-label">PGA</span>
        <span class="metric-value" id="pga">—</span>
        <span class="metric-sub">60s peak</span>
        <span class="metric-unit">Gal (cm/s²)</span>
      </div>
      <div class="metric">
        <span class="metric-label">PGV</span>
        <span class="metric-value" id="pgv">—</span>
        <span class="metric-sub">60s peak</span>
        <span class="metric-unit">cm/s</span>
      </div>
      <div class="metric">
        <span class="metric-label">Shindo</span>
        <span class="metric-value" id="shindo">—</span>
        <span class="metric-sub" id="shindo-jma">JMA —</span>
        <span class="metric-unit">I_inst continuous</span>
      </div>
      <div class="metric">
        <span class="metric-label">MMI</span>
        <span class="metric-value" id="mmi">—</span>
        <span class="metric-sub" id="mmi-sub">—</span>
        <span class="metric-unit">Worden 2012</span>
      </div>
    </div>
  </div>

  <!-- Right: QR + axis bars + events -->
  <div class="side-col">
    <div class="side-section">
      <div class="side-title">Connect Phone</div>
      <div id="qrcode"></div>
      <div class="qr-url" id="qr-url"></div>
    </div>

    <div class="side-section">
      <div class="side-title">Live Axes (Gal)</div>
      <div class="axis-row">
        <span class="axis-name">X</span>
        <div class="axis-track"><div class="axis-bar" id="bar-x" style="background:#44aaff"></div></div>
        <span class="axis-val" id="val-x">0.00</span>
      </div>
      <div class="axis-row">
        <span class="axis-name">Y</span>
        <div class="axis-track"><div class="axis-bar" id="bar-y" style="background:#2ecc71"></div></div>
        <span class="axis-val" id="val-y">0.00</span>
      </div>
      <div class="axis-row">
        <span class="axis-name">Z</span>
        <div class="axis-track"><div class="axis-bar" id="bar-z" style="background:#e74c3c"></div></div>
        <span class="axis-val" id="val-z">0.00</span>
      </div>
      <div style="height:8px"></div>
      <div class="scale-row">
        <div class="scale-seg" style="background:#333;flex:2"></div>
        <div class="scale-seg" style="background:#66aaff"></div>
        <div class="scale-seg" style="background:#44ddaa"></div>
        <div class="scale-seg" style="background:#aadd44"></div>
        <div class="scale-seg" style="background:#ffdd00"></div>
        <div class="scale-seg" style="background:#ffaa00"></div>
        <div class="scale-seg" style="background:#ff6600"></div>
        <div class="scale-seg" style="background:#ff3300"></div>
        <div class="scale-seg" style="background:#cc00cc"></div>
      </div>
      <div class="scale-lbl"><span>0</span><span>1</span><span>2</span><span>3</span><span>4</span><span>5−</span><span>5+</span><span>6−</span><span>6+</span><span>7</span></div>
    </div>

    <div class="side-section" style="flex:0 0 auto">
      <div class="side-title">Events</div>
    </div>
    <div class="evt-log" id="evt-log"></div>
    <button id="dl-btn">⬇ Download CSV</button>
  </div>
</div>

<script>
const CHART_SEC = 60;
const EST_HZ    = 50;
const MAX_PTS   = CHART_SEC * EST_HZ;

const SHINDO_LABELS = ['0','1','2','3','4','5−','5+','6−','6+','7'];
const SHINDO_COLORS = ['#555','#66aaff','#44ddaa','#aadd44','#ffdd00','#ffaa00','#ff6600','#ff3300','#dd0000','#cc00cc'];
const SHINDO_BOUNDS = [0.5,1.5,2.5,3.5,4.5,5.0,5.5,6.0,6.5];
const MMI_ROMAN     = ['','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII'];
const MMI_DESC      = ['','Not felt','Weak','Weak','Light','Moderate','Strong','Very strong','Severe','Violent','Extreme','Extreme','Extreme'];

function shindoIndex(I) {
  if (I === null || I === undefined) return -1;
  for (let i = 0; i < SHINDO_BOUNDS.length; i++) if (I < SHINDO_BOUNDS[i]) return i;
  return 9;
}

// Chart
const dX = new Float32Array(MAX_PTS);
const dY = new Float32Array(MAX_PTS);
const dZ = new Float32Array(MAX_PTS);
const dV = new Float32Array(MAX_PTS);
let ptr = 0;

const chart = new Chart(document.getElementById('chart').getContext('2d'), {
  type: 'line',
  data: {
    labels: new Array(MAX_PTS).fill(''),
    datasets: [
      {label:'X', data:Array.from(dX), borderColor:'#44aaff', borderWidth:1.1, pointRadius:0, tension:.2, fill:false},
      {label:'Y', data:Array.from(dY), borderColor:'#2ecc71', borderWidth:1.1, pointRadius:0, tension:.2, fill:false},
      {label:'Z', data:Array.from(dZ), borderColor:'#e74c3c', borderWidth:1.1, pointRadius:0, tension:.2, fill:false},
      {label:'V', data:Array.from(dV), borderColor:'#f39c12', borderWidth:1.4, pointRadius:0, tension:.2, fill:false},
    ],
  },
  options: {
    animation:false, responsive:true, maintainAspectRatio:false,
    plugins:{legend:{display:false}},
    scales:{
      x:{display:false},
      y:{
        grid:{color:'#13132a'},
        ticks:{color:'#5a607a',font:{family:'monospace',size:9},maxTicksLimit:6,callback:v=>v.toFixed(1)},
        title:{display:true,text:'Gal',color:'#5a607a',font:{size:9}},
      },
    },
  },
});

// QR code
const wsUrl = 'ws://' + location.hostname + ':__WS_PORT__';
document.getElementById('qr-url').textContent = wsUrl;
new QRCode(document.getElementById('qrcode'), {
  text: wsUrl, width: 160, height: 160,
  colorDark: '#ffffff', colorLight: '#080812',
  correctLevel: QRCode.CorrectLevel.M,
});

// CSV log
const csvRows = [['timestamp','ax','ay','az','vec','pga','pgv','shindo_i','shindo','mmi']];

// WebSocket to server
const vws = new WebSocket('ws://' + location.host.replace(':__HTTP_PORT__', '') + ':__WS_PORT__/view');
vws.onopen  = () => console.log('Dashboard connected to server');
vws.onclose = () => { document.getElementById('phone-label').textContent = 'Server disconnected'; };

vws.onmessage = e => {
  const d = JSON.parse(e.data);
  if (d._phone_status !== undefined) {
    const on = d._phone_status === 'connected';
    document.getElementById('phone-dot').className = on ? 'on' : '';
    document.getElementById('phone-label').textContent = 'Phone: ' + (on ? 'connected' : 'disconnected');
    return;
  }

  // Chart
  dX[ptr]=d.ax; dY[ptr]=d.ay; dZ[ptr]=d.az; dV[ptr]=d.vec;
  ptr = (ptr + 1) % MAX_PTS;
  const xs=[],ys=[],zs=[],vs=[];
  for (let i=0;i<MAX_PTS;i++){const idx=(ptr+i)%MAX_PTS;xs.push(dX[idx]);ys.push(dY[idx]);zs.push(dZ[idx]);vs.push(dV[idx]);}
  chart.data.datasets[0].data=xs;
  chart.data.datasets[1].data=ys;
  chart.data.datasets[2].data=zs;
  chart.data.datasets[3].data=vs;
  chart.update('none');

  // Axis bars
  setBar('bar-x','val-x',d.ax,15);
  setBar('bar-y','val-y',d.ay,15);
  setBar('bar-z','val-z',d.az,15);

  // Metrics
  const pgaEl=document.getElementById('pga');
  pgaEl.textContent = d.pga.toFixed(2);
  colorM(pgaEl, d.pga, 2, 20);

  const pgvEl=document.getElementById('pgv');
  pgvEl.textContent = d.pgv.toFixed(3);
  colorM(pgvEl, d.pgv, 0.1, 1);

  const siEl=document.getElementById('shindo');
  const siJma=document.getElementById('shindo-jma');
  if (d.shindo_i !== null && d.shindo_i !== undefined) {
    const idx = shindoIndex(d.shindo_i);
    const dispI = Math.max(0, d.shindo_i);
    siEl.textContent  = dispI.toFixed(2);
    siEl.style.color  = SHINDO_COLORS[Math.max(0,idx)];
    siJma.textContent = 'JMA ' + (d.shindo || '0');
    siJma.style.color = SHINDO_COLORS[Math.max(0,idx)];
  } else {
    siEl.textContent='—'; siEl.style.color='var(--dim)';
    siJma.textContent='JMA —';
  }

  const mmiEl=document.getElementById('mmi');
  const mmiSub=document.getElementById('mmi-sub');
  mmiEl.textContent  = d.mmi.toFixed(2);
  const mmiR = Math.round(d.mmi);
  mmiSub.textContent = (MMI_ROMAN[mmiR]||'I') + ' — ' + (MMI_DESC[mmiR]||'');
  colorM(mmiEl, d.mmi, 4, 6);

  // Event detection (PGA > 3 Gal)
  if (d.pga > 3) maybeLogEvent(d);

  // CSV
  csvRows.push([d.t,d.ax,d.ay,d.az,d.vec,d.pga,d.pgv,d.shindo_i,d.shindo,d.mmi]);
  if (csvRows.length > 100000) csvRows.splice(1, csvRows.length - 100000);
};

let lastEvtPga = 0, lastEvtT = 0;
function maybeLogEvent(d) {
  if (d.pga <= lastEvtPga * 0.5 || Date.now() - lastEvtT < 5000) return;
  lastEvtPga = d.pga; lastEvtT = Date.now();
  const log = document.getElementById('evt-log');
  const ts  = new Date(d.t).toLocaleTimeString('en',{hour12:false});
  const idx = shindoIndex(d.shindo_i);
  const item = document.createElement('div');
  item.className = 'evt-item';
  item.innerHTML = `<span class="evt-ts">${ts}</span>
    <span class="evt-pga">PGA ${d.pga.toFixed(2)} Gal · PGV ${d.pgv.toFixed(3)} · MMI ${d.mmi.toFixed(2)}</span>
    <span class="evt-si" style="color:${SHINDO_COLORS[Math.max(0,idx)]}">S${Math.max(0,d.shindo_i??0).toFixed(2)}</span>`;
  log.prepend(item);
  while (log.children.length > 20) log.removeChild(log.lastChild);
}

function setBar(id,valId,v,max){
  const b=document.getElementById(id);
  const pct=Math.min(Math.abs(v)/max*50,50);
  b.style.left  = v>=0?'50%':(50-pct)+'%';
  b.style.width = pct+'%';
  document.getElementById(valId).textContent=v.toFixed(2);
}
function colorM(el,v,warn,hi){
  el.classList.remove('warn','high');
  if(v>=hi)el.classList.add('high');
  else if(v>=warn)el.classList.add('warn');
}

document.getElementById('dl-btn').addEventListener('click', () => {
  const blob = new Blob([csvRows.map(r=>r.join(',')).join('\n')], {type:'text/csv'});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'seismophone_' + new Date().toISOString().slice(0,19).replace(/:/g,'-') + '.csv';
  a.click();
});
</script>
</body>
</html>"""

# ── WebSocket server ──────────────────────────────────────────────────────
async def ws_handler(websocket):
    global phone_connected

    path = '/'
    try:
        path = websocket.request.path
    except AttributeError:
        path = getattr(websocket, 'path', '/')

    if '/view' in path:
        # Computer dashboard viewer
        viewers.add(websocket)
        # Send current phone status
        try:
            await websocket.send(json.dumps({'_phone_status': 'connected' if phone_connected else 'disconnected'}))
            if latest:
                await websocket.send(json.dumps(latest))
            await websocket.wait_closed()
        except Exception:
            pass
        finally:
            viewers.discard(websocket)
    else:
        # Phone sender
        phone_connected = True
        print(f"  📱 Phone connected")
        status_msg = json.dumps({'_phone_status': 'connected'})
        dead = set()
        for v in viewers:
            try:
                await v.send(status_msg)
            except Exception:
                dead.add(v)
        viewers.difference_update(dead)

        try:
            async for msg in websocket:
                try:
                    data = json.loads(msg)
                    latest.update(data)
                except Exception:
                    continue
                if viewers:
                    dead = set()
                    await asyncio.gather(
                        *[v.send(msg) for v in viewers],
                        return_exceptions=True
                    )
        except Exception:
            pass
        finally:
            phone_connected = False
            print(f"  📵 Phone disconnected")
            status_msg = json.dumps({'_phone_status': 'disconnected'})
            for v in list(viewers):
                try:
                    await v.send(status_msg)
                except Exception:
                    pass

# ── HTTP server (dashboard) ───────────────────────────────────────────────
class DashHandler(BaseHTTPRequestHandler):
    html = (DASHBOARD
            .replace('__WS_PORT__', str(WS_PORT))
            .replace('__HTTP_PORT__', str(HTTP_PORT)))

    def do_GET(self):
        self.send_response(200)
        self.send_header('Content-Type', 'text/html; charset=utf-8')
        self.send_header('Cache-Control', 'no-cache')
        self.end_headers()
        self.wfile.write(self.html.encode())

    def log_message(self, *_):
        pass

def get_local_ips():
    ips = []
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(('8.8.8.8', 80))
        ips.append(s.getsockname()[0])
        s.close()
    except Exception:
        pass
    return ips or ['127.0.0.1']

def run_http():
    srv = HTTPServer(('0.0.0.0', HTTP_PORT), DashHandler)
    srv.serve_forever()

# ── Main ──────────────────────────────────────────────────────────────────
async def main():
    import websockets

    ips = get_local_ips()
    primary_ip = ips[0]

    print()
    print('  ╔══════════════════════════════════════════╗')
    print('  ║         SeismoPhone Server               ║')
    print('  ╚══════════════════════════════════════════╝')
    print()
    print(f'  Dashboard  →  http://localhost:{HTTP_PORT}')
    print()
    print(f'  Phone WebSocket URL:')
    print(f'  ➜  ws://{primary_ip}:{WS_PORT}')
    print()
    print('  Enter the URL above in the phone app (tap ⇆)')
    print('  Phone and computer must be on the same Wi-Fi.')
    print()
    print('  Ctrl+C to stop.')
    print()

    # Start HTTP server in background thread
    t = threading.Thread(target=run_http, daemon=True)
    t.start()

    # Open dashboard
    webbrowser.open(f'http://localhost:{HTTP_PORT}')

    # Start WebSocket server
    async with websockets.serve(ws_handler, '0.0.0.0', WS_PORT):
        await asyncio.Future()

if __name__ == '__main__':
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print('\n  Server stopped.')
