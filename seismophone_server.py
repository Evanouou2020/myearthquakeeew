#!/usr/bin/env python3
"""
SeismoPhone Computer Server
────────────────────────────
Run on your computer:
    pip install websockets
    python3 seismophone_server.py

Two ports — no websockets version hacks needed:
  Dashboard  →  http://localhost:8766   (open in browser)
  Phone WS   →  ws://YOUR-IP:8765      (enter in phone app)
"""

import asyncio
import json
import socket
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

try:
    import websockets
except ImportError:
    print("  ERROR: run  pip install websockets  first")
    raise

HTTP_PORT = 8766
WS_PORT   = 8765

# ── Shared state ───────────────────────────────────────────────────────────
viewers: set = set()
latest: dict = {}
phone_connected = False


# ── Dashboard HTML ─────────────────────────────────────────────────────────
DASHBOARD = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SeismoPhone — Dashboard</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/qrcodejs@1.0.0/qrcode.min.js"></script>
<style>
  *{box-sizing:border-box;margin:0;padding:0}
  :root{--bg:#080812;--surface:#10101e;--border:#1e1e3a;--accent:#44aaff;--green:#2ecc71;--red:#e74c3c;--amber:#f39c12;--text:#c8d0e8;--dim:#5a607a}
  html,body{background:var(--bg);color:var(--text);font-family:'SF Mono','Fira Code','Courier New',monospace;height:100vh;overflow:hidden}
  header{display:flex;align-items:center;justify-content:space-between;padding:10px 20px;background:var(--surface);border-bottom:1px solid var(--border)}
  .logo{font-size:1.1rem;font-weight:700;color:var(--accent);letter-spacing:.04em}
  .logo small{color:var(--dim);font-weight:400;font-size:.68rem;margin-left:8px}
  .hdr-r{display:flex;align-items:center;gap:12px;font-size:.72rem;color:var(--dim)}
  #phone-dot{width:9px;height:9px;border-radius:50%;background:var(--dim);transition:background .3s,box-shadow .3s}
  #phone-dot.on{background:var(--green);box-shadow:0 0 8px var(--green);animation:blink 1.5s infinite}
  @keyframes blink{0%,100%{opacity:1}50%{opacity:.4}}

  .layout{display:grid;grid-template-columns:1fr 290px;height:calc(100vh - 48px)}
  .main{display:flex;flex-direction:column;border-right:1px solid var(--border);overflow:hidden}
  .chart-wrap{flex:1;min-height:0;padding:10px 8px 0;background:var(--surface)}
  .legend{display:flex;gap:14px;padding:4px 14px;font-size:.65rem;color:var(--dim)}
  .legend span{display:flex;align-items:center;gap:5px}
  .dot{width:8px;height:8px;border-radius:50%}

  .mode-bar{display:flex;border-top:1px solid var(--border);flex-shrink:0}
  .mode-btn{flex:1;padding:7px 6px;font-family:inherit;font-size:.7rem;font-weight:600;letter-spacing:.08em;text-transform:uppercase;background:var(--surface);border:none;border-right:1px solid var(--border);color:var(--dim);cursor:pointer;transition:background .15s}
  .mode-btn:last-child{border-right:none}
  .mode-btn.on{background:rgba(68,170,255,.18);color:var(--accent)}
  .mode-btn.max-on{background:rgba(243,156,18,.18);color:var(--amber)}

  .metrics{display:grid;grid-template-columns:1fr 1fr 1fr 1fr;border-top:1px solid var(--border);flex-shrink:0}
  .metric{display:flex;flex-direction:column;align-items:center;justify-content:center;padding:10px 3px 8px;border-right:1px solid var(--border);gap:2px}
  .metric:last-child{border-right:none}
  .metric-label{font-size:.57rem;color:var(--dim);text-transform:uppercase;letter-spacing:.09em}
  .metric-value{font-size:1.45rem;font-weight:700;color:var(--accent);transition:color .2s}
  .metric-value.live-mode{color:var(--accent)}
  .metric-value.max-mode{color:var(--amber)}
  .metric-sub{font-size:.64rem;color:var(--dim)}
  .metric-sub2{font-size:.58rem;color:#3a4060}
  .metric-unit{font-size:.5rem;color:#1e1e38}
  .metric-value.warn{color:var(--amber)!important}
  .metric-value.high{color:var(--red)!important}

  .side{display:flex;flex-direction:column;overflow:hidden}
  .side-sec{padding:12px 14px;border-bottom:1px solid var(--border);flex-shrink:0}
  .side-title{font-size:.63rem;text-transform:uppercase;letter-spacing:.09em;color:var(--dim);margin-bottom:9px}
  #qrcode{display:flex;justify-content:center;margin:6px 0}
  .qr-url{font-size:.6rem;color:var(--dim);text-align:center;word-break:break-all;margin-top:4px}
  .a-row{display:flex;align-items:center;gap:8px;font-size:.64rem;margin-bottom:5px}
  .a-name{width:14px;color:var(--dim)}
  .a-track{flex:1;height:6px;background:var(--border);border-radius:3px;overflow:hidden;position:relative}
  .a-bar{position:absolute;top:0;bottom:0;border-radius:3px;transition:width .07s,left .07s}
  .a-val{width:52px;text-align:right;font-size:.62rem;color:var(--dim)}
  .scale-row{display:flex;gap:2px}
  .scale-seg{flex:1;height:4px;border-radius:2px}
  .scale-lbl{display:flex;justify-content:space-between;font-size:.52rem;color:var(--dim);margin-top:2px}
  .evt-wrap{flex:1;overflow-y:auto}
  .evt-item{padding:5px 14px;border-bottom:1px solid var(--border);font-size:.61rem;animation:fi .3s}
  @keyframes fi{from{opacity:0;background:#1a1400}to{opacity:1}}
  .evt-ts{color:var(--dim)}
  .evt-pga{color:var(--amber);margin:0 5px}
  .evt-si{font-weight:700}
  #dl-btn{display:block;width:calc(100% - 28px);margin:8px 14px;padding:8px;background:var(--surface);border:1px solid var(--border);color:var(--text);font-family:inherit;font-size:.67rem;border-radius:6px;cursor:pointer;text-align:center;text-transform:uppercase;letter-spacing:.06em;flex-shrink:0}
  #dl-btn:hover{border-color:var(--accent);color:var(--accent)}
</style>
</head>
<body>
<header>
  <div class="logo">SeismoPhone <small>computer dashboard</small></div>
  <div class="hdr-r">
    <span id="phone-label">Phone: disconnected</span>
    <div id="phone-dot"></div>
  </div>
</header>
<div class="layout">
  <div class="main">
    <div class="chart-wrap"><canvas id="chart"></canvas></div>
    <div class="legend">
      <span><span class="dot" style="background:#44aaff"></span>X</span>
      <span><span class="dot" style="background:#2ecc71"></span>Y</span>
      <span><span class="dot" style="background:#e74c3c"></span>Z</span>
      <span><span class="dot" style="background:#f39c12"></span>Vector</span>
    </div>
    <div class="mode-bar">
      <button class="mode-btn on" id="btn-live">● LIVE (instant)</button>
      <button class="mode-btn"   id="btn-max">▲ MAX <span id="max-hint" style="font-size:.6rem;font-weight:400;opacity:.7"></span></button>
    </div>
    <div class="metrics">
      <div class="metric"><span class="metric-label">PGA</span><span class="metric-value" id="pga">—</span><span class="metric-sub" id="pga-sub">Gal</span><span class="metric-unit">cm/s²</span></div>
      <div class="metric"><span class="metric-label">PGV</span><span class="metric-value" id="pgv">—</span><span class="metric-sub" id="pgv-sub">cm/s</span><span class="metric-unit">integrated</span></div>
      <div class="metric"><span class="metric-label">Shindo</span><span class="metric-value" id="shindo">—</span><span class="metric-sub" id="shindo-jma" style="font-weight:700">JMA —</span><span class="metric-unit">I = 2log(PGA)+0.94</span></div>
      <div class="metric"><span class="metric-label">MMI</span><span class="metric-value" id="mmi">—</span><span class="metric-sub" id="mmi-sub">—</span><span class="metric-sub2" id="mmi-pgv-sub">PGV→ —</span><span class="metric-unit">Worden 2012</span></div>
    </div>
  </div>
  <div class="side">
    <div class="side-sec">
      <div class="side-title">Phone WebSocket URL</div>
      <div id="qrcode"></div>
      <div class="qr-url" id="qr-url"></div>
    </div>
    <div class="side-sec">
      <div class="side-title">Live Axes (Gal)</div>
      <div class="a-row"><span class="a-name">X</span><div class="a-track"><div class="a-bar" id="bar-x" style="background:#44aaff"></div></div><span class="a-val" id="val-x">0.00</span></div>
      <div class="a-row"><span class="a-name">Y</span><div class="a-track"><div class="a-bar" id="bar-y" style="background:#2ecc71"></div></div><span class="a-val" id="val-y">0.00</span></div>
      <div class="a-row"><span class="a-name">Z</span><div class="a-track"><div class="a-bar" id="bar-z" style="background:#e74c3c"></div></div><span class="a-val" id="val-z">0.00</span></div>
      <div style="height:6px"></div>
      <div class="scale-row">
        <div class="scale-seg" style="background:#999999;flex:2"></div>
        <div class="scale-seg" style="background:#99aaff"></div>
        <div class="scale-seg" style="background:#80ffff"></div>
        <div class="scale-seg" style="background:#7fff00"></div>
        <div class="scale-seg" style="background:#ffff00"></div>
        <div class="scale-seg" style="background:#ffa500"></div>
        <div class="scale-seg" style="background:#ff6600"></div>
        <div class="scale-seg" style="background:#ff4500"></div>
        <div class="scale-seg" style="background:#880000"></div>
      </div>
      <div class="scale-lbl"><span>0</span><span>1</span><span>2</span><span>3</span><span>4</span><span>5−</span><span>5+</span><span>6−</span><span>6+</span><span>7</span></div>
    </div>
    <div class="side-sec" style="padding-bottom:6px"><div class="side-title">Events</div></div>
    <div class="evt-wrap" id="evt-log"></div>
    <button id="dl-btn">⬇ Download CSV</button>
  </div>
</div>
<script>
const SHINDO_BOUNDS=[0.5,1.5,2.5,3.5,4.5,5.0,5.5,6.0,6.5];
const SHINDO_LABELS=['0','1','2','3','4','5−','5+','6−','6+','7'];
const SHINDO_COLORS=['#999999','#99aaff','#80ffff','#7fff00','#ffff00','#ffa500','#ff6600','#ff4500','#ff0000','#880000'];
const MMI_ROMAN=['','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII'];
const MMI_DESC=['','Not felt','Weak','Weak','Light','Moderate','Strong','Very strong','Severe','Violent','Extreme','Extreme','Extreme'];
function shindoFloat(g){return g<0.001?null:2*Math.log10(g)+0.94;}
function shindoIndex(I){if(I===null)return -1;for(let i=0;i<SHINDO_BOUNDS.length;i++)if(I<SHINDO_BOUNDS[i])return i;return 9;}
function mmiFromPGA(g){if(g<0.05)return 1;const l=Math.log10(g);return Math.max(1,g<43.7?1.78*l+1.55:3.70*l-1.60);}
function mmiFromPGV(v){if(v<0.001)return 1;const l=Math.log10(v);return Math.max(1,v<3.36?1.47*l+3.78:3.16*l+2.89);}
function mmiLabel(v){const r=Math.round(v);return r<=12?(MMI_ROMAN[r]||'I'):'XII+';}
function mmiDesc(v){const r=Math.round(v);return r<=12?(MMI_DESC[r]||''):'Beyond scale';}

const CHART_SEC=60,EST_HZ=50,MAX_PTS=CHART_SEC*EST_HZ;
const dX=new Float32Array(MAX_PTS),dY=new Float32Array(MAX_PTS),dZ=new Float32Array(MAX_PTS),dV=new Float32Array(MAX_PTS);
let ptr=0;
const chart=new Chart(document.getElementById('chart').getContext('2d'),{
  type:'line',
  data:{labels:new Array(MAX_PTS).fill(''),datasets:[
    {label:'X',data:Array.from(dX),borderColor:'#44aaff',borderWidth:1.1,pointRadius:0,tension:.2,fill:false},
    {label:'Y',data:Array.from(dY),borderColor:'#2ecc71',borderWidth:1.1,pointRadius:0,tension:.2,fill:false},
    {label:'Z',data:Array.from(dZ),borderColor:'#e74c3c',borderWidth:1.1,pointRadius:0,tension:.2,fill:false},
    {label:'V',data:Array.from(dV),borderColor:'#f39c12',borderWidth:1.4,pointRadius:0,tension:.2,fill:false},
  ]},
  options:{animation:false,responsive:true,maintainAspectRatio:false,plugins:{legend:{display:false}},
    scales:{x:{display:false},y:{grid:{color:'#13132a'},ticks:{color:'#5a607a',font:{family:'monospace',size:9},maxTicksLimit:6,callback:v=>v.toFixed(1)},title:{display:true,text:'Gal',color:'#5a607a',font:{size:9}}}}},
});

// WebSocket connects to WS port 8765 (not same as HTTP port 8766)
const wsUrl='ws://'+location.hostname+':8765';
document.getElementById('qr-url').textContent=wsUrl;
new QRCode(document.getElementById('qrcode'),{text:wsUrl,width:150,height:150,colorDark:'#ffffff',colorLight:'#080812',correctLevel:QRCode.CorrectLevel.M});

let displayMode='live';
let dashMaxPGA=0,dashMaxPGV=0,dashMaxTime=null;
function setMode(m){
  displayMode=m;
  document.getElementById('btn-live').className='mode-btn'+(m==='live'?' on':'');
  document.getElementById('btn-max').className='mode-btn'+(m==='max'?' max-on':'');
}
document.getElementById('btn-live').addEventListener('click',()=>setMode('live'));
document.getElementById('btn-max').addEventListener('click',()=>setMode('max'));

const CSV_HDR='time_iso,time_ms,ax,ay,az,vec,live_pga,live_pgv,live_shindo_i,live_mmi_pga,live_mmi_pgv,max_pga,max_pgv';
const csvRows=[];

const vws=new WebSocket(wsUrl);
vws.onopen=()=>{vws.send(JSON.stringify({role:'view'}));};
vws.onclose=()=>{document.getElementById('phone-label').textContent='Server: disconnected';};
vws.onerror=()=>{document.getElementById('phone-label').textContent='WS error — is server running?';};

vws.onmessage=e=>{
  const d=JSON.parse(e.data);
  if(d._phone_status!==undefined){
    const on=d._phone_status==='connected';
    document.getElementById('phone-dot').className=on?'on':'';
    document.getElementById('phone-label').textContent='Phone: '+(on?'connected':'disconnected');
    return;
  }
  dX[ptr]=d.ax;dY[ptr]=d.ay;dZ[ptr]=d.az;dV[ptr]=d.vec;
  ptr=(ptr+1)%MAX_PTS;
  const xs=[],ys=[],zs=[],vs=[];
  for(let i=0;i<MAX_PTS;i++){const idx=(ptr+i)%MAX_PTS;xs.push(dX[idx]);ys.push(dY[idx]);zs.push(dZ[idx]);vs.push(dV[idx]);}
  chart.data.datasets[0].data=xs;chart.data.datasets[1].data=ys;chart.data.datasets[2].data=zs;chart.data.datasets[3].data=vs;
  chart.update('none');
  setBar('bar-x','val-x',d.ax,15);setBar('bar-y','val-y',d.ay,15);setBar('bar-z','val-z',d.az,15);

  const livePGA=d.live_pga??0,livePGV=d.live_pgv??0;
  if(livePGA>dashMaxPGA){dashMaxPGA=livePGA;dashMaxTime=new Date();}
  if(livePGV>dashMaxPGV)dashMaxPGV=livePGV;

  const isMax=displayMode==='max';
  const pga=isMax?dashMaxPGA:livePGA,pgv=isMax?dashMaxPGV:livePGV;

  const hint=document.getElementById('max-hint');
  hint.textContent=isMax&&dashMaxTime?'peaked '+dashMaxTime.toLocaleTimeString('en',{hour12:false}):'';

  const pgaEl=document.getElementById('pga');
  pgaEl.textContent=pga.toFixed(2);pgaEl.className='metric-value '+(isMax?'max-mode':'live-mode');colorM(pgaEl,pga,2,20);
  document.getElementById('pga-sub').textContent=isMax?'session max':'instant';

  const pgvEl=document.getElementById('pgv');
  pgvEl.textContent=pgv.toFixed(3);pgvEl.className='metric-value '+(isMax?'max-mode':'live-mode');colorM(pgvEl,pgv,.1,1);
  document.getElementById('pgv-sub').textContent=isMax?'session max':'instant';

  const si=shindoFloat(pga),siIdx=shindoIndex(si);
  const siEl=document.getElementById('shindo'),siJma=document.getElementById('shindo-jma');
  if(si===null||pga<0.01){siEl.textContent='—';siEl.style.color='var(--dim)';siJma.textContent='JMA —';siJma.style.color='var(--dim)';}
  else{const dispI=Math.max(0,si);siEl.textContent=dispI.toFixed(2);siEl.style.color=SHINDO_COLORS[Math.max(0,siIdx)];siJma.textContent='JMA '+SHINDO_LABELS[Math.max(0,siIdx)];siJma.style.color=SHINDO_COLORS[Math.max(0,siIdx)];}
  siEl.className='metric-value '+(isMax?'max-mode':'live-mode');

  const mPGA=mmiFromPGA(pga),mPGV=mmiFromPGV(pgv);
  const mmiEl=document.getElementById('mmi');
  mmiEl.textContent=mPGA.toFixed(2);mmiEl.className='metric-value '+(isMax?'max-mode':'live-mode');
  const mR=Math.round(mPGA);
  document.getElementById('mmi-sub').textContent=mmiLabel(mPGA)+' '+mmiDesc(mPGA);
  document.getElementById('mmi-pgv-sub').textContent='PGV→ '+mPGV.toFixed(2);
  colorM(mmiEl,mPGA,4,6);

  if(livePGA>3)maybeLogEvent(livePGA,livePGV);
  const now=new Date();
  csvRows.push([now.toISOString(),now.getTime(),d.ax,d.ay,d.az,d.vec,d.live_pga,d.live_pgv,d.live_shindo_i,d.live_mmi_pga,d.live_mmi_pgv,d.max_pga,d.max_pgv]);
  if(csvRows.length>100000)csvRows.splice(0,csvRows.length-100000);
};

let lastEvtT=0,lastEvtPGA=0;
function maybeLogEvent(pga,pgv){
  if(Date.now()-lastEvtT<5000&&pga<=lastEvtPGA*1.2)return;
  lastEvtT=Date.now();lastEvtPGA=pga;
  const log=document.getElementById('evt-log');
  const ts=new Date().toLocaleTimeString('en',{hour12:false});
  const si=shindoFloat(pga),idx=shindoIndex(si),mi=mmiFromPGA(pga);
  const item=document.createElement('div');item.className='evt-item';
  item.innerHTML=`<span class="evt-ts">${ts}</span><span class="evt-pga">PGA ${pga.toFixed(2)} · PGV ${pgv.toFixed(3)} · MMI ${mi.toFixed(2)}</span><span class="evt-si" style="color:${SHINDO_COLORS[Math.max(0,idx)]}">${si!==null?'S'+Math.max(0,si).toFixed(2):'—'}</span>`;
  log.prepend(item);while(log.children.length>20)log.removeChild(log.lastChild);
}
function setBar(id,valId,v,max){const b=document.getElementById(id);const pct=Math.min(Math.abs(v)/max*50,50);b.style.left=v>=0?'50%':(50-pct)+'%';b.style.width=pct+'%';document.getElementById(valId).textContent=v.toFixed(2);}
function colorM(el,v,warn,hi){el.classList.remove('warn','high');if(v>=hi)el.classList.add('high');else if(v>=warn)el.classList.add('warn');}

document.getElementById('dl-btn').addEventListener('click',()=>{
  const blob=new Blob([CSV_HDR+'\n'+csvRows.map(r=>r.join(',')).join('\n')],{type:'text/csv'});
  const a=document.createElement('a');a.href=URL.createObjectURL(blob);
  a.download='seismophone_'+new Date().toISOString().slice(0,19).replace(/:/g,'-')+'.csv';a.click();
});
</script>
</body>
</html>"""


# ── WebSocket handler — role determined by first message ───────────────────
async def ws_handler(websocket):
    global phone_connected

    try:
        raw = await asyncio.wait_for(websocket.recv(), timeout=10.0)
        first = json.loads(raw)
    except Exception as ex:
        print(f"  ⚠ Handshake error: {ex}")
        return

    role = first.get("role", "phone")

    if role == "view":
        viewers.add(websocket)
        try:
            await websocket.send(json.dumps({
                "_phone_status": "connected" if phone_connected else "disconnected"
            }))
            if latest:
                await websocket.send(json.dumps(latest))
            await websocket.wait_closed()
        except Exception:
            pass
        finally:
            viewers.discard(websocket)

    else:
        phone_connected = True
        print("  📱 Phone connected")
        status = json.dumps({"_phone_status": "connected"})
        for v in list(viewers):
            try:
                await v.send(status)
            except Exception:
                viewers.discard(v)

        try:
            latest.update({k: v for k, v in first.items() if k != "role"})
            payload = json.dumps({k: v for k, v in first.items() if k != "role"})
            dead = set()
            for v in viewers:
                try:
                    await v.send(payload)
                except Exception:
                    dead.add(v)
            viewers.difference_update(dead)
        except Exception:
            pass

        try:
            async for msg in websocket:
                try:
                    data = json.loads(msg)
                    latest.update(data)
                except Exception:
                    continue
                dead = set()
                for v in viewers:
                    try:
                        await v.send(msg)
                    except Exception:
                        dead.add(v)
                viewers.difference_update(dead)
        except Exception:
            pass
        finally:
            phone_connected = False
            print("  📵 Phone disconnected")
            disc = json.dumps({"_phone_status": "disconnected"})
            for v in list(viewers):
                try:
                    await v.send(disc)
                except Exception:
                    pass


# ── HTTP server for dashboard (separate port, no version hacks) ────────────
class DashboardHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        body = DASHBOARD.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass  # suppress per-request logs


def run_http():
    HTTPServer(("0.0.0.0", HTTP_PORT), DashboardHandler).serve_forever()


# ── Helpers ────────────────────────────────────────────────────────────────
def get_local_ip():
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return "127.0.0.1"


# ── Main ───────────────────────────────────────────────────────────────────
async def main():
    ip = get_local_ip()

    threading.Thread(target=run_http, daemon=True).start()

    print()
    print("  ╔══════════════════════════════════════════════╗")
    print("  ║          SeismoPhone Server                  ║")
    print("  ╚══════════════════════════════════════════════╝")
    print()
    print(f"  ┌──────────────────────────────────────────────┐")
    print(f"  │  1. Dashboard — open in browser:             │")
    print(f"  │     http://localhost:{HTTP_PORT}                      │")
    print(f"  │                                              │")
    print(f"  │  2. Phone — enter in app (tap ⇆):            │")
    print(f"  │     ws://{ip}:{WS_PORT}                      │")
    print(f"  └──────────────────────────────────────────────┘")
    print()
    print("  Phone and computer must be on the same Wi-Fi.")
    print("  Ctrl+C to stop.")
    print()

    # Explicit socket binding so macOS doesn't silently restrict to localhost
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    sock.bind(("0.0.0.0", WS_PORT))
    sock.setblocking(False)

    async with websockets.serve(ws_handler, sock=sock):
        await asyncio.Future()


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\n  Server stopped.")
