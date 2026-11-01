"""Line WO, copied from line WG's harness: serve gbdraw's LOSAT browser host (gbdraw 900992caf, exported read-only with git archive)
at /arm/<arm>/, with the threaded command module of each arm. Adapted from gbdraw's
analysis/ov2369/serve.py (probes on thread faults and the main worker's memory size).
usage: serve.py <port>"""
import sys, mimetypes
from functools import partial
from pathlib import Path
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
H = Path(__file__).parent
WGH = Path("/home/kawato/losat-baselines/wg-wasm-guard-20261101/cmdharness")
WEB = WGH.parent / "gbdraw-web-900992caf/gbdraw/web"
ARMS = {
    "devampv1": H / "art/ampv1.wasm",  # dev 74bbd3c7e, widened window, line WG's guard (fill and copy)
    "devampv2": H / "art/ampv2.wasm",  # the same module, the guard with atomic instructions
}
PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 9893
mimetypes.add_type("text/javascript", ".mjs"); mimetypes.add_type("application/wasm", ".wasm")
ISO = {"Cross-Origin-Opener-Policy": "same-origin", "Cross-Origin-Embedder-Policy": "require-corp",
       "Cross-Origin-Resource-Policy": "same-origin", "Cache-Control": "no-store"}
def rep(s, old, new):
    assert old and old in s, old
    return s.replace(old, new)
def thread_worker():
    s = (WEB / "js/workers/losat-wasi-thread-worker.js").read_text()
    s = rep(s, "return { instance, stdoutChunks, stderrChunks };", "return { instance, memory, stdoutChunks, stderrChunks };")
    s = rep(s, "context.instance.exports.wasi_thread_start(tid, startArg);",
            "var __t0 = performance.now(); context.instance.exports.wasi_thread_start(tid, startArg);")
    s = rep(s, "reportThreadFault(tid, error, context);", """
    try {
      const m = context.memory; const before = m.buffer.byteLength; m.grow(0); const after = m.buffer.byteLength;
      const c = new BroadcastChannel('ov2369-probe');
      c.postMessage({ kind: 'thread', tid, startArg, ms: performance.now() - __t0, before, after, name: error?.name, message: error?.message, stack: String(error?.stack || '').slice(0, 3000) });
      c.close();
    } catch (e) {}
    reportThreadFault(tid, error, context);""")
    return s
def declared_initial_pages(path):
    """The initial pages of the module's imported memory (LOSAT/tests/wasi_artifact.js)."""
    import json, subprocess
    out = subprocess.check_output(["node", "/home/kawato/losat-work/.worktrees/wo-guard/LOSAT/tests/wasi_artifact.js", str(path), "threaded-command"], text=True)
    return json.loads(out)["memory"]["initial"]
def job_worker(arm):
    s = (WEB / "js/workers/losat-threaded-worker.js").read_text()
    # gbdraw's host creates the memory with 21 pages, its module's declared minimum; a module that
    # declares more (LOSAT dev: 43) cannot be instantiated on it, so the host gets the module's own.
    s = rep(s, "const DEFAULT_INITIAL_MEMORY_PAGES = 21;", f"const DEFAULT_INITIAL_MEMORY_PAGES = {max(21, declared_initial_pages(ARMS[arm]))};")
    if arm.startswith("dev"):
        # LOSAT dev's command takes the NCBI option names (CLI v2); gbdraw's host writes the v1 names.
        for old, new in (("'--query'", "'-query'"), ("'--subject'", "'-subject'"), ("'--outfmt'", "'-outfmt'"),
                         ("args.push('--num-threads'", "args.push('-num_threads'")):
            s = rep(s, old, new)
    s = rep(s, "    const exitCode = wasi.start(instance);", """    const exitCode = wasi.start(instance);
    try { const c = new BroadcastChannel('ov2369-probe'); c.postMessage({ kind: 'size', bytes: mainMemory.buffer.byteLength, spawnCount }); c.close(); } catch (e) {}""")
    s = rep(s, "  } catch (error) {\n    self.postMessage({\n      id,\n      type: 'run',\n      ok: false,", """  } catch (error) {
    try { const c = new BroadcastChannel('ov2369-probe'); c.postMessage({ kind: 'main', name: error?.name, message: error?.message, stack: String(error?.stack || '').slice(0, 3000) }); c.close(); } catch (e) {}
    self.postMessage({\n      id,\n      type: 'run',\n      ok: false,""")
    return s
class Handler(SimpleHTTPRequestHandler):
    def end_headers(self):
        for k, v in ISO.items(): self.send_header(k, v)
        super().end_headers()
    def log_message(self, *a): pass
    def _send(self, body, ctype):
        b = body.encode() if isinstance(body, str) else body
        self.send_response(200); self.send_header("Content-Type", ctype); self.send_header("Content-Length", str(len(b))); self.end_headers(); self.wfile.write(b)
    def do_GET(self):
        p = self.path.split("?")[0]
        if p.startswith("/arm/"):
            _, _, arm, rest = p.split("/", 3)
            if rest == "harness.html": return self._send((H / "harness.html").read_text(), "text/html")
            if rest.startswith("data/"): return self._send((WGH / rest).read_bytes(), "application/json")
            if rest == "js/workers/losat-wasi-thread-worker.js": return self._send(thread_worker(), "text/javascript")
            if rest == "js/workers/losat-threaded-worker.js": return self._send(job_worker(arm), "text/javascript")
            if rest == "wasm/losat/losat-threaded.wasm": return self._send(ARMS[arm].read_bytes(), "application/wasm")
            self.path = "/" + rest
        return super().do_GET()
ThreadingHTTPServer(("127.0.0.1", PORT), partial(Handler, directory=str(WEB))).serve_forever()
