#!/usr/bin/env python3
"""Line WO: the plain store after an atomic read in Chromium, stack checks off, as in the guard's
tests. On Linux ARM64 the original traps; on x64 it returns 1 (the stale size). usage: run_refresh.py [rounds]"""
import functools, http.server, json, pathlib, sys, threading
from playwright.sync_api import sync_playwright

rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 10

class Handler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cache-Control", "no-store")
        super().end_headers()
    def log_message(self, *args):
        pass

here = pathlib.Path(__file__).parent
server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Handler, directory=str(here)))
threading.Thread(target=server.serve_forever, daemon=True).start()
with sync_playwright() as p:
    for label, args in (("stack checks off", ["--js-flags=--no-wasm-stack-checks"]), ("default flags", [])):
        browser = p.chromium.launch(headless=True, args=args)
        page = browser.new_page()
        page.goto(f"http://127.0.0.1:{server.server_address[1]}/refresh.html")
        page.set_default_timeout(600000)
        print(label, browser.version, json.dumps(page.evaluate("(r) => run(r)", rounds)), flush=True)
        browser.close()
server.shutdown()
