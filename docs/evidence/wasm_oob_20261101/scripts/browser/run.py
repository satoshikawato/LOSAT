#!/usr/bin/env python3
"""Line WO: the V8-only stale-atomic test in browsers. usage: run.py [rounds] [fences]"""
import functools, http.server, json, pathlib, sys, threading
from playwright.sync_api import sync_playwright

rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
fences = int(sys.argv[2]) if len(sys.argv) > 2 else 100000

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
url = f"http://127.0.0.1:{server.server_address[1]}/index.html"
with sync_playwright() as p:
    home = pathlib.Path.home() / ".cache/ms-playwright"
    paths = {"firefox": home / "firefox-1532/firefox/firefox", "webkit": home / "webkit-2311/pw_run.sh"}
    for name in ("chromium", "firefox"):
        try:
            options = {"executable_path": str(paths[name])} if name in paths else {}
            browser = getattr(p, name).launch(headless=True, **options)
            page = browser.new_page()
            page.goto(url)
            page.set_default_timeout(60000)
            result = page.evaluate("([r, f]) => run(r, f)", [rounds, fences])
            print(name, browser.version, json.dumps(result), flush=True)
            browser.close()
        except Exception as err:
            print(name, "FAILED", str(err).splitlines()[0][:200], flush=True)
server.shutdown()
