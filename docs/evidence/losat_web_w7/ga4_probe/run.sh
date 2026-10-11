#!/bin/sh
# Runs the GA4 probe (RESULTS.md) in Chromium, Firefox and WebKit with web/app's
# node_modules (npm ci there first). LOSAT_WEB_WEBKIT_EXECUTABLE may name a WebKit launcher
# (web/app/playwright.config.ts). Extra arguments go to `playwright test`.
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"
APP="$HERE/../../../../web/app"
cd "$HERE"
NODE_PATH="$APP/node_modules" exec "$APP/node_modules/.bin/playwright" test -c "$HERE/playwright.config.ts" "$@"
