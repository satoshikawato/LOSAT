- Web tests: in WebKit, the system memory monitor is off (`WEBKIT_DISABLE_MEMORY_PRESSURE_MONITOR=1`), so a test
  does not depend on how much memory the rest of the machine uses. Under critical memory pressure (faked locally
  as 97 % used), the page's main thread sometimes waited forever in JavaScriptCore's response to the pressure, and
  a search never finished (3 of 18 local runs of the engine tests; the 15-minute test timeout). The application is
  unchanged.
