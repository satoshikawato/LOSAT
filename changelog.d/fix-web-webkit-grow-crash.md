- Web: cancelling a search (or a search that fails) while its threads run no longer crashes the tab in
  WebKit-based browsers. JavaScriptCore crashed the page when a thread grew the shared memory while the
  workers that held instances of it were ending; the main thread now grows the memory to its maximum
  before it ends the Engine worker during a search, so no thread grows it while they end.
