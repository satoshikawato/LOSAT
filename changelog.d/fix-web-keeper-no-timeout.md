- Web (WebKit): when a search on more than one thread is cancelled, or its engine stops with an error, the page
  now waits as long as needed for the step that keeps WebKit from crashing the page (the memory is grown to its
  maximum) before it ends the search's runtime. After a cancel it used to wait at most 1 s, which a busy device
  could exceed, and when the engine stopped it did not wait at all; the runtime then ended without that step.
