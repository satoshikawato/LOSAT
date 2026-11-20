- Web (WebKit): a search on more than one thread that is cancelled in its first moments, or while the page is still
  warming up the engine for it, no longer ends the engine's runtime without first growing its memory, the step that
  keeps WebKit from crashing the page. The page now knows at every moment whether the engine is working on that
  memory, instead of learning it from a message that could arrive too late or, during the warm-up, not at all.
