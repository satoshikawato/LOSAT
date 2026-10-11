// M: the stand-in for the measurement document of variant (a) (server.mjs prefixes
// self.PROBE = {A, G, C}). In a frame it says hello to A, accepts a 'run' only from A's
// origin and from its parent window, then keeps an identifier (localStorage and cookies),
// loads G's gtag stand-in, lets it send a probe hit to C, and answers A. Opened top-level
// with mode=read, it only reads its first-party storage (to detect partitioning).
(() => {
  const { A, G, C } = self.PROBE;
  const q = new URLSearchParams(location.search);
  const run = q.get('run');
  const variant = q.get('variant');
  const top = window.parent === window;
  const out = { variant, top, origin: self.origin, crossOriginIsolated: self.crossOriginIsolated, cspViolations: [], runMessages: [] };
  document.addEventListener('securitypolicyviolation', (e) => {
    out.cspViolations.push({ directive: e.effectiveDirective, blocked: e.blockedURI });
  });

  const cookies = () => Object.fromEntries(document.cookie.split('; ').filter(Boolean).map((c) => c.split('=')));
  function storage(write) {
    const s = {};
    try {
      s.localBefore = localStorage.getItem('probe-id');
      if (write && s.localBefore === null) localStorage.setItem('probe-id', crypto.randomUUID());
      s.localAfter = localStorage.getItem('probe-id');
    } catch (e) {
      s.localError = e.name;
    }
    try {
      s.cookieBefore = cookies();
      if (write) {
        const id = crypto.randomUUID();
        // gtag's _ga is a first-party cookie of the measurement document; without
        // SameSite=None a browser treats it as Lax, which a cross-site frame cannot set.
        if (!s.cookieBefore.probe_none) document.cookie = `probe_none=${id}; SameSite=None; Secure; Max-Age=86400; Path=/`;
        if (!s.cookieBefore.probe_lax) document.cookie = `probe_lax=${id}; Max-Age=86400; Path=/`;
      }
      s.cookieAfter = cookies();
    } catch (e) {
      s.cookieError = e.name;
    }
    return s;
  }
  async function storageAccess() {
    try {
      return typeof document.hasStorageAccess === 'function' ? await document.hasStorageAccess() : 'unsupported';
    } catch (e) {
      return `error:${e.name}`;
    }
  }
  function loadScript() {
    return new Promise((resolve) => {
      const script = document.createElement('script');
      script.src = `${G}/gtag/js?id=G-XXXXXXXXXX&run=${run}`;
      script.onload = () => resolve({ loaded: true, ran: !!self.__probeGtag });
      script.onerror = () => resolve({ loaded: false, ran: !!self.__probeGtag });
      setTimeout(() => resolve({ loaded: 'timeout', ran: !!self.__probeGtag }), 5000);
      document.head.append(script);
    });
  }

  if (top) {
    storageAccess().then((access) => {
      out.storage = storage(false);
      out.hasStorageAccess = access;
      out.done = true;
      window.__probe = out;
    });
    return;
  }

  addEventListener('message', async (e) => {
    const ok = e.origin === A && e.source === window.parent;
    out.runMessages.push({ origin: e.origin, fromParent: e.source === window.parent, ok, target: e.data && e.data.target });
    if (!ok || out.started) return;
    out.started = true;
    out.storage = storage(true);
    out.hasStorageAccess = await storageAccess();
    out.script = await loadScript();
    out.transports = self.__probeGtag ? await self.__probeGtag.send(`${C}/g/collect`, { run, variant }) : null;
    // Lets a second 'run' (the opaque frame's '*' post) arrive before answering.
    await new Promise((resolve) => setTimeout(resolve, 700));
    window.parent.postMessage({ type: 'result', result: out }, A);
  });
  window.parent.postMessage({ type: 'hello' }, A);
})();
