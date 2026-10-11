// A: the stand-in for the app, served with the site's headers. Variant (b) posts the fixed
// event to its own /api/visit; variants a1-a4 embed M's document and talk to it by
// postMessage with origin and source checks. The observations go to window.__probe.
(() => {
  const q = new URLSearchParams(location.search);
  const variant = q.get('variant');
  const run = q.get('run');
  const M = q.get('m');
  const sandbox = q.get('sandbox') ?? 'none';
  const report = { variant, sandbox, crossOriginIsolated: self.crossOriginIsolated, cspViolations: [], messages: [], iframeLoadEvent: false, done: false };
  window.__probe = report;
  document.addEventListener('securitypolicyviolation', (e) => {
    report.cspViolations.push({ directive: e.effectiveDirective, blocked: e.blockedURI });
  });
  const finish = () => {
    report.done = true;
  };

  if (variant === 'b') {
    // The fixed event of variant (b): same origin, under COEP and the unchanged CSP.
    fetch(`/api/visit?run=${run}`, {
      method: 'POST',
      keepalive: true,
      credentials: 'omit',
      cache: 'no-store',
      referrerPolicy: 'no-referrer',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ v: 1, release: 'probe', client: crypto.randomUUID() }),
    }).then(
      (response) => {
        report.visit = { status: response.status };
      },
      (error) => {
        report.visit = { error: String(error) };
      },
    ).finally(finish);
    return;
  }

  const frame = document.createElement('iframe');
  if (sandbox === 'same-origin') frame.setAttribute('sandbox', 'allow-scripts allow-same-origin');
  if (sandbox === 'opaque') frame.setAttribute('sandbox', 'allow-scripts');
  if (variant === 'a4') {
    report.credentiallessAttributeSupported = 'credentialless' in frame;
    frame.setAttribute('credentialless', '');
  }
  frame.addEventListener('load', () => {
    report.iframeLoadEvent = true;
  });
  addEventListener('message', (e) => {
    const fromFrame = e.source === frame.contentWindow;
    const checksPass = e.origin === M && fromFrame;
    const type = e.data && e.data.type;
    report.messages.push({ type, origin: e.origin, fromFrame, checksPass });
    if (type === 'hello') {
      report.hello = { origin: e.origin, fromFrame, checksPass };
      if (checksPass) {
        frame.contentWindow.postMessage({ type: 'run', target: 'M' }, M);
      } else if (fromFrame) {
        // A sandboxed frame without allow-same-origin has an opaque origin: the origin check
        // refuses it, a post addressed to M's origin is dropped, and only '*' reaches it.
        frame.contentWindow.postMessage({ type: 'run', target: 'M' }, M);
        setTimeout(() => frame.contentWindow.postMessage({ type: 'run', target: '*' }, '*'), 300);
      }
    }
    if (type === 'result') {
      report.result = e.data.result;
      report.resultChecksPass = checksPass;
      finish();
    }
  });
  frame.src = `${M}/m.html?variant=${variant}&run=${run}`;
  document.body.append(frame);
  setTimeout(finish, 12000);
})();
