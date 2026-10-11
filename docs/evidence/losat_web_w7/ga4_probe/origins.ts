// The local origins of the GA4 probe (RESULTS.md), on the loopback interface and all
// potentially trustworthy (secure contexts). `localhost` and `127.0.0.1` are different
// sites, so M's frame in A is cross-site as it would be on the real hosts, and G and C are
// cross-site to M, where they are used. (An IPv6 literal would give a third site, but a CSP
// source cannot name one, so G and C share A's site; A never loads them.)
//   A  the app (LOSAT Web's response headers, read from web/app/public/_headers)
//   M  the measurement document's origin (variant (a))
//   G  a stand-in for www.googletagmanager.com (gtag.js with the recorded response headers)
//   C  a stand-in for the collect endpoint (www.google-analytics.com/g/collect); it shares
//      G's server and site, and only records what reaches it
export const PORT_A = Number(process.env['PROBE_PORT_A'] || 4762);
export const PORT_M = Number(process.env['PROBE_PORT_M'] || 4862);
export const PORT_G = Number(process.env['PROBE_PORT_G'] || 4962);

export const A = `http://localhost:${PORT_A}`;
export const M = `http://127.0.0.1:${PORT_M}`;
export const G = `http://localhost:${PORT_G}`;
export const C = G;
