// Stand-in for https://www.googletagmanager.com/gtag/js, served with the response headers
// recorded from it (google-headers.txt). It runs no Google code: it marks that it ran and
// sends one probe hit to C by each transport that gtag.js uses (sendBeacon, fetch with
// keepalive, an image).
self.__probeGtag = {
  ran: true,
  async send(url, { run, variant }) {
    const hit = (transport) => `${url}?v=2&tid=G-XXXXXXXXXX&en=page_view&probe_run=${run}&variant=${variant}&transport=${transport}`;
    const result = {};
    try {
      result.beacon = navigator.sendBeacon(hit('beacon'));
    } catch (e) {
      result.beacon = `error:${e.name}`;
    }
    result.fetch = await fetch(hit('fetch'), { method: 'POST', mode: 'no-cors', keepalive: true, credentials: 'include' }).then(
      () => 'resolved',
      (e) => `rejected:${e.name}`,
    );
    result.image = await new Promise((resolve) => {
      const image = new Image();
      image.onload = () => resolve('load');
      image.onerror = () => resolve('error');
      setTimeout(() => resolve('timeout'), 3000);
      image.src = hit('image');
    });
    return result;
  },
};
