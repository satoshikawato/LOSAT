// A pause that lets the page draw between two steps of long work on its thread.

/** How often a pause also waits for a timer, at least. */
const TIMER_EVERY_MS = 40;
/** How often a pause looks whether the browser has drawn a frame since it last looked. */
const FRAME_EVERY_MS = 50;
/** The longest wait for a frame: a page in the background may draw none. */
const FRAME_WAIT_MS = 100;

let lastTimer = 0;
let lastLook = 0;
let drawn = true;

/**
 * Resolves in a task of its own, after the tasks already queued: a message to itself, which no
 * browser delays as it delays a timer set again and again (to 4 ms or more), so a pause between
 * short steps costs little. At least every TIMER_EVERY_MS the pause also waits for a timer (set
 * from the message's task, so not delayed), and when the browser has drawn no frame for
 * FRAME_EVERY_MS it waits for one: WebKit runs ready messages and timers before its drawing, so a
 * chain of short steps held its page 150-230 ms (fix round 2). Chromium and Firefox draw between
 * the steps, and wait for nothing more.
 */
export async function nextTask(): Promise<void> {
  await new Promise<void>((resolve) => {
    const channel = new MessageChannel();
    channel.port1.onmessage = () => {
      channel.port1.close();
      resolve();
    };
    channel.port2.postMessage(undefined);
  });
  const now = performance.now();
  if (now - lastTimer >= TIMER_EVERY_MS) {
    await new Promise<void>((resolve) => setTimeout(resolve, 0));
    lastTimer = performance.now();
  }
  if (typeof requestAnimationFrame !== 'function' || now - lastLook < FRAME_EVERY_MS) return;
  if (!drawn && document.visibilityState === 'visible') {
    // Go on after the frame is drawn: a timer set in the frame's callback fires after its painting.
    await new Promise<void>((resolve) => {
      const timeout = setTimeout(resolve, FRAME_WAIT_MS);
      requestAnimationFrame(() =>
        setTimeout(() => {
          clearTimeout(timeout);
          resolve();
        }, 0),
      );
    });
  }
  drawn = false;
  lastLook = performance.now();
  requestAnimationFrame(() => (drawn = true));
}
