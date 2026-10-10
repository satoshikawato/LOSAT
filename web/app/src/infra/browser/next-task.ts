// A pause that lets the page draw between two steps of long work on its thread.

/** How often a pause also waits for a timer, at least. */
const TIMER_EVERY_MS = 40;
/** When a pause last waited for a timer. */
let lastTimer = 0;

/**
 * Resolves in a task of its own, after the tasks already queued: a message to itself, which no
 * browser delays as it delays a timer set again and again (to 4 ms or more), so a pause between
 * short steps costs little. At least every TIMER_EVERY_MS the pause also waits for a timer (set
 * from the message's task, so not delayed): WebKit runs posted messages one after another before
 * its timers and its drawing (fix round 2: a chain of message pauses held its page 120-160 ms).
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
  if (performance.now() - lastTimer >= TIMER_EVERY_MS) {
    await new Promise<void>((resolve) => setTimeout(resolve, 0));
    lastTimer = performance.now();
  }
}
