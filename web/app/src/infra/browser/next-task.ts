// A pause that lets the page draw between two steps of long work on its thread.

/**
 * Resolves in a task of its own, after the tasks already queued: a message to itself, which no
 * browser delays as it delays a timer set again and again (to 4 ms or more), so a pause between
 * short steps costs little while the browser can draw between them.
 */
export function nextTask(): Promise<void> {
  return new Promise((resolve) => {
    const channel = new MessageChannel();
    channel.port1.onmessage = () => {
      channel.port1.close();
      resolve();
    };
    channel.port2.postMessage(undefined);
  });
}
