// A minimal request/response channel between the main thread and a worker. The worker
// serves the methods of one object; the main thread gets an object with the same methods,
// each of which returns a promise. Errors keep their name and message.
//
// Transfer rule for results: a MessagePort, and a typed array (a Uint8Array, or a column of
// an HSP table) that owns its whole buffer (at the top level, as a property of the result, or
// as an element of an array that is a property of the result, such as the residues of
// `readResidues`), are transferred instead of copied. A served object must therefore return
// buffers that it does not keep.

export interface RpcEndpoint {
  postMessage(message: unknown, transfer: Transferable[]): void;
  addEventListener(type: string, listener: (event: Event) => void): void;
  start?(): void;
}

interface RpcRequest {
  readonly id: number;
  readonly method: string;
  readonly args: readonly unknown[];
}

type RpcResponse =
  | { readonly id: number; readonly ok: true; readonly value: unknown }
  | { readonly id: number; readonly ok: false; readonly error: { readonly name: string; readonly message: string } };

type Methods<T> = { [K in keyof T]: T[K] extends (...args: never[]) => Promise<unknown> ? K : never }[keyof T] &
  string;

/** Serves `methods` of `service` (which may still be starting) on `endpoint`. */
export function serveRpc<T extends object>(
  endpoint: RpcEndpoint,
  service: T | Promise<T>,
  methods: readonly Methods<T>[],
): void {
  const allowed = new Set<string>(methods);
  const ready = Promise.resolve(service);
  // A service that fails to start is reported to every call, not as an unhandled rejection.
  ready.catch(() => undefined);
  endpoint.addEventListener('message', (event) => {
    const request = (event as MessageEvent<RpcRequest>).data;
    void (async () => {
      try {
        if (!allowed.has(request.method)) throw new Error(`unknown method ${request.method}`);
        const target = await ready;
        const method = (target as Record<string, (...args: unknown[]) => unknown>)[request.method]!;
        const value = await method.apply(target, [...request.args]);
        const response: RpcResponse = { id: request.id, ok: true, value };
        endpoint.postMessage(response, transferables(value));
      } catch (error) {
        const response: RpcResponse = {
          id: request.id,
          ok: false,
          error: {
            name: error instanceof Error ? error.name : 'Error',
            message: error instanceof Error ? error.message : String(error),
          },
        };
        endpoint.postMessage(response, []);
      }
    })();
  });
  endpoint.start?.();
}

/**
 * Returns an object whose `methods` call the served object behind `endpoint`. If the
 * endpoint reports an error before its first answer (a worker that failed to start),
 * every pending and later call rejects with it.
 */
export function rpcClient<T extends object>(endpoint: RpcEndpoint, methods: readonly Methods<T>[]): T {
  let nextId = 0;
  let broken: Error | undefined;
  const pending = new Map<number, { resolve(value: unknown): void; reject(error: Error): void }>();
  let answered = false;
  endpoint.addEventListener('message', (event) => {
    answered = true;
    const response = (event as MessageEvent<RpcResponse>).data;
    const call = pending.get(response.id);
    if (call === undefined) return;
    pending.delete(response.id);
    if (response.ok) {
      call.resolve(response.value);
    } else {
      const error = new Error(response.error.message);
      error.name = response.error.name;
      call.reject(error);
    }
  });
  endpoint.addEventListener('error', (event) => {
    // An error before the first answer means the worker failed to start. Later errors are
    // uncaught exceptions of a running worker; the served calls report their own errors.
    if (answered) return;
    const detail = (event as ErrorEvent).message;
    broken = new Error(`the worker stopped${detail ? `: ${detail}` : ''}`);
    for (const call of pending.values()) call.reject(broken);
    pending.clear();
  });
  endpoint.start?.();
  const client: Record<string, (...args: unknown[]) => Promise<unknown>> = {};
  for (const method of methods) {
    client[method] = (...args) =>
      new Promise((resolve, reject) => {
        if (broken !== undefined) {
          reject(broken);
          return;
        }
        const id = ++nextId;
        pending.set(id, { resolve, reject });
        const request: RpcRequest = { id, method, args };
        endpoint.postMessage(request, []);
      });
  }
  return client as T;
}

function transferables(value: unknown): Transferable[] {
  const found: Transferable[] = [];
  const add = (item: unknown) => {
    if (item instanceof MessagePort) found.push(item);
    else if (
      ArrayBuffer.isView(item) &&
      !(item instanceof DataView) &&
      item.buffer instanceof ArrayBuffer &&
      item.byteOffset === 0 &&
      item.byteLength === item.buffer.byteLength &&
      !found.includes(item.buffer)
    ) {
      found.push(item.buffer);
    }
  };
  add(value);
  if (typeof value === 'object' && value !== null && !ArrayBuffer.isView(value)) {
    for (const item of Object.values(value)) {
      add(item);
      if (Array.isArray(item)) for (const element of item) add(element);
    }
  }
  return found;
}
