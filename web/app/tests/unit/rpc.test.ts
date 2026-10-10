import { afterEach, describe, expect, it } from 'vitest';
import { rpcClient, serveRpc, type RpcEndpoint } from '../../src/infra/data-worker/rpc';

interface Service {
  add(a: number, b: number): Promise<number>;
  fail(): Promise<never>;
  channel(): Promise<MessagePort>;
  bytes(): Promise<{ readonly data: Uint8Array }>;
  pieces(): Promise<{ readonly id: string; readonly residues: readonly Uint8Array[] }>;
  hidden(): Promise<string>;
}

class Impl implements Service {
  async add(a: number, b: number) {
    return a + b;
  }
  async fail(): Promise<never> {
    const error = new Error('Not enough temporary storage');
    error.name = 'StorageFullError';
    throw error;
  }
  async channel() {
    const channel = new MessageChannel();
    channel.port2.onmessage = (event) => channel.port2.postMessage(`echo ${String(event.data)}`);
    return channel.port1;
  }
  /** The arrays that the last `bytes` and `pieces` returned, to see whether they were transferred. */
  readonly returned: Uint8Array[] = [];
  async bytes() {
    const data = new Uint8Array([1, 2, 3]);
    this.returned.push(data);
    return { data };
  }
  async pieces() {
    const residues = [new Uint8Array([65, 67]), new Uint8Array([71]), new Uint8Array(new ArrayBuffer(8), 2, 3)];
    this.returned.push(...residues);
    return { id: 'x', residues };
  }
  async hidden() {
    return 'secret';
  }
}

const open: MessagePort[] = [];
afterEach(() => {
  for (const port of open.splice(0)) port.close();
});

function connect(service: Service | Promise<Service> = new Impl()) {
  const channel = new MessageChannel();
  open.push(channel.port1, channel.port2);
  serveRpc(channel.port2, service, ['add', 'fail', 'channel', 'bytes', 'pieces']);
  return rpcClient<Service>(channel.port1, ['add', 'fail', 'channel', 'bytes', 'pieces', 'hidden']);
}

describe('rpc', () => {
  it('calls the served methods and returns their results', async () => {
    const client = connect();
    expect(await client.add(2, 3)).toBe(5);
    expect(await client.bytes()).toEqual({ data: new Uint8Array([1, 2, 3]) });
  });

  it('transfers the typed arrays of a result and of an array in it when they own their buffers', async () => {
    const impl = new Impl();
    const client = connect(impl);
    expect(await client.bytes()).toEqual({ data: new Uint8Array([1, 2, 3]) });
    const got = await client.pieces();
    expect(got.id).toBe('x');
    expect(got.residues.map((bytes) => [...bytes])).toEqual([[65, 67], [71], [0, 0, 0]]);
    const [data, first, second, view] = impl.returned;
    expect([data!.byteLength, first!.byteLength, second!.byteLength]).toEqual([0, 0, 0]);
    // A view into a larger buffer is copied, so the served object's buffer stays usable.
    expect(view!.byteLength).toBe(3);
  });

  it('keeps the name and message of an error', async () => {
    const error = await connect()
      .fail()
      .catch((e: unknown) => e as Error);
    expect(error.name).toBe('StorageFullError');
    expect(error.message).toBe('Not enough temporary storage');
  });

  it('serves only the listed methods', async () => {
    await expect(connect().hidden()).rejects.toThrow('unknown method hidden');
  });

  it('transfers a returned MessagePort', async () => {
    const port = await connect().channel();
    open.push(port);
    const reply = new Promise((resolve) => {
      port.onmessage = (event) => resolve(event.data);
    });
    port.postMessage('hello');
    expect(await reply).toBe('echo hello');
  });

  it('rejects every call when the service fails to start', async () => {
    const client = connect(Promise.reject(new Error('no storage')));
    await expect(client.add(1, 1)).rejects.toThrow('no storage');
  });

  it('keeps working after an error event of a worker that has answered', async () => {
    const target = new EventTarget();
    const endpoint: RpcEndpoint = {
      postMessage: (message) => {
        const { id } = message as { id: number };
        queueMicrotask(() => target.dispatchEvent(new MessageEvent('message', { data: { id, ok: true, value: 3 } })));
      },
      addEventListener: (type, listener) => target.addEventListener(type, listener),
    };
    const client = rpcClient<Service>(endpoint, ['add']);
    expect(await client.add(1, 2)).toBe(3);
    target.dispatchEvent(new Event('error'));
    expect(await client.add(1, 2)).toBe(3);
  });

  it('rejects pending and later calls when the worker fails to start', async () => {
    const target = new EventTarget();
    const endpoint: RpcEndpoint = {
      postMessage: () => undefined,
      addEventListener: (type, listener) => target.addEventListener(type, listener),
    };
    const client = rpcClient<Service>(endpoint, ['add']);
    const pending = client.add(1, 2);
    target.dispatchEvent(new Event('error'));
    await expect(pending).rejects.toThrow('the worker stopped');
    await expect(client.add(1, 2)).rejects.toThrow('the worker stopped');
  });
});
