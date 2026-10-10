import { describe, expect, it } from 'vitest';
import { ExportWriter, writeFile } from '../../src/application/export-writer';
import type { ExportSink } from '../../src/ports/download';
import { memoryDownloader, type SavedFile } from './support/memory-downloader';

const decoder = new TextDecoder();

function recordingSink(): { sink: ExportSink; blocks: string[] } {
  const blocks: string[] = [];
  return {
    blocks,
    sink: {
      write: async (bytes) => {
        blocks.push(decoder.decode(bytes));
      },
      close: async () => undefined,
      abort: () => undefined,
    },
  };
}

describe('ExportWriter (design §12.1)', () => {
  it('writes text in blocks of about blockChars characters, in order', async () => {
    const { sink, blocks } = recordingSink();
    const writer = new ExportWriter(sink, 10);
    for (let i = 0; i < 7; i++) await writer.text(`row ${i}\n`);
    await writer.flush();
    expect(blocks.join('')).toBe(Array.from({ length: 7 }, (_, i) => `row ${i}\n`).join(''));
    expect(blocks.length).toBeGreaterThan(1);
    for (const block of blocks.slice(0, -1)) expect(block.length).toBeGreaterThanOrEqual(10);
    expect(writer.written).toBe(blocks.join('').length);
  });

  it('writes bytes after the text before them, as they are', async () => {
    const { sink, blocks } = recordingSink();
    const writer = new ExportWriter(sink);
    await writer.text('a');
    await writer.bytes(new TextEncoder().encode('b'));
    await writer.text('c');
    await writer.flush();
    expect(blocks).toEqual(['a', 'b', 'c']);
  });

  it('counts UTF-8 bytes, not characters', async () => {
    const { sink } = recordingSink();
    const writer = new ExportWriter(sink);
    await writer.text('é');
    await writer.flush();
    expect(writer.written).toBe(2);
  });
});

describe('writeFile', () => {
  it('saves the file once it is written', async () => {
    const saved: SavedFile[] = [];
    const bytes = await writeFile(memoryDownloader((file) => saved.push(file)), 'x.txt', 'text/plain', async (writer) => {
      await writer.text('hello\n');
    });
    expect(bytes).toBe(6);
    expect(saved.map((file) => [file.name, file.mime, decoder.decode(file.bytes)])).toEqual([['x.txt', 'text/plain', 'hello\n']]);
  });

  it('saves nothing when writing fails, and throws the error again', async () => {
    const saved: SavedFile[] = [];
    const failure = writeFile(memoryDownloader((file) => saved.push(file)), 'x.txt', 'text/plain', async (writer) => {
      await writer.text('partial');
      await writer.flush();
      throw new Error('read failed');
    });
    await expect(failure).rejects.toThrow('read failed');
    expect(saved).toEqual([]);
  });
});
