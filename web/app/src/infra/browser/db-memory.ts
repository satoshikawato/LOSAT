// sessionStorage-backed implementation of DatabaseMemory. Wrapped in try/catch because
// browsers may block storage access in private windows or restricted iframes.
import type { DatabaseMemory, DatabaseSlot, RememberedDatabase } from '../../ports/db-memory';

const PREFIX = 'losat.web.db.';

export function browserDatabaseMemory(): DatabaseMemory {
  return {
    read(slot: DatabaseSlot): RememberedDatabase | undefined {
      try {
        const raw = window.sessionStorage.getItem(`${PREFIX}${slot}`);
        if (raw === null) return undefined;
        const parsed = JSON.parse(raw) as unknown;
        if (typeof parsed !== 'object' || parsed === null) return undefined;
        const name = (parsed as { name?: unknown }).name;
        const folder = (parsed as { folder?: unknown }).folder;
        if (typeof name !== 'string' || name.trim() === '') return undefined;
        return typeof folder === 'string' && folder !== '' ? { name, folder } : { name };
      } catch {
        return undefined;
      }
    },
    write(slot: DatabaseSlot, value: RememberedDatabase | undefined): void {
      try {
        const key = `${PREFIX}${slot}`;
        if (value === undefined) window.sessionStorage.removeItem(key);
        else window.sessionStorage.setItem(key, JSON.stringify(value));
      } catch {
        // Storage blocked: degrade silently without affecting searches.
      }
    },
  };
}
