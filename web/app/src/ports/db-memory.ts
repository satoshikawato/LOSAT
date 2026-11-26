// Persistence port for remembering that a previous run/session referenced a local BLAST
// database folder, so after a page reload the UI can prompt "Choose the folder again."
// (browsers do not retain directory handles across reloads). Only metadata (name/folder)
// is remembered — never sequence data.

export interface RememberedDatabase {
  readonly name: string;
  readonly folder?: string;
}

export type DatabaseSlot = 'subject' | 'filtering';

export interface DatabaseMemory {
  read(slot: DatabaseSlot): RememberedDatabase | undefined;
  write(slot: DatabaseSlot, value: RememberedDatabase | undefined): void;
}
