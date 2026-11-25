// Local BLAST database discovery, molecule-type matching, filter option metadata and
// run-directory assembly (WQA D12b + RG-D4 W2). Single owner reused by Subject "Database" (-db)
// and BLASTN "Repeats database" (-filtering_db).
import type { DirectoryFile, SequenceKind } from 'losat/domain';

export const DB_FLAG = '-db';
export const FILTERING_DB_FLAG = '-filtering_db';
export const TAXDB_INDEX_FILE = 'taxdb.bti';
export const TAXDB_DATA_FILE = 'taxdb.btd';
export const TAXONOMY_SQLITE_FILE = 'taxonomy4blast.sqlite3';
const AUX_TAXONOMY_FILES = [TAXDB_INDEX_FILE, TAXDB_DATA_FILE, TAXONOMY_SQLITE_FILE] as const;

/** Safety cap for reading small `.nal`/`.pal` alias text files during folder discovery. */
const MAX_ALIAS_BYTES = 256 * 1024;

export interface DatabaseFilterField {
  readonly flag: string;
  readonly label: string;
  readonly kind: 'file' | 'text' | 'flag';
  readonly placeholder?: string;
}

/** Database filter options whose value is a list file mounted by name into the run directory. */
export const DB_FILE_FILTER_FIELDS: readonly DatabaseFilterField[] = Object.freeze([
  { flag: '-seqidlist', label: 'Include SeqIDs list', kind: 'file' },
  { flag: '-negative_seqidlist', label: 'Exclude SeqIDs list', kind: 'file' },
  { flag: '-taxidlist', label: 'Include taxonomy IDs list', kind: 'file' },
  { flag: '-negative_taxidlist', label: 'Exclude taxonomy IDs list', kind: 'file' },
  { flag: '-gilist', label: 'Include GI list', kind: 'file' },
  { flag: '-negative_gilist', label: 'Exclude GI list', kind: 'file' },
  { flag: '-ipglist', label: 'Include IPG list', kind: 'file' },
  { flag: '-negative_ipglist', label: 'Exclude IPG list', kind: 'file' },
]);

/** Database filter options that take text values (comma-separated taxids or masking algorithm IDs). */
export const DB_TEXT_FILTER_FIELDS: readonly DatabaseFilterField[] = Object.freeze([
  { flag: '-taxids', label: 'Include taxonomy IDs', kind: 'text', placeholder: '9606,10090' },
  { flag: '-negative_taxids', label: 'Exclude taxonomy IDs', kind: 'text', placeholder: '9606' },
  { flag: '-db_soft_mask', label: 'Database soft mask ID', kind: 'text', placeholder: '11' },
  { flag: '-db_hard_mask', label: 'Database hard mask ID', kind: 'text', placeholder: '11' },
]);

/** Bare boolean flags for database filtering. */
export const DB_BOOL_FILTER_FIELDS: readonly DatabaseFilterField[] = Object.freeze([
  { flag: '-no_taxid_expansion', label: 'Do not expand taxonomy IDs to descendants', kind: 'flag' },
]);

export const ALL_DB_FILTER_FIELDS: readonly DatabaseFilterField[] = Object.freeze([
  ...DB_TEXT_FILTER_FIELDS,
  ...DB_BOOL_FILTER_FIELDS,
  ...DB_FILE_FILTER_FIELDS,
]);

export const ALL_DB_FILTER_FLAGS: ReadonlySet<string> = new Set(ALL_DB_FILTER_FIELDS.map((field) => field.flag));
export const DB_FILE_FILTER_FLAGS: ReadonlySet<string> = new Set(DB_FILE_FILTER_FIELDS.map((field) => field.flag));

export interface BlastDbCandidate {
  /** Exact relative database path passed to `-db` or `-filtering_db` (without `.nal`/`.pal`/`.nin`/`.pin`). */
  readonly name: string;
  readonly kind: SequenceKind;
  readonly entry: 'alias' | 'volume';
  /** Optional `TITLE` line from `.nal`/`.pal`. */
  readonly title?: string;
  /** Component volumes listed in `DBLIST` (for alias files). */
  readonly volumes?: readonly string[];
}

export interface BlastDbDiscovery {
  /** Top-level folder name when chosen via folder picker (`webkitRelativePath`), else undefined. */
  readonly folder?: string;
  /** Valid BLAST databases discovered (with alias-component volumes hidden), sorted by name. */
  readonly databases: readonly BlastDbCandidate[];
  /** All safe picked files, ready for zero-copy mounting into `RunDirectory`. */
  readonly files: readonly DirectoryFile[];
  /** Total bytes across all mounted files. */
  readonly totalBytes: number;
  /** Auxiliary taxonomy files detected in the picked folder/files. */
  readonly taxonomy: {
    readonly taxdb: boolean;
    readonly sqlite: boolean;
  };
}

/**
 * Safe relative path for a picked File, using `/` separators and rejecting empty, `.`, `..`,
 * or NUL-containing segments (matching `RunDirectory`'s path rules in `packages/losat`).
 */
export function relativeFilePath(file: File): string | undefined {
  const rel = (file as File & { readonly webkitRelativePath?: string }).webkitRelativePath;
  const raw = (rel && rel.trim() !== '' ? rel : file.name).replaceAll('\\', '/').replace(/^\/+/, '');
  if (raw === '' || raw.includes('\0')) return undefined;
  const parts = raw.split('/').filter((segment) => segment !== '');
  if (parts.length === 0 || parts.some((segment) => segment === '.' || segment === '..')) return undefined;
  return parts.join('/');
}

/** Parses `TITLE` and `DBLIST` tokens from a BLAST `.nal` or `.pal` alias file. */
export function parseAliasMetadata(text: string): { readonly title?: string; readonly dblist: readonly string[] } {
  let title: string | undefined;
  const dblist: string[] = [];
  for (const rawLine of text.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (line === '' || line.startsWith('#')) continue;
    if (line.startsWith('TITLE')) {
      const value = line.slice('TITLE'.length).trim();
      if (value !== '') title = value;
      continue;
    }
    if (!line.startsWith('DBLIST')) continue;
    const rest = line.slice('DBLIST'.length).trim();
    let i = 0;
    while (i < rest.length) {
      while (i < rest.length && /\s/.test(rest[i]!)) i++;
      if (i >= rest.length) break;
      let token: string;
      if (rest[i] === '"') {
        const end = rest.indexOf('"', i + 1);
        token = end < 0 ? rest.slice(i + 1) : rest.slice(i + 1, end);
        i = end < 0 ? rest.length : end + 1;
      } else {
        const start = i;
        while (i < rest.length && !/\s/.test(rest[i]!)) i++;
        token = rest.slice(start, i);
      }
      if (token !== '') dblist.push(token.replaceAll('\\', '/'));
    }
  }
  return title === undefined ? { dblist } : { title, dblist };
}

function parentDir(path: string): string {
  const slash = path.lastIndexOf('/');
  return slash < 0 ? '' : path.slice(0, slash);
}

function baseName(path: string): string {
  const slash = path.lastIndexOf('/');
  return slash < 0 ? path : path.slice(slash + 1);
}

function joinPath(dir: string, rel: string): string {
  const clean = rel.replace(/^\/+/, '');
  return dir === '' ? clean : `${dir}/${clean}`;
}

/**
 * Discovers valid BLAST databases (`.nal`/`.pal`, `.nin`/`.pin`, multi-volume sets) and
 * auxiliary taxonomy files (`taxdb.bti`/`taxdb.btd`, `taxonomy4blast.sqlite3`) from picked files.
 * Reads only tiny `.nal`/`.pal` alias text files; binary volume files are never read into memory.
 */
export async function discoverBlastDatabases(picked: readonly File[]): Promise<BlastDbDiscovery> {
  const seenPaths = new Set<string>();
  const files: DirectoryFile[] = [];
  let totalBytes = 0;
  let folder: string | undefined;
  let allShareFolder = picked.length > 0;

  for (const file of picked) {
    const path = relativeFilePath(file);
    if (path === undefined || seenPaths.has(path)) continue;
    seenPaths.add(path);
    files.push({ name: path, file });
    totalBytes += file.size;
    const slash = path.indexOf('/');
    if (slash > 0) {
      const top = path.slice(0, slash);
      if (folder === undefined) folder = top;
      else if (folder !== top) allShareFolder = false;
    } else {
      allShareFolder = false;
    }
  }
  if (!allShareFolder) folder = undefined;

  // Detect auxiliary taxonomy files anywhere in the pick.
  const basenames = new Set(files.map((entry) => baseName(entry.name)));
  const taxdb = basenames.has(TAXDB_INDEX_FILE) && basenames.has(TAXDB_DATA_FILE);
  const sqlite = basenames.has(TAXONOMY_SQLITE_FILE);

  // Collect `.nal`, `.nin`, `.pal`, `.pin` entries (case-sensitive, matching NCBI FindBlastDBs).
  interface RawEntry {
    readonly name: string;
    readonly kind: SequenceKind;
    readonly entry: 'alias' | 'volume';
    readonly file: Blob;
  }
  const byKey = new Map<string, RawEntry>();
  for (const entry of files) {
    const ext = entry.name.length > 4 ? entry.name.slice(-4) : '';
    let kind: SequenceKind | undefined;
    let kindEntry: 'alias' | 'volume' | undefined;
    if (ext === '.nal') {
      kind = 'nucleotide';
      kindEntry = 'alias';
    } else if (ext === '.nin') {
      kind = 'nucleotide';
      kindEntry = 'volume';
    } else if (ext === '.pal') {
      kind = 'protein';
      kindEntry = 'alias';
    } else if (ext === '.pin') {
      kind = 'protein';
      kindEntry = 'volume';
    }
    if (kind === undefined || kindEntry === undefined) continue;
    const stem = entry.name.slice(0, -4);
    const key = `${kind}:${stem}`;
    const existing = byKey.get(key);
    // Prefer alias over index when both exist for the same stem.
    if (existing === undefined || (existing.entry === 'volume' && kindEntry === 'alias')) {
      byKey.set(key, { name: stem, kind, entry: kindEntry, file: entry.file });
    }
  }

  // Read alias metadata (`TITLE`, `DBLIST`) to hide component volumes of multi-volume aliases.
  const aliasMeta = new Map<string, { readonly title?: string; readonly dblist: readonly string[] }>();
  const hiddenComponents = new Set<string>();
  const aliasStems = new Set<string>();

  for (const [key, item] of byKey) {
    if (item.entry !== 'alias') continue;
    aliasStems.add(key);
    if (item.file.size > MAX_ALIAS_BYTES) continue;
    try {
      const text = await item.file.text();
      const meta = parseAliasMetadata(text);
      aliasMeta.set(key, meta);
      const dir = parentDir(item.name);
      for (const token of meta.dblist) {
        hiddenComponents.add(`${item.kind}:${token}`);
        hiddenComponents.add(`${item.kind}:${joinPath(dir, token)}`);
      }
    } catch {
      // Unreadable alias file is still surfaced so the engine can report its own error if chosen.
    }
  }

  const databases: BlastDbCandidate[] = [];
  for (const [key, item] of byKey) {
    if (item.entry === 'volume') {
      if (hiddenComponents.has(key)) continue;
      // Also hide numbered volume stems (`foo.00`) when an alias `foo` of the same molecule type exists.
      const numbered = /^(.*)\.\d{2,}$/.exec(item.name);
      if (numbered !== null && aliasStems.has(`${item.kind}:${numbered[1]!}`)) continue;
    }
    const meta = aliasMeta.get(key);
    databases.push(
      Object.freeze({
        name: item.name,
        kind: item.kind,
        entry: item.entry,
        ...(meta?.title === undefined ? {} : { title: meta.title }),
        ...(meta !== undefined && meta.dblist.length > 0 ? { volumes: Object.freeze([...meta.dblist]) } : {}),
      }),
    );
  }

  databases.sort((a, b) => (a.name < b.name ? -1 : a.name > b.name ? 1 : a.kind < b.kind ? -1 : a.kind > b.kind ? 1 : 0));

  return Object.freeze({
    ...(folder === undefined ? {} : { folder }),
    databases: Object.freeze(databases),
    files: Object.freeze(files),
    totalBytes,
    taxonomy: Object.freeze({ taxdb, sqlite }),
  });
}

/** Databases from `discovery` matching the required molecule type. */
export function matchingDatabases(discovery: BlastDbDiscovery | undefined, kind: SequenceKind): readonly BlastDbCandidate[] {
  if (discovery === undefined) return [];
  return discovery.databases.filter((db) => db.kind === kind);
}

/** Clear user-facing message when picked files hold no BLAST database of the required molecule type. */
export function noMatchingDatabaseMessage(discovery: BlastDbDiscovery, kind: SequenceKind, contextLabel: string): string {
  const otherKind: SequenceKind = kind === 'nucleotide' ? 'protein' : 'nucleotide';
  const others = discovery.databases.filter((db) => db.kind === otherKind);
  if (others.length > 0) {
    const shown = others
      .slice(0, 3)
      .map((db) => db.name)
      .join(', ');
    const more = others.length > 3 ? ` and ${others.length - 3} more` : '';
    return `Found ${otherKind} database (${shown}${more}), not a ${kind} database for ${contextLabel}.`;
  }
  const exts = kind === 'nucleotide' ? '.nal or .nin' : '.pal or .pin';
  return `No ${kind} BLAST database (${exts}) was found in the chosen files.`;
}

/**
 * Builds the deduplicated `DirectoryFile[]` mounted into `RunDirectory` for a search run:
 * - all picked files of the subject DB (with `taxdb.bti`, `taxdb.btd`, `taxonomy4blast.sqlite3`
 *   also mounted at root so NCBI's `SeqDB_ResolveDbPath` finds them in `/work`)
 * - all picked files of BLASTN's `-filtering_db` repeats DB
 * - all picked DB filter list files (`-seqidlist`, `-taxidlist`, etc.), mounted at root by `file.name`
 *
 * Prevents duplicate paths and file-vs-directory collisions that `RunDirectory.open` rejects.
 */
export function assembleSearchDirectory(parts: {
  readonly subjectDb?: { readonly discovery: BlastDbDiscovery; readonly selectedName: string };
  readonly filterFiles?: Readonly<Record<string, File>>;
  readonly filteringDb?: { readonly discovery: BlastDbDiscovery; readonly selectedName: string };
}): readonly DirectoryFile[] {
  const entries = new Map<string, Blob>();
  const directories = new Set<string>();

  const parentPrefixes = (name: string): string[] => {
    const segments = name.split('/');
    const prefixes: string[] = [];
    let prefix = '';
    for (let i = 0; i < segments.length - 1; i++) {
      prefix = prefix === '' ? segments[i]! : `${prefix}/${segments[i]!}`;
      prefixes.push(prefix);
    }
    return prefixes;
  };

  const canAdd = (name: string): boolean => !directories.has(name) && !parentPrefixes(name).some((prefix) => entries.has(prefix));

  const add = (name: string, file: Blob, overwrite = false): void => {
    if (!overwrite && entries.has(name)) return;
    if (!canAdd(name)) return;
    entries.set(name, file);
    for (const prefix of parentPrefixes(name)) directories.add(prefix);
  };

  if (parts.subjectDb !== undefined) {
    for (const entry of parts.subjectDb.discovery.files) add(entry.name, entry.file);
    // Alias auxiliary taxonomy files to root (`/work`) so `SeqDB_ResolveDbPath` resolves them
    // without needing BLASTDB set. Prefer files in the selected DB's directory, then anywhere in the pick.
    const dbDir = parentDir(parts.subjectDb.selectedName);
    for (const aux of AUX_TAXONOMY_FILES) {
      if (entries.has(aux)) continue;
      const inDbDir = dbDir === '' ? undefined : parts.subjectDb.discovery.files.find((entry) => entry.name === `${dbDir}/${aux}`);
      const anywhere = inDbDir ?? parts.subjectDb.discovery.files.find((entry) => baseName(entry.name) === aux);
      if (anywhere !== undefined) add(aux, anywhere.file);
    }
  }

  if (parts.filteringDb !== undefined) {
    for (const entry of parts.filteringDb.discovery.files) add(entry.name, entry.file);
  }

  if (parts.filterFiles !== undefined) {
    for (const file of Object.values(parts.filterFiles)) {
      const path = relativeFilePath(file);
      if (path !== undefined) add(file.name, file, true);
    }
  }

  return Object.freeze([...entries.entries()].map(([name, file]) => Object.freeze({ name, file })));
}
