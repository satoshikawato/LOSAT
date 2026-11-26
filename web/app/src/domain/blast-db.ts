// Local BLAST database discovery, molecule-type matching, filter option metadata,
// run-directory path normalization, duplicate checks and run-directory assembly
// (WQA D12b + RG-D4 W2). Single owner reused by Subject "Database" (-db),
// BLASTN "Repeats database" (-filtering_db), and WindowMasker (-infmt blastdb).
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

export type DirectoryEntries = readonly DirectoryFile[];

export interface BlastDbTaxonomy {
  readonly taxdb: boolean;
  readonly sqlite: boolean;
  /** Present auxiliary taxonomy files in the database's own folder. */
  readonly entries: DirectoryEntries;
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
}

export interface AssembledSearchDirectory {
  readonly files: DirectoryEntries;
  /**
   * Relative folder of the selected `-db` whose auxiliary taxonomy files were mounted at the
   * run directory root (`/work`); CLI reproduction commands prepend `BLASTDB=<blastdbEnv>`.
   */
  readonly blastdbEnv?: string;
}

/** Standard user-facing error when two different files claim the same run-directory path. */
export function duplicateRunFileMessage(name: string): string {
  return `A file named "${name}" is already in the run directory.`;
}

/**
 * Normalizes a picked or folder-relative entry name for a run directory: uses `/` separators,
 * trims whitespace, strips leading `./` and `/` segments, and rejects empty, interior `.`, `..`,
 * or NUL-containing paths (matching `RunDirectory`'s path rules in `packages/losat`).
 */
export function normalizeRunPath(raw: string): string | undefined {
  const cleaned = raw.replaceAll('\\', '/').trim();
  if (cleaned === '' || cleaned.includes('\0')) return undefined;
  const parts: string[] = [];
  for (const segment of cleaned.split('/')) {
    if (segment === '' || (parts.length === 0 && segment === '.')) continue;
    if (segment === '.' || segment === '..') return undefined;
    parts.push(segment);
  }
  return parts.length === 0 ? undefined : parts.join('/');
}

/** Safe relative run-directory path for a picked File (preferring `relativeName`, then `webkitRelativePath`, then `file.name`). */
export function pickedFilePath(file: File, relativeName?: string): string | undefined {
  const rel = relativeName ?? (file as File & { readonly webkitRelativePath?: string }).webkitRelativePath;
  return normalizeRunPath(rel && rel.trim() !== '' ? rel : file.name);
}

function parentPrefixes(name: string): string[] {
  const segments = name.split('/');
  const prefixes: string[] = [];
  let prefix = '';
  for (let i = 0; i < segments.length - 1; i++) {
    prefix = prefix === '' ? segments[i]! : `${prefix}/${segments[i]!}`;
    prefixes.push(prefix);
  }
  return prefixes;
}

/** Two picks of one disk file: the browser creates a new File each time; name, size, and lastModified are what it exposes. */
function sameDiskFile(a: Blob, b: Blob): boolean {
  return a instanceof File && b instanceof File && a.name === b.name && a.size === b.size && a.lastModified === b.lastModified;
}

/**
 * Single owner for assembling unique run-directory entries: allows the exact same Blob/File
 * reference or a second pick of the same disk File (`name`, `size`, and `lastModified` equal)
 * to be added twice under the same path (mounting it once), and throws `duplicateRunFileMessage`
 * whenever two different files claim the same path or collide as file vs. directory.
 */
export class RunDirectoryBuilder {
  private readonly entries = new Map<string, Blob>();
  private readonly directories = new Set<string>();

  add(name: string, file: Blob): void {
    const existing = this.entries.get(name);
    if (existing === file || (existing !== undefined && sameDiskFile(existing, file))) return;
    if (existing !== undefined || this.directories.has(name)) {
      throw new Error(duplicateRunFileMessage(name));
    }
    for (const prefix of parentPrefixes(name)) {
      if (this.entries.has(prefix)) {
        throw new Error(duplicateRunFileMessage(prefix));
      }
    }
    this.entries.set(name, file);
    for (const prefix of parentPrefixes(name)) this.directories.add(prefix);
  }

  has(name: string): boolean {
    return this.entries.has(name);
  }

  build(): readonly DirectoryFile[] {
    return Object.freeze([...this.entries.entries()].map(([name, file]) => Object.freeze({ name, file })));
  }
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

function joinPath(dir: string, rel: string): string {
  const clean = rel.replace(/^\/+/, '');
  return dir === '' ? clean : `${dir}/${clean}`;
}

function classifyBlastDbFile(
  name: string,
): { readonly stem: string; readonly kind: SequenceKind; readonly entry: 'alias' | 'volume' } | undefined {
  // Case-sensitive extensions matching NCBI `FindBlastDBs` (`seqdb.cpp:1438-1449`: *.pin, *.pal, *.nin, *.nal).
  if (name.length <= 4) return undefined;
  const ext = name.slice(-4);
  const stem = name.slice(0, -4);
  if (ext === '.nal') return { stem, kind: 'nucleotide', entry: 'alias' };
  if (ext === '.nin') return { stem, kind: 'nucleotide', entry: 'volume' };
  if (ext === '.pal') return { stem, kind: 'protein', entry: 'alias' };
  if (ext === '.pin') return { stem, kind: 'protein', entry: 'volume' };
  return undefined;
}

/**
 * Synchronous, name-only BLAST database discovery shared by Subject Database (`-db`),
 * BLASTN Repeats database (`-filtering_db`), and WindowMasker (`-infmt blastdb`).
 * Matches NCBI `FindBlastDBs` case-sensitive `.nal`/`.nin`/`.pal`/`.pin` extensions,
 * prefers aliases over index volumes of the same stem, and hides numbered volume stems
 * (`foo.00`) when an alias `foo` of the same molecule type exists (plus any explicit
 * `hiddenComponents` keys `${kind}:${stem}` learned from alias `DBLIST` lines).
 */
export function discoverBlastDatabaseNames(
  fileNames: readonly string[],
  hiddenComponents?: ReadonlySet<string>,
): readonly BlastDbCandidate[] {
  const byKey = new Map<string, { readonly name: string; readonly kind: SequenceKind; readonly entry: 'alias' | 'volume' }>();
  const aliasStems = new Set<string>();

  for (const name of fileNames) {
    const classified = classifyBlastDbFile(name);
    if (classified === undefined) continue;
    const key = `${classified.kind}:${classified.stem}`;
    if (classified.entry === 'alias') aliasStems.add(key);
    const existing = byKey.get(key);
    if (existing === undefined || (existing.entry === 'volume' && classified.entry === 'alias')) {
      byKey.set(key, { name: classified.stem, kind: classified.kind, entry: classified.entry });
    }
  }

  const databases: BlastDbCandidate[] = [];
  for (const [key, item] of byKey) {
    if (item.entry === 'volume') {
      if (hiddenComponents?.has(key)) continue;
      const numbered = /^(.*)\.\d{2,}$/.exec(item.name);
      if (numbered !== null && aliasStems.has(`${item.kind}:${numbered[1]!}`)) continue;
    }
    databases.push(Object.freeze({ name: item.name, kind: item.kind, entry: item.entry }));
  }

  databases.sort((a, b) => (a.name < b.name ? -1 : a.name > b.name ? 1 : a.kind < b.kind ? -1 : a.kind > b.kind ? 1 : 0));
  return Object.freeze(databases);
}

const byShortestName = (a: BlastDbCandidate, b: BlastDbCandidate) => a.name.length - b.name.length || a.name.localeCompare(b.name);

/**
 * Picks the primary BLAST database name from a list of run-directory paths: an alias file
 * (`.nal`/`.pal`) first, otherwise a volume (`.nin`/`.pin`), shortest name first.
 */
export function primaryBlastDbName(fileNames: readonly string[]): string | undefined {
  const candidates = discoverBlastDatabaseNames(fileNames);
  if (candidates.length === 0) return undefined;
  const aliases = candidates.filter((db) => db.entry === 'alias').sort(byShortestName);
  if (aliases.length > 0) return aliases[0]!.name;
  return [...candidates].sort(byShortestName)[0]!.name;
}

/**
 * Discovers valid BLAST databases (`.nal`/`.pal`, `.nin`/`.pin`, multi-volume sets) from picked
 * files. Reads only tiny `.nal`/`.pal` alias text files; binary volume files are never read into memory.
 */
export async function discoverBlastDatabases(picked: readonly File[]): Promise<BlastDbDiscovery> {
  const builder = new RunDirectoryBuilder();
  const filesByName = new Map<string, File>();
  let totalBytes = 0;
  let folder: string | undefined;
  let allShareFolder = picked.length > 0;

  for (const file of picked) {
    const path = pickedFilePath(file);
    if (path === undefined) continue;
    if (builder.has(path)) {
      builder.add(path, file);
      continue;
    }
    builder.add(path, file);
    filesByName.set(path, file);
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

  const files = builder.build();
  const names = files.map((entry) => entry.name);
  const initial = discoverBlastDatabaseNames(names);

  // Read alias metadata (`TITLE`, `DBLIST`) to hide component volumes of multi-volume aliases.
  const aliasMeta = new Map<string, { readonly title?: string; readonly dblist: readonly string[] }>();
  const hiddenComponents = new Set<string>();

  for (const item of initial) {
    if (item.entry !== 'alias') continue;
    const ext = item.kind === 'nucleotide' ? '.nal' : '.pal';
    const aliasFile = filesByName.get(`${item.name}${ext}`);
    if (aliasFile === undefined || aliasFile.size > MAX_ALIAS_BYTES) continue;
    try {
      const text = await aliasFile.text();
      const meta = parseAliasMetadata(text);
      aliasMeta.set(`${item.kind}:${item.name}`, meta);
      const dir = parentDir(item.name);
      for (const token of meta.dblist) {
        hiddenComponents.add(`${item.kind}:${token}`);
        hiddenComponents.add(`${item.kind}:${joinPath(dir, token)}`);
      }
    } catch {
      // Unreadable alias file is still surfaced so the engine can report its own error if chosen.
    }
  }

  const databases = discoverBlastDatabaseNames(names, hiddenComponents).map((item) => {
    const meta = aliasMeta.get(`${item.kind}:${item.name}`);
    if (meta === undefined) return item;
    return Object.freeze({
      ...item,
      ...(meta.title === undefined ? {} : { title: meta.title }),
      ...(meta.dblist.length > 0 ? { volumes: Object.freeze([...meta.dblist]) } : {}),
    });
  });

  return Object.freeze({
    ...(folder === undefined ? {} : { folder }),
    databases: Object.freeze(databases),
    files,
    totalBytes,
  });
}

/**
 * Auxiliary taxonomy files (`taxdb.bti`, `taxdb.btd`, `taxonomy4blast.sqlite3`) present in the
 * selected database's own folder (never from unrelated sibling subfolders).
 */
export function databaseTaxonomy(discovery: BlastDbDiscovery, selectedName: string): BlastDbTaxonomy {
  const dbDir = parentDir(selectedName);
  const byName = new Map(discovery.files.map((entry) => [entry.name, entry]));
  const present: DirectoryFile[] = [];
  let hasBti = false;
  let hasBtd = false;
  let hasSqlite = false;

  for (const aux of AUX_TAXONOMY_FILES) {
    const entry = byName.get(joinPath(dbDir, aux));
    if (entry === undefined) continue;
    present.push(entry);
    if (aux === TAXDB_INDEX_FILE) hasBti = true;
    else if (aux === TAXDB_DATA_FILE) hasBtd = true;
    else if (aux === TAXONOMY_SQLITE_FILE) hasSqlite = true;
  }

  return Object.freeze({
    taxdb: hasBti && hasBtd,
    sqlite: hasSqlite,
    entries: Object.freeze(present),
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
 * Builds the `DirectoryFile[]` mounted into `RunDirectory` for a search run:
 * - all picked files of the subject DB (with `taxdb.bti`, `taxdb.btd`, `taxonomy4blast.sqlite3`
 *   from the selected DB's own folder also mounted at root so NCBI's `SeqDB_ResolveDbPath` finds
 *   them in `/work`, recording `blastdbEnv` so CLI reproduction commands prepend `BLASTDB=<dir>`)
 * - all picked files of BLASTN's `-filtering_db` repeats DB
 * - all picked DB filter list files (`-seqidlist`, `-taxidlist`, etc.), mounted at root by `file.name`
 *
 * Refuses with `duplicateRunFileMessage` whenever two different files claim the same mount path
 * or collide as file vs. directory. The same `File` reference or a second pick of the same disk
 * file (`name`, `size`, and `lastModified` equal) is mounted once.
 */
export function assembleSearchDirectory(parts: {
  readonly subjectDb?: { readonly discovery: BlastDbDiscovery; readonly selectedName: string };
  readonly filterFiles?: Readonly<Record<string, File>>;
  readonly filteringDb?: { readonly discovery: BlastDbDiscovery; readonly selectedName: string };
}): AssembledSearchDirectory {
  const dir = new RunDirectoryBuilder();
  let blastdbEnv: string | undefined;

  if (parts.subjectDb !== undefined) {
    for (const entry of parts.subjectDb.discovery.files) dir.add(entry.name, entry.file);
    // Alias auxiliary taxonomy files from the selected DB's own folder to root (`/work`)
    // so `SeqDB_ResolveDbPath` resolves them as `BLASTDB=<dbDir>` does on the CLI.
    const dbDir = parentDir(parts.subjectDb.selectedName);
    if (dbDir !== '') {
      const tax = databaseTaxonomy(parts.subjectDb.discovery, parts.subjectDb.selectedName);
      for (const entry of tax.entries) {
        const auxName = entry.name.slice(dbDir.length + 1);
        dir.add(auxName, entry.file);
        blastdbEnv = dbDir;
      }
    }
  }

  if (parts.filteringDb !== undefined) {
    for (const entry of parts.filteringDb.discovery.files) dir.add(entry.name, entry.file);
  }

  if (parts.filterFiles !== undefined) {
    for (const file of Object.values(parts.filterFiles)) {
      const name = normalizeRunPath(file.name);
      if (name === undefined) throw new Error(`Invalid file name: "${file.name}".`);
      dir.add(name, file);
    }
  }

  const files = dir.build();
  return Object.freeze(blastdbEnv === undefined ? { files } : { files, blastdbEnv });
}
