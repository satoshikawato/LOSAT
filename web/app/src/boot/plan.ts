// The decisions of the boot page (plan §5.9, design §14, docs/web/release_manifest.md): which
// release to start, whether to pin it, and what to say when it cannot start. Pure: boot.ts
// does the requests and the DOM. A pinned release is never switched silently; only the user's
// "Use release <latest>", "Use release <id> (stored)" (or an update the user chose in the app)
// changes the pin.
import { parseReleaseManifest, releaseFileUrl, type ReleaseManifest } from '../domain/release-manifest';
import { RELEASE_ID_PATTERN, type SiteIndex } from '../domain/site-index';
import { RELEASE_CACHE_PREFIX } from '../sw/protocol';

/** The pinned release in localStorage (S16 decision 3); an app setting, not research data. */
export const PIN_KEY = 'losat-web:release';

/** The pinned release, or undefined when there is none, it is not a release ID, or storage cannot be read. */
export function readPin(storage: () => Pick<Storage, 'getItem'>): string | undefined {
  try {
    const value = storage().getItem(PIN_KEY);
    return value !== null && RELEASE_ID_PATTERN.test(value) ? value : undefined;
  } catch {
    return undefined;
  }
}

/** Pins a release; false when storage cannot be written (the page starts anyway). */
export function writePin(storage: () => Pick<Storage, 'setItem'>, release: string): boolean {
  try {
    storage().setItem(PIN_KEY, release);
    return true;
  } catch {
    return false;
  }
}

export type MessageKind = 'first-visit-offline' | 'no-index' | 'not-served' | 'offline-not-stored' | 'failed' | 'pin-not-saved';

/** What the boot page shows instead of the app; `useLatest` offers the site's latest release. */
export interface BootMessage {
  readonly kind: MessageKind;
  readonly text: string;
  readonly useLatest?: string;
  /** "Try again": reloads the page (the network or the site may be back). */
  readonly retry?: boolean;
  /** Releases that this browser stores whole, each offered as "Use release <id> (stored)": they start without the network. */
  readonly useStored?: readonly string[];
  /** A stored release that the user chose and that was no longer stored at the click (withoutStored). */
  readonly noLongerStored?: string;
}

/** The messages that a reload may resolve, or a stored release may stand in for: no network, or the site did not answer. */
const RETRIED: ReadonlySet<MessageKind> = new Set<MessageKind>(['first-visit-offline', 'no-index', 'offline-not-stored', 'failed']);

/**
 * Whether a message gets "Try again" and the stored releases: where the network or the site
 * failed, and for a release that the site does not serve when no other release is offered (the
 * site index was not read: offline, timed out, a server error; or it names that same release).
 * A reload may then read the index, and a stored release starts without it; without them the
 * message would be a dead end (code review 3 L3). With "Use release <latest>" it needs neither.
 */
function offersChoices(message: BootMessage): boolean {
  return RETRIED.has(message.kind) || (message.kind === 'not-served' && message.useLatest === undefined);
}

/**
 * Adds the choices of a message that the network may resolve (screen review 1 I6; offersChoices):
 * "Try again", and each release in `stored` other than `release` (the one that could not start)
 * and the latest already offered. A stored release is only offered: the pin changes when the user
 * chooses it, never by itself.
 */
export function withChoices(message: BootMessage, stored: readonly string[], release?: string): BootMessage {
  if (!offersChoices(message)) return message;
  const others = stored.filter((id) => id !== release && id !== message.useLatest);
  return { ...message, retry: true, ...(others.length > 0 ? { useStored: others } : {}) };
}

/**
 * The message again after the user chose "Use release <release> (stored)" and the release was no
 * longer stored whole at the click (`stored`, read then): another tab's cleanup may have removed it
 * since the message was shown, for the boot page holds no release lock (code review 3 I1). It says
 * so and offers the stored releases of `message` that remain, with its other choices; the pin is
 * not changed.
 */
export function withoutStored(message: BootMessage, release: string, stored: readonly string[]): BootMessage {
  const { useStored, ...rest } = message;
  const remaining = (useStored ?? []).filter((id) => id !== release && stored.includes(id));
  return { ...rest, noLongerStored: release, ...(remaining.length > 0 ? { useStored: remaining } : {}) };
}

/** The sentence for a stored release that was no longer stored when the user chose it. */
export function noLongerStoredText(release: string): string {
  return `Release ${release} is no longer stored in this browser.`;
}

/** The release whose cache the site's service worker names `name` (src/sw/protocol.ts); undefined for another cache. */
export function releaseOfCache(name: string): string | undefined {
  if (!name.startsWith(RELEASE_CACHE_PREFIX)) return undefined;
  const release = name.slice(RELEASE_CACHE_PREFIX.length);
  return RELEASE_ID_PATTERN.test(release) ? release : undefined;
}

/** The sentence before the stored releases' buttons. */
export function storedText(stored: readonly string[]): string {
  return stored.length === 1
    ? `Release ${stored[0]!} is stored in this browser and starts without the network.`
    : 'These releases are stored in this browser and start without the network.';
}

/** The release to start; `pin` says whether the boot pins it (a first visit). */
export interface Start {
  readonly kind: 'start';
  readonly release: string;
  readonly pin: boolean;
}

/**
 * The release to start: the pinned one if there is one, else the site's latest (pinned now).
 * `index` is the site index, undefined when it was not read (offline, or the request failed).
 */
export function chooseRelease(pin: string | undefined, index: SiteIndex | undefined, online: boolean): Start | BootMessage {
  if (pin !== undefined) return { kind: 'start', release: pin, pin: false };
  if (index !== undefined) return { kind: 'start', release: index.latest, pin: true };
  if (!online) return { kind: 'first-visit-offline', text: 'LOSAT Web needs a network connection the first time it is opened.' };
  return { kind: 'no-index', text: "LOSAT Web could not read this site's list of releases. Reload the page to try again." };
}

/** The answer to the request for a release manifest; undefined when there was none (network error). */
export interface FetchedManifest {
  readonly status: number;
  /** The body read as JSON; undefined with `json: false` when it is not JSON. */
  readonly body: unknown;
  readonly json: boolean;
}

export type ManifestOutcome =
  | { readonly kind: 'ok'; readonly manifest: ReleaseManifest }
  /** The site does not serve the release: 404 or 410, or an answer that is not its manifest. */
  | { readonly kind: 'not-served'; readonly reason: string }
  /** No answer: the network failed, or the page is offline and the release is not stored. */
  | { readonly kind: 'unreachable' }
  /** Another status (a server error). */
  | { readonly kind: 'failed'; readonly reason: string };

/**
 * Whether an answer that would mean "not served" (404, 410, not JSON) must be read once more
 * past the HTTP cache first: under /r/ the site's rules make every answer immutable, a 404
 * included, so a release that was missing once and is served again (a rollback undone) would
 * otherwise stay "not served" for as long as the browser keeps that 404.
 */
export function recheckManifest(fetched: FetchedManifest | undefined): boolean {
  return fetched !== undefined && (fetched.status === 404 || fetched.status === 410 || (fetched.status >= 200 && fetched.status <= 299 && !fetched.json));
}

export function judgeManifest(release: string, fetched: FetchedManifest | undefined): ManifestOutcome {
  if (fetched === undefined) return { kind: 'unreachable' };
  if (fetched.status === 404 || fetched.status === 410) return { kind: 'not-served', reason: `the site answered ${fetched.status}` };
  if (fetched.status < 200 || fetched.status > 299) return { kind: 'failed', reason: `the site answered ${fetched.status}` };
  if (!fetched.json) return { kind: 'not-served', reason: 'the answer is not JSON' };
  try {
    return { kind: 'ok', manifest: parseReleaseManifest(fetched.body, release) };
  } catch (error) {
    return { kind: 'not-served', reason: error instanceof Error ? error.message : String(error) };
  }
}

/** The message for a release that cannot start. */
export function failureMessage(release: string, outcome: Exclude<ManifestOutcome, { kind: 'ok' }>, index: SiteIndex | undefined, online: boolean): BootMessage {
  switch (outcome.kind) {
    case 'not-served':
      return {
        kind: 'not-served',
        text: `This browser uses LOSAT Web release ${release}, which this site no longer serves.`,
        ...(index !== undefined && index.latest !== release ? { useLatest: index.latest } : {}),
      };
    case 'unreachable':
      return online
        ? { kind: 'failed', text: `LOSAT Web could not load release ${release}. Check the network connection and reload the page.` }
        : {
            kind: 'offline-not-stored',
            text: `This browser uses LOSAT Web release ${release}, which is not stored for offline use. Connect to the network and reload the page.`,
          };
    case 'failed':
      return { kind: 'failed', text: `LOSAT Web could not load release ${release} (${outcome.reason}). Reload the page to try again.` };
  }
}

/** The message when the user chose another release and it could not be pinned. */
export function pinNotSavedMessage(release: string): BootMessage {
  return { kind: 'pin-not-saved', text: `This browser did not let LOSAT Web save the choice of release ${release}.` };
}

/** The absolute URLs of what the boot page adds to the document: the style sheets, then the module script. */
export function entryUrls(manifest: ReleaseManifest): { readonly script: string; readonly styles: readonly string[] } {
  return {
    script: releaseFileUrl(manifest.release, manifest.entry.script),
    styles: manifest.entry.styles.map((path) => releaseFileUrl(manifest.release, path)),
  };
}
