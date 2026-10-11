// The boot script of the site (`/boot.js`, built by build/site.ts; docs/web/release_manifest.md):
// it starts the pinned release, or pins the site's latest on a first visit, by adding the
// release's style sheets and entry script to this page. It reads no research data and sends
// nothing but two GETs: `/releases.json` and the release's manifest. When a release cannot
// start, it says why and offers "Try again" and the releases that the service worker stores
// (it reads the names of the site's caches, nothing else). The decisions are in plan.ts.
import { releaseManifestUrl } from '../domain/release-manifest';
import { parseSiteIndex, SITE_INDEX_PATH, type SiteIndex } from '../domain/site-index';
import { MANIFEST_NAME, releasePath, SW_SCOPE, SW_URL } from '../sw/protocol';
import {
  chooseRelease,
  entryUrls,
  failureMessage,
  judgeManifest,
  pinNotSavedMessage,
  readPin,
  recheckManifest,
  releaseOfCache,
  storedText,
  withChoices,
  writePin,
  type BootMessage,
  type FetchedManifest,
} from './plan';

const INDEX_TIMEOUT_MS = 5_000;
const MANIFEST_TIMEOUT_MS = 20_000;

const online = (): boolean => navigator.onLine !== false;
// Reading `window.localStorage` itself can throw (storage blocked); plan.ts calls this inside try/catch.
const storage = (): Storage => window.localStorage;

/** The site index from the network, never from a cache; undefined when it cannot be read. */
async function fetchIndex(): Promise<SiteIndex | undefined> {
  try {
    const response = await fetch(SITE_INDEX_PATH, { cache: 'no-store', credentials: 'same-origin', signal: AbortSignal.timeout(INDEX_TIMEOUT_MS) });
    return response.ok ? parseSiteIndex(await response.json()) : undefined;
  } catch {
    return undefined;
  }
}

/**
 * The release's manifest, with the default cache mode: a service worker answers it when the
 * release is stored. An answer that means "not served" is read once more past the HTTP cache
 * (plan.ts recheckManifest); without an answer to that, the first one stands.
 */
async function fetchManifest(release: string): Promise<FetchedManifest | undefined> {
  const first = await requestManifest(release, 'default');
  return recheckManifest(first) ? ((await requestManifest(release, 'reload')) ?? first) : first;
}

async function requestManifest(release: string, cache: RequestCache): Promise<FetchedManifest | undefined> {
  let response: Response;
  try {
    response = await fetch(releaseManifestUrl(release), { cache, credentials: 'same-origin', signal: AbortSignal.timeout(MANIFEST_TIMEOUT_MS) });
  } catch {
    return undefined;
  }
  try {
    return { status: response.status, body: await response.json(), json: true };
  } catch {
    return { status: response.status, body: undefined, json: false };
  }
}

/**
 * The releases that the site's service worker stores whole: a release cache that holds its
 * manifest, which is stored last (src/sw/prepare.ts). None where the page has no CacheStorage
 * or cannot read it.
 */
async function storedReleases(): Promise<string[]> {
  try {
    if (typeof caches === 'undefined') return [];
    const stored: string[] = [];
    for (const name of await caches.keys()) {
      const release = releaseOfCache(name);
      if (release === undefined) continue;
      if ((await caches.match(window.location.origin + releasePath(release, MANIFEST_NAME), { cacheName: name })) !== undefined) {
        stored.push(release);
      }
    }
    return stored.sort();
  } catch {
    return [];
  }
}

function load(element: HTMLLinkElement | HTMLScriptElement): Promise<void> {
  return new Promise((resolve, reject) => {
    element.addEventListener('load', () => resolve(), { once: true });
    element.addEventListener('error', () => reject(new Error('not loaded')), { once: true });
    document.head.append(element);
  });
}

/** Adds the style sheets, waits for them (so that the app never shows unstyled), then the entry script. */
async function start(script: string, styles: readonly string[]): Promise<void> {
  await Promise.all(
    styles.map((href) => {
      const link = document.createElement('link');
      link.rel = 'stylesheet';
      link.href = href;
      return load(link);
    }),
  );
  const element = document.createElement('script');
  element.type = 'module';
  element.src = script;
  await load(element);
}

function button(testid: string, text: string, onClick: () => void): HTMLButtonElement {
  const element = document.createElement('button');
  element.type = 'button';
  element.setAttribute('data-testid', testid);
  element.textContent = text;
  element.addEventListener('click', onClick);
  return element;
}

/** Pins the release that the user chose, and starts it. */
function use(release: string): void {
  if (writePin(storage, release)) window.location.reload();
  else show(pinNotSavedMessage(release));
}

function show(message: BootMessage): void {
  const root = document.getElementById('boot');
  if (root === null) return;
  const text = document.createElement('p');
  text.setAttribute('data-testid', 'boot-message');
  text.setAttribute('data-kind', message.kind);
  text.textContent = message.text;
  root.replaceChildren(text);
  const stored = message.useStored ?? [];
  if (stored.length > 0) {
    const note = document.createElement('p');
    note.setAttribute('data-testid', 'boot-stored');
    note.textContent = storedText(stored);
    root.append(note);
  }
  const actions = document.createElement('p');
  actions.className = 'boot-actions';
  const latest = message.useLatest;
  if (latest !== undefined) actions.append(button('boot-use-latest', `Use release ${latest}`, () => use(latest)));
  for (const release of stored) {
    const element = button('boot-use-stored', `Use release ${release} (stored)`, () => use(release));
    element.setAttribute('data-release', release);
    actions.append(element);
  }
  if (message.retry === true) {
    const retry = button('boot-retry', 'Try again', () => window.location.reload());
    retry.className = 'secondary';
    actions.append(retry);
  }
  if (actions.childElementCount > 0) root.append(actions);
  root.hidden = false;
}

async function boot(): Promise<void> {
  const pin = readPin(storage);
  const index = online() ? fetchIndex() : Promise.resolve(undefined);
  // With a pin the index is needed only when the pinned release cannot start.
  const choice = chooseRelease(pin, pin === undefined ? await index : undefined, online());
  if (choice.kind !== 'start') {
    show(withChoices(choice, await storedReleases()));
    return;
  }
  // A pin that cannot be written is not an error: this visit starts the release anyway.
  if (choice.pin) writePin(storage, choice.release);
  const outcome = judgeManifest(choice.release, await fetchManifest(choice.release));
  if (outcome.kind !== 'ok') {
    show(withChoices(failureMessage(choice.release, outcome, await index, online()), await storedReleases(), choice.release));
    return;
  }
  const { script, styles } = entryUrls(outcome.manifest);
  try {
    await start(script, styles);
  } catch {
    show(withChoices(failureMessage(choice.release, { kind: 'unreachable' }, undefined, online()), await storedReleases(), choice.release));
    return;
  }
  registerServiceWorker();
}

/**
 * Registers the site's service worker (src/sw/worker.ts) once the release has started, so that
 * it never delays the start. It stores a release only when the app asks it to. A browser
 * without service workers, or a registration that fails, is not an error: the app works online.
 */
function registerServiceWorker(): void {
  if (!('serviceWorker' in navigator)) return;
  navigator.serviceWorker.register(SW_URL, { scope: SW_SCOPE, updateViaCache: 'none' }).catch(() => undefined);
}

boot().catch(() => show(withChoices({ kind: 'failed', text: 'LOSAT Web could not start. Reload the page to try again.' }, [])));
