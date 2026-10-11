// The boot script of the site (`/boot.js`, built by build/site.ts; docs/web/release_manifest.md):
// it starts the pinned release, or pins the site's latest on a first visit, by adding the
// release's style sheets and entry script to this page. It reads no research data and sends
// nothing but two GETs: `/releases.json` and the release's manifest. The decisions are in plan.ts.
import { releaseManifestUrl } from '../domain/release-manifest';
import { parseSiteIndex, SITE_INDEX_PATH, type SiteIndex } from '../domain/site-index';
import {
  chooseRelease,
  entryUrls,
  failureMessage,
  judgeManifest,
  pinNotSavedMessage,
  readPin,
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

/** The release's manifest, with the default cache mode: a service worker answers it when the release is stored. */
async function fetchManifest(release: string): Promise<FetchedManifest | undefined> {
  let response: Response;
  try {
    response = await fetch(releaseManifestUrl(release), { credentials: 'same-origin', signal: AbortSignal.timeout(MANIFEST_TIMEOUT_MS) });
  } catch {
    return undefined;
  }
  try {
    return { status: response.status, body: await response.json(), json: true };
  } catch {
    return { status: response.status, body: undefined, json: false };
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

function show(message: BootMessage): void {
  const root = document.getElementById('boot');
  if (root === null) return;
  const text = document.createElement('p');
  text.setAttribute('data-testid', 'boot-message');
  text.setAttribute('data-kind', message.kind);
  text.textContent = message.text;
  root.replaceChildren(text);
  const latest = message.useLatest;
  if (latest !== undefined) {
    const button = document.createElement('button');
    button.type = 'button';
    button.setAttribute('data-testid', 'boot-use-latest');
    button.textContent = `Use release ${latest}`;
    button.addEventListener('click', () => {
      if (writePin(storage, latest)) window.location.reload();
      else show(pinNotSavedMessage(latest));
    });
    const actions = document.createElement('p');
    actions.append(button);
    root.append(actions);
  }
  root.hidden = false;
}

async function boot(): Promise<void> {
  const pin = readPin(storage);
  const index = online() ? fetchIndex() : Promise.resolve(undefined);
  // With a pin the index is needed only when the pinned release cannot start.
  const choice = chooseRelease(pin, pin === undefined ? await index : undefined, online());
  if (choice.kind !== 'start') {
    show(choice);
    return;
  }
  // A pin that cannot be written is not an error: this visit starts the release anyway.
  if (choice.pin) writePin(storage, choice.release);
  const outcome = judgeManifest(choice.release, await fetchManifest(choice.release));
  if (outcome.kind !== 'ok') {
    show(failureMessage(choice.release, outcome, await index, online()));
    return;
  }
  const { script, styles } = entryUrls(outcome.manifest);
  try {
    await start(script, styles);
  } catch {
    show(failureMessage(choice.release, { kind: 'unreachable' }, undefined, online()));
  }
}

boot().catch(() => show({ kind: 'failed', text: 'LOSAT Web could not start. Reload the page to try again.' }));
