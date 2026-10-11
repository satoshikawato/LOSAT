// Small settings of the application kept in this browser (localStorage), never research
// data: whether "Getting started" was closed, for one. Every key starts with `losat-web:`
// (S16 decision 3). Storage can be missing or refuse (a private window, blocked site data):
// a read then gives undefined and a write false, and the page works without it. The browser
// implementation is src/infra/browser/local-settings.ts.

export const LOCAL_SETTINGS_PREFIX = 'losat-web:';

export interface LocalSettings {
  /** The value kept under `key`, or undefined when there is none or storage cannot be read. */
  read(key: string): string | undefined;
  /** Keeps `value` under `key`; false when storage refused it. */
  write(key: string, value: string): boolean;
}
