// The application's small settings in localStorage (ports/local-settings.ts). The accessor
// itself can throw (blocked site data), so every read and write is wrapped.
import { LOCAL_SETTINGS_PREFIX, type LocalSettings } from '../../ports/local-settings';

function checked(key: string): string {
  if (!key.startsWith(LOCAL_SETTINGS_PREFIX)) throw new Error(`a setting's key starts with ${LOCAL_SETTINGS_PREFIX}: ${key}`);
  return key;
}

export const browserLocalSettings: LocalSettings = {
  read(key) {
    checked(key);
    try {
      return window.localStorage.getItem(key) ?? undefined;
    } catch {
      return undefined;
    }
  },
  write(key, value) {
    checked(key);
    try {
      window.localStorage.setItem(key, value);
      return true;
    } catch {
      return false;
    }
  },
};
