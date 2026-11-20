// The registry file of `-conffile` in the application (line RG item E1w-3b): the name under which the
// argv reads it, as the engine's start gives it (one rule for the form, the draft and the runs), and
// the sentence for a chosen file that can no longer be read.
import type { StartInfo } from 'losat/ports';

/**
 * The name of the registry file that the argv's `-conffile` gives, as the engine's start read it: the
 * form asks for a file under it, and the validation and the runs read the one chosen. Undefined without
 * `-conffile`, and for an empty name (`-conffile=`), which names no file that the engine can hold: the
 * engine's message then stands.
 */
export function registryName(start: StartInfo | undefined): string | undefined {
  return start?.conffile || undefined;
}

/**
 * Why the file chosen for `-conffile NAME` cannot be given to the engine: it could not be read (a
 * browser refuses a file from the disk that changed, or went, after it was chosen), with the browser's
 * reason, and what to do.
 */
export function unreadableConffile(name: string, fileName: string, error: unknown): string {
  const reason = (error instanceof Error ? error.message : String(error)).trim().replace(/\.$/, '');
  return `The file chosen for -conffile ${name} (${fileName}) could not be read: ${reason}. Choose it again.`;
}
