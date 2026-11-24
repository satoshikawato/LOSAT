// Where messages of a run are, in the words of the results (line WQA's sweep, D7, 2026-11-20): the error of an output
// that the run could not write (`ResultSetRef.failed`), and the run's warnings (`LoadedRun.messages`), each at the place
// that the engine gave, never one guessed from the program or the argv. The run's diagnostics are its Warnings in Run
// details; the file of its -logfile is a tab of the Outputs (TBLASTX writes the messages of its run there); a -logfile
// that names no file writes them nowhere.
import { runFileOption, type MessagePlace, type OutputPlace, type RunFile } from 'losat/domain';
import type { RunView } from '../application/coordinator';

/**
 * A piece of the words of a place: an option name, which the page shows unbroken (`option-name`, as other option names
 * in the results: the narrow screen broke "-logfile" after its hyphen), or plain words (`PlaceWords.vue` shows them).
 */
export interface PlacePart {
  readonly text: string;
  readonly option?: true;
}

/** What follows "its error is" (or "their error is"): where an error of `run` is, in parts (`placeText` joins them). */
export function placeParts(place: OutputPlace, run: RunView): readonly PlacePart[] {
  switch (place) {
    // An output that a failed run does not keep (line WQA's sweep, D15, 2026-11-20): the Outputs of a run do not list it.
    case 'run':
      return [{ text: 'the error that ended the run' }];
    case 'diagnostics':
      return [{ text: 'under Warnings in Run details' }];
    case 'logfile': {
      const names = logFiles(run)
        .map((file) => file.name)
        .join(', ');
      return [{ text: "in the file of the run's " }, LOGFILE, { text: `${names === '' ? '' : `, ${names},`} under Outputs` }];
    }
    case 'nowhere':
      return [{ text: "not kept, as the run's " }, LOGFILE, { text: ' writes its messages nowhere' }];
  }
}

const LOGFILE: PlacePart = { text: '-logfile', option: true };

/** The words of `placeParts` as one text. */
export function placeText(place: OutputPlace, run: RunView): string {
  return joined(placeParts(place, run));
}

/** What Run details says under Warnings where the engine wrote none to the diagnostics and they are where its messages go. */
export const NO_WARNINGS = 'The engine wrote no warnings.';

/**
 * What Run details says under Warnings besides the run's diagnostics (`diagnostics`, stream 3), from where the run's
 * messages went (`messages`, the engine's word): where its warnings are when they did not go to the diagnostics (the
 * file of its -logfile, or nowhere), said also when the diagnostics hold lines; that the engine wrote none where the
 * place holds nothing (the diagnostics, or the file of the -logfile, empty) and the diagnostics are empty; else nothing.
 * A place of nowhere may have taken warnings or none: its words say only that none is kept.
 */
export function warningsPlaceParts(messages: MessagePlace, diagnostics: string, run: RunView): readonly PlacePart[] | undefined {
  if (messages === 'diagnostics' || (messages === 'logfile' && emptyFile(run))) {
    return diagnostics === '' ? [{ text: NO_WARNINGS }] : undefined;
  }
  return [{ text: `${diagnostics === '' ? "The run's warnings are" : "The run's other warnings are"} ` }, ...placeParts(messages, run), { text: '.' }];
}

/**
 * Whether the file of the run's -logfile is kept and holds nothing: a TBLASTX search without warnings leaves the file
 * empty, as the CLI's (line WQA's sweep, D7, 2026-11-20). A file that is not kept says nothing; the run's search
 * strategy file (-export_search_strategy) holds no messages.
 */
function emptyFile(run: RunView): boolean {
  const files = logFiles(run);
  return files.length > 0 && files.every((file) => file.bytes.length === 0);
}

/** The kept file of the run's -logfile (none or one), of the files that the run made (`RunRecord.runFiles`). */
function logFiles(run: RunView): readonly RunFile[] {
  return (run.record.runFiles ?? []).filter((file) => runFileOption(file) === 'logfile');
}

/** The words of `warningsPlaceParts` as one text, or undefined where Warnings says nothing besides the diagnostics. */
export function warningsPlaceText(messages: MessagePlace, diagnostics: string, run: RunView): string | undefined {
  const parts = warningsPlaceParts(messages, diagnostics, run);
  return parts === undefined ? undefined : joined(parts);
}

function joined(parts: readonly PlacePart[]): string {
  return parts.map((part) => part.text).join('');
}
