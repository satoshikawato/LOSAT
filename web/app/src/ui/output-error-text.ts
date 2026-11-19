// Where the error of an output that a run could not write is, in the words of the results (line WQA's sweep, D7,
// 2026-11-20): the place that the engine gave (`ResultSetRef.failed`), never one guessed from the program or the
// argv. The run's diagnostics are its Warnings in Run details; the file of its -logfile is a tab of the Outputs
// (TBLASTX writes the messages of its run there, as the CLI does); a -logfile that names no file writes them nowhere.
import type { OutputErrorPlace } from 'losat/ports';
import type { RunView } from '../application/coordinator';

/** What follows "its error is" (or "their error is"): where the error of an output of `run` that could not be written is. */
export function outputErrorPlace(error: OutputErrorPlace, run: RunView): string {
  switch (error) {
    case 'diagnostics':
      return 'under Warnings in Run details';
    case 'logfile': {
      const names = (run.record.runFiles ?? []).map((file) => file.name).join(', ');
      return `in the file of the run's -logfile${names === '' ? '' : `, ${names},`} under Outputs`;
    }
    case 'nowhere':
      return "not kept, as the run's -logfile writes its messages nowhere";
  }
}
