// `describe` and `validate` over a serial reactor (docs/web/abi_v2.md §4, §9). The Data
// worker answers them with its own serial instance, so that the application can check a
// search while the Engine worker is inside a long `run` (the Engine worker's thread is
// then busy). Validation there gives the same result as in the module that runs the
// search: the argv of the application has no -num_threads, the only option whose
// validation depends on the build; the Engine worker validates the argv it runs, with its
// -num_threads, again before the search.
import { OUTPUT_FORMATS, type OutputFormat } from '../../domain/output-format';
import type { ProgramId } from '../../domain/programs';
import type { EngineGateway, ParameterDescription, ProgramDescription, ValidationResult } from '../../ports/engine';
import type { ReactorAbi } from './abi';

export type EngineControl = Pick<EngineGateway, 'describe' | 'validate'>;

interface DescribeJson {
  readonly program: ProgramId;
  readonly formats: readonly number[];
  readonly parameters: ReadonlyArray<{
    readonly flag: string;
    readonly help: string;
    readonly takes_value: boolean;
    readonly default?: string;
    readonly choices?: readonly string[];
  }>;
  readonly query_gencodes?: readonly number[];
  readonly subject_gencodes?: readonly number[];
}

/** Maps the *describe* JSON to the port's field names. */
export function toProgramDescription(json: unknown): ProgramDescription {
  const value = json as DescribeJson;
  const formats = value.formats.filter((format): format is OutputFormat => OUTPUT_FORMATS.includes(format as OutputFormat));
  return {
    program: value.program,
    formats,
    parameters: value.parameters.map(
      (parameter): ParameterDescription => ({
        flag: parameter.flag,
        help: parameter.help,
        takesValue: parameter.takes_value,
        ...(parameter.default === undefined ? {} : { defaultValue: parameter.default }),
        ...(parameter.choices === undefined ? {} : { choices: parameter.choices }),
      }),
    ),
    ...(value.query_gencodes === undefined ? {} : { queryGencodes: value.query_gencodes }),
    ...(value.subject_gencodes === undefined ? {} : { subjectGencodes: value.subject_gencodes }),
  };
}

/** EngineControl over a serial reactor that `reactor` opens (and reopens after a stop). */
export function reactorControl(reactor: () => Promise<ReactorAbi>): EngineControl {
  return {
    async describe(program) {
      return toProgramDescription((await reactor()).describe(program));
    },
    async validate(argv): Promise<ValidationResult> {
      const message = (await reactor()).validate(argv);
      return message === undefined ? { ok: true } : { ok: false, message };
    },
  };
}
