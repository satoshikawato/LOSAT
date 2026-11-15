// `describe`, `start` and `validate` over a serial reactor (docs/web/abi_v2.md §4, §9). The Data
// worker answers them with its own serial instance, so that the application can check a
// search while the Engine worker is inside a long `run` (the Engine worker's thread is
// then busy). Validation there gives the same result as in the module that runs the
// search: the argv of the application has no -num_threads, the only option whose
// validation depends on the build; the Engine worker validates the argv it runs, with its
// -num_threads, again before the search. A validation with the registry file of `-conffile` runs on
// an instance of its own (line RG E1w-3b, DECISIONS E1w3b-4), opened for it from the compiled module
// and dropped after it, as each inflation is (inflater.ts): the shared instance never holds the
// file, so no other call of it (a scan, a check, another validation) meets the file, whose start
// would read a file named `.ncbirc` or `<program>.ini` as the CLI reads one in its directory, and
// the memory that the file takes goes with that instance, however large the file chosen.
import { OUTPUT_FORMATS, type OutputFormat } from '../../domain/output-format.js';
import type { ResolvedDefault } from '../../domain/parameters.js';
import type { ProgramId } from '../../domain/programs.js';
import type { EngineGateway, ParameterDescription, ProgramDescription, StartResult, ValidationResult } from '../../ports/engine.js';
import { FILE_REGISTRY, loadFile, type ReactorAbi } from './abi.js';

export type EngineControl = Pick<EngineGateway, 'describe' | 'start' | 'validate'>;

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
  readonly defaults?: DefaultsJson;
}

type ValuesJson = Readonly<Record<string, string>>;

/** The *describe* JSON's `defaults` (abi_v2.md §9): tables by task, of all options at once. */
interface DefaultsJson {
  readonly tasks: Readonly<Record<string, ValuesJson>>;
  readonly by: Readonly<Record<string, Readonly<Record<string, Readonly<Record<string, ValuesJson>>>>>>;
  readonly unless: Readonly<Record<string, readonly string[]>>;
  /** The tasks where an option is told given from not given, by flag. */
  readonly given: Readonly<Record<string, readonly string[]>>;
}

/** One option's part of the `defaults` tables; undefined when they do not name it. */
function resolvedDefault(defaults: DefaultsJson | undefined, flag: string): ResolvedDefault | undefined {
  if (defaults === undefined) return undefined;
  const byTask: Record<string, string> = {};
  for (const [task, values] of Object.entries(defaults.tasks)) {
    const value = values[flag];
    if (value !== undefined) byTask[task] = value;
  }
  const by: Record<string, Record<string, Record<string, string>>> = {};
  for (const [other, table] of Object.entries(defaults.by)) {
    for (const [otherValue, tasks] of Object.entries(table)) {
      for (const [task, values] of Object.entries(tasks)) {
        const value = values[flag];
        if (value === undefined) continue;
        ((by[other] ??= {})[otherValue] ??= {})[task] = value;
      }
    }
  }
  const unless = defaults.unless[flag];
  const given = defaults.given[flag];
  if (Object.keys(byTask).length === 0 && unless === undefined && given === undefined) return undefined;
  return {
    byTask,
    ...(Object.keys(by).length === 0 ? {} : { by }),
    ...(unless === undefined ? {} : { unless }),
    ...(given === undefined ? {} : { given }),
  };
}

/** Maps the *describe* JSON to the port's field names. */
export function toProgramDescription(json: unknown): ProgramDescription {
  const value = json as DescribeJson;
  const formats = value.formats.filter((format): format is OutputFormat => OUTPUT_FORMATS.includes(format as OutputFormat));
  return {
    program: value.program,
    formats,
    parameters: value.parameters.map((parameter): ParameterDescription => {
      const resolved = resolvedDefault(value.defaults, parameter.flag);
      return {
        flag: parameter.flag,
        help: parameter.help,
        takesValue: parameter.takes_value,
        ...(parameter.default === undefined ? {} : { defaultValue: parameter.default }),
        ...(parameter.choices === undefined ? {} : { choices: parameter.choices }),
        ...(resolved === undefined ? {} : { resolvedDefault: resolved }),
      };
    }),
    ...(value.query_gencodes === undefined ? {} : { queryGencodes: value.query_gencodes }),
    ...(value.subject_gencodes === undefined ? {} : { subjectGencodes: value.subject_gencodes }),
  };
}

/**
 * EngineControl over a serial reactor that `reactor` opens (and reopens after a stop), the shared
 * instance; `instantiate` opens a new instance of the serial module, which a validation with a
 * registry file takes for itself.
 */
export function reactorControl(reactor: () => Promise<ReactorAbi>, instantiate: () => Promise<ReactorAbi>): EngineControl {
  return {
    async describe(program) {
      return toProgramDescription((await reactor()).describe(program));
    },
    async start(argv): Promise<StartResult> {
      return (await reactor()).start(argv);
    },
    async validate(argv, registry): Promise<ValidationResult> {
      let abi: ReactorAbi;
      if (registry === undefined) abi = await reactor();
      else {
        // The validation's own instance, never the shared one; one that cannot be opened fails the
        // validation with its reason. The file is read in chunks from the Blob (`loadFile`); a file
        // that the engine cannot hold is refused with its message (`not enough memory for the registry
        // file ...`). Nothing is released: the instance goes after the validation, with the file.
        try {
          abi = await instantiate();
        } catch (error) {
          throw new Error(`the engine could not open an instance for the file of -conffile: ${error instanceof Error ? error.message : String(error)}`, {
            cause: error,
          });
        }
        await loadFile(abi, FILE_REGISTRY, registry.name, registry.file);
      }
      const message = abi.validate(argv);
      return message === undefined ? { ok: true } : { ok: false, message };
    },
  };
}
