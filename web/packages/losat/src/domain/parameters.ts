// The defaults that the engine resolves (`describe`, docs/web/abi_v2.md §9).

/**
 * The default that the engine resolves for an option that a search leaves out (`describe`'s
 * `defaults`, docs/web/abi_v2.md §9), by task and from other options.
 */
export interface ResolvedDefault {
  /** The value of each task with no other option given (the key '' for a program without -task). */
  readonly byTask: Readonly<Record<string, string>>;
  /**
   * For each option that changes this one (-matrix, -word_size): the value by task for each
   * value of that option that the engine accepts.
   */
  readonly by?: Readonly<Record<string, Readonly<Record<string, Readonly<Record<string, string>>>>>>;
  /** The value holds only while none of these options is given (the hit list options). */
  readonly unless?: readonly string[];
  /**
   * The tasks where BLAST+ tells the option given from not given: a value set is written even
   * when it equals the default.
   */
  readonly given?: readonly string[];
}
