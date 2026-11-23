// Reader for NCBI's PrintUsageXml output (`<tool> -xmlhelp`), used by `ToolEngine.describeTool`.
// Pure TypeScript with no DOMParser so that it runs identically in browsers, workers and Node.

export type ToolParameterKind = 'Boolean' | 'File_In' | 'File_Out' | 'Int8' | 'Integer' | 'Real' | 'String' | 'Flag';

const KEY_KINDS: ReadonlySet<string> = new Set<ToolParameterKind>([
  'Boolean',
  'File_In',
  'File_Out',
  'Int8',
  'Integer',
  'Real',
  'String',
]);

export interface ToolParameterRange {
  readonly min: number;
  readonly max: number;
  readonly kind: 'integer' | 'double';
}

export interface ToolParameterDescription {
  /** Option name with leading `-` (for example `-unit`, `-mk_counts`). */
  readonly flag: string;
  readonly help: string;
  readonly takesValue: boolean;
  readonly kind: ToolParameterKind;
  readonly mandatory: boolean;
  readonly group?: string;
  readonly synopsis?: string;
  readonly defaultValue?: string;
  readonly choices?: readonly string[];
  readonly caseSensitive?: boolean;
  readonly range?: ToolParameterRange;
  readonly constraintDescription?: string;
  readonly alias?: string;
  readonly hidden?: boolean;
  readonly append?: boolean;
}

export interface ToolDescription {
  readonly tool: string;
  readonly version: string;
  readonly description: string;
  readonly detailedDescription?: string;
  readonly parameters: readonly ToolParameterDescription[];
  /** Pairs `[first, second]` of `-flag` names from `<first_excludes_second>`. */
  readonly excludes: readonly (readonly [string, string])[];
  /** Pairs `[first, second]` of `-flag` names from `<first_requires_second>`. */
  readonly requires: readonly (readonly [string, string])[];
}

interface XmlElement {
  readonly tag: string;
  readonly attrs: Readonly<Record<string, string>>;
  readonly children: readonly XmlElement[];
  readonly text: string;
}

const NAMED_ENTITIES: Readonly<Record<string, string>> = {
  amp: '&',
  apos: "'",
  gt: '>',
  lt: '<',
  quot: '"',
};

function decodeNumericEntity(whole: string, raw: string, radix: number): string {
  const code = Number.parseInt(raw, radix);
  if (Number.isInteger(code) && code >= 0 && code <= 0x10ffff) return String.fromCodePoint(code);
  throw new Error(`invalid numeric XML entity ${whole}`);
}

/** Decodes XML character and named references in a single pass. */
export function decodeXmlEntities(text: string): string {
  if (!text.includes('&')) return text;
  return text.replace(/&(?:#x([0-9a-fA-F]+)|#(\d+)|([a-zA-Z]+));/g, (whole, hex: string | undefined, dec: string | undefined, named: string | undefined) => {
    if (hex !== undefined) return decodeNumericEntity(whole, hex, 16);
    if (dec !== undefined) return decodeNumericEntity(whole, dec, 10);
    const decoded = named === undefined ? undefined : NAMED_ENTITIES[named];
    if (decoded === undefined) throw new Error(`unknown XML entity ${whole}`);
    return decoded;
  });
}

function isWhitespace(code: number): boolean {
  return code === 0x20 || code === 0x09 || code === 0x0a || code === 0x0d;
}

function isNameChar(code: number): boolean {
  return (
    (code >= 0x41 && code <= 0x5a) ||
    (code >= 0x61 && code <= 0x7a) ||
    (code >= 0x30 && code <= 0x39) ||
    code === 0x5f || // _
    code === 0x2d || // -
    code === 0x3a || // :
    code === 0x2e // .
  );
}

/** Parses an XML document of the PrintUsageXml subset into its root element. */
function parseXml(input: string): XmlElement {
  let pos = 0;
  const len = input.length;

  const skipWs = (): void => {
    while (pos < len && isWhitespace(input.charCodeAt(pos))) pos++;
  };

  const skipComment = (): void => {
    const end = input.indexOf('-->', pos + 4);
    if (end < 0) throw new Error('unterminated XML comment');
    pos = end + 3;
  };

  const skipProlog = (): void => {
    while (pos < len) {
      skipWs();
      if (input.startsWith('<?', pos)) {
        const end = input.indexOf('?>', pos + 2);
        if (end < 0) throw new Error('unterminated XML processing instruction');
        pos = end + 2;
      } else if (input.startsWith('<!--', pos)) {
        skipComment();
      } else {
        break;
      }
    }
  };

  const readName = (): string => {
    const start = pos;
    while (pos < len && isNameChar(input.charCodeAt(pos))) pos++;
    if (pos === start) throw new Error(`expected XML name at offset ${pos}`);
    return input.slice(start, pos);
  };

  const readElement = (): XmlElement => {
    if (input.charCodeAt(pos) !== 0x3c) throw new Error(`expected '<' at offset ${pos}`);
    pos++;
    const tag = readName();
    const attrs: Record<string, string> = {};
    while (pos < len) {
      skipWs();
      if (input.startsWith('/>', pos)) {
        pos += 2;
        return { tag, attrs, children: [], text: '' };
      }
      if (input.charCodeAt(pos) === 0x3e) {
        pos++;
        break;
      }
      const attrName = readName();
      skipWs();
      if (input.charCodeAt(pos) !== 0x3d) throw new Error(`expected '=' after attribute ${attrName} at offset ${pos}`);
      pos++;
      skipWs();
      const quote = input.charCodeAt(pos);
      if (quote !== 0x22 && quote !== 0x27) throw new Error(`expected quoted attribute value for ${attrName} at offset ${pos}`);
      pos++;
      const valueStart = pos;
      const valueEnd = input.indexOf(quote === 0x22 ? '"' : "'", valueStart);
      if (valueEnd < 0) throw new Error(`unterminated attribute value for ${attrName}`);
      attrs[attrName] = decodeXmlEntities(input.slice(valueStart, valueEnd));
      pos = valueEnd + 1;
    }

    const children: XmlElement[] = [];
    let text = '';
    while (pos < len) {
      if (input.startsWith('</', pos)) {
        pos += 2;
        const closeTag = readName();
        skipWs();
        if (input.charCodeAt(pos) !== 0x3e) throw new Error(`expected '>' after closing tag </${closeTag}> at offset ${pos}`);
        pos++;
        if (closeTag !== tag) throw new Error(`mismatched closing tag </${closeTag}> for <${tag}>`);
        return { tag, attrs, children, text };
      }
      if (input.startsWith('<!--', pos)) {
        skipComment();
      } else if (input.startsWith('<![CDATA[', pos)) {
        const end = input.indexOf(']]>', pos + 9);
        if (end < 0) throw new Error('unterminated CDATA section');
        text += input.slice(pos + 9, end);
        pos = end + 3;
      } else if (input.charCodeAt(pos) === 0x3c) {
        children.push(readElement());
      } else {
        const next = input.indexOf('<', pos);
        if (next < 0) throw new Error(`unterminated element <${tag}>`);
        text += decodeXmlEntities(input.slice(pos, next));
        pos = next;
      }
    }
    throw new Error(`unterminated element <${tag}>`);
  };

  skipProlog();
  if (pos >= len) throw new Error('empty XML document');
  const root = readElement();
  skipProlog();
  if (pos < len) throw new Error(`unexpected trailing content at offset ${pos}`);
  return root;
}

function child(element: XmlElement, tag: string): XmlElement | undefined {
  return element.children.find((entry) => entry.tag === tag);
}

function childText(element: XmlElement, tag: string): string | undefined {
  const entry = child(element, tag);
  return entry === undefined ? undefined : entry.text;
}

function parseNumberBound(raw: string | undefined, label: string): number {
  if (raw === undefined || raw.trim() === '') throw new Error(`missing ${label} bound in XML constraint`);
  const parsed = Number(raw.trim());
  if (!Number.isFinite(parsed)) throw new Error(`invalid ${label} bound ${JSON.stringify(raw)} in XML constraint`);
  return parsed;
}

function parseConstraint(
  constraint: XmlElement,
): Pick<ToolParameterDescription, 'choices' | 'caseSensitive' | 'range' | 'constraintDescription'> {
  const constraintDescription = childText(constraint, 'description');
  const desc = constraintDescription === undefined ? {} : { constraintDescription };
  const strings = child(constraint, 'Strings');
  if (strings !== undefined) {
    const choices = strings.children.filter((entry) => entry.tag === 'value').map((entry) => entry.text);
    const caseSensitive = strings.attrs['case_sensitive'] === undefined ? undefined : strings.attrs['case_sensitive'] === 'true';
    return {
      choices,
      ...(caseSensitive === undefined ? {} : { caseSensitive }),
      ...desc,
    };
  }
  const integers = child(constraint, 'Integers');
  if (integers !== undefined) {
    return {
      range: {
        min: parseNumberBound(childText(integers, 'min'), 'min'),
        max: parseNumberBound(childText(integers, 'max'), 'max'),
        kind: 'integer',
      },
      ...desc,
    };
  }
  const doubles = child(constraint, 'Doubles');
  if (doubles !== undefined) {
    return {
      range: {
        min: parseNumberBound(childText(doubles, 'min'), 'min'),
        max: parseNumberBound(childText(doubles, 'max'), 'max'),
        kind: 'double',
      },
      ...desc,
    };
  }
  return desc;
}

function parseParameter(element: XmlElement): ToolParameterDescription {
  const name = element.attrs['name'];
  if (!name) throw new Error(`<${element.tag}> is missing a name attribute`);
  const takesValue = element.tag === 'key';
  let kind: ToolParameterKind;
  if (takesValue) {
    const rawType = element.attrs['type'];
    if (rawType === undefined || !KEY_KINDS.has(rawType)) {
      throw new Error(`<key name="${name}"> has unsupported type ${JSON.stringify(rawType)}`);
    }
    kind = rawType as ToolParameterKind;
  } else {
    kind = 'Flag';
  }
  const mandatory = element.attrs['optional'] !== 'true';
  const help = childText(element, 'description') ?? '';
  const group = childText(element, 'group');
  const synopsis = childText(element, 'synopsis');
  const defaultElement = child(element, 'default');
  const alias = childText(element, 'alias');
  const constraint = child(element, 'constraint');
  const parsedConstraint = constraint === undefined ? {} : parseConstraint(constraint);
  const flags = child(element, 'flags');
  const hidden = flags !== undefined && child(flags, 'hidden') !== undefined;
  const append = flags !== undefined && child(flags, 'append') !== undefined;
  return {
    flag: `-${name}`,
    help,
    takesValue,
    kind,
    mandatory,
    ...(group === undefined ? {} : { group }),
    ...(synopsis === undefined ? {} : { synopsis }),
    ...(defaultElement === undefined ? {} : { defaultValue: defaultElement.text }),
    ...parsedConstraint,
    ...(alias === undefined ? {} : { alias }),
    ...(hidden ? { hidden: true } : {}),
    ...(append ? { append: true } : {}),
  };
}

function parseDependencyPair(element: XmlElement): readonly [string, string] {
  const arg1 = childText(element, 'arg1');
  const arg2 = childText(element, 'arg2');
  if (!arg1 || !arg2) throw new Error(`<${element.tag}> is missing arg1 or arg2`);
  return [`-${arg1}`, `-${arg2}`];
}

/** Parses NCBI PrintUsageXml (`-xmlhelp`) output into a typed `ToolDescription`. */
export function parseToolDescription(xml: string): ToolDescription {
  const root = parseXml(xml);
  if (root.tag !== 'ncbi_application') throw new Error(`expected <ncbi_application> root element, got <${root.tag}>`);
  const program = child(root, 'program');
  if (program === undefined) throw new Error('missing <program> element in -xmlhelp output');
  const tool = childText(program, 'name');
  if (!tool) throw new Error('missing <program><name> in -xmlhelp output');
  const version = childText(program, 'version') ?? '';
  const description = childText(program, 'description') ?? '';
  const detailed = childText(program, 'detailed_description');
  const argumentsElement = child(root, 'arguments');
  if (argumentsElement === undefined) throw new Error('missing <arguments> element in -xmlhelp output');

  const parameters: ToolParameterDescription[] = [];
  const excludes: Array<readonly [string, string]> = [];
  const requires: Array<readonly [string, string]> = [];

  for (const entry of argumentsElement.children) {
    if (entry.tag === 'key' || entry.tag === 'flag') {
      parameters.push(parseParameter(entry));
    } else if (entry.tag === 'dependencies') {
      for (const dep of entry.children) {
        if (dep.tag === 'first_excludes_second') excludes.push(parseDependencyPair(dep));
        else if (dep.tag === 'first_requires_second') requires.push(parseDependencyPair(dep));
      }
    }
  }

  return {
    tool,
    version,
    description,
    ...(detailed ? { detailedDescription: detailed } : {}),
    parameters,
    excludes,
    requires,
  };
}
