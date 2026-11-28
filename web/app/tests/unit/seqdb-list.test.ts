// The one writer, reader and refusal of the SeqDB name lists that LOSAT Web writes (line WW, WW4d): what NCBI BLAST+
// 2.17.0 does with each form was measured for every option (WW/ww4d/probe/probe-compare.tsv, 2026-11-28).
import { describe, expect, it } from 'vitest';
import { isSeqdbList, listNameRefusal, SEQDB_LIST_FLAGS, seqdbNameList, seqdbNames } from '../../src/domain/list-flags';

describe('SeqDB name lists (line WW, WW4d)', () => {
  it('writes a list as the command line reads it: names joined by spaces, a name with a space between double quotes', () => {
    expect(seqdbNameList('-db', ['ndb'])).toBe('ndb');
    expect(seqdbNameList('-db', ['my dbs/ndb'])).toBe('"my dbs/ndb"');
    expect(seqdbNameList('-db', ['ndb', 'tdb'])).toBe('ndb tdb');
    expect(seqdbNameList('-filtering_db', ['my dbs/ndb', 'tdb'])).toBe('"my dbs/ndb" tdb');
    expect(seqdbNameList('-in', ['input.fa', 'my genome.fa.gz'], 'file')).toBe('input.fa "my genome.fa.gz"');
    // A tab is no separator of the reader, so the writer leaves it as it is.
    expect(seqdbNameList('-db', ['a\tb'])).toBe('a\tb');
    // Only an empty name is left out (it reads back as nothing); a name of spaces is a name, quoted as
    // SeqDB_CombineAndQuote quotes any name with a space (line WW, WW4d R1, review finding 4).
    expect(seqdbNameList('-db', [''])).toBe('');
    expect(seqdbNameList('-db', ['a', '', 'b'])).toBe('a b');
    expect(seqdbNameList('-db', [' '])).toBe('" "');
    expect(seqdbNameList('-db', ['', ' ', 'b'])).toBe('" " b');
  });

  it('throws for a name that the list cannot carry, with the reason', () => {
    expect(() => seqdbNameList('-db', ['q"d/ndb'])).toThrow('A database named "q"d/ndb" cannot be in -db: a double quote there starts or ends a name.');
    expect(() => seqdbNameList('-in', ['a"b.fa'], 'file')).toThrow('A file named "a"b.fa" cannot be in -in: a double quote there starts or ends a name.');
  });

  it("reads a list back as the command line does, for each form of the probe", () => {
    expect(seqdbNames('ndb')).toEqual(['ndb']);
    expect(seqdbNames('my dbs/ndb')).toEqual(['my', 'dbs/ndb']);
    expect(seqdbNames('"my dbs/ndb"')).toEqual(['my dbs/ndb']);
    expect(seqdbNames('"ndb"')).toEqual(['ndb']);
    expect(seqdbNames('q"d/ndb')).toEqual(['q', 'd/ndb']);
    expect(seqdbNames('ndb tdb')).toEqual(['ndb', 'tdb']);
    expect(seqdbNames('"my dbs/ndb" tdb')).toEqual(['my dbs/ndb', 'tdb']);
    expect(seqdbNames('  a   b ')).toEqual(['a', 'b']);
    expect(seqdbNames('"my db')).toEqual(['my db']);
    expect(seqdbNames('')).toEqual([]);
    for (const names of [['ndb'], ['my dbs/ndb'], ['ndb', 'my dbs/tdb', 'x y z/w']]) {
      expect(seqdbNames(seqdbNameList('-db', names))).toEqual(names);
    }
  });

  it('knows which options are SeqDB lists, per tool where the flag means other things', () => {
    expect([...SEQDB_LIST_FLAGS].sort()).toEqual(['-db', '-filtering_db', '-in']);
    expect(isSeqdbList({ flag: '-db' })).toBe(true);
    expect(isSeqdbList({ flag: '-filtering_db' })).toBe(true);
    expect(isSeqdbList({ command: 'makeblastdb', flag: '-in' })).toBe(true);
    expect(isSeqdbList({ command: 'windowmasker', flag: '-in', inputFormat: 'blastdb' })).toBe(true);
    expect(isSeqdbList({ command: 'windowmasker', flag: '-in', inputFormat: 'fasta' })).toBe(false);
    expect(isSeqdbList({ command: 'windowmasker', flag: '-in' })).toBe(false);
    expect(isSeqdbList({ command: 'blast_formatter', flag: '-in' })).toBe(false);
    expect(isSeqdbList({ command: 'makeblastdb', flag: '-out' })).toBe(false);
    expect(isSeqdbList({ command: 'makeblastdb', flag: '-mask_data' })).toBe(false);
  });

  it('refuses a name with a double quote in every SeqDB list, and nothing in a list of one file name', () => {
    const reason = 'a double quote there starts or ends a name';
    expect(listNameRefusal({ flag: '-db' }, 'q"d/ndb', 'database')).toBe(`A database named "q"d/ndb" cannot be in -db: ${reason}.`);
    expect(listNameRefusal({ flag: '-filtering_db' }, 'q"d/rep', 'database')).toBe(`A database named "q"d/rep" cannot be in -filtering_db: ${reason}.`);
    expect(listNameRefusal({ command: 'makeblastdb', flag: '-in' }, 'a"b.fa')).toBe(`A file named "a"b.fa" cannot be in -in: ${reason}.`);
    expect(listNameRefusal({ command: 'windowmasker', flag: '-in', inputFormat: 'blastdb' }, 'q"d/wm', 'database')).toBe(
      `A database named "q"d/wm" cannot be in -in: ${reason}.`,
    );
    // A space is no refusal: the writer quotes the name.
    expect(listNameRefusal({ flag: '-db' }, 'my dbs/ndb', 'database')).toBeUndefined();
    expect(listNameRefusal({ flag: '-filtering_db' }, 'my dbs/rep', 'database')).toBeUndefined();
    // WindowMasker's -in of FASTA input is one file's name.
    expect(listNameRefusal({ command: 'windowmasker', flag: '-in', inputFormat: 'fasta' }, 'a"b.fa')).toBeUndefined();
    expect(listNameRefusal({ command: 'windowmasker', flag: '-in' }, 'a"b.fa')).toBeUndefined();
    // makeblastdb's comma lists.
    expect(listNameRefusal({ command: 'makeblastdb', flag: '-mask_data' }, 'a,b.asnb')).toBe('A file named "a,b.asnb" cannot be in -mask_data: a comma separates its files.');
  });
});
