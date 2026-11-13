// Example sets (line WF-2, docs/web/fetch_design.md §3, §9): for each mode, a public subject set
// that LOSAT downloads only when the user asks ("Example"), and a query that is part of the
// table. One click sends one request, for the file's URL, to the host that the in-place sentence
// names; nothing of the user's is sent (web/AGENTS.md rule 4). The files change with the
// sources' releases, so the table holds no hash: the saved copy records what came
// (ports/example.ts). The queries tell one story: Bacillus subtilis 168 against RefSeq and
// Escherichia coli K-12 (fetched once from NCBI and UniProt on 2026-11-03; the recA CDS
// translates to P16971).
import type { ProgramId, SequenceKind } from 'losat/domain';
import { programById, sequenceKind } from './programs';

export type ExampleId = 'refseq-16s-bacteria' | 'refseq-its-fungi' | 'uniprot-ecoli-k12' | 'pdbaa' | 'ecoli-k12-genome';

export interface ExampleSet {
  readonly id: ExampleId;
  /** The words of the set's button. */
  readonly label: string;
  /** The modes that the set serves, in the order of the mode tabs. */
  readonly programs: readonly ProgramId[];
  /** The subject set: a public file (gzip), downloaded only when the user asks. */
  readonly url: string;
  /** The name of the subject source: the file's name without ".gz". */
  readonly fileName: string;
  /** Bytes downloaded, as measured on 2026-11-03, for the in-place sentence (a release may change it). */
  readonly downloadBytes: number;
  /** Where the set comes from and on what terms, for the help. */
  readonly provenance: string;
  /** The query of each kind that the set's modes search with, as FASTA. */
  readonly queries: Readonly<Partial<Record<SequenceKind, string>>>;
}

/** Bytes that one Example keeps at most once decompressed (pdbaa, the largest, is 101 MB). */
export const EXAMPLE_BYTE_LIMIT = 256 * 1024 * 1024;

const BSU_16S = [
  '>NR_102783.2 Bacillus subtilis subsp. subtilis strain 168 16S ribosomal RNA, complete sequence',
  'TTATCGGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGA',
  'CAGATGGGAGCTTGCTCCCTGATGTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTGTAAGA',
  'CTGGGATAACTCCGGGAAACCGGGGCTAATACCGGATGGTTGTTTGAACCGCATGGTTCAAACATAAAAG',
  'GTGGCTTCGGCTACCACTTACAGATGGACCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAA',
  'GGCGACGATGCGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTA',
  'CGGGAGGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGAGTGATGA',
  'AGGTTTTCGGATCGTAAAGCTCTGTTGTTAGGGAAGAACAAGTGCCGTTCGAATAGGGCGGTACCTTGAC',
  'GGTACCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTG',
  'TCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTTCTTAAGTCTGATGTGAAAGCCCCCGGCTCAAC',
  'CGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTG',
  'AAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGAG',
  'CGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGT',
  'TAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAG',
  'ACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGC',
  'GAAGAACCTTACCAGGTCTTGACATCCTCTGACAATCCTAGAGATAGGACGTCCCCTTCGGGGGCAGAGT',
  'GACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAAC',
  'CCTTGATCTTAGTTGCCAGCATTCAGTTGGGCACTCTAAGGTGACTGCCGGTGACAAACCGGAGGAAGGT',
  'GGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACAGAACAAAG',
  'GGCAGCGAAACCGCGAGGTTAAGCCAATCCCACAAATCTGTTCTCAGTTCGGATCGCAGTCTGCAACTCG',
  'ACTGCGTGAAGCTGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGT',
  'ACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTTAGGAGCCAGC',
  'CGCCGAAGGTGGGACAGATGATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGAT',
  'CACCTCCTTT',
].join('\n') + '\n';

const SCE_ITS = [
  '>NR_111007.1 Saccharomyces cerevisiae CBS 1171 ITS region; from TYPE material',
  'AAAGAAATTTAATAATTTTGAAAATGGATTTTTTTGTTTTGGCAAGAGCATGAGAGCTTTTACTGGGCAA',
  'GAAGACAAGAGATGGAGAGTCCAGCCGGGCCTGCGCTTAAGTGCGCGGTCTTGCTAGGCTTGTAAGTTTC',
  'TTTCTTGCTATTCCAAACGGTGAGAGATTTCTGTGCTTTTGTTATAGGACAATTAAAACCGTTTCAATAC',
  'AACACACTGTGGAGTTTTCATATCTTTGCAACTTTTTCTTTGGGCATTCGAGCAATCGGGGCCCAGAGGT',
  'AACAAACACAAACAATTTTATCTATTCATTAAATTTTTGTCAAAAACAAGAATTTTCGTAACTGGAAATT',
  'TTAAAAATATTAAAAACTTTCAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGCGAAATGC',
  'GATACGTAATGTGAATTGCAGAATTCCGTGAATCATCGAATCTTTGAACGCACATTGCCCCCTTGGTATT',
  'CCAGGGGGCATGCCTGTTTGAGCGTCATTTCCTTCTCAAACATTCTGTTTGGTAGTGAGTGATACTCTTT',
  'GGAGTTAACTTGAAATTGCTGGCCTTTTCATTGGATGTTTTTTTTCCAAAGAGAGGTTTCTCTGCGTGCT',
  'TGAGGTATAATGCAAGTACGGTCGTTTTAGGTTTTACCAACTGCGGCTAATCTTTTTTATACTGAGCGTA',
  'TTGGAACGTTATCGATAAGAAGAGAGCGTCTAGGCAACAATGTTCTTAAAGT',
].join('\n') + '\n';

const BSU_RECA_PROTEIN = [
  '>sp|P16971|RECA_BACSU Protein RecA OS=Bacillus subtilis (strain 168) OX=224308 GN=recA PE=1 SV=2',
  'MSDRQAALDMALKQIEKQFGKGSIMKLGEKTDTRISTVPSGSLALDTALGIGGYPRGRII',
  'EVYGPESSGKTTVALHAIAEVQQQGGQAAFIDAEHALDPVYAQKLGVNIEELLLSQPDTG',
  'EQALEIAEALVRSGAVDIVVVDSVAALVPKAEIEGDMGDSHVGLQARLMSQALRKLSGAI',
  'NKSKTIAIFINQIREKVGVMFGNPETTPGGRALKFYSSVRLEVRRAEQLKQGNDVMGNKT',
  'KIKVVKNKVAPPFRTAEVDIMYGEGISKEGEIIDLGTELDIVQKSGSWYSYEEERLGQGR',
  'ENAKQFLKENKDIMLMIQEQIREHYGLDNNGVVQQQAEETQEELEFEE',
].join('\n') + '\n';

const BSU_RECA_CDS = [
  '>AL009126.3:1764645-1765691 Bacillus subtilis subsp. subtilis str. 168 genome assembly, chromosome: 1',
  'ATGAGTGATCGTCAGGCAGCCTTAGATATGGCTCTTAAACAAATAGAAAAACAGTTCGGCAAAGGTTCCA',
  'TTATGAAACTGGGAGAAAAGACAGATACAAGAATTTCTACTGTACCAAGCGGCTCCCTCGCTCTTGATAC',
  'AGCACTGGGAATTGGCGGATATCCTCGCGGACGGATTATTGAAGTATACGGTCCTGAAAGCTCAGGTAAA',
  'ACAACTGTGGCGCTTCATGCGATTGCTGAAGTTCAGCAGCAGGGCGGACAAGCCGCGTTTATCGATGCGG',
  'AGCATGCGTTAGATCCGGTATACGCGCAAAAGCTCGGTGTTAACATCGAAGAGCTTTTACTGTCTCAGCC',
  'TGACACAGGCGAGCAGGCGCTTGAAATTGCGGAAGCATTGGTTCGAAGCGGGGCAGTTGACATTGTCGTT',
  'GTCGACTCTGTAGCCGCTCTCGTTCCGAAAGCGGAAATTGAAGGCGACATGGGAGATTCGCATGTCGGTT',
  'TACAAGCACGCTTAATGTCTCAAGCGCTTCGTAAGCTTTCAGGGGCCATTAACAAATCGAAGACAATCGC',
  'GATTTTCATTAACCAAATTCGTGAAAAAGTCGGTGTTATGTTCGGGAACCCGGAAACAACACCTGGCGGC',
  'CGTGCGTTGAAATTCTATTCTTCCGTGCGTCTTGAAGTGCGCCGTGCTGAACAGCTGAAACAAGGCAACG',
  'ACGTAATGGGGAACAAAACGAAAATCAAAGTCGTGAAAAACAAGGTGGCTCCGCCGTTCCGTACAGCCGA',
  'GGTTGACATTATGTACGGAGAAGGCATTTCAAAAGAAGGCGAAATCATTGATCTAGGAACTGAACTTGAT',
  'ATCGTGCAAAAAAGCGGTTCATGGTACTCTTATGAAGAAGAGCGTCTTGGCCAAGGCCGTGAAAATGCAA',
  'AACAATTCCTGAAAGAAAATAAAGATATCATGCTGATGATCCAGGAGCAAATTCGCGAACATTACGGCTT',
  'GGATAATAACGGAGTAGTGCAGCAGCAAGCTGAAGAGACACAAGAAGAACTCGAATTTGAAGAATAA',
].join('\n') + '\n';

const NCBI_FTP = 'https://ftp.ncbi.nlm.nih.gov';
const REFSEQ_TERMS = 'NCBI places no restrictions on the use of RefSeq (ncbi.nlm.nih.gov/home/about/policies).';

/** The sets, in the order of their buttons within a mode. */
export const EXAMPLE_SETS: readonly ExampleSet[] = Object.freeze([
  {
    id: 'refseq-16s-bacteria',
    label: '16S rRNA (bacteria)',
    programs: ['blastn'],
    url: `${NCBI_FTP}/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz`,
    fileName: 'bacteria.16SrRNA.fna',
    downloadBytes: 8_038_469,
    provenance: `RefSeq Targeted Loci, bacterial 16S rRNA (NCBI). ${REFSEQ_TERMS}`,
    queries: { nucleotide: BSU_16S },
  },
  {
    id: 'refseq-its-fungi',
    label: 'ITS (fungi)',
    programs: ['blastn'],
    url: `${NCBI_FTP}/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz`,
    fileName: 'fungi.ITS.fna',
    downloadBytes: 2_956_620,
    provenance: `RefSeq Targeted Loci, fungal ITS (NCBI). ${REFSEQ_TERMS}`,
    queries: { nucleotide: SCE_ITS },
  },
  {
    id: 'uniprot-ecoli-k12',
    label: 'E. coli K-12 proteins (UniProt reviewed)',
    programs: ['blastp', 'blastx'],
    url: 'https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:83333&format=fasta&compressed=true',
    fileName: 'uniprot_reviewed_ecoli_k12.fasta',
    downloadBytes: 981_194,
    provenance: 'UniProtKB/Swiss-Prot entries of Escherichia coli K-12 (taxon 83333), UniProt, CC BY 4.0.',
    queries: { protein: BSU_RECA_PROTEIN, nucleotide: BSU_RECA_CDS },
  },
  {
    id: 'pdbaa',
    label: 'PDB proteins (pdbaa)',
    programs: ['blastp', 'blastx'],
    url: `${NCBI_FTP}/blast/db/FASTA/pdbaa.gz`,
    fileName: 'pdbaa',
    downloadBytes: 36_755_745,
    provenance: 'Protein chains of the Protein Data Bank, as NCBI BLAST\'s pdbaa (NCBI). PDB data are CC0.',
    queries: { protein: BSU_RECA_PROTEIN, nucleotide: BSU_RECA_CDS },
  },
  {
    id: 'ecoli-k12-genome',
    label: 'E. coli K-12 genome',
    programs: ['tblastn', 'tblastx'],
    url: `${NCBI_FTP}/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz`,
    fileName: 'GCF_000005845.2_ASM584v2_genomic.fna',
    downloadBytes: 1_379_902,
    provenance: `Escherichia coli K-12 MG1655, RefSeq assembly GCF_000005845.2 (NCBI). ${REFSEQ_TERMS}`,
    queries: { protein: BSU_RECA_PROTEIN, nucleotide: BSU_RECA_CDS },
  },
] satisfies ExampleSet[]);

export function exampleById(id: ExampleId): ExampleSet {
  const set = EXAMPLE_SETS.find((candidate) => candidate.id === id);
  if (set === undefined) throw new Error(`unknown example ${id}`);
  return set;
}

/** The sets of a mode. */
export function examplesFor(program: ProgramId): readonly ExampleSet[] {
  return EXAMPLE_SETS.filter((set) => set.programs.includes(program));
}

/** The query of a set for a mode: the FASTA of the mode's query kind. */
export function exampleQuery(set: ExampleSet, program: ProgramId): string {
  const query = set.queries[sequenceKind(programById(program), 'query')];
  if (query === undefined) throw new Error(`example ${set.id} has no query for ${program}`);
  return query;
}

/** The host that a set is downloaded from, as the in-place sentence names it. */
export function exampleHost(set: ExampleSet): string {
  return new URL(set.url).host;
}
