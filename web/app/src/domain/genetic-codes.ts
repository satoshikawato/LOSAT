// Names of NCBI's genetic codes, for the labels of the genetic-code fields only. They are
// the `name` entries of NCBI's c++/src/objects/seqfeat/gc.prt. Which codes a program
// accepts is the engine's list (`describe` `query_gencodes` / `subject_gencodes`).
const NAMES: Readonly<Record<number, string>> = {
  1: 'Standard',
  2: 'Vertebrate Mitochondrial',
  3: 'Yeast Mitochondrial',
  4: 'Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma',
  5: 'Invertebrate Mitochondrial',
  6: 'Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear',
  9: 'Echinoderm Mitochondrial; Flatworm Mitochondrial',
  10: 'Euplotid Nuclear',
  11: 'Bacterial, Archaeal and Plant Plastid',
  12: 'Alternative Yeast Nuclear',
  13: 'Ascidian Mitochondrial',
  14: 'Alternative Flatworm Mitochondrial',
  15: 'Blepharisma Macronuclear',
  16: 'Chlorophycean Mitochondrial',
  21: 'Trematode Mitochondrial',
  22: 'Scenedesmus obliquus Mitochondrial',
  23: 'Thraustochytrium Mitochondrial',
  24: 'Rhabdopleuridae Mitochondrial',
  25: 'Candidate Division SR1 and Gracilibacteria',
  26: 'Pachysolen tannophilus Nuclear',
  27: 'Karyorelict Nuclear',
  28: 'Condylostoma Nuclear',
  29: 'Mesodinium Nuclear',
  30: 'Peritrich Nuclear',
  31: 'Blastocrithidia Nuclear',
  32: 'Balanophoraceae Plastid',
  33: 'Cephalodiscidae Mitochondrial',
};

/** "11. Bacterial, Archaeal and Plant Plastid", or the number alone for a code without a name here. */
export function geneticCodeLabel(id: number): string {
  const name = NAMES[id];
  return name === undefined ? String(id) : `${id}. ${name}`;
}
