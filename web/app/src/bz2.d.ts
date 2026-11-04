// The `bz2` package (1.0.1) ships no types. It runs only in a Worker or Node (see infra/data/inflate.ts).
declare module 'bz2' {
  const bz2: { decompress(bytes: Uint8Array, checkCRC?: boolean): Uint8Array };
  export default bz2;
}
