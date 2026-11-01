"use strict";
// Writes guarded.wasm: the page's probe module through the fixed guard. usage: node mk_guarded.js <wasi_shared_memory.js> <fences>
const fs = require("fs");
const html = fs.readFileSync(`${__dirname}/index.html`, "utf8");
const code = html.slice(html.indexOf("<script>") + 8, html.indexOf("async function run"));
const bytes = new Function(`${code}; return moduleBytes;`)()(Number(process.argv[3]));
const guarded = require(process.argv[2]).guardSharedMemory(Buffer.from(bytes));
fs.writeFileSync(`${__dirname}/guarded.wasm`, guarded.bytes);
console.log(guarded.record);
