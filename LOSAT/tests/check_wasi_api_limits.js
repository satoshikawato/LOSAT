"use strict";
const fs = require("node:fs");
const assert = require("node:assert/strict");
const { WASI } = require("node:wasi");
const { prepareArtifact } = require("./wasi_artifact");
const { createThreadHost } = require("./wasi_thread_host");
const { runPair } = require("./check_wasi_reactor");

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3152-3187
// arg_desc.SetConstraint(kArgNumThreads, new CArgAllowValuesGreaterThanOrEqual(1));
// Explicit unsupported requests and host limits fail before any worker starts.
(async () => {
  const [artifact, kind, fixtures, output] = process.argv.slice(2);
  fs.mkdirSync(output,{recursive:true});
  let host;
  if (kind === "serial-reactor") {
    // NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
    // CATCH_ALL(status); ... return status;
    // Reuse the validated module; reactor initialization and ownership are unchanged.
    const bytes=fs.readFileSync(artifact), {module}=prepareArtifact(bytes,kind);
    const wasi=new WASI({version:"preview1",env:process.env,preopens:{"/":"/"},returnOnExit:true});
    const instance=await WebAssembly.instantiate(module,{wasi_snapshot_preview1:wasi.wasiImport});
    wasi.initialize(instance);
    host={instance,memory:instance.exports.memory,events:[],waitForWorkers:async()=>{},close:async()=>{}};
  } else host=await createThreadHost(artifact,kind);
  const records=[];
  try {
    for (const program of ["blastn","tblastx","blastp"]) {
      const seq=fs.readFileSync(`${fixtures}/${program === "blastp"?"aa1":"nuc1"}.fasta`,"utf8");
      const cap=process.env.LOSAT_WASI_THREAD_CAP;
      for (const n of [1,2,4,1]) {
        const response=runPair(host,program,seq,seq,"6",["-num_threads",String(n)]);
        const events=await host.waitForWorkers() || [];
        const expectedSuccess=cap === undefined ? n===1 : /^[1-9]\d*$/.test(cap) && n<=Number(cap);
        assert.equal(response.status,expectedSuccess?0:-1,response.error);
        assert.equal(events.length,0);
        if (!expectedSuccess) {
          assert.equal(response.result.length,0);
          assert.match(response.error,cap===undefined?/unsupported.*parallel/:/LOSAT_WASI_THREAD_CAP/);
        }
        const id=`${program}-${records.length}-n${n}`;
        fs.writeFileSync(`${output}/${id}.out`,response.result);
        fs.writeFileSync(`${output}/${id}.error`,response.error);
        records.push({id,status:response.status,events});
        fs.writeFileSync(`${output}/runs.json`,JSON.stringify(records,null,2));
      }
    }
  } finally {await host.close();}
})().catch(error=>{console.error(error);process.exitCode=1;});
