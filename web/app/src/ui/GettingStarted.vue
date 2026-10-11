<script setup lang="ts">
// "Getting started" (S16): the shortest tutorial, six one-line steps from choosing a program
// to saving a session, what stays in the browser, and the optional browser self-test. The
// header's "Getting started" opens and closes it (application/onboarding.ts); it is open on
// the very first visit.
import type { GettingStarted } from '../application/onboarding';
import type { SelfTest } from '../application/self-test';
import { PROGRAMS } from '../domain/programs';
import SelfTestPanel from './SelfTestPanel.vue';
import './onboarding.css';

defineProps<{ gettingStarted: GettingStarted; selfTest: SelfTest }>();
/** The programs that can be searched in this release ("BLASTN, BLASTP, TBLASTN or TBLASTX"). */
const labels = PROGRAMS.filter((program) => program.unavailable === undefined).map((program) => program.label);
const programs = labels.length <= 1 ? (labels[0] ?? '') : `${labels.slice(0, -1).join(', ')} or ${labels.at(-1)}`;
</script>

<template>
  <section id="getting-started" class="getting-started" aria-labelledby="getting-started-title" data-testid="getting-started">
    <div class="getting-started-head">
      <h2 id="getting-started-title">Getting started</h2>
      <button type="button" class="link" data-testid="getting-started-close" @click="gettingStarted.close()">Close</button>
    </div>
    <ol class="getting-started-steps">
      <li>Choose a program on the Search tab: {{ programs }}.</li>
      <li>Give a query and a subject: paste FASTA, open FASTA files, or choose “Try an example”.</li>
      <li>Choose Run LOSAT. The search runs in this browser; the queue shows how far it is.</li>
      <li>Read the results: the Graphic Summary, the Descriptions and the Alignments, on one page.</li>
      <li>Collect HSPs in Candidates and extract their sequences as FASTA.</li>
      <li>Save a session file to keep the runs, and open it later without searching again.</li>
    </ol>
    <p class="local-note">
      Your sequences, results and files stay in this browser: LOSAT Web sends none of them anywhere, and saves a file only
      when you choose to save it.
    </p>
    <SelfTestPanel :self-test="selfTest" />
  </section>
</template>
