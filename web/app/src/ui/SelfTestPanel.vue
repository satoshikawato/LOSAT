<script setup lang="ts">
// "Test this browser (optional)" in the Getting started panel (S16): the user may start the
// self-test (application/self-test.ts), which runs the examples with this browser's engine and
// compares the SHA-256 of each output with the release's. Each case says pass, fail (the
// expected and the actual SHA-256, their first characters) or error, with how the engine ran
// it; the browser's facts that the engine depends on come first.
import { computed } from 'vue';
import type { SelfTest, SelfTestCaseResult } from '../application/self-test';
import { digestPrefix } from '../domain/examples';
import { programById } from '../domain/programs';
import { formatCount } from './format';
import { useStore } from './useStore';

const props = defineProps<{ selfTest: SelfTest }>();
const state = useStore(props.selfTest.state);
const unavailable = props.selfTest.unavailable();
const plan = props.selfTest.plan();
const programs = plan.programs.map((program) => programById(program).label);
const threadsText = plan.threads.map((n) => `${n} ${n === 1 ? 'thread' : 'threads'}`).join(' and ');
const running = computed(() => state.value.status === 'running');

function joinWords(words: readonly string[]): string {
  return words.length <= 1 ? (words[0] ?? '') : `${words.slice(0, -1).join(', ')} and ${words.at(-1)}`;
}

function caseName(result: SelfTestCaseResult): string {
  return `${programById(result.program).label} at ${result.threads} ${result.threads === 1 ? 'thread' : 'threads'}`;
}

function engineText(result: SelfTestCaseResult): string | undefined {
  const runtime = result.runtime;
  if (runtime === undefined) return undefined;
  const module = runtime.path === 'threaded' ? 'threaded engine' : 'serial engine';
  const ran = `${module}, ${runtime.threads} ${runtime.threads === 1 ? 'thread' : 'threads'}`;
  const time = result.ms === undefined ? '' : `, ${formatCount(Math.round(result.ms))} ms`;
  return `${ran}${time}${runtime.fallbackReason === undefined ? '' : ` (serial because ${runtime.fallbackReason})`}`;
}
</script>

<template>
  <details class="self-test" data-testid="self-test">
    <summary>Test this browser (optional)</summary>
    <p class="hint">
      Runs the examples of {{ joinWords(programs) }} in this browser at {{ threadsText }} and compares the SHA-256 of
      their outfmt 0, 6 and 7 with the values that came with this release. Nothing is sent, and the test runs do not appear
      in the queue, the results, Candidates or session files.
    </p>
    <p v-if="unavailable" class="notice" data-testid="self-test-unavailable">{{ unavailable }}</p>
    <p v-for="skip in state.skipped" :key="skip.program" class="hint">
      {{ programById(skip.program).label }} is not tested: {{ skip.reason }}
    </p>
    <button type="button" data-testid="self-test-run" :disabled="unavailable !== undefined || running" @click="selfTest.run()">
      {{ running ? 'Testing…' : 'Run the self-test' }}
    </button>
    <p v-if="state.refusal && !running" class="error" role="alert" data-testid="self-test-refusal">{{ state.refusal }}</p>
    <div v-if="state.cases.length > 0" class="self-test-result" data-testid="self-test-result" :data-status="state.status">
      <p v-if="state.facts" class="hint" data-testid="self-test-facts">
        Cross-origin isolated: {{ state.facts.crossOriginIsolated ? 'yes' : 'no (searches run on one thread)' }}.
        Logical processors: {{ state.facts.logicalProcessors ?? 'not reported' }}.
        <template v-if="state.facts.storage">
          Outputs kept in: {{ state.facts.storage.backend === 'opfs' ? "the browser's private storage (OPFS)" : 'memory' }}<template
            v-if="state.facts.storage.fallbackReason"
          > ({{ state.facts.storage.fallbackReason }})</template>.
        </template>
      </p>
      <ul class="self-test-cases">
        <li
          v-for="result in state.cases"
          :key="`${result.program}-${result.threads}`"
          :data-testid="`self-test-case-${result.program}-${result.threads}`"
          :data-outcome="result.outcome?.kind ?? result.state"
        >
          <span class="self-test-case">{{ caseName(result) }}:</span>
          <template v-if="result.outcome === undefined">
            <span class="hint">{{ result.state === 'running' ? 'running…' : 'waiting' }}</span>
          </template>
          <template v-else-if="result.outcome.kind === 'pass'">
            <span class="ok">pass</span>
          </template>
          <template v-else-if="result.outcome.kind === 'fail'">
            <span class="error">fail</span>
            <span v-for="mismatch in result.outcome.mismatches" :key="mismatch.format" class="self-test-mismatch">
              outfmt {{ mismatch.format }}: expected {{ digestPrefix(mismatch.expected) }}…, got {{ digestPrefix(mismatch.got) }}…
            </span>
          </template>
          <template v-else>
            <span class="error">error</span> <span class="self-test-message">{{ result.outcome.message }}</span>
          </template>
          <span v-if="engineText(result)" class="hint self-test-engine">{{ engineText(result) }}</span>
        </li>
      </ul>
      <p v-if="state.summary" role="status" data-testid="self-test-summary">{{ state.summary }}</p>
    </div>
  </details>
</template>
