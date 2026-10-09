<script setup lang="ts">
import type { Attention, AttentionState } from '../application/attention';

defineProps<{ attention: Attention; state: AttentionState }>();
</script>

<template>
  <section class="attention" data-testid="attention">
    <h2>While searching</h2>
    <label v-if="state.wakeLock !== 'unsupported'" class="flag-label">
      <input
        type="checkbox"
        :checked="state.keepAwake"
        data-testid="keep-awake"
        @change="attention.setKeepAwake(($event.target as HTMLInputElement).checked)"
      />
      Keep the screen on while a search runs
    </label>
    <p v-else class="muted" data-testid="keep-awake-unsupported">This browser cannot keep the screen on.</p>
    <p v-if="state.keepAwake && state.wakeLock === 'held'" class="muted" data-testid="wake-lock-held">
      The screen stays on until the search ends.
    </p>
    <p v-if="state.wakeLock === 'failed'" class="warning" data-testid="wake-lock-failed">
      The browser did not keep the screen on: {{ state.wakeLockError }}
    </p>
    <p class="muted" :class="{ warning: state.busy }" data-testid="attention-note">
      Searches run in this tab. Keep it open and in the foreground: browsers slow down or stop the work of hidden tabs,
      and closing or reloading the tab ends the searches and discards their results.
    </p>
  </section>
</template>
