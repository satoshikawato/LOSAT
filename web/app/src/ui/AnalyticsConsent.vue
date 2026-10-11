<script setup lang="ts">
// The consent notice of usage statistics (docs/web/analytics.md): shown only when the site
// measures and no choice is stored, or when "Change" reopens it. It sits in the page's flow,
// so it never covers the search form, and it takes no focus.
import type { Analytics } from '../application/analytics';
import { useStore } from './useStore';

const props = defineProps<{ analytics: Analytics }>();
const state = useStore(props.analytics.state);
</script>

<template>
  <section
    v-if="state.available && state.asking"
    class="analytics-consent"
    aria-labelledby="analytics-consent-title"
    data-testid="analytics-consent"
  >
    <h2 id="analytics-consent-title">Usage statistics</h2>
    <p>
      May LOSAT Web count your visits? If you allow it, each time this page opens it sends one event through this site to
      Google Analytics: that LOSAT Web was opened, and its release ID<template v-if="analytics.release"> ({{ analytics.release }})</template>,
      with a random ID kept in this browser so that a returning visitor is counted once.
    </p>
    <p><strong>Your sequences, file names, results, notes and settings are never sent.</strong></p>
    <p class="analytics-actions">
      <button type="button" data-testid="analytics-allow" @click="analytics.allow()">Allow</button>
      <button type="button" data-testid="analytics-deny" @click="analytics.deny()">Don't allow</button>
      <span class="hint">You can change this at any time under “Usage statistics” at the bottom of the page.</span>
    </p>
  </section>
</template>
