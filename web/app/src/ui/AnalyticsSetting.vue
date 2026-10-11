<script setup lang="ts">
// The footer line of usage statistics (docs/web/analytics.md): the choice in effect, and
// "Change", which reopens the consent notice. Shown only when the site measures.
import type { Analytics } from '../application/analytics';
import { useStore } from './useStore';

const props = defineProps<{ analytics: Analytics }>();
const state = useStore(props.analytics.state);
</script>

<template>
  <footer v-if="state.available" class="app-footer">
    <span data-testid="analytics-setting" :data-on="String(state.on)">Usage statistics: {{ state.on ? 'on' : 'off' }}</span>
    <span v-if="state.unsaved" class="muted" data-testid="analytics-unsaved">
      (this browser could not keep your choice, so nothing is sent and you will be asked again on your next visit)
    </span>
    <button v-if="!state.asking" type="button" data-testid="analytics-change" @click="analytics.change()">Change</button>
  </footer>
</template>
