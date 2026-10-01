<script setup lang="ts">
import { ref } from 'vue';
import type { Coordinator } from '../application/coordinator';
import { useStore } from './useStore';
import SearchPanel from './SearchPanel.vue';
import QueuePanel from './QueuePanel.vue';
import ResultsPanel from './ResultsPanel.vue';
import StorageStatus from './StorageStatus.vue';

const props = defineProps<{ coordinator: Coordinator; usesFakeEngine: boolean }>();
const state = useStore(props.coordinator.state);
const tab = ref<'search' | 'results'>('search');
</script>

<template>
  <header class="app-header">
    <h1>LOSAT Web</h1>
    <span class="local-note">Local processing: your sequences stay in this browser.</span>
  </header>
  <p v-if="usesFakeEngine" class="banner" data-testid="fake-engine-banner">
    Development build: a fake engine is in use. Outputs are not search results.
  </p>
  <nav class="tabs" aria-label="Main">
    <button :aria-pressed="tab === 'search'" data-testid="tab-search" @click="tab = 'search'">Search</button>
    <button :aria-pressed="tab === 'results'" data-testid="tab-results" @click="tab = 'results'">Results</button>
  </nav>
  <main class="layout">
    <section class="primary">
      <SearchPanel v-if="tab === 'search'" :coordinator="coordinator" />
      <ResultsPanel v-else :coordinator="coordinator" :runs="state.runs" />
    </section>
    <aside class="secondary">
      <QueuePanel :coordinator="coordinator" :runs="state.runs" />
      <StorageStatus :storage="state.storage" />
    </aside>
  </main>
</template>
