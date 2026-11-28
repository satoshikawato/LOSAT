import { createApp as createVueApp } from 'vue';
import AppView from './ui/AppView.vue';
import { createApp } from './composition';
import './ui/styles.css';

const { coordinator, draft, toolDraft, maskingDraft, blastFormatterDraft, makeblastdbDraft, examples, counts, results, attention, availability, usesFakeEngine } =
  createApp();
createVueApp(AppView, {
  coordinator,
  draft,
  toolDraft,
  maskingDraft,
  blastFormatterDraft,
  makeblastdbDraft,
  examples,
  counts,
  results,
  attention,
  availability,
  usesFakeEngine,
}).mount('#app');
