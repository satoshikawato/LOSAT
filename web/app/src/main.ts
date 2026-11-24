import { createApp as createVueApp } from 'vue';
import AppView from './ui/AppView.vue';
import { createApp } from './composition';
import './ui/styles.css';

const { coordinator, draft, toolDraft, examples, counts, results, attention, availability, usesFakeEngine } = createApp();
createVueApp(AppView, { coordinator, draft, toolDraft, examples, counts, results, attention, availability, usesFakeEngine }).mount('#app');
