import { createApp as createVueApp } from 'vue';
import AppView from './ui/AppView.vue';
import { createApp } from './composition';
import './ui/styles.css';

const { coordinator, draft, results, candidates, exporter, attention, usesFakeEngine, session, runFiles, examples, selfTest, gettingStarted } = createApp();
createVueApp(AppView, { coordinator, draft, results, candidates, exporter, attention, usesFakeEngine, session, runFiles, examples, selfTest, gettingStarted }).mount('#app');
