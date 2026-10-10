import { createApp as createVueApp } from 'vue';
import AppView from './ui/AppView.vue';
import { createApp } from './composition';
import './ui/styles.css';

const { coordinator, draft, results, candidates, attention, usesFakeEngine, runFiles } = createApp();
createVueApp(AppView, { coordinator, draft, results, candidates, attention, usesFakeEngine, runFiles }).mount('#app');
