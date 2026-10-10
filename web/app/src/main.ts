import { createApp as createVueApp } from 'vue';
import AppView from './ui/AppView.vue';
import { createApp } from './composition';
import './ui/styles.css';

const { coordinator, draft, results, candidates, attention, usesFakeEngine, session } = createApp();
createVueApp(AppView, { coordinator, draft, results, candidates, attention, usesFakeEngine, session }).mount('#app');
