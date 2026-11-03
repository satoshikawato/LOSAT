import { createApp as createVueApp } from 'vue';
import AppView from './ui/AppView.vue';
import { createApp } from './composition';
import './ui/styles.css';

const { coordinator, draft, examples, results, attention, usesFakeEngine } = createApp();
createVueApp(AppView, { coordinator, draft, examples, results, attention, usesFakeEngine }).mount('#app');
