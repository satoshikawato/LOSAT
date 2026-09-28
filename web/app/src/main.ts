import { createApp as createVueApp } from 'vue';
import AppView from './ui/AppView.vue';
import { createApp } from './composition';
import './ui/styles.css';

const { coordinator, usesFakeEngine } = createApp();
createVueApp(AppView, { coordinator, usesFakeEngine }).mount('#app');
