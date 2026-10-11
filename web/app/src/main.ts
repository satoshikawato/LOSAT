import { createApp as createVueApp } from 'vue';
import AppView from './ui/AppView.vue';
import { createApp } from './composition';
import './ui/styles.css';

const { coordinator, draft, results, candidates, exporter, attention, usesFakeEngine, session, runFiles, release } = createApp();
createVueApp(AppView, { coordinator, draft, results, candidates, exporter, attention, usesFakeEngine, session, runFiles, release }).mount('#app');
// After the mount: the release manager waits until the page is idle before it does anything.
release.start();
