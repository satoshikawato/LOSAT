import { createApp as createVueApp } from 'vue';
import AppView from './ui/AppView.vue';
import { createApp } from './composition';
import './ui/styles.css';

const { coordinator, draft, results, candidates, exporter, attention, usesFakeEngine, session, runFiles, examples, selfTest, gettingStarted, release, analytics } = createApp();
createVueApp(AppView, { coordinator, draft, results, candidates, exporter, attention, usesFakeEngine, session, runFiles, examples, selfTest, gettingStarted, release, analytics }).mount('#app');
// After the mount: the release manager waits until the page is idle before it does anything.
release.start();
