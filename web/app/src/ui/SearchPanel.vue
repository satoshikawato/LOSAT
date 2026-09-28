<script setup lang="ts">
import { ref } from 'vue';
import type { Coordinator } from '../application/coordinator';
import { PROGRAMS, type ProgramId } from '../domain/programs';

const props = defineProps<{ coordinator: Coordinator }>();
const program = ref<ProgramId>('blastn');
const queryText = ref('');
const subjectText = ref('');
const message = ref('');
const busy = ref(false);

async function addToQueue(): Promise<void> {
  busy.value = true;
  message.value = '';
  try {
    const result = await props.coordinator.enqueue({
      program: program.value,
      query: { text: queryText.value },
      subject: { text: subjectText.value },
      parameters: [],
      requestedThreads: 'auto',
    });
    message.value = result.ok ? 'Added to the queue.' : result.message;
  } finally {
    busy.value = false;
  }
}
</script>

<template>
  <h2>Search</h2>
  <fieldset class="programs">
    <legend>Program</legend>
    <label v-for="p in PROGRAMS" :key="p.id">
      <input v-model="program" type="radio" name="program" :value="p.id" :data-testid="`program-${p.id}`" />
      {{ p.label }}
    </label>
  </fieldset>
  <div class="inputs">
    <label>
      Query
      <textarea v-model="queryText" rows="6" spellcheck="false" data-testid="query-input" />
    </label>
    <label>
      Subject
      <textarea v-model="subjectText" rows="6" spellcheck="false" data-testid="subject-input" />
    </label>
  </div>
  <button :disabled="busy" data-testid="add-to-queue" @click="addToQueue">Add to queue</button>
  <p v-if="message" role="status" data-testid="search-message">{{ message }}</p>
</template>
