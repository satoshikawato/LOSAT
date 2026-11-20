<script setup lang="ts">
// The words of a place of messages of a run (message-place-text.ts), with each option name unbroken (`option-name`), as
// the other option names in the results (line WQA's sweep, D7, 2026-11-20).
import type { PlacePart } from './message-place-text';

defineProps<{ parts: readonly PlacePart[] }>();
</script>

<template>
  <template v-for="(part, i) in parts" :key="i"><span v-if="part.option" class="option-name">{{ part.text }}</span><template v-else>{{ part.text }}</template></template>
</template>
