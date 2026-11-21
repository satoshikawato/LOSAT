<script setup lang="ts">
// The words of a place of messages of a run (message-place-text.ts), with each option name written by OptionName, as
// the other option names in the results (line WQA's sweep, D7, 2026-11-20).
import type { PlacePart } from './message-place-text';
import OptionName from './OptionName.vue';

defineProps<{ parts: readonly PlacePart[] }>();
</script>

<template>
  <template v-for="(part, i) in parts" :key="i"><OptionName v-if="part.option" :name="part.text" /><template v-else>{{ part.text }}</template></template>
</template>
