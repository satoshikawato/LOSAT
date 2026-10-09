<script setup lang="ts">
// The words of an argv or a command line, joined by single spaces as given: a line may break
// between words and inside a long value, but not after an option's leading hyphen (S13 screen
// review L2). The text itself is unchanged.
import { computed } from 'vue';

const props = defineProps<{ text: string }>();
const words = computed(() => props.text.split(' '));
</script>

<template>
  <template v-for="(word, i) in words" :key="i">{{ i === 0 ? '' : ' ' }}<span v-if="word.startsWith('-')" class="option-name">{{ word }}</span><template v-else>{{ word }}</template></template>
</template>
