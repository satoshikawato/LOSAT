<script setup lang="ts">
// The words of an argv or a command line, quoted as on a command line (`shellWords`) and joined
// by single spaces: a line may break between words and inside a long value, and inside an option's
// name only where the name is longer than the line, never after its leading hyphen (OptionName; S13
// screen review L2, line WQA's sweep 2 S2-1). A value that starts with '-' (`-penalty -3`) is a
// value, as `describe` says which options take one (WA-FM-07). The text is `shellWords`'s.
import { computed } from 'vue';
import { optionNames } from '../domain/argv';
import OptionName from './OptionName.vue';

const props = defineProps<{ words: readonly string[]; takesValue?: ((flag: string) => boolean) | undefined }>();
const parts = computed(() => optionNames(props.words, props.takesValue));
</script>

<template>
  <template v-for="(part, i) in parts" :key="i">{{ i === 0 ? '' : ' ' }}<OptionName v-if="part.name" :name="part.text" /><template v-else>{{ part.text }}</template></template>
</template>
