<script setup lang="ts">
// The words of an argv or a command line, quoted as on a command line (`shellWords`) and joined
// by single spaces: a line may break between words and inside a long value, but not after an
// option's leading hyphen (S13 screen review L2). A value that starts with '-' (`-penalty -3`)
// is a value, as `describe` says which options take one (WA-FM-07). The text is `shellWords`'s.
import { computed } from 'vue';
import { optionNames } from '../domain/argv';

const props = defineProps<{ words: readonly string[]; takesValue?: ((flag: string) => boolean) | undefined }>();
const parts = computed(() => optionNames(props.words, props.takesValue));
</script>

<template>
  <template v-for="(part, i) in parts" :key="i">{{ i === 0 ? '' : ' ' }}<span v-if="part.name" class="option-name">{{ part.text }}</span><template v-else>{{ part.text }}</template></template>
</template>
