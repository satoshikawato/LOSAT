<script setup lang="ts">
// An option's name (`-max_target_seqs`), as a command line or a sentence shows it. A line takes the name whole, or
// leaves it whole to the next line; only a name longer than its line breaks, at the line's end. Its lead (the hyphen and
// the letter after it) is never broken, so that no line ends with a lone hyphen (S13 screen review L2). The name used to
// be one piece that could not break, which made Run details' Arguments and the results' Options line wider than a 320 or
// 360 px window (-min_raw_gapped_score; line WQA's sweep 2, S2-1, 2026-11-20). The text is the name's, unchanged.
import { computed } from 'vue';

const props = defineProps<{ name: string }>();
/** The end of the lead: the letter after the first hyphen (a name that the command line quotes starts with its quote). */
const lead = computed(() => props.name.indexOf('-') + 2);
</script>

<template>
  <span class="option-name"><span class="option-lead">{{ name.slice(0, lead) }}</span>{{ name.slice(lead) }}</span>
</template>
