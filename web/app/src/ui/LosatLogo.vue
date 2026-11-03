<script setup lang="ts">
// The LOSAT logo in the header (DW-43), filled in the text colour. The heading's name is the text
// beside it; the drawing is hidden from assistive technology.
import { LOGO_A, LOGO_LETTERS, LOGO_VIEWBOX } from './logo';
</script>

<template>
  <svg class="losat-logo" :viewBox="LOGO_VIEWBOX" aria-hidden="true" focusable="false" data-testid="losat-logo">
    <path v-for="(d, i) in LOGO_LETTERS" :key="i" :d="d" />
    <path :d="LOGO_A" />
  </svg>
</template>
