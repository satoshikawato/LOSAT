<script setup lang="ts">
// The loading indicator (DW-43): the logo's A drawn in one stroke from the left, as gbdraw's loader
// draws its logo. A white pen in a mask reveals the filled A; the motion is `.losat-loader-pen` in
// styles.css. It sits before the words of a state that is in progress, which say what it is.
import { useId } from 'vue';
import { LOADER_PEN, LOADER_VIEWBOX, LOGO_A } from './logo';

const mask = `losat-loader-${useId()}`;
</script>

<template>
  <svg class="losat-loader" :viewBox="LOADER_VIEWBOX" aria-hidden="true" focusable="false" data-testid="losat-loader">
    <mask :id="mask" maskUnits="userSpaceOnUse" x="70" y="140" width="20" height="20">
      <path class="losat-loader-pen" pathLength="100" :d="LOADER_PEN" />
    </mask>
    <path :d="LOGO_A" :mask="`url(#${mask})`" />
  </svg>
</template>
