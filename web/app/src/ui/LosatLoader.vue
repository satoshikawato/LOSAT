<script setup lang="ts">
// The loading indicator (DW-43): the logo's A drawn in one stroke from the left foot over the apex,
// which branches left and right along the base, as gbdraw's loader draws its logo. White pens in a
// mask reveal the filled A; the motion is `.losat-loader` in styles.css. It sits before the words of
// a state that is in progress, which say what it is.
import { useId } from 'vue';
import { LOADER_PENS, LOADER_VIEWBOX, LOGO_A } from './logo';

const mask = `losat-loader-${useId()}`;
</script>

<template>
  <svg class="losat-loader" :viewBox="LOADER_VIEWBOX" aria-hidden="true" focusable="false" data-testid="losat-loader">
    <mask :id="mask" maskUnits="userSpaceOnUse" x="70" y="140" width="20" height="20">
      <path
        v-for="pen in LOADER_PENS"
        :key="pen.id"
        class="losat-loader-pen"
        :class="`losat-loader-pen-${pen.id}`"
        :pathLength="pen.length"
        :d="pen.d"
      />
    </mask>
    <g :mask="`url(#${mask})`">
      <path v-for="(d, i) in LOGO_A" :key="i" :d="d" />
    </g>
  </svg>
</template>
