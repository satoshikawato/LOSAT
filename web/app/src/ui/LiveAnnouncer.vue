<script setup lang="ts">
// The two live regions of announcer.ts, in the page from its start and not seen. Each telling
// replaces the region's text with a new node, so that the same text told again is read again.
import { told } from './announcer';
</script>

<template>
  <div class="visually-hidden">
    <p role="status" data-testid="announce-status"><span :key="told.status.value.key">{{ told.status.value.text }}</span></p>
    <p role="alert" data-testid="announce-alert"><span :key="told.alert.value.key">{{ told.alert.value.text }}</span></p>
  </div>
</template>
