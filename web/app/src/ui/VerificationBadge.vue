<script setup lang="ts">
// The verification badge of a run (plan §6.1, domain/verification.ts). The compact form
// sits next to the run; the full form, in the run's details, says why.
import type { Badge } from '../domain/verification';

defineProps<{ badge: Badge; compact?: boolean }>();
defineEmits<{ details: [] }>();
</script>

<template>
  <button
    v-if="compact"
    type="button"
    class="verification-badge"
    :data-level="badge.level"
    data-testid="verification-badge"
    :title="badge.details.join(' ')"
    @click="$emit('details')"
  >
    {{ badge.label }}<template v-if="badge.exceptions.length > 0"> · approved exception</template>
  </button>
  <div v-else class="verification" data-testid="verification-details" :data-level="badge.level">
    <p><span class="verification-badge" :data-level="badge.level">{{ badge.label }}</span></p>
    <ul>
      <li v-for="(line, i) in badge.details" :key="i">{{ line }}</li>
    </ul>
    <p v-for="(line, i) in badge.exceptions" :key="`e${i}`" class="notice" data-testid="verification-exception">{{ line }}</p>
  </div>
</template>
