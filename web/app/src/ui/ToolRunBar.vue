<script setup lang="ts">
// Primary Run button + live command preview row for WindowMasker (shown above and below options).
import CommandText from './CommandText.vue';
import LosatLoader from './LosatLoader.vue';

defineProps<{
  submitting: boolean;
  buttonLabel: string;
  buttonTestId: string;
  commandWords: readonly string[];
  takesValue: (flag: string) => boolean;
  commandTestId?: string;
}>();

const emit = defineEmits<{ run: [] }>();
</script>

<template>
  <div class="blast-row">
    <button
      type="button"
      class="primary-action"
      :data-testid="buttonTestId"
      :aria-disabled="submitting || undefined"
      @click="emit('run')"
    >
      <LosatLoader v-if="submitting" />{{ buttonLabel }}
    </button>
    <p class="wm-command-line" :data-testid="commandTestId">
      <code><CommandText :words="commandWords" :takes-value="takesValue" /></code>
    </p>
  </div>
</template>
