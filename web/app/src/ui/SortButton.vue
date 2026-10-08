<script setup lang="ts" generic="K extends string">
// A column header that sorts its list by the column's engine value: the first click sorts
// in the column's natural direction, the next reverses it.
import type { SortSpec } from '../domain/result-index';

const props = withDefaults(
  defineProps<{ label: string; sortKey: K; sort: SortSpec<K>; scope: string; firstDescending?: boolean }>(),
  { firstDescending: false },
);
const emit = defineEmits<{ sort: [key: K, descending: boolean] }>();

function click(): void {
  if (props.sort.key === props.sortKey) emit('sort', props.sortKey, !props.sort.descending);
  else emit('sort', props.sortKey, props.firstDescending);
}
</script>

<template>
  <span role="columnheader" :aria-sort="sort.key === sortKey ? (sort.descending ? 'descending' : 'ascending') : 'none'">
    <button type="button" class="sort-button" :data-testid="`${scope}-sort-${sortKey}`" @click="click">
      {{ label }}<span v-if="sort.key === sortKey" aria-hidden="true">{{ sort.descending ? ' ▼' : ' ▲' }}</span>
    </button>
  </span>
</template>
