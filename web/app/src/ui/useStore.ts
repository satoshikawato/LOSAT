// Bridges an application Store to Vue reactivity. The application layer stays free of Vue.
import { onUnmounted, shallowRef, type ShallowRef } from 'vue';
import type { Store } from '../application/store';

export function useStore<T>(store: Store<T>): Readonly<ShallowRef<T>> {
  const state = shallowRef(store.get());
  const unsubscribe = store.subscribe((value) => {
    state.value = value;
  });
  onUnmounted(unsubscribe);
  return state;
}
