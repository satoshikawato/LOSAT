// Whether the screen is narrow (the same query as the phone rules of styles.css), kept up to date.
// The lists of fixed row height give a row two lines there and need its height in script.
import { onMounted, onUnmounted, ref, type Ref } from 'vue';

const NARROW = '(max-width: 800px)';

export function useNarrow(): Ref<boolean> {
  const narrow = ref(false);
  let query: MediaQueryList | undefined;
  const onChange = (event: MediaQueryListEvent) => (narrow.value = event.matches);
  onMounted(() => {
    query = window.matchMedia(NARROW);
    narrow.value = query.matches;
    query.addEventListener('change', onChange);
  });
  onUnmounted(() => query?.removeEventListener('change', onChange));
  return narrow;
}
