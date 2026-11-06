<script setup lang="ts" generic="T extends string | number">
// A row of tabs in ARIA's tabs pattern (WA2-AX-12): a tablist of tabs, the chosen one
// aria-selected, each naming the panel it shows (aria-controls). One Tab stop, the chosen tab
// (roving tabindex); Left and Right move to the tab before and after it (around the ends), Home
// and End to the first and the last. A tab reached so is chosen at once, as its panel shows at
// once; with `manual` (a panel that takes a while to show, as an output read from storage), the
// keys move only the focus, and Enter or Space chooses (ARIA APG). The look stays NCBI's tab
// strip (styles.css .tabs).
import { nextTick } from 'vue';

const props = defineProps<{
  /**
   * Each tab: its value, its words, its test ID, its ID (which its panel is labelled by), the ID of
   * the panel that it shows, and a class.
   */
  tabs: readonly {
    readonly value: T;
    readonly label: string;
    readonly testid: string;
    readonly id: string;
    readonly panel: string;
    readonly class?: string;
  }[];
  label: string;
  /** The arrow keys move the focus only; Enter or Space chooses the tab. */
  manual?: boolean;
}>();
const chosen = defineModel<T>({ required: true });
const buttons: HTMLButtonElement[] = [];
function keep(index: number, element: unknown): void {
  if (element instanceof HTMLButtonElement) buttons[index] = element;
}
/** Moves the focus to the chosen tab. */
defineExpose({ focus: () => focusChosen() });

function focusChosen(): void {
  const index = props.tabs.findIndex((tab) => tab.value === chosen.value);
  buttons[index]?.focus();
}

async function onKey(event: KeyboardEvent): Promise<void> {
  if (event.ctrlKey || event.metaKey || event.altKey) return;
  const count = props.tabs.length;
  const focused = buttons.findIndex((button, index) => index < count && button === document.activeElement);
  const at = focused >= 0 ? focused : props.tabs.findIndex((tab) => tab.value === chosen.value);
  const next: Record<string, number> = { ArrowRight: (at + 1) % count, ArrowLeft: (at - 1 + count) % count, Home: 0, End: count - 1 };
  const index = next[event.key];
  if (index === undefined || count === 0) return;
  event.preventDefault();
  if (props.manual) {
    buttons[index]?.focus();
    return;
  }
  chosen.value = props.tabs[index]!.value;
  await nextTick();
  focusChosen();
}
</script>

<template>
  <div role="tablist" :aria-label="label" class="tabs" @keydown="onKey">
    <button
      v-for="(tab, index) in tabs"
      :id="tab.id"
      :key="String(tab.value)"
      :ref="(element) => keep(index, element)"
      type="button"
      role="tab"
      :class="tab.class"
      :aria-selected="chosen === tab.value"
      :aria-controls="tab.panel"
      :tabindex="chosen === tab.value ? 0 : -1"
      :data-testid="tab.testid"
      @click="chosen = tab.value"
    >
      {{ tab.label }}
    </button>
  </div>
</template>
