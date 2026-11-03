<script setup lang="ts">
// The "?" beside a control and the help it opens (stop-describing-ui: the control shows, the
// help waits behind "?"). As NCBI's helpbox, the help shows while the "?" is pointed at with a
// mouse or pen or reached with the keyboard, and a click keeps it open until the next click,
// Escape or a click elsewhere (a touch screen has no hover: a tap opens, a tap closes). The
// help floats under the control (styles.css .help-tip), shifted left when it would leave the
// viewport, so the page does not move while the pointer passes over the "?".
import { computed, nextTick, onUnmounted, ref, watch } from 'vue';

const props = defineProps<{ id: string; label: string; text: string }>();

const pinned = ref(false);
const pointed = ref(false);
const shown = computed(() => pinned.value || pointed.value);
const root = ref<HTMLElement | null>(null);
const box = ref<HTMLElement | null>(null);
/** Pixels the help is moved left so that it ends inside the viewport. */
const shift = ref(0);

const point = (event: PointerEvent, on: boolean) => {
  if (event.pointerType !== 'touch') pointed.value = on;
};
/** Keyboard focus shows the help; the focus a click or a tap gives the button does not. */
const focused = (event: FocusEvent) => {
  const target = event.target as Element;
  try {
    pointed.value = target.matches(':focus-visible');
  } catch {
    pointed.value = false;
  }
};
/** A click pins the help; the next click closes it, although the pointer is still over the "?". */
const toggle = () => {
  if (pinned.value) close();
  else pinned.value = true;
};
const close = () => {
  pinned.value = false;
  pointed.value = false;
};
const outside = (event: PointerEvent) => {
  if (root.value !== null && !root.value.contains(event.target as Node)) pinned.value = false;
};

watch(pinned, (on) => {
  if (on) document.addEventListener('pointerdown', outside, true);
  else document.removeEventListener('pointerdown', outside, true);
});
watch(shown, async (on) => {
  shift.value = 0;
  if (!on) return;
  await nextTick();
  const rect = box.value?.getBoundingClientRect();
  if (rect !== undefined) {
    const over = Math.ceil(rect.right - (window.innerWidth - 16));
    shift.value = Math.min(Math.max(0, over), Math.max(0, Math.floor(rect.left - 16)));
  }
});
onUnmounted(() => document.removeEventListener('pointerdown', outside, true));
</script>

<template>
  <span ref="root" class="help-tip">
    <button
      type="button"
      class="help-toggle"
      :aria-expanded="shown"
      :aria-controls="shown ? `${props.id}-help` : undefined"
      :aria-describedby="shown ? `${props.id}-help` : undefined"
      :aria-label="props.label"
      :data-testid="`${props.id}-help-toggle`"
      @click="toggle"
      @pointerenter="point($event, true)"
      @pointerleave="point($event, false)"
      @focus="focused"
      @blur="pointed = false"
      @keydown.escape="close"
    >
      ?
    </button>
    <p
      v-if="shown"
      :id="`${props.id}-help`"
      ref="box"
      class="help"
      role="tooltip"
      :data-testid="`${props.id}-help`"
      :data-pinned="pinned"
      :style="shift > 0 ? { transform: `translateX(-${shift}px)` } : undefined"
    >
      {{ props.text }}
    </p>
  </span>
</template>
