<script lang="ts">
/** Closes the help that the pointer opened last, which is shown while the pointer is not over it any more. */
let pointing: (() => void) | undefined;
</script>

<script setup lang="ts">
// The "?" beside a control and the help it opens (stop-describing-ui: the control shows, the
// help waits behind "?"). As NCBI's helpbox, the help shows while the "?" is pointed at with a
// mouse or pen or reached with the keyboard, and a click keeps it open until the next click,
// Escape or a click elsewhere (a touch screen has no hover: a tap opens, a tap closes). The
// help floats under the control (styles.css .help-tip), shifted left when it would leave the
// viewport, so the page does not move while the pointer passes over the "?". As WCAG 1.4.13
// asks (WA2-AX-07), the pointer can move from the "?" onto the help without closing it, and
// Escape closes the help wherever the focus is, which stays where it is.
import { computed, nextTick, onUnmounted, ref, watch } from 'vue';

const props = defineProps<{ id: string; label: string; text: string }>();

const pinned = ref(false);
/** Pointed at with a mouse or pen: the "?", the help, or the way between them. */
const pointed = ref(false);
/** Reached with the keyboard. */
const keyed = ref(false);
const shown = computed(() => pinned.value || pointed.value || keyed.value);
const root = ref<HTMLElement | null>(null);
const box = ref<HTMLElement | null>(null);
/** Pixels the help is moved left so that it ends inside the viewport. */
const shift = ref(0);

/**
 * Whether a point is in the box around the "?" and the help, the way from one to the other: the
 * help floats a few pixels under the control, not always under the "?".
 */
function around(x: number, y: number): boolean {
  const rects = [root.value, box.value].flatMap((element) => (element === null ? [] : [element.getBoundingClientRect()]));
  if (rects.length === 0) return false;
  return (
    x >= Math.min(...rects.map((r) => r.left)) &&
    x <= Math.max(...rects.map((r) => r.right)) &&
    y >= Math.min(...rects.map((r) => r.top)) &&
    y <= Math.max(...rects.map((r) => r.bottom))
  );
}
/** After the pointer left the "?" and the help on its way between them: closes once it leaves that way. */
const follow = (event: PointerEvent) => {
  if (around(event.clientX, event.clientY)) return;
  pointed.value = false;
  unfollow();
};
const unfollow = () => document.removeEventListener('pointermove', follow, true);
/** The pointer's part of the help closes (a pinned help, or one the keyboard shows, stays). */
const unpoint = () => {
  pointed.value = false;
  unfollow();
};
/**
 * Only one help is open by the pointer: a "?" inside the box of another help (the way from its
 * "?" to the help) closes that one as the pointer comes to it.
 */
const enter = (event: PointerEvent) => {
  if (event.pointerType === 'touch') return;
  if (pointing !== undefined && pointing !== unpoint) pointing();
  pointing = unpoint;
  unfollow();
  pointed.value = true;
};
const leave = (event: PointerEvent) => {
  if (event.pointerType === 'touch' || !pointed.value) return;
  if (around(event.clientX, event.clientY)) document.addEventListener('pointermove', follow, true);
  else pointed.value = false;
};
/** Keyboard focus shows the help; the focus a click or a tap gives the button does not. */
const focused = (event: FocusEvent) => {
  const target = event.target as Element;
  try {
    keyed.value = target.matches(':focus-visible');
  } catch {
    keyed.value = false;
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
  keyed.value = false;
  unfollow();
};
const outside = (event: PointerEvent) => {
  if (root.value !== null && !root.value.contains(event.target as Node)) pinned.value = false;
};
/** Escape closes the shown help, wherever the focus is; the focus stays (WCAG 1.4.13: dismissible). */
const escape = (event: KeyboardEvent) => {
  // The Escape that cancels an IME's composition is not for the page.
  if (event.isComposing) return;
  if (event.key === 'Escape') close();
};

watch(pinned, (on) => {
  if (on) document.addEventListener('pointerdown', outside, true);
  else document.removeEventListener('pointerdown', outside, true);
});
watch(shown, async (on) => {
  if (on) document.addEventListener('keydown', escape);
  else document.removeEventListener('keydown', escape);
  shift.value = 0;
  if (!on) return;
  await nextTick();
  const rect = box.value?.getBoundingClientRect();
  if (rect !== undefined) {
    const over = Math.ceil(rect.right - (window.innerWidth - 16));
    shift.value = Math.min(Math.max(0, over), Math.max(0, Math.floor(rect.left - 16)));
  }
});
onUnmounted(() => {
  if (pointing === unpoint) pointing = undefined;
  document.removeEventListener('pointerdown', outside, true);
  document.removeEventListener('keydown', escape);
  unfollow();
});
</script>

<template>
  <!-- The pointer over the "?" or over the help (a child, although it floats elsewhere) shows it. -->
  <span ref="root" class="help-tip" @pointerenter="enter" @pointerleave="leave">
    <button
      type="button"
      class="help-toggle"
      :aria-expanded="shown"
      :aria-controls="shown ? `${props.id}-help` : undefined"
      :aria-describedby="shown ? `${props.id}-help` : undefined"
      :aria-label="props.label"
      :data-testid="`${props.id}-help-toggle`"
      @click="toggle"
      @focus="focused"
      @blur="keyed = false"
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
