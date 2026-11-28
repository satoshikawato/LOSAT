<script setup lang="ts">
// Input and Counts blocks of the WindowMasker tab (WW3b DESIGN §2.1):
// - Make counts & Mask: -infmt (fasta, blastdb, seqids), FASTA paste/files/dropzone with
//   Combined (combined_input.fa) or Separate runs, or an editable -in name with Open files /
//   Open folder (added into the run directory) for blastdb, seqids or -fa_list (ToolFastaInputs,
//   which the Masking tab shows too).
// - Convert counts: counts file picker (session-made counts or Choose file…).
// - Mask: Counts (-ustat) picker (session-made counts or Choose file…).
import { computed, ref } from 'vue';
import type { SessionCountsEntry } from '../application/coordinator';
import type { ToolDraft, ToolDraftState } from '../application/tool-draft';
import { stageBlockParameter, toToolParameterField } from '../domain/tools';
import type { SavedCountsEntry } from '../ports/counts';
import type { FormModel } from './form-model';
import OptionName from './OptionName.vue';
import ParameterField from './ParameterField.vue';
import ToolCountsPicker from './ToolCountsPicker.vue';
import ToolFastaInputs from './ToolFastaInputs.vue';

const props = defineProps<{
  draft: ToolDraft;
  state: ToolDraftState;
  sessionCounts: readonly SessionCountsEntry[];
  savedCounts: readonly SavedCountsEntry[];
  model: FormModel;
}>();

const dragging = ref(false);

const infmtField = computed(() => {
  const param = stageBlockParameter(props.state.description, props.state.stage, '-infmt');
  return param === undefined ? undefined : toToolParameterField(param);
});

function onDrop(event: DragEvent): void {
  dragging.value = false;
  const files = event.dataTransfer?.files;
  if (files !== undefined && files.length > 0) props.draft.dropInputFiles([...files]);
}
</script>

<template>
  <fieldset
    class="search-block input-panel"
    :class="{ dragging }"
    data-testid="wm-input-block"
    @dragover.prevent="dragging = true"
    @dragleave="dragging = false"
    @drop.prevent="onDrop"
  >
    <legend>Input</legend>

    <template v-if="state.stage === 'convert'">
      <div class="param-row">
        <label class="param-label" for="wm-convert-select">
          Counts file <OptionName class="flag" name="-in" />
        </label>
        <ToolCountsPicker
          prefix="convert"
          :source="state.convertSource"
          :session-counts="sessionCounts"
          :saved-counts="savedCounts"
          file-aria-label="Choose counts file to convert"
          @select-session="(e) => draft.setConvertSession(e.runId, e.runNumber, e.key, e.name, e.size)"
          @select-saved="(e) => draft.setConvertSaved(e.taxid, e.size)"
          @choose-file="(file) => draft.setConvertFile(file)"
          @clear="draft.clearConvert()"
        />
      </div>
    </template>

    <template v-else>
      <ParameterField v-if="infmtField" :model="model" :field="infmtField" />

      <ToolFastaInputs :draft="draft" :state="state" prefix="wm" />
    </template>
  </fieldset>

  <fieldset v-if="state.stage === 'mask'" class="search-block" data-testid="wm-counts-block">
    <legend>Counts</legend>
    <div class="param-row">
      <label class="param-label" for="wm-ustat-select">
        Counts file <OptionName class="flag" name="-ustat" />
      </label>
      <ToolCountsPicker
        prefix="ustat"
        :source="state.ustatSource"
        :session-counts="sessionCounts"
        :saved-counts="savedCounts"
        file-aria-label="Choose counts file for -ustat"
        @select-session="(e) => draft.setUstatSession(e.runId, e.runNumber, e.key, e.name, e.size)"
        @select-saved="(e) => draft.setUstatSaved(e.taxid, e.size)"
        @choose-file="(file) => draft.setUstatFile(file)"
        @clear="draft.clearUstat()"
      />
    </div>
  </fieldset>
</template>
