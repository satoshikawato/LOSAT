- LOSAT Web: grouped local BLAST database (`-db`) controls under `Database masking and filters`
  with an exclusive database masking selector (`None`, `Soft mask (-db_soft_mask)`, `Hard mask
  (-db_hard_mask)`) that preserves the algorithm ID or name when switching between soft and hard
  masking, validates missing algorithm values before running, and updates placeholder guidance for
  nucleotide (`11 or dust`) and protein (`21 or seg`) databases.
- LOSAT Web: wired BLASTN WindowMasker (`-window_masker_taxid`, `-window_masker_db`) interactions
  with `-db` and `-no_taxid_expansion` so enabling WindowMasker while `-no_taxid_expansion` is
  checked defaults to a counts file (`-window_masker_db`), disables conflicting checkboxes and
  select options in place, auto-opens `Database masking and filters` on refusal or conflict, offers
  one-click resolution (`Use Counts file` / `Clear -no_taxid_expansion`) when both flags become
  active across FASTA/database switches, and keeps `?` help toggles inline with file-filter buttons
  at narrow viewports (320-390 px).
