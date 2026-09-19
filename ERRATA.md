# Errata

## Version 3 (doi:10.5281/zenodo.22845658)

An independent recomputation of every structure-prediction number from the raw Boltz-2 outputs (Codex, OpenAI, 20 September 2026) reproduced all ipTM values, ranges, medians and interface counts, and found these errors in version 3. A corrected version 4 is being prepared.

- **SMC2 contact residues (Results, N-terminal paragraph).** The native residue numbers of the SMC2 neck segments were off by the construct's 12-residue linker. Correct, counting SMC2 residues touched by confident kleisin residues in each protein's own numbering: human 169–210 and 986–1032 (stated 168–232 and 992–1048); *H. exemplaris* b 168–212 and 994–1040; *R. varieornatus* b 171–211 and 1000–1043 (stated "about 168–250 and 990–1055"). `caph2_structure_analysis.py --n2` now writes these in column `smc_spans_native`.
- **"Copy b closes the ring at both ends"** holds only for *H. exemplaris* b and *R. varieornatus* b, and for separately predicted interfaces, not a whole ring.
- **First N-terminal construct:** *H. exemplaris* b (0.44) did score above its wrong-partner and scrambled controls (0.29, 0.26); the others did not.
- **C-terminal interfaces:** *Paramacrobiotus* a has a higher local interface error (9.3 Å) than the other C-terminal predictions (2.6–3.8 Å).
- **Interface counts** (49, 39, 39) are for the selected sample; the other samples give 40–44, 37–38 and, for human, 38–42.
- **Controls:** scrambled regions also lose their sequence alignments; wrong-partner controls exist for human and *H. exemplaris* only.

Code fixes in this repository:
- the best sample's model number is now read from its file name, not its list position (no effect on the published numbers: every best sample was model 0);
- ipTM ranges are kept at full precision, with medians and sample counts;
- the comment on the CAP-H control length is corrected (200 residues).
