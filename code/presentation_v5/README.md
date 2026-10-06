# Submission artwork terminology patch

The three included baseline SVG files preserve the existing frozen plotted geometry. The helper changes only the reviewed prespecification labels (and the Figure 2 TMD-signal-reveal wording). It checks numerical text and all non-text SVG elements for equality. It does not read statistical tables, fit models or recalculate plotted values. The original scientific plotting scripts remain unchanged.

From the code-package root, using a Python environment with lxml:

```sh
python presentation_v5/update_frozen_svg_labels.py --input presentation_v5/frozen_svg_inputs/Figure_2.svg --figure Figure_2 --output Figure_2_v5.svg
```

The other accepted figure identifiers are Figure_3 and Supplementary_Figure_2, with correspondingly named frozen inputs. Output paths must not already exist. These administrative commands were tested against the exact v5 SVG output hashes. Baseline labels in frozen_svg_inputs are historical source inputs, not the current journal artwork.

The submission exports were rendered with rsvg-convert, an opaque white background and 400-dpi PNG output. PDF canvas sizes were retained from v4. Figures 1 and 4 and the other four supplementary figures were copied unchanged from v4. The new helper's lxml dependency belongs to document production, not to the historical scientific software stack.
