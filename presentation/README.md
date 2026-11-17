# QUBO Machine for Protein Sampling

Open `index.html` in Chrome or Edge. The presentation has 14 separate HTML slides with a planned speaking time of 20 minutes, followed by questions. Fonts, artwork, and charts are local, so presenting does not require internet access.

## Present

- Press Right Arrow, Page Down, or Space for the next slide.
- Press Left Arrow or Page Up for the previous slide.
- Press Home or End to jump to the opening or closing slide.
- Press F11 for browser fullscreen, which persists across page navigation.
- Press N to open `notes.html` at the current slide in a separate tab or window.
- Swipe horizontally on a touch device to navigate.

The F key also requests page fullscreen. Browsers may exit page fullscreen when navigating between HTML files, so F11 is recommended for a talk. No navigation controls, slide numbers, frames, or progress indicators are shown on the slides.

## Presenter scripts

- `script-neil.html` contains Neil Huang’s complete read-aloud script for slides 1–7, from 00:00 to 10:00.
- `script-daniel.html` contains Daniel Sun’s complete read-aloud script for slides 8–14, from 10:00 to 20:00.
- `notes.html` combines both scripts and opens at the current slide when you press N.

The scripts include the handoff between presenters and follow the presentation clock. Rehearse at a conversational pace, leaving space to point to diagrams. Source references are tucked into expandable sections beneath the spoken text.

## Content and evidence

The opening page retains the original title, authors, and ASCII protein artwork from the supplied single-slide PowerPoint, with a much larger, angled, edge-cropped composition. Subsequent slides use a connected read-aloud narrative to explain the model, encoding, corrected sampler, software implementation, measured results, and next experiments. Speaker notes include timing and source paths for every slide.

Five explanatory diagrams show the two sampling routes, interacting residue choices, one-hot encoding, the corrected proposal loop, and repeated states in a sampling trace. These are explicitly illustrative models, not experimental results or predicted protein structures. The original results charts remain grounded in the saved measurements. On phones, wide conceptual diagrams scroll horizontally within the slide so their labels stay readable.

Charts show the saved three-repeat small-case timing comparison, the reported residue-frequency uncertainty, and the reference replay storage reduction. Timing bars use medians calculated from `output/week07/method_repeats_verified/comparison.json`. The uncertainty and storage values come from `docs/week08/README.md`. The deck distinguishes software simulation from physical hardware execution and exploratory candidates from biological validation.

## Visual reference

The design follows the supplied Claude Code documentation screenshot: near-black background, muted body text, warm accent, serif headings, and unboxed content. The fonts are the actual Anthropic Sans Text and Anthropic Serif Display files referenced by the documentation site, downloaded from `assets.claude.ai/Fonts/`. Their original licensing remains applicable.

Body text uses 16px type with a 26.4px line height. Headings use 40px type. List indentation and 6px bullet geometry follow the documentation CSS. Large displays scale the slide uniformly; narrow screens reflow it for reading. The screenshot is the visual authority where the live website has changed.

## Edit and rebuild

Edit slide copy in `build.py`, read-aloud scripts in `speaker_scripts.py`, conceptual diagrams in `diagrams.py`, styling in `styles.css`, and navigation in `slides.js`. Rebuild from the project root with Python's standard library:

```sh
python presentation/build.py
```

To re-extract the original ASCII artwork as SVG, supply the PowerPoint path:

```sh
python presentation/build.py --source-pptx "C:/path/to/QUBO Machine for Protein Sampling.pptx"
```

The generated HTML files can also be edited directly, but rebuilding replaces those edits. Keep the complete `presentation` directory together when moving or sharing the deck. The original PowerPoint is unchanged.
