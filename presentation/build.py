from pathlib import Path
from html import escape
import json
import statistics
import xml.etree.ElementTree as ET
import zipfile
import argparse
import base64
from diagrams import build_diagrams
from speaker_scripts import SCRIPTS

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
ASSETS = ROOT / 'assets'
ACCENT = '#d4a27f'
WHITE = '#faf9f5'
GRAY = '#a3a3a3'


def svg_file(name, width, height, content):
    font = base64.b64encode((ASSETS / 'fonts/AnthropicSans-Regular.otf').read_bytes()).decode('ascii')
    style = f'<style>@font-face{{font-family:Anthropic;src:url(data:font/otf;base64,{font}) format("opentype")}}</style>'
    markup = f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">{style}<g font-family="Anthropic, sans-serif" font-size="16" fill="{GRAY}">{content}</g></svg>'
    (ASSETS / name).write_text(markup, encoding='utf-8')


def label(horizontal, vertical, text, color=GRAY, size=16):
    return f'<text x="{horizontal}" y="{vertical}" fill="{color}" font-size="{size}">{escape(text)}</text>'


def extract_art(source):
    namespaces = {'a': 'http://schemas.openxmlformats.org/drawingml/2006/main', 'p': 'http://schemas.openxmlformats.org/presentationml/2006/main'}
    with zipfile.ZipFile(source) as archive:
        root = ET.fromstring(archive.read('ppt/slides/slide1.xml'))
    shapes = root.findall('.//p:sp', namespaces)
    art = max(shapes, key=lambda shape: sum(len(text.text or '') for text in shape.findall('.//a:t', namespaces)))
    lines = [''.join(text.text or '' for text in paragraph.findall('.//a:t', namespaces)) for paragraph in art.findall('.//a:p', namespaces)]
    if len(lines) == 1:
        lines = lines[0].splitlines()
    lines = [line.rstrip() for line in lines if line.strip()]
    left = min(len(line) - len(line.lstrip()) for line in lines)
    lines = [line[left:] for line in lines]
    width = max(map(len, lines)) * 6.02
    content = ''.join(f'<text x="0" y="{(index + 1) * 10}" xml:space="preserve">{escape(line)}</text>' for index, line in enumerate(lines))
    svg_file('protein-ascii.svg', width, len(lines) * 10 + 5, f'<g font-family="Consolas, monospace" font-size="10" fill="{ACCENT}">{content}</g>')


def charts():
    comparison = json.loads((PROJECT / 'output/week07/method_repeats_verified/comparison.json').read_text())
    content = label(0, 20, 'Elapsed time to the first qualified checkpoint, in seconds.')
    for index, method in enumerate(['direct', 'encoded']):
        values = [run['elapsed_seconds'] for run in comparison['runs'] if run['method'] == method]
        vertical = 75 + index * 92
        median = statistics.median(values)
        color = WHITE if method == 'direct' else ACCENT
        content += label(0, vertical, 'Direct JAX' if method == 'direct' else 'Encoded THRML', WHITE)
        content += f'<rect x="185" y="{vertical-20}" width="{median*19}" height="27" fill="{color}" opacity=".85"/>'
        content += label(200 + median * 19, vertical, f'{median:.3f} s', color)
        content += label(185, vertical + 34, f'Three runs ranged from {min(values):.3f} to {max(values):.3f} seconds.', GRAY, 14)
    svg_file('timing.svg', 920, 220, content)
    content = label(0, 20, 'Frequency uncertainty remains above the declared limit.')
    for index, (name, value, color) in enumerate([('Required limit', .0025, WHITE), ('Position 65 · F', .00683, ACCENT)]):
        vertical = 78 + index * 84
        content += label(0, vertical, name, WHITE)
        content += f'<rect x="185" y="{vertical-20}" width="{value*80000}" height="27" fill="{color}" opacity=".85"/>'
        content += label(200 + value * 80000, vertical, f'{value:.5f}', color)
    svg_file('uncertainty.svg', 920, 215, content)
    content = label(0, 20, 'Saved artifacts for the identical reference replay, in decimal megabytes.')
    for index, (name, value, color) in enumerate([('Before', 94465955, WHITE), ('After', 12638863, ACCENT)]):
        vertical = 78 + index * 84
        content += label(0, vertical, name, WHITE)
        content += f'<rect x="120" y="{vertical-20}" width="{value/160000}" height="27" fill="{color}" opacity=".85"/>'
        content += label(135 + value / 160000, vertical, f'{value/1000000:.2f} MB', color)
    svg_file('storage.svg', 920, 215, content)


def paragraphs(*sentences):
    return ''.join(f'<p>{sentence}</p>' for sentence in sentences)


def bullets(*sentences):
    return '<ul>' + ''.join(f'<li>{sentence}</li>' for sentence in sentences) + '</ul>'


def figure(name, description, caption=''):
    caption_html = f'<figcaption>{caption}</figcaption>' if caption else ''
    return f'<figure><img src="assets/{name}" alt="{escape(description)}">{caption_html}</figure>'


slides = [
    ('QUBO Machine for Protein Sampling', 45, '', '', 'Introduce Neil Huang and Daniel Sun. Frame the talk as an investigation of protein sequence sampling with a binary energy model. The opening artwork is extracted from the supplied PowerPoint, not a new structure prediction.', 'Supplied PowerPoint, slide 1.'),
    ('We want better samples from the same model.', 90, 'The research question',
     paragraphs('We want to explore many plausible protein sequences, rather than find just one low-energy sequence.') + figure('model-routes.svg', 'One protein model feeds direct sampling and corrected binary sampling, both targeting the same sequence probabilities.', '') + paragraphs('Our question is whether the binary route can give us reliable samples in less time.', 'To answer it fairly, we keep the protein model fixed and check accuracy before comparing speed.'),
     'Explain the difference between finding one low-energy sequence and sampling a distribution of sequences. Sampling is useful for exploring alternative candidates and estimating residue frequencies. A faster sampler would improve access to this fixed model, not automatically improve the model itself. Introduce the talk in three stages: target, encoded method, and measured evidence.', 'docs/implementation_roadmap.md; docs/biology/technical_reference.md'),
    ('We start with a model of chorismate mutase.', 90, 'The biological target',
     paragraphs('Our model has <strong>96 aligned positions</strong>, with an allowed set of amino acids at each position.') + figure('potts-model.svg', 'A five-position amino-acid example: alanine, cysteine, glycine, threonine and alanine. A curve connects cysteine and threonine to represent a pairwise interaction.', 'Here is a five-position example of how residue choices interact.') + paragraphs('Choosing cysteine at one position can change which amino acid the model prefers at another.', 'The model combines these interactions to score the whole sequence.'),
     'Chorismate mutase is the protein family used for this computational benchmark. Explain an alignment column and why an alignment can include gaps. The coordinate anchor maps to Q5SLA5 residues 2–90, which is distinct from the number of aligned model positions. The anchor is outside model support at two positions and cannot initialize a chain. Restrictions specify the distribution being sampled rather than being repaired after sampling.', 'docs/week04/benchmark_contract.md; data/benchmark/position_map.csv; data/benchmark/restriction_cases.json'),
    ('The energy defines a distribution over sequences.', 90, 'The sampling target',
     paragraphs('We describe those preferences and interactions with a Potts model, which assigns each sequence an energy.', 'We then turn that energy into a probability using the relationship below.') + '<p class="equation">π(x) ∝ exp(−β E(x)).</p>' + paragraphs('In words, lower-energy sequences are more likely, and the inverse temperature β controls that preference.', 'We still want a range of sequences, so our goal is to sample this distribution rather than minimize energy.', '<strong>This distribution is the common target for both sampling methods.</strong>'),
     'Introduce x as a complete supported sequence and beta as inverse temperature. The source energy is the negative sum of fields and pair couplings. Emphasize that sampling preserves variation: a sample is not necessarily the minimum-energy sequence. The direct sampler provides a baseline in the original categorical state space. Fixed positions are removed computationally while their energy contributions are retained.', 'docs/encoded_method.md, Target and support'),
    ('We turn each residue choice into bits.', 90, 'The QUBO encoding',
     paragraphs('To use a binary sampler, we first translate each residue choice into a group of on-or-off variables.') + figure('one-hot.svg', 'At an illustrative site allowing A, C and G, choosing C maps to the bits 0, 1, 0; selecting 1, 1, 0 is invalid.', 'This one-hot representation gives each allowed residue its own bit.') + paragraphs('We add an energy penalty when a group selects no residue or more than one residue.', 'For every valid sequence, the binary energy remains exactly the same as the original protein energy.'),
     'Define QUBO as quadratic unconstrained binary optimization, used here as an energy representation for sampling. Give a verbal example: if a site allows A, C, and G, then 100, 010, and 001 are valid choices, while 000 and 110 are invalid. The penalty is finite, so invalid words can still occur during the inner trajectory. Substituting z=(s+1)/2 produces the Ising representation. Exact energy agreement on valid states is the essential compiler property.', 'docs/encoded_method.md, Exact compiler'),
    ('We propose a change, then check it.', 120, 'The sampling method',
     paragraphs('Once the sequence is encoded, we run a short Ising trajectory to propose a new binary state.') + figure('proposal-model.svg', 'The current sequence passes through a finite Ising trajectory; valid endpoints receive correction and invalid endpoints retain the current state.', 'A fair coin chooses forward or reverse update order for a fixed number of sweeps.') + paragraphs('We use a Metropolis–Hastings correction to keep valid proposals consistent with the protein distribution.', 'If the endpoint is invalid, we keep the old sequence and record it again instead of drawing a replacement.'),
     'Spend time on why this is a proposal mechanism rather than an assumption that the Ising sampler has equilibrated. A mixture of a forward trajectory and its reverse is reversible under the proposal distribution. This gives the proposal ratio used in the correction. Invalid endpoints have zero target probability and are rejected once, without redrawing until valid. Redrawing would change the proposal law. Rejections must be recorded to preserve the intended distribution.', 'docs/encoded_method.md, Finite proposal and proof'),
    ('We need the sampler to keep exploring.', 75, 'Following the samples',
     paragraphs('After checking the target distribution, we ask how much new information the sampler produces.') + figure('sampling-trace.svg', 'An example sequence trace repeats ACG before moving to ATG, then returns to ACG.') + paragraphs('Repeated sequences belong in the trace, but long repeats give us less new information.', 'We therefore measure exploration and uncertainty alongside the time each run takes.'),
     'Explain the cancellation of target and proposal energy differences when the temperatures match. Acceptance of a valid endpoint does not imply fast exploration, especially if most endpoints are invalid or unchanged. The implementation can model coefficient scaling and quantization in software, but those are explicit surrogate models, not calibrated hardware behavior. Transition checks validate small exact systems while empirical diagnostics assess larger runs.', 'docs/encoded_method.md, Finite proposal and proof; Scaling, rounding, and numerical limits'),
    ('The current implementation runs in software.', 75, 'Implementation status',
     paragraphs('We have implemented both routes, so we can now test that distinction between correctness and speed.', 'The encoded route runs in NumPy, JAX, or the official THRML software package.') + bullets('We first check small systems where we can calculate the exact probabilities.', 'We then retain the sampling traces, seeds, timings, and diagnostics for each run.') + paragraphs('For the timing comparison, we run the direct method in JAX and the encoded method in THRML.'),
     'THRML is the official software integration in this project and is not a physical chip execution. Mention that the adapter uses singleton blocks and checks update ordering. A 64-spin cap bounds official THRML compilation in this adapter. Dense logical interactions still need a verified mapping to hardware. Do not confuse software support with access to a physical device.', 'docs/encoded_method.md, Official THRML integration; README.md'),
    ('The direct sampler is faster on the small case.', 90, 'Measured timing',
     paragraphs('We began with a four-state target so that both methods could be checked against a known answer.') + figure('timing.svg', 'Direct median 17.885 seconds; encoded THRML median 26.754 seconds, from three runs each.', 'Each bar shows the median of three repetitions on the same small target.') + paragraphs('All six runs passed at <strong>512 observations per chain</strong>, which gives us a common accuracy checkpoint.', 'At that checkpoint, the encoded route took about <strong>1.50 times as long</strong>, so direct sampling won this comparison.'),
     'Describe the qualification gate before discussing speed. Each method has three repeats, and all six runs qualified at the first declared checkpoint with four chains and 2,048 total observations per run. The plotted values are recorded elapsed times, including the costs captured by the workflow, rather than isolated spin update times. The ranges show limited repeated measurements. The practical conclusion is that direct sampling wins on this particular small case. Larger-target behavior needs separate evidence.', 'output/week07/method_repeats_verified/comparison.json; docs/following_weeks.md'),
    ('The larger run needs more precise estimates.', 90, 'Broad-reference results',
     paragraphs('The small case tells us the comparison works, but we also need a trustworthy reference for the broader target.', 'Our larger direct run produced <strong>131,072 observations</strong> in about <strong>202 seconds</strong>.') + bullets('Its energy effective sample size reached 5,833, which cleared the required minimum of 4,000.', 'However, the mean uncertainties for energy and distance were still too large.', 'The frequency checks also left 121 uncertainties above the limit and 11 unknown.') + paragraphs('The next step is to bring those uncertainties within the limits required for the reference.'),
     'Explain effective sample size as a way to account for correlation between observations. The larger run used 32,768 observations per chain and new seeds. It improved several diagnostics, but qualification requires every declared criterion, not one attractive number. No candidates were exported from the reference run. This reference is intended to support trustworthy evaluation of the broader target, so its failure blocks a broad accuracy claim.', 'docs/week08/README.md; output/week08/reference_scaling/reference_32768/manifest.json'),
    ('We can see where exploration falls short.', 90, 'Where exploration is limited',
     paragraphs('To understand what is holding the reference back, we looked at individual residue frequencies.') + figure('uncertainty.svg', 'Position 65 phenylalanine uncertainty is 0.00683, versus the required 0.0025 limit.', 'Phenylalanine at position 65 has the largest observed residue uncertainty.') + paragraphs('Its uncertainty is more than twice the limit, so this estimate is still too imprecise.', 'We also never observed 11 residue pairs involving positions 11 and 15, which leaves their uncertainties unknown.'),
     'Define frequency uncertainty as uncertainty in a model-derived occurrence estimate, not uncertainty from an activity assay. The largest observed residue uncertainty is about 2.73 times its threshold. Positions 82 and 44 are the next largest contributors in the report. The 11 unknown pair uncertainties are not known zeros. They motivate better exploration of coupled residue choices rather than relaxing criteria or inventing biological restrictions.', 'docs/week08/README.md; output/week08/reference_scaling/reference_32768_frequency_blockers.csv'),
    ('We made the saved results much smaller.', 75, 'An engineering result',
     paragraphs('While investigating those limits, we also removed work that the reference run did not need.') + figure('storage.svg', 'Saved artifacts decreased from 94.47 MB to 12.64 MB.', 'Replaying the earlier reference produced identical trace arrays and diagnostic values.') + paragraphs('By skipping unused candidate construction, we reduced saved artifacts by <strong>86.6%</strong>.', 'We kept the same samples and diagnostics, so the smaller files are easier to store and inspect.'),
     'The reference workflow previously materialized a candidate pool even though it exported no candidates. Removing that unnecessary work greatly reduces serialized artifacts while retaining the full trace and diagnostics. The replay used the earlier 16,384-observation-per-chain setup, not the larger run just discussed. Exact replay agreement supports the engineering change. Timing from different measurement boundaries should not be used to manufacture a speed claim.', 'docs/week08/README.md; output/week08/reference_scaling/replay_equivalence.json'),
    ('The next experiment should improve exploration.', 120, 'Next steps',
     paragraphs('The remaining problem is therefore exploration, and that is where the next experiment should focus.') + bullets('First, we need to understand why the sampler misses particular combinations of residues.', 'We can then test changes while keeping the target and accuracy criteria explicit.', 'Once the reference passes, we can make a meaningful timing comparison on the broader target.') + paragraphs('We will then take selected sequences through biological review and experimental validation.'),
     'Separate three future decisions. First, improve the broad reference enough to meet uncertainty requirements. Second, evaluate encoded settings against the direct baseline using comparable end-to-end costs and accuracy. Third, review model-derived candidates biologically and obtain assay evidence. Changing model support changes the scientific question and must be justified rather than used to pass diagnostics. These are next steps, not completed experiments.', 'docs/week08/README.md; docs/following_weeks.md'),
    ('We now have a clear experiment to build on.', 60, 'Takeaways',
     paragraphs('We started by asking whether a binary sampler could explore the same protein model more efficiently.', 'We now have an exact encoding and a corrected sampling method that let us test that question.') + paragraphs('So far, direct sampling is faster on the small case, and the broader reference still needs better exploration.', '<strong>Our next milestone is to make that broader comparison reliable enough to tell us where the binary route is useful.</strong>'),
     'Close by returning to the original question. We now have an implemented and carefully specified route from a restricted protein energy model to binary proposals with a correction. The recorded timing favors the direct method on the small benchmark, and the broad reference is not yet qualified. The value of the current work is a reproducible basis for deciding whether the encoded route becomes useful. Invite questions after the twenty-minute presentation.', 'docs/encoded_method.md; docs/week08/README.md; output/week07/method_repeats_verified/comparison.json'),
]


def document(title, body, script=True):
    return f'<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1">\n<title>{escape(title)} | Protein sampling</title>\n<link rel="stylesheet" href="styles.css">\n' + ('<script src="slides.js" defer></script>\n' if script else '') + f'</head>\n<body>\n{body}\n</body>\n</html>\n'


def build():
    notes = '<main class="notes"><h1>Presenter scripts</h1><p>Neil introduces the model and method. Daniel presents the results and next experiment.</p><p><a href="script-neil.html">Read Neil’s script</a> · <a href="script-daniel.html">Read Daniel’s script</a></p><p>The full talk is paced for 20 minutes, with the handoff at 10 minutes.</p>'
    presenter_sections = {'Neil Huang': [], 'Daniel Sun': []}
    elapsed = 0
    for index, (title, seconds, eyebrow, content, narration, source) in enumerate(slides):
        filename = 'index.html' if index == 0 else f'slide-{index+1:02}.html'
        if index == 0:
            body = '<main class="slide cover" aria-label="Title slide"><img class="protein" src="assets/protein-ascii.svg" alt="Original ASCII protein artwork extracted from the supplied PowerPoint."><div class="title"><h1>QUBO Machine<br>for Protein Sampling</h1><p class="authors">Neil Huang &amp; Daniel Sun</p></div></main>'
        else:
            body = f'<main class="slide"><p class="eyebrow">{eyebrow}</p><h1>{title}</h1>{content}</main>'
        (ROOT / filename).write_text(document(title, body), encoding='utf-8')
        finish = elapsed + seconds
        timing = f'{elapsed//60:02}:{elapsed%60:02}–{finish//60:02}:{finish%60:02}'
        speech = SCRIPTS[index + 1]
        speaker = speech['speaker']
        spoken_html = ''.join(f'<p>{escape(paragraph)}</p>' for paragraph in speech['paragraphs'])
        section = f'<section id="slide-{index+1}"><p class="eyebrow">{timing} · {speaker}</p><h2><a href="{filename}">{title}</a></h2>{spoken_html}<details><summary>Sources</summary><p>{escape(source)}</p></details></section>'
        notes += section
        presenter_sections[speaker].append(section)
        elapsed = finish
    assert elapsed == 1200, elapsed
    (ROOT / 'notes.html').write_text(document('Presenter scripts', notes + '</main>', False), encoding='utf-8')
    for speaker, sections in presenter_sections.items():
        short_name = speaker.split()[0].lower()
        role = 'The model and sampling method · 00:00–10:00' if short_name == 'neil' else 'The results and next experiment · 10:00–20:00'
        intro = f'<main class="notes"><p class="eyebrow">{role}</p><h1>{speaker}’s script</h1><p>This is your ten-minute speaking part. The times below follow the full presentation clock.</p><p><a href="notes.html">Open the complete talk</a></p>'
        (ROOT / f'script-{short_name}.html').write_text(document(f'{speaker}’s script', intro + ''.join(sections) + '</main>', False), encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-pptx', type=Path)
    options = parser.parse_args()
    ASSETS.mkdir(exist_ok=True)
    if options.source_pptx:
        extract_art(options.source_pptx)
    charts()
    build_diagrams(ASSETS)
    build()
    print(f'Built {len(slides)} slides with 20 minutes of speaker notes.')
