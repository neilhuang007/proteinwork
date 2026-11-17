from html import escape
from pathlib import Path
import base64

ACCENT = '#d4a27f'
WHITE = '#faf9f5'
MUTED = '#a3a3a3'
FAINT = '#51463e'


def text(horizontal, vertical, value, color=MUTED, size=16, anchor='start'):
    return f'<text x="{horizontal}" y="{vertical}" fill="{color}" font-size="{size}" text-anchor="{anchor}">{escape(value)}</text>'


def path(points, color=FAINT, width=1.5):
    return f'<path d="{points}" fill="none" stroke="{color}" stroke-width="{width}"/>'


def arrow(start, end, vertical):
    return path(f'M {start} {vertical} H {end} m -7 -5 7 5 -7 5', ACCENT)


def residue(horizontal, vertical, symbol, active=False, radius=21):
    color = ACCENT if active else '#252321'
    foreground = '#090909' if active else WHITE
    return f'<circle cx="{horizontal}" cy="{vertical}" r="{radius}" fill="{color}"/>' + text(horizontal, vertical + 6, symbol, foreground, 18, 'middle')


def save(folder, name, content, height=240):
    font = base64.b64encode((folder / 'fonts/AnthropicSans-Regular.otf').read_bytes()).decode('ascii')
    style = f'<style>@font-face{{font-family:Anthropic;src:url(data:font/otf;base64,{font}) format("opentype")}}text{{font-family:Anthropic,sans-serif}}</style>'
    svg = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 980 {height}" width="980" height="{height}">{style}{content}</svg>'
    (folder / name).write_text(svg, encoding='utf-8')


def build_diagrams(folder: Path):
    content = text(0, 26, 'We begin with one protein model.', WHITE)
    content += text(0, 63, 'Residue preferences') + text(0, 92, '+ pairwise interactions')
    content += arrow(235, 305, 83)
    content += text(335, 36, 'We take two routes.', WHITE)
    content += text(335, 80, 'Direct residue sampling', WHITE)
    content += text(335, 128, 'Binary proposals + correction', ACCENT)
    content += path('M 620 74 C 680 74 665 96 710 96 M 620 123 C 680 123 665 96 710 96', ACCENT)
    content += arrow(710, 755, 96)
    for row, sequence in enumerate(['ACGTA', 'AGGTA', 'ACGAA']):
        for column, symbol in enumerate(sequence):
            content += text(782 + column * 33, 57 + row * 37, symbol, ACCENT if (row, column) in [(1, 1), (2, 3)] else WHITE, 22)
    content += text(0, 207, 'Both routes should reproduce the same sequence probabilities.', ACCENT)
    save(folder, 'model-routes.svg', content)

    content = text(490, 30, 'The model scores these two residue choices together.', ACCENT, 18, 'middle')
    content += path('M 285 139 C 285 64 695 64 695 139', ACCENT, 1.8)
    names = ['Alanine', 'Cysteine', 'Glycine', 'Threonine', 'Alanine']
    for column, symbol in enumerate('ACGTA'):
        horizontal = 80 + column * 205
        if column < 4:
            content += path(f'M {horizontal + 38} 177 H {horizontal + 167}', FAINT, 1)
        content += residue(horizontal, 177, symbol, column in [1, 3], radius=30)
        content += text(horizontal, 239, names[column], WHITE, 16, 'middle')
        content += text(horizontal, 267, f'Position {column + 1}', MUTED, 13, 'middle')
    save(folder, 'potts-model.svg', content, 295)

    content = text(0, 24, 'Suppose one position allows alanine, cysteine, or glycine.', WHITE)
    for index, symbol in enumerate('ACG'):
        horizontal = 330 + index * 110
        content += text(horizontal, 65, symbol, MUTED, 16, 'middle')
        content += residue(horizontal, 108, '1' if index == 1 else '0', index == 1, 26)
    content += residue(62, 108, 'C', True, 30) + arrow(120, 255, 108)
    content += arrow(610, 700, 108) + text(730, 114, 'Exactly one bit is active.', WHITE)
    content += text(0, 195, 'C becomes 0 1 0.', ACCENT)
    content += text(330, 195, '1 1 0 is invalid because it selects two residues.', MUTED)
    content += text(0, 232, 'The compiler preserves the protein energy for every valid encoding.', MUTED, 14)
    save(folder, 'one-hot.svg', content, 250)

    content = text(0, 24, 'We propose a binary state, then decide what to record.', WHITE)
    content += text(0, 99, 'Current sequence', WHITE) + text(0, 136, 'A C G', ACCENT, 24)
    content += arrow(155, 215, 110)
    content += text(242, 99, 'Finite Ising trajectory', WHITE) + text(242, 133, 'Forward or reverse sweeps', MUTED, 14)
    content += arrow(425, 485, 110)
    content += text(512, 99, 'Check the endpoint.', WHITE)
    content += path('M 670 110 C 709 110 709 64 748 64 M 670 110 C 709 110 709 177 748 177', ACCENT)
    content += text(770, 60, 'If valid, apply correction.', ACCENT, 16)
    content += text(770, 87, 'Accept it or keep the old state.', MUTED, 14)
    content += text(770, 175, 'If invalid, keep the old state.', WHITE, 14)
    content += path('M 748 184 C 670 244 95 244 65 160', FAINT)
    content += text(285, 250, 'We record every attempt, including repeats.', MUTED, 14)
    save(folder, 'proposal-model.svg', content, 270)

    content = text(0, 24, 'A sampling run includes both moves and repeated sequences.', WHITE)
    for index, sequence in enumerate(['ACG', 'ACG', 'ACG', 'ATG', 'ATG', 'ACG']):
        horizontal = 60 + index * 158
        if index < 5:
            content += arrow(horizontal + 47, horizontal + 105, 109)
        content += text(horizontal, 117, sequence, ACCENT if index in [3, 5] else WHITE, 25, 'middle')
        content += text(horizontal, 163, 'Start' if index == 0 else ('Move' if index in [3, 5] else 'Repeat'), MUTED, 14, 'middle')
    content += text(0, 226, 'In this example, six observations contain two distinct sequences.', MUTED, 14)
    save(folder, 'sampling-trace.svg', content, 245)
