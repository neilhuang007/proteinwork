const pages = ['index.html', ...Array.from({length:13}, (_, index) => `slide-${String(index + 2).padStart(2, '0')}.html`)];
const filename = location.pathname.split('/').pop() || 'index.html';
const current = pages.indexOf(filename);
function navigate(offset) {
  const destination = pages[current + offset];
  if (destination) location.href = destination;
}
document.addEventListener('keydown', async event => {
  if (event.altKey || event.ctrlKey || event.metaKey || event.target.closest('a, button, input, textarea')) return;
  if (['ArrowRight', 'ArrowDown', 'PageDown', ' '].includes(event.key)) {
    event.preventDefault();
    navigate(1);
  } else if (['ArrowLeft', 'ArrowUp', 'PageUp'].includes(event.key)) {
    event.preventDefault();
    navigate(-1);
  } else if (event.key === 'Home') {
    location.href = pages[0];
  } else if (event.key === 'End') {
    location.href = pages[pages.length - 1];
  } else if (event.key.toLowerCase() === 'f') {
    try {
      if (document.fullscreenElement) await document.exitFullscreen();
      else await document.documentElement.requestFullscreen();
    } catch {}
  } else if (event.key.toLowerCase() === 'n') {
    window.open(`notes.html#slide-${current + 1}`, 'speaker-notes');
  }
});
let touchStart = null;
document.addEventListener('touchstart', event => { touchStart = event.changedTouches[0].clientX; }, {passive:true});
document.addEventListener('touchend', event => {
  if (touchStart === null) return;
  const distance = event.changedTouches[0].clientX - touchStart;
  if (Math.abs(distance) > 65) navigate(distance < 0 ? 1 : -1);
  touchStart = null;
}, {passive:true});
