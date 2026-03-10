const dnaColors = {
  A: '#16A34A',
  T: '#DC2626',
  C: '#2563EB',
  G: '#D97706',
  N: '#64748B'
};

function createSvgEl(name, attrs = {}) {
  const el = document.createElementNS('http://www.w3.org/2000/svg', name);
  for (const [key, value] of Object.entries(attrs)) {
    el.setAttribute(key, String(value));
  }
  return el;
}

function clearNode(node) {
  while (node.firstChild) node.removeChild(node.firstChild);
}

function renderSequenceRow(parent, sequence, x, y, boxWidth, opts = {}) {
  const {
    highlightStart = -1,
    highlightLength = 0,
    mismatchIndex = -1,
    pamStart = -1,
    label = '',
    fontSize = 20
  } = opts;

  if (label) {
    parent.appendChild(createSvgEl('text', {
      x,
      y: y - 16,
      fill: '#50607B',
      'font-family': 'Sora, sans-serif',
      'font-size': 13,
      'font-weight': 600
    })).textContent = label;
  }

  sequence.split('').forEach((base, index) => {
    const boxX = x + index * (boxWidth + 4);
    const inGuideWindow = index >= highlightStart && index < highlightStart + highlightLength;
    const inPamWindow = pamStart >= 0 && index >= pamStart && index < pamStart + 3;
    const isMismatch = index === mismatchIndex;

    let fill = '#FFFFFF';
    let stroke = '#D8E1EC';

    if (inGuideWindow) {
      fill = 'rgba(29, 78, 216, 0.10)';
      stroke = '#1D4ED8';
    }
    if (inPamWindow) {
      fill = 'rgba(180, 83, 9, 0.14)';
      stroke = '#B45309';
    }
    if (isMismatch) {
      fill = 'rgba(220, 38, 38, 0.14)';
      stroke = '#DC2626';
    }

    parent.appendChild(createSvgEl('rect', {
      x: boxX,
      y,
      width: boxWidth,
      height: 36,
      rx: 10,
      fill,
      stroke,
      'stroke-width': isMismatch ? 2.5 : 1.5
    }));

    const text = createSvgEl('text', {
      x: boxX + boxWidth / 2,
      y: y + 24,
      'text-anchor': 'middle',
      fill: isMismatch ? '#B91C1C' : (dnaColors[base] || '#162033'),
      'font-family': 'JetBrains Mono, monospace',
      'font-size': fontSize,
      'font-weight': 700
    });
    text.textContent = base;
    parent.appendChild(text);
  });
}

(function initGuideSimulation() {
  const genomeGroup = document.getElementById('guideGenomeTrack');
  const queryGroup = document.getElementById('guideQueryTrack');
  const pamGroup = document.getElementById('guidePamBadge');
  const statusEl = document.getElementById('guideStatus');
  const nextBtn = document.getElementById('guideNextBtn');
  const mismatchBtn = document.getElementById('guideMismatchBtn');
  const resetBtn = document.getElementById('guideResetBtn');

  if (!genomeGroup || !queryGroup || !pamGroup) return;

  const genome = 'CCATTAGATTGTTTGTAGATTTGACACGGATAC';
  const guide = 'ATTGTTTGTAGATTTGACAC';
  const pam = 'CGG';
  const genomePamStart = 28;
  const guideStart = 7;
  let step = 0;
  let mismatch = false;

  const messages = [
    'Step 1: The guide scans for a sequence match in the genome. At this point, matching letters are necessary but not yet sufficient.',
    'Step 2: A strong candidate appears. The guide aligns with a 20-base target window in the genome.',
    'Step 3: Now check the PAM. Here the neighboring motif is compatible, so the site stays in play.',
    'Step 4: With both guide match and PAM satisfied, this location becomes a plausible CRISPR target.'
  ];

  function draw() {
    clearNode(genomeGroup);
    clearNode(queryGroup);
    clearNode(pamGroup);

    const highlightLength = step >= 1 ? guide.length : 0;
    const pamStart = step >= 2 ? genomePamStart : -1;
    const mismatchIndex = mismatch ? guideStart + 11 : -1;

    renderSequenceRow(genomeGroup, genome, 46, 70, 20, {
      highlightStart: guideStart,
      highlightLength,
      mismatchIndex,
      pamStart,
      label: 'Genome letters'
    });

    const querySequence = mismatch
      ? guide.slice(0, 11) + 'C' + guide.slice(12) + pam
      : guide + pam;

    renderSequenceRow(queryGroup, querySequence, 122, 196, 20, {
      highlightStart: 0,
      highlightLength: guide.length,
      mismatchIndex: mismatch ? 11 : -1,
      pamStart: guide.length,
      label: 'Guide + PAM expectation'
    });

    const bracket = createSvgEl('path', {
      d: 'M197 160 C197 148 197 148 209 148 L612 148 C624 148 624 148 624 160',
      fill: 'none',
      stroke: '#1D4ED8',
      'stroke-width': step >= 1 ? 2.5 : 0,
      'stroke-linecap': 'round'
    });
    genomeGroup.appendChild(bracket);

    if (step >= 2) {
      const pamBadge = createSvgEl('g');
      pamBadge.appendChild(createSvgEl('rect', {
        x: 700,
        y: 168,
        width: 118,
        height: 36,
        rx: 18,
        fill: mismatch ? 'rgba(220,38,38,0.12)' : 'rgba(180,83,9,0.14)',
        stroke: mismatch ? '#DC2626' : '#B45309',
        'stroke-width': 1.6
      }));
      const text = createSvgEl('text', {
        x: 759,
        y: 191,
        'text-anchor': 'middle',
        fill: mismatch ? '#B91C1C' : '#9A3412',
        'font-family': 'Sora, sans-serif',
        'font-size': 14,
        'font-weight': 700
      });
      text.textContent = mismatch ? 'Mismatch detected' : 'PAM = CGG';
      pamBadge.appendChild(text);
      pamGroup.appendChild(pamBadge);
    }

    if (mismatch) {
      statusEl.textContent = 'A single mismatch in the guide window can be enough to lower performance or break the target entirely, depending on position and assay design.';
    } else {
      statusEl.textContent = messages[step];
    }
  }

  nextBtn.addEventListener('click', () => {
    mismatch = false;
    step = (step + 1) % messages.length;
    draw();
  });

  mismatchBtn.addEventListener('click', () => {
    mismatch = !mismatch;
    if (step < 1) step = 1;
    draw();
  });

  resetBtn.addEventListener('click', () => {
    step = 0;
    mismatch = false;
    draw();
  });

  draw();
})();

(function initConservationSimulation() {
  const spikeGroup = document.getElementById('conservationSpikeGrid');
  const polyGroup = document.getElementById('conservationPolyGrid');
  const barsGroup = document.getElementById('conservationBars');
  const slider = document.getElementById('yearSlider');
  const yearLabel = document.getElementById('yearLabel');
  const statusEl = document.getElementById('conservationStatus');

  if (!spikeGroup || !polyGroup || !barsGroup || !slider || !yearLabel || !statusEl) return;

  const years = [2020, 2021, 2022, 2023, 2024, 2025, 2026];
  const spikeConservation = [98, 88, 76, 65, 58, 49, 42];
  const polyConservation = [99, 98, 98, 97, 97, 96, 95];

  function drawGrid(group, x, y, activeCells, color, fadedColor) {
    clearNode(group);
    for (let row = 0; row < 4; row += 1) {
      for (let col = 0; col < 7; col += 1) {
        const index = row * 7 + col;
        group.appendChild(createSvgEl('rect', {
          x: x + col * 42,
          y: y + row * 42,
          width: 32,
          height: 32,
          rx: 9,
          fill: index < activeCells ? color : fadedColor,
          stroke: 'rgba(22,32,51,0.06)'
        }));
      }
    }
  }

  function drawBars(yearIndex) {
    clearNode(barsGroup);

    const specs = [
      {
        x: 68,
        y: 310,
        width: 230,
        label: 'Exposed region retained',
        fill: '#F97316',
        pct: spikeConservation[yearIndex]
      },
      {
        x: 520,
        y: 310,
        width: 230,
        label: 'Constrained region retained',
        fill: '#0F766E',
        pct: polyConservation[yearIndex]
      }
    ];

    specs.forEach((spec) => {
      const panel = createSvgEl('rect', {
        x: spec.x - 10,
        y: spec.y - 22,
        width: spec.width + 54,
        height: 42,
        rx: 14,
        fill: 'rgba(255,255,255,0.74)',
        stroke: 'rgba(22,32,51,0.06)'
      });
      barsGroup.appendChild(panel);

      const label = createSvgEl('text', {
        x: spec.x,
        y: spec.y - 2,
        fill: '#50607B',
        'font-family': 'Karla, sans-serif',
        'font-size': 12,
        'font-weight': 600
      });
      label.textContent = spec.label;
      barsGroup.appendChild(label);

      barsGroup.appendChild(createSvgEl('rect', {
        x: spec.x,
        y: spec.y + 8,
        width: spec.width,
        height: 12,
        rx: 6,
        fill: 'rgba(22,32,51,0.08)'
      }));

      barsGroup.appendChild(createSvgEl('rect', {
        x: spec.x,
        y: spec.y + 8,
        width: spec.width * (spec.pct / 100),
        height: 12,
        rx: 6,
        fill: spec.fill
      }));

      const pct = createSvgEl('text', {
        x: spec.x + spec.width + 10,
        y: spec.y + 19,
        fill: spec.fill,
        'font-family': 'Sora, sans-serif',
        'font-size': 12,
        'font-weight': 700
      });
      pct.textContent = `${spec.pct}%`;
      barsGroup.appendChild(pct);
    });
  }

  function update() {
    const yearIndex = Number(slider.value);
    yearLabel.textContent = String(years[yearIndex]);

    const spikeCells = Math.round((spikeConservation[yearIndex] / 100) * 28);
    const polyCells = Math.round((polyConservation[yearIndex] / 100) * 28);
    drawGrid(spikeGroup, 68, 118, spikeCells, '#FB7185', '#FFE4E6');
    drawGrid(polyGroup, 520, 118, polyCells, '#2DD4BF', '#CCFBF1');
    drawBars(yearIndex);

    statusEl.textContent = `${years[yearIndex]}: the exposed region keeps ${spikeConservation[yearIndex]}% of exact targets, while the constrained replication region keeps ${polyConservation[yearIndex]}%. This is why conserved internal targets often age better than highly visible surface targets.`;
  }

  slider.addEventListener('input', update);
  update();
})();

(function initIndexSimulation() {
  const barsGroup = document.getElementById('indexSimBars');
  const playBtn = document.getElementById('indexPlayBtn');
  const resetBtn = document.getElementById('indexResetBtn');
  const input = document.getElementById('indexQueryInput');
  const statusEl = document.getElementById('indexStatus');

  if (!barsGroup || !playBtn || !resetBtn || !input || !statusEl) return;

  let timer = null;

  function buildIntervals(query) {
    const clean = query.toUpperCase().replace(/[^ACGT]/g, '').slice(0, 8) || 'ATTGTT';
    const base = 3200000;
    return clean.split('').map((char, index) => {
      const shrink = Math.max(1, Math.round(base / Math.pow(5.2, index + 1)));
      return {
        char,
        label: `Step ${index + 1}`,
        count: index === clean.length - 1 ? Math.max(1, Math.round(shrink / 4)) : shrink
      };
    });
  }

  function renderIntervals(intervals, activeCount = 0) {
    clearNode(barsGroup);

    const maxCount = intervals.length ? intervals[0].count : 1;

    intervals.forEach((interval, index) => {
      const y = 120 + index * 24;
      const width = 410 * (interval.count / maxCount);
      const active = index < activeCount;

      const label = createSvgEl('text', {
        x: 56,
        y: y + 12,
        fill: '#50607B',
        'font-family': 'Sora, sans-serif',
        'font-size': 12,
        'font-weight': 600
      });
      label.textContent = `${interval.label}: add ${interval.char}`;
      barsGroup.appendChild(label);

      barsGroup.appendChild(createSvgEl('rect', {
        x: 196,
        y,
        width: 430,
        height: 14,
        rx: 7,
        fill: 'rgba(22,32,51,0.06)'
      }));

      barsGroup.appendChild(createSvgEl('rect', {
        x: 196,
        y,
        width,
        height: 14,
        rx: 7,
        fill: active ? '#1D4ED8' : '#93C5FD'
      }));

      const count = createSvgEl('text', {
        x: 632,
        y: y + 12,
        fill: active ? '#1E3A8A' : '#64748B',
        'font-family': 'JetBrains Mono, monospace',
        'font-size': 12,
        'font-weight': 700,
        'text-anchor': 'start'
      });
      count.textContent = `${interval.count.toLocaleString()} candidates`;
      barsGroup.appendChild(count);
    });
  }

  function stopTimer() {
    if (timer) {
      window.clearInterval(timer);
      timer = null;
    }
  }

  function play() {
    stopTimer();
    const query = input.value.toUpperCase().replace(/[^ACGT]/g, '').slice(0, 8) || 'ATTGTT';
    input.value = query;
    const intervals = buildIntervals(query);
    let active = 0;
    renderIntervals(intervals, 0);
    statusEl.textContent = `Running a conceptual exact search for ${query}. Watch the candidate interval shrink as each base removes impossible matches.`;

    timer = window.setInterval(() => {
      active += 1;
      renderIntervals(intervals, active);
      if (active >= intervals.length) {
        stopTimer();
        const last = intervals[intervals.length - 1];
        statusEl.textContent = `After ${query.length} steps, the candidate set collapses to about ${last.count.toLocaleString()} exact matches in this teaching example. Indexed search wins by eliminating impossibilities early.`;
      }
    }, 520);
  }

  playBtn.addEventListener('click', play);
  resetBtn.addEventListener('click', () => {
    stopTimer();
    input.value = 'ATTGTT';
    const intervals = buildIntervals(input.value);
    renderIntervals(intervals, 0);
    statusEl.textContent = 'Use a short DNA pattern. The bars show a teaching example of how indexed exact search collapses possible matches into a tiny interval.';
  });

  input.addEventListener('keydown', (event) => {
    if (event.key === 'Enter') {
      event.preventDefault();
      play();
    }
  });

  renderIntervals(buildIntervals(input.value), 0);
})();

(function initSectionTracking() {
  const railLinks = [...document.querySelectorAll('.learn-rail a')];
  const sections = railLinks
    .map((link) => document.querySelector(link.getAttribute('href')))
    .filter(Boolean);

  if (!railLinks.length || !sections.length || !('IntersectionObserver' in window)) return;

  const linkById = new Map(railLinks.map((link) => [link.getAttribute('href').slice(1), link]));

  const observer = new IntersectionObserver((entries) => {
    entries.forEach((entry) => {
      if (!entry.isIntersecting) return;
      railLinks.forEach((link) => link.classList.remove('active'));
      const activeLink = linkById.get(entry.target.id);
      if (activeLink) activeLink.classList.add('active');
    });
  }, {
    rootMargin: '-30% 0px -55% 0px',
    threshold: 0.1
  });

  sections.forEach((section) => observer.observe(section));
})();
