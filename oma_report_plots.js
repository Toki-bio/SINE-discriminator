/* Plot controls for oma_report.html — hide/show, fixed axes, legend colour. */
(function () {
  'use strict';

  var axisLocks = {};
  var colorPicker = null;

  function allTraces(divId) {
    var gd = document.getElementById(divId);
    if (!gd || !gd.data) return [];
    return gd.data.map(function (_, i) { return i; });
  }

  function setAllVisible(divId, show) {
    var gd = document.getElementById(divId);
    if (!gd || !gd.data || !gd.data.length) return;
    var vis = gd.data.map(function () { return show ? true : 'legendonly'; });
    Plotly.restyle(gd, { visible: vis }, allTraces(divId));
    setTimeout(function () { restoreAxisLock(divId); }, 0);
  }

  function addControls(divId) {
    var plot = document.getElementById(divId);
    if (!plot) return;
    var bar = document.createElement('div');
    bar.className = 'plot-controls';
    bar.style.cssText = 'margin:4px 0 6px 0;font-size:.72rem;';
    bar.innerHTML =
      '<button type="button" class="plot-btn" data-act="none">Hide all</button> ' +
      '<button type="button" class="plot-btn" data-act="all">Show all</button>';
    plot.parentNode.insertBefore(bar, plot);
    bar.querySelectorAll('.plot-btn').forEach(function (btn) {
      btn.style.cssText =
        'margin-right:4px;padding:1px 6px;font-size:.72rem;line-height:1.2;cursor:pointer;';
      btn.addEventListener('click', function () {
        setAllVisible(divId, btn.getAttribute('data-act') === 'all');
      });
    });
  }

  function minMax(vals) {
    var lo = Infinity;
    var hi = -Infinity;
    vals.forEach(function (v) {
      if (v == null || v !== v) return;
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    });
    return [lo, hi];
  }

  function padRange(lo, hi, frac) {
    if (!isFinite(lo) || !isFinite(hi)) return null;
    if (lo === hi) {
      lo -= 0.5;
      hi += 0.5;
    }
    var m = (hi - lo) * (frac || 0.06);
    return [lo - m, hi + m];
  }

  function collectXY(gd) {
    var xs = [];
    var ys = [];
    gd.data.forEach(function (tr) {
      (tr.x || []).forEach(function (v) { xs.push(v); });
      (tr.y || []).forEach(function (v) { ys.push(v); });
    });
    return { xs: xs, ys: ys };
  }

  function captureKdeAxes(divId) {
    var gd = document.getElementById(divId);
    if (!gd || !gd.data || !gd.data.length) return;
    var xy = collectXY(gd);
    var xr = minMax(xy.xs);
    var yr = minMax(xy.ys);
    var xPad = padRange(Math.max(0, xr[0]), xr[1], 0.04);
    axisLocks[divId] = {
      'xaxis.range': xPad ? [0, xPad[1]] : [0, null],
      'yaxis.range': padRange(0, yr[1], 0.08),
      'xaxis.autorange': false,
      'yaxis.autorange': false,
      uirevision: divId + '-fixed',
    };
    applyAxisLock(divId);
  }

  function captureViolinAxes(divId) {
    var gd = document.getElementById(divId);
    if (!gd || !gd.data || !gd.data.length) return;
    var ys = [];
    gd.data.forEach(function (tr) {
      (tr.y || []).forEach(function (v) { ys.push(v); });
    });
    var yr = minMax(ys);
    axisLocks[divId] = {
      'yaxis.range': [0, padRange(0, yr[1], 0.06)[1]],
      'yaxis.autorange': false,
      uirevision: divId + '-fixed',
    };
    applyAxisLock(divId);
  }

  function capturePcaAxes(divId) {
    var gd = document.getElementById(divId);
    if (!gd || !gd.data || !gd.data.length) return;
    var xs = [];
    var ys = [];
    gd.data.forEach(function (tr) {
      (tr.x || []).forEach(function (v) { xs.push(v); });
      (tr.y || []).forEach(function (v) { ys.push(v); });
    });
    var xr = minMax(xs);
    var yr = minMax(ys);
    axisLocks[divId] = {
      'xaxis.range': padRange(xr[0], xr[1], 0.06),
      'yaxis.range': padRange(yr[0], yr[1], 0.06),
      'xaxis.autorange': false,
      'yaxis.autorange': false,
      uirevision: divId + '-fixed',
    };
    applyAxisLock(divId);
  }

  function applyAxisLock(divId) {
    var gd = document.getElementById(divId);
    var lock = axisLocks[divId];
    if (!gd || !lock) return;
    Plotly.relayout(gd, lock);
  }

  function restoreAxisLock(divId) {
    applyAxisLock(divId);
  }

  function traceColor(tr) {
    if (tr.line && tr.line.color) return tr.line.color;
    if (tr.marker && tr.marker.color) return tr.marker.color;
    if (tr.fillcolor) return tr.fillcolor;
    return '#3366cc';
  }

  function ensureColorPicker() {
    if (colorPicker) return colorPicker;
    var wrap = document.createElement('div');
    wrap.id = 'plot-color-picker';
    wrap.style.cssText =
      'display:none;position:fixed;z-index:10000;background:#fff;' +
      'border:1px solid #ccc;border-radius:6px;padding:8px 10px;' +
      'box-shadow:0 4px 14px rgba(0,0,0,.18);font-size:.78rem;';
    wrap.innerHTML =
      '<div style="margin-bottom:6px;font-weight:600" id="plot-color-label"></div>' +
      '<input type="color" id="plot-color-input" style="width:100%;height:32px;border:0;padding:0;cursor:pointer">' +
      '<div style="margin-top:6px;text-align:right">' +
      '<button type="button" id="plot-color-close" style="font-size:.72rem;cursor:pointer">Close</button>' +
      '</div>';
    document.body.appendChild(wrap);
    wrap.querySelector('#plot-color-close').addEventListener('click', hideColorPicker);
    document.addEventListener('click', function (ev) {
      if (wrap.style.display !== 'none' &&
          !wrap.contains(ev.target) &&
          !ev.target.closest('.legend')) {
        hideColorPicker();
      }
    });
    colorPicker = wrap;
    return wrap;
  }

  function hideColorPicker() {
    if (!colorPicker) return;
    colorPicker.style.display = 'none';
    colorPicker._gd = null;
    colorPicker._idx = null;
  }

  function openColorPicker(gd, idx, clientX, clientY) {
    var tr = gd.data[idx];
    if (!tr) return;
    var picker = ensureColorPicker();
    var input = picker.querySelector('#plot-color-input');
    var label = picker.querySelector('#plot-color-label');
    label.textContent = 'Line colour: ' + (tr.name || ('trace ' + idx));
    input.value = rgbToHex(traceColor(tr));
    picker._gd = gd;
    picker._idx = idx;
    picker.style.display = 'block';
    picker.style.left = Math.min(clientX, window.innerWidth - 220) + 'px';
    picker.style.top = Math.min(clientY, window.innerHeight - 120) + 'px';
    input.oninput = function () {
      applyTraceColor(picker._gd, picker._idx, input.value);
    };
  }

  function rgbToHex(c) {
    if (!c) return '#3366cc';
    if (c.charAt(0) === '#') {
      if (c.length === 4) {
        return '#' + c[1] + c[1] + c[2] + c[2] + c[3] + c[3];
      }
      return c.slice(0, 7);
    }
    var m = c.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)/);
    if (!m) return '#3366cc';
    function h(n) {
      var s = parseInt(n, 10).toString(16);
      return s.length === 1 ? '0' + s : s;
    }
    return '#' + h(m[1]) + h(m[2]) + h(m[3]);
  }

  function applyTraceColor(gd, idx, hex) {
    if (!gd || idx == null) return;
    var tr = gd.data[idx];
    var patch = {};
    if (tr.type === 'violin') {
      patch.fillcolor = hex;
      patch['line.color'] = hex;
    } else if (tr.type === 'scatter' || tr.mode === 'lines') {
      patch['line.color'] = hex;
    } else if (tr.marker) {
      patch['marker.color'] = hex;
    } else {
      patch['line.color'] = hex;
    }
    Plotly.restyle(gd, patch, idx);
    setTimeout(function () { restoreAxisLock(gd.id); }, 0);
  }

  function wireAxisLock(divId) {
    var gd = document.getElementById(divId);
    if (!gd) return;
    gd.on('plotly_legendclick', function () {
      setTimeout(function () { restoreAxisLock(divId); }, 0);
      return true;
    });
    gd.on('plotly_legenddoubleclick', function (ev) {
      openColorPicker(gd, ev.curveNumber, ev.event.clientX, ev.event.clientY);
      return false;
    });
    gd.on('plotly_restyle', function () {
      setTimeout(function () { restoreAxisLock(divId); }, 0);
    });
  }

  function init() {
    addControls('plot_pctid_kde');
    addControls('plot_sim_violins');
    addControls('plot_pca');
    setTimeout(function () {
      captureKdeAxes('plot_pctid_kde');
      wireAxisLock('plot_pctid_kde');
      captureViolinAxes('plot_sim_violins');
      wireAxisLock('plot_sim_violins');
      capturePcaAxes('plot_pca');
      wireAxisLock('plot_pca');
    }, 300);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
