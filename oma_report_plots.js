/* Plot controls for oma_report.html — hide/show + fixed axis scales. */
(function () {
  'use strict';

  var axisLocks = {};

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
    axisLocks[divId] = {
      'xaxis.range': padRange(xr[0], xr[1], 0.06),
      'yaxis.range': padRange(Math.max(0, yr[0]), yr[1], 0.08),
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
    var names = [];
    gd.data.forEach(function (tr) {
      if (tr.name) names.push(tr.name);
      (tr.y || []).forEach(function (v) { ys.push(v); });
    });
    var yr = minMax(ys);
    axisLocks[divId] = {
      'yaxis.range': padRange(Math.max(0, yr[0]), yr[1], 0.06),
      'yaxis.autorange': false,
      'xaxis.autorange': false,
      'xaxis.categoryorder': 'array',
      'xaxis.categoryarray': names,
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

  function wireAxisLock(divId) {
    var gd = document.getElementById(divId);
    if (!gd) return;
    gd.on('plotly_legendclick', function () {
      setTimeout(function () { restoreAxisLock(divId); }, 0);
      return true;
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
