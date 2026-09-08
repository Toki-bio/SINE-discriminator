/* PCA plot controls for oma_report.html */
(function () {
  'use strict';

  var axisLock = null;

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

  function capturePcaAxes() {
    var gd = document.getElementById('plot_pca');
    if (!gd || !gd.data || !gd.data.length) return;
    var xs = [];
    var ys = [];
    gd.data.forEach(function (tr) {
      (tr.x || []).forEach(function (v) { xs.push(v); });
      (tr.y || []).forEach(function (v) { ys.push(v); });
    });
    var xr = minMax(xs);
    var yr = minMax(ys);
    axisLock = {
      'xaxis.range': padRange(xr[0], xr[1], 0.06),
      'yaxis.range': padRange(yr[0], yr[1], 0.06),
      'xaxis.autorange': false,
      'yaxis.autorange': false,
    };
    Plotly.relayout(gd, axisLock);
  }

  function restorePcaAxes() {
    var gd = document.getElementById('plot_pca');
    if (gd && axisLock) Plotly.relayout(gd, axisLock);
  }

  function init() {
    var plot = document.getElementById('plot_pca');
    if (!plot) return;
    var bar = document.createElement('div');
    bar.style.cssText = 'margin:4px 0 6px 0;font-size:.72rem;';
    bar.innerHTML =
      '<button type="button" class="plot-btn" data-act="none">Hide all</button> ' +
      '<button type="button" class="plot-btn" data-act="all">Show all</button>';
    plot.parentNode.insertBefore(bar, plot);
    bar.querySelectorAll('.plot-btn').forEach(function (btn) {
      btn.style.cssText =
        'margin-right:4px;padding:1px 6px;font-size:.72rem;cursor:pointer;';
      btn.addEventListener('click', function () {
        var gd = document.getElementById('plot_pca');
        if (!gd || !gd.data) return;
        var show = btn.getAttribute('data-act') === 'all';
        var vis = gd.data.map(function () { return show ? true : 'legendonly'; });
        Plotly.restyle(gd, { visible: vis });
        setTimeout(restorePcaAxes, 0);
      });
    });
    setTimeout(function () {
      capturePcaAxes();
      plot.on('plotly_legendclick', function () {
        setTimeout(restorePcaAxes, 0);
        return true;
      });
    }, 300);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
