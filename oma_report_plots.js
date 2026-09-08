/* Plot controls for oma_report.html — uncheck-all + fixed PCA axis scale. */
(function () {
  'use strict';

  function plotIds() {
    return ['plot_div_kde', 'plot_div_kde_filtered', 'plot_pctid_kde',
            'plot_sim_violins', 'plot_pca'];
  }

  function allTraces(divId) {
    var gd = document.getElementById(divId);
    if (!gd || !gd.data) return [];
    return gd.data.map(function (_, i) { return i; });
  }

  function setAllVisible(divId, visible) {
    var gd = document.getElementById(divId);
    if (!gd || !gd.data || !gd.data.length) return;
    var vis = gd.data.map(function () { return visible; });
    Plotly.restyle(gd, { visible: vis }, allTraces(divId));
  }

  function addControls(divId, label) {
    var plot = document.getElementById(divId);
    if (!plot) return;
    var bar = document.createElement('div');
    bar.className = 'plot-controls';
    bar.style.cssText = 'margin:6px 0 8px 0;font-size:.85rem;';
    bar.innerHTML =
      '<button type="button" class="plot-btn" data-act="none">Uncheck all</button> ' +
      '<button type="button" class="plot-btn" data-act="all">Check all</button>' +
      (label ? ' <span class="muted small">' + label + '</span>' : '');
    plot.parentNode.insertBefore(bar, plot);
    bar.querySelectorAll('.plot-btn').forEach(function (btn) {
      btn.style.cssText = 'margin-right:6px;padding:3px 10px;cursor:pointer;';
      btn.addEventListener('click', function () {
        setAllVisible(divId, btn.getAttribute('data-act') === 'all');
        if (divId === 'plot_pca') restorePcaRange();
      });
    });
  }

  var pcaRange = null;

  function capturePcaRange() {
    var gd = document.getElementById('plot_pca');
    if (!gd || !gd.data || !gd.data.length) return;
    var xs = [], ys = [];
    gd.data.forEach(function (tr) {
      if (!tr.x || tr.visible === 'legendonly') return;
      tr.x.forEach(function (v) { xs.push(v); });
      tr.y.forEach(function (v) { ys.push(v); });
    });
    if (!xs.length) {
      gd.data.forEach(function (tr) {
        (tr.x || []).forEach(function (v) { xs.push(v); });
        (tr.y || []).forEach(function (v) { ys.push(v); });
      });
    }
    if (!xs.length) return;
    function pad(min, max) {
      var m = (max - min) * 0.06 || 0.5;
      return [min - m, max + m];
    }
    pcaRange = {
      x: pad(Math.min.apply(null, xs), Math.max.apply(null, xs)),
      y: pad(Math.min.apply(null, ys), Math.max.apply(null, ys)),
    };
    Plotly.relayout(gd, {
      'xaxis.range': pcaRange.x,
      'yaxis.range': pcaRange.y,
      uirevision: 'oma-pca-fixed',
    });
  }

  function restorePcaRange() {
    if (!pcaRange) return;
    var gd = document.getElementById('plot_pca');
    if (!gd) return;
    Plotly.relayout(gd, {
      'xaxis.range': pcaRange.x,
      'yaxis.range': pcaRange.y,
    });
  }

  function wirePcaLegendLock() {
    var gd = document.getElementById('plot_pca');
    if (!gd) return;
    gd.on('plotly_legendclick', function () {
      setTimeout(restorePcaRange, 0);
      return true;
    });
    gd.on('plotly_restyle', function () {
      setTimeout(restorePcaRange, 0);
    });
  }

  function init() {
    addControls('plot_div_kde', 'bitscore divergence (all assigned copies)');
    addControls('plot_div_kde_filtered', 'same metric, copies above step2 similarity floor');
    addControls('plot_pctid_kde', 'ssearch36 %identity (step4, Tal gallery metric)');
    addControls('plot_sim_violins', '');
    addControls('plot_pca', 'axes stay fixed when toggling legend');
    setTimeout(function () {
      capturePcaRange();
      wirePcaLegendLock();
    }, 300);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
