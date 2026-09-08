/* Plot controls for oma_report.html — uncheck-all + fixed PCA axis scale. */
(function () {
  'use strict';

  function plotIds() {
    return ['plot_pctid_kde', 'plot_sim_violins', 'plot_pca'];
  }

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
    addControls('plot_pctid_kde');
    addControls('plot_sim_violins');
    addControls('plot_pca');
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
