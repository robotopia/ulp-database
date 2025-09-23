// Fixing the "JS modulo bug"
// - https://web.archive.org/web/20090717035140if_/javascript.about.com/od/problemsolving/a/modulobug.htm
// - https://stackoverflow.com/questions/4467539/javascript-modulo-gives-a-negative-result-for-negative-numbers
Number.prototype.mod = function (n) {
  "use strict";
  return ((this % n) + n) % n;
};

function calc_tspin(pulse_phase, ephemeris) {
  return (2*pulse_phase*ephemeris.p0 / (1 + Math.sqrt(1 - 2*ephemeris.p1*pulse_phase)))
}

// Function to calculate the pulse numbers and phases of TOAs
function calc_pulse_phase(mjd, ephemeris) {
  // mjd and pepoch should be in days,
  // P (the p0) in seconds
  // Pdot (the p1) is dimensionless
  var t = 86400*(mjd - ephemeris.pepoch); // in seconds
  var t_p0 = t/ephemeris.p0;
  var pulse_phase = t_p0 - 0.5*ephemeris.p1*t_p0*t_p0; // Dimensionless
  var pulse = Math.round(pulse_phase)
  var phase = (pulse_phase + 0.5).mod(1) - 0.5;
  return {pulse: pulse, phase: phase, pulse_phase: pulse_phase};
}

function calc_mjd(pulse_phase, ephemeris) {
  // This is the inverse of calc_pulse_phase
  const p0 = ephemeris.p0;
  const p1 = ephemeris.p1
  const pepoch = ephemeris.pepoch;
  const retval = calc_tspin(pulse_phase, ephemeris)/86400 + pepoch;
  return retval;
}

function calc_jacobian(pulse_phase, ephemeris) {
  const p0 = ephemeris.p0;
  const p1 = ephemeris.p1;
  const pepoch = ephemeris.pepoch;
  const tspin = calc_tspin(pulse_phase, ephemeris);

  // The elements of the Jacobian
  const dt_dpepoch = 86400;  // Because of units
  const dt_dp0 = tspin/p0;
  const dt_dp1 = tspin*tspin*tspin/(2*p0*(2*pulse_phase*p0 - tspin));
  const dt_dm = 0; // TODO: add DM into the mix

  // Return the Jacobian as an array
  return [dt_dpepoch, dt_dp0, dt_dp1, dt_dm];
}

function residual_err(pulse_phase, ephemeris, covariance_matrix) {
  const J = math.matrix(calc_jacobian(pulse_phase, ephemeris));
  const S = math.matrix(covariance_matrix);
  const JSJ = math.multiply(math.transpose(J), S, J);
  return math.sqrt(math.squeeze(JSJ));
}

function generate_toas(mjd_start, mjd_end, ephemeris) {
  pp_start = calc_pulse_phase(mjd_start, ephemeris);
  pp_end = calc_pulse_phase(mjd_end, ephemeris);

  // Get the pulse numbers
  mjds = [];
  for (let n = Math.ceil(pp_start.pulse_phase); n <= pp_end.pulse_phase; n++) {
    mjds.push(calc_mjd(n, ephemeris));
  }

  return mjds;
}

// Function to create the plot elements
function create_residual_plot_elements(parentDivId) {
  // parentDivId should be a string

  var svg = d3.select("#" + parentDivId).append("svg")
  var parentDiv = document.getElementById(parentDivId);

  svg.attr("width", parentDiv.clientWidth)
    .attr("height", parentDiv.clientHeight)

  // The main group that uses graph coordinates
  var g = svg.append('g')

  // Add axis labels
  var xlabel = g.append("text")
    .attr("text-anchor", "middle")
    .text("MJD")
    .attr("fill", "var(--bs-body-color)")
    .classed("axis-label", true);
  var x2label = g.append("text")
    .attr("text-anchor", "middle")
    .text("Pulse number")
    .attr("fill", "var(--bs-body-color)")
    .classed("axis-label", true);
  var ylabel = g.append("text")
    .attr("transform", "rotate(-90)")
    .attr("text-anchor", "middle")
    .text("Residual (phase)")
    .attr("fill", "var(--bs-body-color)")
    .classed("axis-label", true);
  var y2label = g.append("text")
    .attr("transform", "rotate(-90)")
    .attr("text-anchor", "middle")
    .text("Residual (s)")
    .attr("fill", "var(--bs-body-color)")
    .classed("axis-label", true);

  // Add the axes themselves (as groups)
  var xaxis = g.append("g").attr("pointer-events", "bounding-box").attr("cursor", "ew-resize");
  var yaxis = g.append("g").attr("pointer-events", "bounding-box").attr("cursor", "ns-resize");
  var x2axis = g.append("g").attr("pointer-events", "bounding-box").attr("cursor", "ew-resize");
  var y2axis = g.append("g").attr("pointer-events", "bounding-box").attr("cursor", "ns-resize");

  var dashed_line_color = (parentDiv.classList.contains("dark") ? "#aaee1180" : "#55770880");
  // Add a dashed vertical line to mark PEPOCH
  var pepoch_path = g.append("path")
    .style("stroke", dashed_line_color)
    .style("stroke-width", "2")
    .style("stroke-dasharray", "7");

  // Add a dashed horizontal line to mark the zero residual line
  var zero_residual_path = g.append("path")
    .style("stroke", dashed_line_color)
    .style("stroke-width", "2")
    .style("stroke-dasharray", "7");

  // Add a similar dashed line to be used for changing the period
  var period_path = g.append("path")
      .style("stroke", dashed_line_color)
      .style("stroke-width", "2")
      .style("stroke-dasharray", "7")
      .style("stroke-opacity", "0");

  // Construct axes
  var x = d3.scaleLinear();
  var y = d3.scaleLinear();
  var x2 = d3.scaleLinear();
  var y2 = d3.scaleLinear();

  // Return all these new elements in an object
  return {
    parentDiv: parentDiv,
    svg: svg,
    g: g,
    xlabel: xlabel,
    ylabel: ylabel,
    x2label: x2label,
    y2label: y2label,
    xaxis: xaxis,
    yaxis: yaxis,
    x2axis: x2axis,
    y2axis: y2axis,
    x: x,
    y: y,
    x2: x2,
    y2: y2,
    pepoch_path: pepoch_path,
    period_path: period_path,
    zero_residual_path: zero_residual_path
  }
}

function set_residual_plot_dimensions(plot, xlim, ylim, margins, ephemeris) {
  // margins should be of the form:
  //     {top: 50, right: 80, bottom: 50, left: 80}
  // xlim and ylim should be of the form:
  //     [-5, 5]
  // ephemeris should be an object of the form:
  //     {p0: 1000.0, pepoch: 60000.0, dm: 100.0}
  // plot should be an object returned by create_residual_plot_elements()

  // Attach margins amd lims to plot object
  plot.margins = margins;
  plot.xlim = [...xlim]; // Shallow copy
  plot.ylim = [...ylim];

  // Set up graph coordinates
  plot.g.attr("transform", "translate(" + margins.left + "," + margins.top + ")");

  // Calculate and set the width and height of the graph area
  plot.width = plot.parentDiv.clientWidth - margins.left - margins.right;
  plot.height = plot.parentDiv.clientHeight - margins.top - margins.bottom;

  // Set axis limits
  plot.x.domain(xlim).range([0, plot.width]);
  plot.x2
    .domain([calc_pulse_phase(xlim[0], ephemeris).pulse_phase,
             calc_pulse_phase(xlim[1], ephemeris).pulse_phase])
    .range([0, plot.width]);

  plot.y.domain(ylim).range([plot.height, 0]);
  plot.y2
    .domain([ylim[0]*ephemeris.p0, ylim[1]*ephemeris.p0])
    .range([plot.height, 0]);

  // Construct X axes
  plot.xaxis.attr("transform", "translate(0," + plot.height + ")");

  plot.xlabel.attr("x", plot.width/2)
    .attr("y", plot.height + 0.75*margins.bottom)
  plot.x2label.attr("x", plot.width/2)
    .attr("y", -0.75*margins.top)

  plot.xaxis.call(d3.axisBottom(plot.x));
  plot.xaxis.selectAll("text").style("fill", "var(--bs-body-color)");
  plot.xaxis.selectAll("path").style("stroke", "var(--bs-body-color)");
  plot.xaxis.selectAll("line").style("stroke", "var(--bs-body-color)");

  plot.x2axis.call(d3.axisTop(plot.x2));
  plot.x2axis.selectAll("text").style("fill", "var(--bs-body-color)");
  plot.x2axis.selectAll("path").style("stroke", "var(--bs-body-color)");
  plot.x2axis.selectAll("line").style("stroke", "var(--bs-body-color)");

  // Construct Y axis
  plot.y2axis.attr("transform", "translate(" + plot.width + ",0)");

  plot.ylabel.attr("y", -0.5*margins.left)
    .attr("x", -plot.height/2)
  plot.y2label.attr("y", plot.width + 0.75*margins.right)
    .attr("x", -plot.height/2)

  plot.yaxis.call(d3.axisLeft(plot.y));
  plot.yaxis.selectAll("text").style("fill", "var(--bs-body-color)");
  plot.yaxis.selectAll("path").style("stroke", "var(--bs-body-color)");
  plot.yaxis.selectAll("line").style("stroke", "var(--bs-body-color)");

  plot.y2axis.call(d3.axisRight(plot.y2));
  plot.y2axis.selectAll("text").style("fill", "var(--bs-body-color)");
  plot.y2axis.selectAll("path").style("stroke", "var(--bs-body-color)");
  plot.y2axis.selectAll("line").style("stroke", "var(--bs-body-color)");

  // Plot the PEPOCH line
  plot.pepoch_path.attr("d", "M " + plot.x2(0) + ",0 l 0," + plot.height);

  // Plot the zero-residual line
  plot.zero_residual_path.attr("d", "M 0," + plot.y(0) + " l " + plot.width + ",0");
}

function add_residual_data(plot, json_url, appearance, ephemeris, barycentre, label, loader_element) {
  // plot should be an object returned by create_residual_plot_elements()
  // json_url should return a JSON object of the form:
  //     [{mjd: 60001.0, mjd_err: 1e-4, bc_correction: 0.001, freq_MHz: 200.0}, ... ]
  // appearance should be ... TODO
  // ephemeris should be an object of the form:
  //     {p0: 1000.0, pepoch: 60000.0, dm: 100.0, p1: 1e-12}

  if (typeof loader_element !== "undefined") {
    loader_element.style.display = "block";
  }

  d3.json(json_url, function(data) {

    const datapoints = plot.g.append("g")
    plot[label] = {};

    plot[label].err = datapoints
      .selectAll(".err")
      .data(data)
      .enter()
      .append("a")
      .attr("xlink:href", function(datum) { return datum.detail_link; })
      .append("path")
      .style("stroke", appearance.color)
      .style("stroke-width", appearance.stroke_width);

    plot[label].points = datapoints
      .selectAll(".data")
      .data(data)
      .enter()
      .append("a")
      .attr("xlink:href", function(datum) { return datum.detail_link; })
      .append("circle")
      .attr("r", appearance.circle_radius)
      .style("fill", appearance.color);

    plot[label].appearance = appearance;

    position_residual_data(plot, ephemeris, barycentre, label);

    if (typeof loader_element !== "undefined") {
      loader_element.style.display = "none";
    }
  });

}

function add_covariance_err_data(plot, ephemeris, covariance) {
  // plot should be an object returned by create_residual_plot_elements()
  // ephemeris should be an object of the form:
  //     {p0: 1000.0, pepoch: 60000.0, dm: 100.0, p1: 1e-12}
  // Result goes into plot.covariance_err_data

  [xmin, xmax] = plot.xlim;
  dx = (xmax - xmin)/1000.0;
  plot['covariance_err_data'] = Array.from(Array(1000), (_,i) => {
    pulse_phase = calc_pulse_phase(i*dx + xmin, ephemeris).pulse_phase;
    residual = 3 * residual_err(pulse_phase, ephemeris, covariance); // 3σ
    return {
      pulse_phase: pulse_phase,
      residual: residual
    }
  });

  plot['covariance_err'] = plot.g.append("path");
}

function position_covariance_err(plot) {
  const covariance_err_func = d3.area()
    .x(function(d) { return plot.x2(d.pulse_phase) })
    .y1(function(d) { return plot.y2(d.residual) })
    .y0(function(d) { return plot.y2(-d.residual) });

  plot.covariance_err
    .attr('d', covariance_err_func(plot.covariance_err_data))
    .attr('fill', '#f004');
}

function calc_dmdelay(dm, freq_MHz) {
  // Returns answer in units of days
  dmdelay = (4.148808e3 / 86400.0) * dm / (freq_MHz*freq_MHz);
  return dmdelay;
}

function apply_corrections(toa, ephemeris, barycentre) {
  // Calculate MJD of ToA with dedispersion and barycentring corrections
  // (the latter only if requested; 'barycentre' is boolean)
  var adjusted_mjd = toa.mjd - calc_dmdelay(ephemeris.dm, toa.freq_MHz);
  if (barycentre === true) {
    adjusted_mjd += toa.bc_correction;
  }
  return adjusted_mjd;
}

function position_residual_data(plot, ephemeris, barycentre, label) {
  // plot should be an object returned by create_residual_plot_elements()
  // ephemeris should be an object of the form:
  //     {p0: 1000.0, pepoch: 60000.0, dm: 100.0, p1: 1e-12}

  if (plot[label].appearance.display_points === true) {
    plot[label].points.attr("cx", function (toa) {
      return plot.x2(calc_pulse_phase(apply_corrections(toa, ephemeris, barycentre), ephemeris).pulse);
    }).attr("cy", function (toa) {
      return plot.y(calc_pulse_phase(apply_corrections(toa, ephemeris, barycentre), ephemeris).phase);
    }).attr("visibility", "visible");
  } else {
    plot[label].points.attr("visibility", "hidden");
  }

  if (plot[label].appearance.display_err === true) {
    plot[label].err.attr("d", function (toa) {
      const pulse_phase_lo = calc_pulse_phase(apply_corrections(toa, ephemeris, barycentre) - toa.mjd_err, ephemeris);
      const pulse_phase_hi = calc_pulse_phase(apply_corrections(toa, ephemeris, barycentre) + toa.mjd_err, ephemeris);

      const pulse_lo = pulse_phase_lo.pulse;
      const pulse_hi = pulse_phase_hi.pulse;

      const phase_lo = pulse_phase_lo.phase;
      const phase_hi = pulse_phase_hi.phase;

      const xpos_lo = plot.x2(pulse_lo);
      const xpos_hi = plot.x2(pulse_hi);

      const ypos_lo = plot.y(phase_lo);
      const ypos_hi = plot.y(phase_hi);

      var path;

      if (phase_lo < phase_hi) {
        // The errors do not straddle the phase = 0.5 boundary between adjacent pulses,
        // so draw a single line
        path = " M " + xpos_lo + "," + ypos_lo +
          " L " + xpos_hi + "," + ypos_hi;
      } else {
        // The errors *do* straddle the phase = 0.5 boundary between adjacent pulses,
        // so draw two line segments
        path = " M " + xpos_lo + "," + ypos_lo +
          " L " + xpos_lo + "," + plot.y(0.5) +
          " M " + xpos_hi + "," + plot.y(-0.5) +
          " L " + xpos_hi + "," + ypos_hi;
      }

      return path;
    }).attr("visibility", "visible");
  } else {
    plot[label].err.attr("visibility", "hidden");
  }
}

// A handy function for plotting the period path
function plot_period_line(plot, pos) {

  // Draw the svg line from the origin to the clicked point
  // Everything here in "g" coords
  let xpos = pos[0] - plot.margins.left; // (Now converted to "g" coords)
  let ypos = pos[1] - plot.margins.top;
  let xorig = plot.x2(0);
  let yorig = plot.y(0);
  let slope = (ypos - yorig) / (xpos - xorig);
  let yintercept = ypos - slope*xpos;

  // Draw the path connecting the clicked point and the origin
  plot.period_path
    .attr("d", "M 0," + yintercept + " L " + plot.width + "," + (slope*plot.width + yintercept))
    .style("stroke-opacity", "0.5");

}

/*********************\
* Zooming with scroll *
**********************/

function plot_zoom(plot, pos, ephemeris, barycentre) {
  xpos = plot.x.invert(pos[0] - plot.margins.left);
  ypos = plot.y.invert(pos[1] - plot.margins.top);
  [xmin, xmax] = plot.xlim;
  [ymin, ymax] = plot.ylim;

  const scale_factor = d3.event.wheelDelta < 0 ? 0.05 : -0.05;

  if (xpos > xmin && xpos < xmax) { // if mouse is not over the y-axis itself (if it were, we leave x-zoom untouched)
    xmin -= scale_factor*(xpos - xmin);
    xmax += scale_factor*(xmax - xpos);
  }

  if (ypos > ymin && ypos < ymax) { // if mouse is not over the x-axis itself (if it were, we leave y-zoom untouched)
    ymin -= scale_factor*(ypos - ymin);
    ymax += scale_factor*(ymax - ypos);
  }

  set_residual_plot_dimensions(plot, [xmin, xmax], [ymin, ymax], plot.margins, ephemeris);
}

/*****************************************************************\
*  Mouse interaction functions for editing the period and/or axes *
\*****************************************************************/

// Global state
svg_mousedown = false;
axis_mousedown = false;
period_ref_factor = null;
prev_pos = {x: 0, y: 0}; // Dummy values. Set on mousedown

function edit_period_mousedown() {
  svg_mousedown = true;
  period_ref_factor = null;
}

function edit_axis_mousedown(pos) {
  d3.event.stopPropagation();
  axis_mousedown = true;
  prev_pos.x = plot.x.invert(pos[0] - plot.margins.left);
  prev_pos.y = plot.y.invert(pos[1] - plot.margins.top);
}

function edit_period_mousemove(plot, ephemeris, pos, barycentre, label) {

  var pepoch = pepoch_element.value;
  var p0 = p0_element.value;

  // Record for later use the world coords of the clicked event, as well as the p0
  let mjd = plot.x.invert(pos[0] - plot.margins.left) - ephemeris.pepoch; // world coords of click pos
  let res = plot.y2.invert(pos[1] - plot.margins.top);          // world coords of click pos
  let factor = res / (mjd*86400);

  if (period_ref_factor !== null) {
    ephemeris.p0 *= (1 - (factor - period_ref_factor));

    // Update plot
    position_residual_data(plot, ephemeris, barycentre, label);
  }

  period_ref_factor = factor;

  plot_period_line(plot, pos);
}

function edit_axis_mousemove(plot, pos, axis) {
  if (axis_mousedown === true) {
    d3.event.stopPropagation();
    if (axis == "x") {
      xpos = plot.x.invert(pos[0] - plot.margins.left);
      xdiff = xpos - prev_pos.x;
      plot.xlim[0] -= xdiff;
      plot.xlim[1] -= xdiff;
    } else if (axis == "y") {
      ypos = plot.y.invert(pos[1] - plot.margins.top);
      ydiff = ypos - prev_pos.y;
      plot.ylim[0] -= ydiff;
      plot.ylim[1] -= ydiff;
    }
  }
}

function edit_mouseup(plot) {

  if (svg_mousedown === true) {
    svg_mousedown = false;

    plot.period_path
      .style("stroke-opacity", "0.0");
  }

  if (axis_mousedown === true) {
    axis_mousedown = false;
  }
}

