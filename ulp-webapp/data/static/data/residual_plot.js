// Function to calculate the pulse numbers and phases of TOAs
function calc_pulse_phase(mjd, pepoch, P) {
  // mjd and pepoch should be in days,
  // P (the folding period) in seconds
  var pulse_phase = 86400*(mjd - pepoch)/P; // Dimensionless
  var pulse = Math.round(pulse_phase)
  var phase = (pulse_phase + 0.5) % 1 - 0.5;
  return {pulse: pulse, phase: phase, pulse_phase: pulse_phase};
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
    .attr("fill", "var(--bs-body-color)");
  var x2label = g.append("text")
    .attr("text-anchor", "middle")
    .text("Pulse number")
    .attr("fill", "var(--bs-body-color)");
  var ylabel = g.append("text")
    .attr("transform", "rotate(-90)")
    .attr("text-anchor", "middle")
    .text("Residual (phase)")
    .attr("fill", "var(--bs-body-color)");
  var y2label = g.append("text")
    .attr("transform", "rotate(-90)")
    .attr("text-anchor", "middle")
    .text("Residual (s)")
    .attr("fill", "var(--bs-body-color)");

  // Add the axes themselves (as groups)
  var xaxis = g.append("g");
  var yaxis = g.append("g");
  var x2axis = g.append("g");
  var y2axis = g.append("g");

  // Add a dashed vertical line to mark PEPOCH
  var pepoch_path = g.append("path")
      .style("stroke", "#aaee1180")
      .style("stroke-width", "2")
      .style("stroke-dasharray", "7");

  // Add a similar dashed line to be used for changing the period
  var period_path = g.append("path")
      .style("stroke", "#aaee1180")
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
    period_path: period_path
  }
}

function set_residual_plot_dimensions(plot, xlim, ylim, margins, ephemeris) {
  // margins should be of the form:
  //     {top: 50, right: 80, bottom: 50, left: 80}
  // xlim and ylim should be of the form:
  //     [-5, 5]
  // ephemeris should be an object of the form:
  //     {folding_period: 1000.0, pepoch: 60000.0}
  // plot should be an object returned by create_residual_plot_elements()

  // Attach margins amd lims to plot object
  plot.margins = margins;
  plot.xlim = xlim;
  plot.ylim = ylim;

  // Set up graph coordinates
  plot.g.attr("transform", "translate(" + margins.left + "," + margins.top + ")");

  // Calculate and set the width and height of the graph area
  plot.width = plot.parentDiv.clientWidth - margins.left - margins.right;
  plot.height = plot.parentDiv.clientHeight - margins.top - margins.bottom;

  // Set axis limits
  plot.x.domain(xlim).range([0, plot.width]);
  plot.x2
    .domain([calc_pulse_phase(xlim[0], ephemeris.pepoch, ephemeris.folding_period).pulse_phase,
             calc_pulse_phase(xlim[1], ephemeris.pepoch, ephemeris.folding_period).pulse_phase])
    .range([0, plot.width]);

  plot.y.domain(ylim).range([plot.height, 0]);
  plot.y2
    .domain([ylim[0]*ephemeris.folding_period, ylim[1]*ephemeris.folding_period])
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
}

function add_residual_data(plot, toa_url, color, ephemeris) {
  // plot should be an object returned by create_residual_plot_elements()
  // toas should be an array of objects of the form:
  //     [{'mjd': 60001.0, 'mjd_err': 1e-4}, ... ]
  // color can be any string representing a valid color

  d3.json(toa_url, function(toas) {

    const datapoints = plot.g.append("g")

    plot.toa_err_points = datapoints
      .selectAll(".err")
      .data(toas)
      .enter()
      .append("path")
      .style("stroke", color)
      .style("stroke-width", "2");

    plot.toa_points = datapoints
      .selectAll(".data")
      .data(toas)
      .enter()
      .append("circle")
      .attr("r", 3)
      .style("fill", color);

    position_residual_data(plot, ephemeris);
  });
}

function position_residual_data(plot, ephemeris) {
  // plot should be an object returned by create_residual_plot_elements()
  // ephemeris should be an object of the form:
  //     {folding_period: 1000.0, pepoch: 60000.0}

  plot.toa_points
    .attr("cx", function (toa) { return plot.x2(calc_pulse_phase(toa.mjd, ephemeris.pepoch, ephemeris.folding_period).pulse); })
    .attr("cy", function (toa) { return plot.y(calc_pulse_phase(toa.mjd, ephemeris.pepoch, ephemeris.folding_period).phase); })

  plot.toa_err_points
    .attr("d", function (toa) {
      pulse_phase_lo = calc_pulse_phase(toa.mjd - toa.mjd_err, ephemeris.pepoch, ephemeris.folding_period);
      pulse_phase_hi = calc_pulse_phase(toa.mjd + toa.mjd_err, ephemeris.pepoch, ephemeris.folding_period);

      pulse = pulse_phase_lo.pulse;
      phase_lo = pulse_phase_lo.phase;
      phase_hi = pulse_phase_hi.phase;

      let xpos = plot.x2(pulse);

      return " M " + xpos + "," + plot.y(phase_lo) +
             " L " + xpos + "," + plot.y(phase_hi);
    })

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

function plot_zoom(plot, pos, ephemeris) {
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
  position_residual_data(plot, ephemeris);
}

/******************************************************\
*  Mouse interaction functions for editing the period  *
\******************************************************/

// Global state
svg_mousedown = false;
period_ref_factor = null;

function edit_period_mousedown() {
  svg_mousedown = true;
  period_ref_factor = null;
}

function edit_period_mousemove(plot, ephemeris, pos) {

  var pepoch = pepoch_element.value;
  var folding_period = folding_period_element.value;

  // Record for later use the world coords of the clicked event, as well as the folding period
  let mjd = plot.x.invert(pos[0] - plot.margins.left) - ephemeris.pepoch; // world coords of click pos
  let res = plot.y2.invert(pos[1] - plot.margins.top);          // world coords of click pos
  let factor = res / (mjd*86400);

  if (period_ref_factor !== null) {
    ephemeris.folding_period *= (1 - (factor - period_ref_factor));

    // Update plot
    position_residual_data(plot, ephemeris);
  }

  period_ref_factor = factor;

  plot_period_line(plot, pos);
}

function edit_period_mouseup(plot) {

  svg_mousedown = false;

  plot.period_path
    .style("stroke-opacity", "0.0");
}
