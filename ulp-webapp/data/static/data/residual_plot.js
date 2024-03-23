// Function to calculate the pulse numbers and phases of TOAs
function calc_pulse_phase(mjd, pepoch, P) {
  // mjd and pepoch should be in days,
  // P (the folding period) in seconds
  var pulse_phase = 86400*(mjd - pepoch)/P; // Dimensionless
  var pulse = Math.floor(pulse_phase)
  var phase = (pulse_phase + 0.5) % 1 - 0.5;
  return {pulse: pulse, phase: phase, pulse_phase: pulse_phase};
}

// Function to create the plot elements
function create_residual_plot_elements(parentDivId) {
  // parentDivId should be a string

  var svg = d3.select("#" + parentDivId).append("svg");

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

  // Return all these new elements in an object
  return {
    svg: svg,
    g: g,
    xlabel: xlabel,
    ylabel: ylabel,
    x2label: x2label,
    y2label: y2label,
    xaxis: xaxis,
    yaxis: yaxis,
    x2axis: x2axis,
    y2axis: y2axis
  }
}

