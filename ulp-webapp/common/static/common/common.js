function update_permissions(url, csrf_token, app, model, pk, group_or_user, name, permission_type, checkbox) {

  const xhr = new XMLHttpRequest();
  xhr.open('PUT', url);
  xhr.setRequestHeader("X-CSRFToken", csrf_token);
  xhr.setRequestHeader("Content-type", "application/json");

  const data = {
    app: app,
    model: model,
    pk: pk,
    group_or_user: group_or_user,
    name: name,
    permission_type: permission_type,
    permit: checkbox.checked
  }

  xhr.onload = () => {
    if (xhr.status != 200) {
      checkbox.style.backgroundColor = "#600";
      alert(xhr.responseText);
      return;
    }

    checkbox.style.backgroundColor = "field";
  }

  xhr.send(JSON.stringify(data));
  checkbox.style.backgroundColor = "#880";
}

/***********************
 * Fold:
 * Convert timestamps to pulse numbers and phases given ephemeris information
 *
 * Inputs:
 * - mjds (array of floats): The times (in MJD) to be converted
 * - pepoch (float): The reference timestamp (MJD) when pulse = phase = 0
 * - period (P0): The folding period (in seconds)
 *
 * Phase is defined to go from -0.5 to 0.5. The resulting phase and pulse must
 * satisfy:
 *   mjds = (pulses + phase)*period/86400 + pepoch
 ***********************/
function fold(mjds, pepoch, period) {
  //console.log("mjds: ", mjds);
  //console.log("pepoch: ", pepoch);
  //console.log("period: ", period);
  let pulse_phases = mjds.map((mjd) => (mjd - pepoch)/(period/86400.0));
  let pulses = pulse_phases.map((pulse_phase) => Math.floor(pulse_phase + 0.5));
  let phases = [];
  for (let i = 0; i < mjds.length; i++) {
    phases[i] = pulse_phases[i] - pulses[i];
  }

  return {pulses: pulses, phases: phases};
}

function get_plotly_default_layout() {
  return {
    showlegend: true,
    font: {color: "rgb(222,226, 230)"},
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    xaxis: {tickformat: "f"},
    hoverlabel: {font: {color: 'black'}}
  };
}

function on_value_change(element) {
  element.style.backgroundColor = "#880";
}

function freq_color(freq, min_freq, max_freq) {
  // Also see this very nice tool to build a custom colormap:
  //   https://leonardocolor.io/scales.html

  // However, I actually want a function to return a specific value, and
  // interpolating a colormap manually is too hard.
  min_h = 301;
  max_h = 235;
  // Map frequency in range to (0,1)
  f = (freq - min_freq) / (max_freq - min_freq);
  // ...then map to desired hue space
  h = f*(max_h - min_h) + min_h;
  return 'hsl(' + h + ',1,0.65)';
}
