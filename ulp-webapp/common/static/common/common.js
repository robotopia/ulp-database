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
  let pulse_phases = mjds.map((mjd) => (mjd - pepoch)/(period/86400.0));
  let pulses = pulse_phases.map((pulse_phase) => Math.floor(pulse_phase + 0.5));
  let phases = [];
  for (let i = 0; i < mjds.length; i++) {
    phases[i] = pulse_phases[i] - pulses[i] - 0.5;
  }

  return {pulses: pulses, phases: phases}
}

