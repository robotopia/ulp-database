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

