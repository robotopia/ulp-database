PermissionFieldset = (
    "Permissions", {
        "description": "Who owns, can see, and can edit this instance.",
        "classes": ["collapse"],
        "fields": ["owner", ("can_view_groups", "can_edit_groups"), ("can_view_users", "can_edit_users")],
    },
)
