from django import template

register = template.Library()

@register.filter
def err(obj, arg):
    return obj.err(arg)

@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)

