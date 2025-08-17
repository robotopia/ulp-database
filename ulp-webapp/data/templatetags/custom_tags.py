from django import template

register = template.Library()

@register.filter
def err(obj, arg):
    return obj.err(arg)

@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)

@register.filter
def toa_format(obj, format):
    if not obj.raw_mjd:
        return None
    return Time(float(obj.raw_mjd), scale='utc', format='mjd').to_value(format)
