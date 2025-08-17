from django import template
from astropy.time import Time
import astropy.units as u

register = template.Library()

@register.filter
def err(obj, arg):
    return obj.err(arg)

@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)

@register.filter
def toa_format(obj, format):
    if not hasattr(obj, "raw_mjd"):
        return None
    t = Time(float(obj.raw_mjd), scale='utc', format='mjd')
    #print(f'{t.to_value(format) = }')
    return t.to_value(format)

@register.filter
def convert_unit(obj, arg_string):
    from_unit, to_unit = arg_string.split(',')
    return str((float(obj) * u.Unit(from_unit)).to(to_unit).value)
