from django.apps import apps
from django.core.serializers.json import Serializer as JSONSerializer
from django.core.serializers.python import Deserializer as PythonDeserializer
from django.core.serializers.base import DeserializationError

import json
import itertools

class WorkingEphemerisCovarianceSerializer(JSONSerializer):

    def get_dump_object(self, obj):
        dumped_object = super().get_dump_object(obj)
        model_fields = ['pepoch', 'p0', 'p1', 'dm']  # This should be defined somehow in models
        N = len(model_fields)

        matrix = []
        for i in range(N):
            row = []
            for j in range(N):
                field = f'{model_fields[i]}_{model_fields[j]}' if i <= j else f'{model_fields[j]}_{model_fields[i]}'
                row.append(dumped_object["fields"][field] or 0)
            matrix.append(row)

        dumped_object['fields'] = model_fields
        dumped_object['matrix'] = matrix

        return dumped_object


def WorkingEphemerisCovarianceDeserializer(stream_or_string, **options):
    """Deserialize a stream or string of JSON data."""
    if not isinstance(stream_or_string, (bytes, str)):
        stream_or_string = stream_or_string.read()
    if isinstance(stream_or_string, bytes):
        stream_or_string = stream_or_string.decode()

    def transform_object(obj):

        model_fields  = obj.pop('fields')
        matrix        = obj.pop('matrix')

        N = len(model_fields)
        obj['fields'] = {
            f'{model_fields[i]}_{model_fields[j]}': matrix[i][j]
            for i, j in itertools.combinations_with_replacement(range(N), 2)
        }
        return obj

    try:
        objs = [transform_object(obj) for obj in json.loads(stream_or_string)]
        yield from PythonDeserializer(objs, **options)
    except (GeneratorExit, DeserializationError):
        raise
    except Exception as exc:
        raise DeserializationError() from exc

