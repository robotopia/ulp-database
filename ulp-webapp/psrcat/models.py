# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey and OneToOneField has `on_delete` set to the desired behavior
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from django.db import models


class Ancillary(models.Model):
    ancillary_id = models.AutoField(primary_key=True)
    parameter = models.ForeignKey('Parameter', models.DO_NOTHING, blank=True, null=True)
    value = models.TextField(blank=True, null=True)
    description = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'ancillary'


class Association(models.Model):
    association_id = models.AutoField(primary_key=True)
    pulsar = models.ForeignKey('Pulsar', models.DO_NOTHING, blank=True, null=True)
    associationtype = models.ForeignKey('Associationtype', models.DO_NOTHING, db_column='associationType_id', blank=True, null=True)  # Field name made lowercase.
    citation = models.ForeignKey('Citation', models.DO_NOTHING, blank=True, null=True)
    confidence = models.IntegerField(blank=True, null=True)
    name = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'association'


class Associationtype(models.Model):
    associationtype_id = models.AutoField(db_column='associationType_id', primary_key=True)  # Field name made lowercase.
    label = models.TextField(blank=True, null=True)
    description = models.TextField(blank=True, null=True)
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'associationType'


class Binarytype(models.Model):
    binarytype_id = models.AutoField(db_column='binaryType_id', primary_key=True)  # Field name made lowercase.
    pulsar = models.ForeignKey('Pulsar', models.DO_NOTHING, blank=True, null=True)
    binarytypeoptions = models.ForeignKey('Binarytypeoptions', models.DO_NOTHING, db_column='binaryTypeOptions_id', blank=True, null=True)  # Field name made lowercase.
    confidence = models.IntegerField(blank=True, null=True)
    citation = models.ForeignKey('Citation', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'binaryType'


class Binarytypeoptions(models.Model):
    binarytypeoptions_id = models.AutoField(db_column='binaryTypeOptions_id', primary_key=True)  # Field name made lowercase.
    label = models.TextField(blank=True, null=True)
    description = models.TextField(blank=True, null=True)
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'binaryTypeOptions'


class Citation(models.Model):
    citation_id = models.AutoField(primary_key=True)
    v1label = models.TextField(blank=True, null=True)
    label = models.TextField(blank=True, null=True)
    title = models.TextField(blank=True, null=True)
    author = models.TextField(blank=True, null=True)
    journal = models.TextField(blank=True, null=True)
    year = models.TextField(blank=True, null=True)
    month = models.TextField(blank=True, null=True)
    volume = models.TextField(blank=True, null=True)
    number = models.TextField(blank=True, null=True)
    pages = models.TextField(blank=True, null=True)
    doi = models.TextField(blank=True, null=True)
    url = models.TextField(blank=True, null=True)
    version = models.ForeignKey('Version', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'citation'


class Derived(models.Model):
    derived_id = models.AutoField(primary_key=True)
    parameter = models.ForeignKey('Parameter', models.DO_NOTHING, blank=True, null=True)
    method = models.TextField(blank=True, null=True)
    methodversion = models.TextField(db_column='methodVersion', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'derived'


class Derivedfromparameter(models.Model):
    derivedfromparameter_id = models.AutoField(db_column='derivedFromParameter_id', primary_key=True)  # Field name made lowercase.
    derived = models.ForeignKey(Derived, models.DO_NOTHING, blank=True, null=True)
    parameter = models.ForeignKey('Parameter', models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'derivedFromParameter'


class Distance(models.Model):
    distance_id = models.AutoField(primary_key=True)
    pulsar = models.ForeignKey('Pulsar', models.DO_NOTHING, blank=True, null=True)
    citation = models.ForeignKey(Citation, models.DO_NOTHING, blank=True, null=True)
    value = models.TextField(blank=True, null=True)
    uncertainty = models.TextField(blank=True, null=True)
    label = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'distance'


class Fitparameters(models.Model):
    fitparameters_id = models.AutoField(db_column='fitParameters_id', primary_key=True)  # Field name made lowercase.
    units = models.TextField(blank=True, null=True)
    ephemeris = models.TextField(blank=True, null=True)
    clock = models.TextField(blank=True, null=True)
    citation_id = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'fitParameters'


class Linkedset(models.Model):
    linkedset_id = models.AutoField(db_column='linkedSet_id', primary_key=True)  # Field name made lowercase.
    citation = models.ForeignKey(Citation, models.DO_NOTHING, blank=True, null=True)
    description = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'linkedSet'


class Name(models.Model):
    name_id = models.AutoField(primary_key=True)
    pulsar = models.ForeignKey('Pulsar', models.DO_NOTHING, blank=True, null=True)
    citation = models.ForeignKey(Citation, models.DO_NOTHING, blank=True, null=True)
    name = models.TextField(blank=True, null=True)
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'name'


class Observingsystem(models.Model):
    observingsystem_id = models.AutoField(db_column='observingSystem_id', primary_key=True)  # Field name made lowercase.
    systemlabel = models.TextField(db_column='systemLabel', blank=True, null=True)  # Field name made lowercase.
    centralfrequency = models.FloatField(db_column='centralFrequency', blank=True, null=True)  # Field name made lowercase.
    bandwidth = models.FloatField(blank=True, null=True)
    telescope = models.TextField(blank=True, null=True)
    approximate = models.IntegerField(blank=True, null=True)
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'observingSystem'


class Parameter(models.Model):
    parameter_id = models.AutoField(primary_key=True)
    pulsar = models.ForeignKey('Pulsar', models.DO_NOTHING, blank=True, null=True)
    citation = models.ForeignKey(Citation, models.DO_NOTHING, blank=True, null=True)
    linkedset = models.ForeignKey(Linkedset, models.DO_NOTHING, db_column='linkedSet_id', blank=True, null=True)  # Field name made lowercase.
    fitparameters = models.ForeignKey(Fitparameters, models.DO_NOTHING, db_column='fitParameters_id', blank=True, null=True)  # Field name made lowercase.
    observingsystem = models.ForeignKey(Observingsystem, models.DO_NOTHING, db_column='observingSystem_id', blank=True, null=True)  # Field name made lowercase.
    parametertype = models.ForeignKey('Parametertype', models.DO_NOTHING, db_column='parameterType_id', blank=True, null=True)  # Field name made lowercase.
    timederivative = models.IntegerField(db_column='timeDerivative', blank=True, null=True)  # Field name made lowercase.
    companionnumber = models.IntegerField(db_column='companionNumber', blank=True, null=True)  # Field name made lowercase.
    value = models.TextField(blank=True, null=True)
    uncertainty = models.TextField(blank=True, null=True)
    referencetime = models.TextField(db_column='referenceTime', blank=True, null=True)  # Field name made lowercase.
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'parameter'


class Parametertype(models.Model):
    parametertype_id = models.AutoField(db_column='parameterType_id', primary_key=True)  # Field name made lowercase.
    label = models.TextField(blank=True, null=True)
    unit = models.TextField(blank=True, null=True)
    description = models.TextField(blank=True, null=True)
    timingflag = models.IntegerField(db_column='timingFlag', blank=True, null=True)  # Field name made lowercase.
    datatype = models.TextField(db_column='dataType', blank=True, null=True)  # Field name made lowercase.
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'parameterType'


class Pulsar(models.Model):
    pulsar_id = models.AutoField(primary_key=True)
    jname = models.TextField(blank=True, null=True)
    survey = models.ForeignKey('Survey', models.DO_NOTHING, blank=True, null=True)
    citation = models.ForeignKey(Citation, models.DO_NOTHING, blank=True, null=True)
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'pulsar'


class Sourcetype(models.Model):
    sourcetype_id = models.AutoField(db_column='sourceType_id', primary_key=True)  # Field name made lowercase.
    pulsar = models.ForeignKey(Pulsar, models.DO_NOTHING, blank=True, null=True)
    sourcetypeoptions = models.ForeignKey('Sourcetypeoptions', models.DO_NOTHING, db_column='sourceTypeOptions_id', blank=True, null=True)  # Field name made lowercase.
    confidence = models.IntegerField(blank=True, null=True)
    citation = models.ForeignKey(Citation, models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'sourceType'


class Sourcetypeoptions(models.Model):
    sourcetypeoptions_id = models.AutoField(db_column='sourceTypeOptions_id', primary_key=True)  # Field name made lowercase.
    label = models.TextField(blank=True, null=True)
    description = models.TextField(blank=True, null=True)
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'sourceTypeOptions'


class Survey(models.Model):
    survey_id = models.AutoField(primary_key=True)
    label = models.TextField(blank=True, null=True)
    shortlabel = models.TextField(db_column='shortLabel', blank=True, null=True)  # Field name made lowercase.
    telescope = models.TextField(blank=True, null=True)
    receiver = models.TextField(blank=True, null=True)
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'survey'


class Surveytopulsar(models.Model):
    surveytopulsar_id = models.AutoField(db_column='surveyToPulsar_id', primary_key=True)  # Field name made lowercase.
    survey = models.ForeignKey(Survey, models.DO_NOTHING, blank=True, null=True)
    pulsar = models.ForeignKey(Pulsar, models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'surveyToPulsar'


class Tag(models.Model):
    tag_id = models.AutoField(primary_key=True)
    taglabel = models.TextField(db_column='tagLabel', blank=True, null=True)  # Field name made lowercase.
    tagstring = models.TextField(db_column='tagString', blank=True, null=True)  # Field name made lowercase.
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'tag'


class Tagtocitation(models.Model):
    tagtocitation_id = models.AutoField(db_column='tagToCitation_id', primary_key=True)  # Field name made lowercase.
    citation = models.ForeignKey(Citation, models.DO_NOTHING, blank=True, null=True)
    tag = models.ForeignKey(Tag, models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'tagToCitation'


class Tagtolinkedset(models.Model):
    tagtolinkedset_id = models.AutoField(db_column='tagToLinkedSet_id', primary_key=True)  # Field name made lowercase.
    linkedset = models.ForeignKey(Linkedset, models.DO_NOTHING, db_column='linkedSet_id', blank=True, null=True)  # Field name made lowercase.
    tag = models.ForeignKey(Tag, models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'tagToLinkedSet'


class Tagtopulsar(models.Model):
    tagtopulsar_id = models.AutoField(db_column='tagToPulsar_id', primary_key=True)  # Field name made lowercase.
    pulsar = models.ForeignKey(Pulsar, models.DO_NOTHING, blank=True, null=True)
    tag = models.ForeignKey(Tag, models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'tagToPulsar'


class Version(models.Model):
    version_id = models.AutoField(primary_key=True)
    version = models.TextField(blank=True, null=True)
    entrydate = models.TextField(db_column='entryDate', blank=True, null=True)  # Field name made lowercase.
    notes = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'version'
