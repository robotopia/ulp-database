# Generated by Django 5.0.2 on 2024-07-03 02:38

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('data', '0007_ephemerismeasurement'),
    ]

    operations = [
        migrations.AddField(
            model_name='timeofarrival',
            name='telescope',
            field=models.CharField(blank=True, choices=[('0', 'ALMA'), ('1', 'AO'), ('2', 'ARCA'), ('3', 'ASKAP'), ('4', 'ATA'), ('5', 'ATST'), ('6', 'Allen Telescope Array'), ('7', 'Anderson Mesa'), ('8', 'Anglo-Australian Observatory'), ('9', 'Apache Point'), ('10', 'Apache Point Observatory'), ('11', 'Arecibo'), ('12', 'Arecibo Observatory'), ('13', 'Astroparticle Research with Cosmics in the Abyss'), ('14', 'Atacama Large Millimeter Array'), ('15', 'Australian Square Kilometre Array Pathfinder'), ('16', 'BAO'), ('17', 'BBSO'), ('18', 'Beijing XingLong Observatory'), ('19', 'Big Bear Solar Observatory'), ('20', 'Black Moshannon Observatory'), ('21', 'CAHA'), ('22', 'CAHA'), ('23', 'CHARA'), ('24', 'CHIME'), ('25', 'Canada-France-Hawaii Telescope'), ('26', 'Canadian Hydrogen Intensity Mapping Experiment'), ('27', 'Catalina Observatory'), ('28', 'Catalina Observatory: 61 inch telescope'), ('29', 'Centro Astronomico Hispano-Aleman, Almeria'), ('30', 'Cerro Armazones Observatory'), ('31', 'Cerro Pachon'), ('32', 'Cerro Paranal'), ('33', 'Cerro Tololo'), ('34', 'Cerro Tololo Interamerican Observatory'), ('35', 'Cima Ekar 182 cm Telescope'), ('36', 'Cima Ekar Observing Station'), ('37', 'Ckoirama'), ('38', 'Ckoirama Observatory'), ('39', 'DCT'), ('40', 'DKIST'), ('41', 'DRAO'), ('42', 'DRAO 26m Telescope'), ('43', 'Daniel K. Inouye Solar Telescope'), ('44', 'Discovery Channel Telescope'), ('45', 'Dominion Astrophysical Observatory'), ('46', 'Dominion Radio Astrophysical Observatory'), ('47', 'Effelsberg'), ('48', 'Effelsberg 100-m Radio Telescope'), ('49', 'FAST'), ('50', 'Five-hundred-meter Aperture Spherical radio Telescope'), ('51', 'G1'), ('52', 'GBT'), ('53', 'GEO'), ('54', 'GEO600 Gravitational Wave Detector'), ('55', 'GEO_600'), ('56', 'GMRT'), ('57', 'Gemini North'), ('58', 'Gemini South'), ('59', 'Giant Metrewave Radio Telescope'), ('60', 'Green Bank Observatory'), ('61', 'Green Bank Telescope'), ('62', 'H1'), ('63', 'HALO'), ('64', 'HET'), ('65', 'Hale Telescope'), ('66', 'Haleakala Observatories'), ('67', 'Happy Jack'), ('68', 'Hat Creek'), ('69', 'Hat Creek Radio Observatory'), ('70', 'Helium And Lead Observatory'), ('71', 'Hobby Eberly Telescope'), ('72', 'Hyper-Kamiokande'), ('73', 'HyperK'), ('74', 'Hyperk'), ('75', 'IAO'), ('76', 'ICECUBE'), ('77', 'IceCube'), ('78', 'IceCube Neutrino Observatory'), ('79', 'Indian Astronomical Observatory'), ('80', 'JCMT'), ('81', 'James Clerk Maxwell Telescope'), ('82', 'Jansky Very Large Array'), ('83', 'John Galt Telescope'), ('84', 'K1'), ('85', 'KAGRA'), ('86', 'KM3NeT ARCA'), ('87', 'KM3NeT ORCA'), ('88', 'KM3NeT arca'), ('89', 'KM3NeT orca'), ('90', 'Kamioka Gravitational Wave Detector'), ('91', 'Keck Observatory'), ('92', 'Kitt Peak'), ('93', 'Kitt Peak National Observatory'), ('94', 'L1'), ('95', 'LDT'), ('96', 'LHO'), ('97', 'LHO_4k'), ('98', 'LIGO Hanford Observatory'), ('99', 'LIGO Livingston Observatory'), ('100', 'LLO'), ('101', 'LLO_4k'), ('102', 'LOFAR'), ('103', 'LSST'), ('104', 'LSST 1.4m'), ('105', 'LSST 8.4m'), ('106', 'LSST AuxTel'), ('107', 'LWA1'), ('108', 'La Silla Observatory'), ('109', 'La Silla Observatory (ESO)'), ('110', 'Large Binocular Telescope'), ('111', 'Las Campanas Observatory'), ('112', 'Lick Observatory'), ('113', 'Long Wavelength Array 1'), ('114', 'Low-Frequency Array'), ('115', 'Lowell Discovery Telescope'), ('116', 'Lowell Observatory'), ('117', 'Lowell Observatory - Anderson Mesa'), ('118', 'Lowell Observatory - Mars Hill'), ('119', 'MDM Observatory'), ('120', 'MEERKAT'), ('121', 'MJO'), ('122', 'MOA'), ('123', 'MWA'), ('124', 'Manastash Ridge Observatory'), ('125', 'Mars Hill'), ('126', 'McDonald Observatory'), ('127', 'Medicina'), ('128', 'Medicina Dish'), ('129', 'Medicina Radio Telescope'), ('130', 'MeerKAT'), ('131', 'Michigan-Dartmouth-MIT Observatory'), ('132', 'Mont Mégantic Observatory'), ('133', 'Mount Graham International Observatory'), ('134', 'Mount Wilson Observatory'), ('135', 'Mt Graham'), ('136', 'Mt John'), ('137', 'Mt. Ekar 182 cm Telescope'), ('138', 'Mt. Stromlo Observatory'), ('139', 'Multiple Mirror Telescope'), ('140', 'Murchison Widefield Array'), ('141', 'Murriyang'), ('142', 'NANCAY'), ('143', 'NASA Infrared Telescope Facility'), ('144', 'NOV'), ('145', 'NOVA'), ('146', 'NOvA'), ('147', 'NPOI'), ('148', 'NST'), ('149', 'Nancay'), ('150', 'Nancay Radio Telescope'), ('151', 'National Observatory of Venezuela'), ('152', 'Navy Precision Optical Interferometer'), ('153', 'Noto'), ('154', 'Noto Radio Telescope'), ('155', 'NuMI Off-axis νe Appearance'), ('156', 'OAJ'), ('157', 'OAJ'), ('158', 'OAO'), ('159', 'OAO'), ('160', 'OARMA'), ('161', 'OARMA'), ('162', 'OCA'), ('163', 'OMM'), ('164', 'ORCA'), ('165', 'OT'), ('166', 'Observatoire SIRENE'), ('167', 'Observatoire de Haute Provence'), ('168', 'Observatoire du Mont Mégantic'), ('169', 'Observatorio Astrofisico de Javalambre'), ('170', 'Observatorio Astronomico Nacional, San Pedro Martir'), ('171', 'Observatorio Astronomico Nacional, Tonantzintla'), ('172', 'Observatorio Astronomico Ramon Maria Aller, Santiago de Compostela'), ('173', 'Observatorio Cerro Armazones'), ('174', 'Observatorio Ckoirama'), ('175', 'Observatorio Ramon Maria Aller'), ('176', 'Observatorio de Calar Alto'), ('177', 'Observatorio del Teide'), ('178', 'Observatorio del Teide, Tenerife'), ('179', 'Okayama Astrophysical Observatory'), ('180', 'Oscillation Research with Cosmics in the Abyss'), ('181', 'Otehiwai'), ('182', 'Otehiwai Observatory'), ('183', 'Owens Valley Radio Observatory'), ('184', 'PTO'), ('185', 'Palomar'), ('186', 'Paranal Observatory'), ('187', 'Paranal Observatory (ESO)'), ('188', 'Parkes'), ('189', 'Perkins'), ('190', 'Roque de los Muchachos'), ('191', 'Roque de los Muchachos, La Palma'), ('192', 'Royal Observatory Greenwich'), ('193', 'Rubin'), ('194', 'Rubin AuxTel'), ('195', 'Rubin Observatory'), ('196', 'SAAO'), ('197', 'SALT'), ('198', 'SNO+'), ('199', 'SPO'), ('200', 'SRT'), ('201', 'Sac Peak'), ('202', 'Sacramento Peak'), ('203', 'Sacramento Peak Observatory'), ('204', 'Sardinia Radio Telescope'), ('205', 'Siding Spring Observatory'), ('206', 'Southern African Large Telescope'), ('207', 'Subaru'), ('208', 'Subaru Telescope'), ('209', 'Sudbury Neutrino Observatory +'), ('210', 'Sunspot'), ('211', 'Super-Kamiokande'), ('212', 'SuperK'), ('213', 'Superk'), ('214', 'Sutherland'), ('215', 'TNO'), ('216', 'TNO'), ('217', 'TUBITAK National Observatory'), ('218', 'TUG'), ('219', 'Thai National Observatory'), ('220', 'The Hale Telescope'), ('221', 'UKIRT'), ('222', 'United Kingdom Infrared Telescope'), ('223', 'V1'), ('224', 'VIRGO'), ('225', 'Vainu Bappu Observatory'), ('226', 'Very Large Array'), ('227', 'Virgo'), ('228', 'Virgo Observatory'), ('229', 'W. M. Keck Observatory'), ('230', 'WIYN'), ('231', 'WIYN 3.5 m'), ('232', 'WIYN Observatory'), ('233', 'Whipple'), ('234', 'Whipple Observatory'), ('235', 'Winer'), ('236', 'Winer Observatory'), ('237', 'Wise Observatory'), ('238', 'aao'), ('239', 'alma'), ('240', 'ao'), ('241', 'apo'), ('242', 'arca'), ('243', 'arecibo'), ('244', 'askap'), ('245', 'bbso'), ('246', 'bmo'), ('247', 'cfht'), ('248', 'chime'), ('249', 'ckoirama'), ('250', 'ctio'), ('251', 'dao'), ('252', 'dct'), ('253', 'dkist'), ('254', 'drao'), ('255', 'effelsberg'), ('256', 'ekar'), ('257', 'example_site'), ('258', 'fast'), ('259', 'flwo'), ('260', 'gbt'), ('261', 'gemini_north'), ('262', 'gemini_south'), ('263', 'gemn'), ('264', 'gems'), ('265', 'geo_600'), ('266', 'gmrt'), ('267', 'greenwich'), ('268', 'haleakala'), ('269', 'halo'), ('270', 'hcro'), ('271', 'het'), ('272', 'hyperK'), ('273', 'hyperk'), ('274', 'iao'), ('275', 'icecube'), ('276', 'irtf'), ('277', 'jcmt'), ('278', 'kagra'), ('279', 'keck'), ('280', 'km3net arca'), ('281', 'km3net orca'), ('282', 'kpno'), ('283', 'lapalma'), ('284', 'lasilla'), ('285', 'lbt'), ('286', 'lco'), ('287', 'ldt'), ('288', 'lho_4k'), ('289', 'lick'), ('290', 'llo_4k'), ('291', 'lo-am'), ('292', 'lo-mh'), ('293', 'lofar'), ('294', 'lowell'), ('295', 'lwa1'), ('296', 'mars_hill'), ('297', 'mcdonald'), ('298', 'mdm'), ('299', 'medicina'), ('300', 'meerkat'), ('301', 'mh'), ('302', 'mma'), ('303', 'mmt'), ('304', 'mro'), ('305', 'mso'), ('306', 'mtbigelow'), ('307', 'mwa'), ('308', 'mwo'), ('309', 'nancay'), ('310', 'noto'), ('311', 'nova'), ('312', 'oca'), ('313', 'ohp'), ('314', 'omm'), ('315', 'orca'), ('316', 'ovro'), ('317', 'paranal'), ('318', 'parkes'), ('319', 'pto'), ('320', 'rubin'), ('321', 'rubin_aux'), ('322', 'salt'), ('323', 'sirene'), ('324', 'sno+'), ('325', 'spm'), ('326', 'spo'), ('327', 'srt'), ('328', 'sso'), ('329', 'superK'), ('330', 'superk'), ('331', 'teide'), ('332', 'tona'), ('333', 'tug'), ('334', 'ukirt'), ('335', 'vbo'), ('336', 'virgo'), ('337', 'vla'), ('338', 'winer'), ('339', 'wise'), ('340', 'wise'), ('341', 'wiyn')], help_text='The telescope at which the detection of this TOA was made', max_length=4, null=True),
        ),
    ]
