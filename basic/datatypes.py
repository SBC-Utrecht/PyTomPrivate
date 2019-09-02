#DATATYPES

datatype = [('DefocusU', 'f4'),
            ('DefocusV', 'f4'),
            ('DefocusAngle', 'f4'),
            ('Voltage', 'i4'),
            ('SphericalAberration', 'f4'),
            ('AmplitudeContrast', 'f4'),
            ('PhaseShift', 'f4'),
            ('PixelSpacing', 'f4'),
            ('MarkerDiameter', 'i4'),
            ('TiltAngle', 'f4'),
            ('RotationTheta', 'f4'),
            ('InPlaneRotation', 'f4'),
            ('TranslationX', 'f4'),
            ('TranslationY', 'f4'),
            ('TranslationZ', 'f4'),
            ('Magnification', 'f4'),
            ('Intensity', 'f4'),
            ('ImageSize', 'i4'),
            ('AcquisitionOrder', 'i4'),
            ('FileName', 'U1000')]

headerText = ''
units = ['um', 'um', 'deg', 'kV', 'mm', '', 'deg', 'A', 'A', 'deg', 'deg', 'deg', 'px', 'px', 'px', '', '','px', '', '' ]
fmt='%11.6f %11.6f %6.2f %4d %6.2f %4.2f %11.6f %11.6f %4d %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %5.3f %5.3f %4d %3d %s'

for n, h in enumerate(datatype):
    headerText += '{} {}\n'.format(h[0], '({})'.format(units[n])*(units[n]!=''))

for n, h in enumerate(datatype):
    headerText0 += '{} {}\n'.format(h[0], '({})'.format(units[n])*(units[n]!=''))

datatypeMR = [('MarkerIndex', 'i4'),
              ('OffsetX',     'f4'),
              ('OffsetY',     'f4'),
              ('OffsetZ',     'f4'),
              ('PositionX',   'f4'),
              ('PositionY',   'f4'),
              ('PositionZ',   'f4')]

headerMarkerResults = ''
unitsMR = ['', 'px', 'px', 'px', 'px', 'px', 'px']
fmtMR='%3d %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f'
for n, h in enumerate(datatypeMR):
    headerMarkerResults += '{} {}\n'.format(h[0], '({})'.format(unitsMR[n])*(unitsMR[n]!=''))


datatypeAR = [('AlignmentTransX', 'f4'),
              ('AlignmentTransY', 'f4'),
              ('TiltAngle',       'f4'),
              ('InPlaneRotation', 'f4'),
              ('Magnification',   'f4'),
              ('FileName', 'U1000')]

headerAlignmentResults = ''
unitsAR = ['px', 'px', 'degrees', 'degrees', '', '']
fmtAR='%15.10f %15.10f %15.10f %15.10f %15.10f %s'
for n, h in enumerate(datatypeAR):
    headerAlignmentResults += '{} {}\n'.format(h[0], '({})'.format(unitsAR[n])*(unitsAR[n]!=''))
