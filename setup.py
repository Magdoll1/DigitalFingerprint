from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext

#ext_modules = [Extension("hello", ["ReadDF.pyx"])]
#
#setup(
##		name = 'ReadDF',
#		cmdclass = {'build_ext': build_ext},
#		ext_modules = ext_modules
#)

setup(name='DigitalFingerprint',
		version='0.0.1',
		packages=['DigitalFingerprint', 'DigitalFingerprint.utils', 'DigitalFingerprint.newick'],
		package_dir={'DigitalFingerprint':'src/', 'DigitalFingerprint.utils':'src/utils/',\
				'DigitalFingerprint.newick':'src/newick/'},
		)
