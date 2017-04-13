from distutils.core import setup
setup(name='tractmapper',
      version='0.1.0',
      description="Map DTI tract atlas to subject space and extract biber coordinates",
      author="Tom Wright",
      author_email="tom@maladmin.com",
      py_modules=['get_subject_tract_coordinates', 'parse_mrml'])
