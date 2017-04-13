from distutils.core import setup
import os

basedir = os.path.dirname('__file__')
datadir = os.path.join(basedir, 'data')
clusterdir = os.path.join(datadir, 'clusters')
cluster_f = os.listdir(clusterdir)
cluster_f = [os.path.join(clusterdir, f) for f in cluster_f]

setup(name='tractmapper',
      version='0.1.0',
      description="Map DTI tract atlas to subject space and extract fiber coordinates",
      author="Tom Wright",
      author_email="tom@maladmin.com",
      py_modules=['get_subject_tract_coordinates', 'parse_mrml', 'tempdir'],
      scripts=['get_subject_tract_coordinates.py', 'parse_mrml.py'],
      data_files=[('data', ['data/clustered_whole_brain.vtp',
                            'data/clustered_tracts_display_100_percent_aem.mrml']),
                  ('data/clusters/', cluster_f)])
