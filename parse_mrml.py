"""
Parse a mrml (slicer hierarchy file) created by manually adding fiber clusters
to a tracts hierarchy.
Returns results in json format

Usage:
    parse_mrml.py [options] <filename>

Options:
    -h --help   Show this screen
    --quiet     Minimise logging
    --verbose   Maximise logging
"""
from __future__ import absolute_import
import xml.etree.ElementTree as ET
import json
from docopt import docopt
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)


class MapTracts(object):
    """
    Object to extract tract membership for clusers defined in an mrml file
    """
    def __init__(self, fname):
        root = self.load_mrml(fname)
        tracts = self.find_tract_names(root)
        for tract_name, tract_id in tracts.iteritems():
            clusters = self.find_clusters(tract_id, root)
            tracts[tract_name] = clusters

        self.tract_map = tracts

    def load_mrml(self, fname):
        '''
        Loads a mrml file, returns an xml root object
        '''
        tree = ET.parse(fname)
        root = tree.getroot()
        return(root)

    def find_tract_names(self, root):
        '''
        Parses and XML root object to find tract names
        returns a list of tracts
        '''
        items = root.findall('ModelHierarchy')
        tracts = {el.get('name'): el.get('id')
                  for el in items
                  if not el.get('name').startswith('ModelHierarchy')}
        logger.info('Found {} tracts:{}'.format(len(tracts),
                                                ' : '.join(tracts.keys())))
        return(tracts)

    def find_clusters(self, nodeid, root):
        """
        Parses the xml root object to find clusters assigned to a single tract
        identified by nodeid
        returns a list of tuplets (cluster, cluster_file)
        """
        cluster_files = []
        items = root.findall('.//*[@parentNodeRef="{}"]'.format(nodeid))
        fiberBundles = [i.get('associatedNodeRef') for i in items]
        fiberBundles = set(fiberBundles)

        for bundle_name in fiberBundles:
            bundle = root.findall('.//*FiberBundle[@id="{}"]'.format(bundle_name))
            storage_node_name = bundle[0].get('storageNodeRef')
            cluster_name = bundle[0].get('name')
            storage_node = root.findall('.//*FiberBundleStorage[@id="{}"]'
                                        .format(storage_node_name))
            cluster_files.append((cluster_name, storage_node[0].get('fileName')))

        return(cluster_files)


def main():
    args = docopt(__doc__)
    if args['--quiet']:
        logger.setLevel(logging.ERROR)
    if args['--verbose']:
        logger.setLevel(logging.DEBUG)

    tract_map = tractmap(args['<filename>'])

    return(json.dumps(tract_map.tract_map))


if __name__ == '__main__':
    print main()
