#!/usr/bin/env python
"""
Calculates tract coordinates in subject space.
Transforms a DTI atlas to subject space and returns streamlines coordinates
    associated with tracts.

Usage:
    get_subject_tract_coordinates.py [options] <subjectFile>
    get_subject_tract_coordinates.py [options] <subjectFile> <anatFile>

Arguments:
    <subjectFile>   Full path to a tractography file
    <anatFile>      Full path to the subject nifti format DTI file

Options:
    --cluster-pattern=<pattern>     A regular expression used to limit files
                                    in <clusterDir>
                                    [default: ^.*cluster_\d{5}]
    --cleanup                       Delete temporary files.
    --debug                         Extra logging information
    --quiet                         Only log errors
    --mirtk_file=<file>             Path to the MIRTK singularity container
                                    [default: MIRTK.img]
    --work_dir=<dir>                Where to create intermediate files.
    --output=<output>               Path to the output file
    --atlas_file=<atlas_file>       Path to a tractography atlas file (vtp or vtk)
                                    [default: ../data/clustered_whole_brain.vtp]
    --cluster_dir=<cluster_dir>     Path to a folder containing the atlas tract clusters
                                    [default: ../data/clusters/]
    --mrml_file=<mrml_file>         Path to the atlas mrml (Slicer) file mapping clusters to tracts
                                    [default: ../data/clustered_tracts_display_100_percent_aem.mrml]

Returns:
    A json object with the start and end coordinates of fibers organised
    by tract
    Coordinates are in the same units as anatFile (mm)
    {'tract1': {start: [(x, y, z), (x, y, z)],
                end: [(x, y, x), (x, y, z)]}}}

Dependencies:
    These need to be on your PATH
        wm_register_to_atlas_new.py -   https://github.com/SlicerDMRI/whitematteranalysis
        tractconverter.py           -   https://github.com/MarcCote/tractconverter
    Path can be specified
        MIRTK.img                   -   singularity container

Details:
    --atlas_file, --cluster_dir or --mrml_file can be specified. If a relative
    path is provided it is interpreted relative to __file__
"""
import os
import subprocess
import tempdir
import logging
import sys
import shutil
import re
import glob
import json
from docopt import docopt
import numpy as np
import numpy.linalg as npl
from parse_mrml import MapTracts
import nibabel as nib
from nibabel import trackvis as tv

logging.basicConfig()
logger = logging.getLogger(__name__)

def __run_cmd(command):
    '''
    Wrapper for subprocess.call_check
    exits if return code is not 0
    '''
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError as e:
        msg = 'Error: {} failed with exit code:{}'.format(e.cmd, e.returncode)
        logger.error(msg)
        sys.exit(msg)
    except OSError as e:
        if e.errno == 2:
            msg = ('Failed to find command: {}. Did you load required modules?'
                   .format(command))
        else:
            msg = ('Cmd:{} failed. With excuse:{}'.format(command, str(e)))
        logger.error(msg)
        sys.exit(msg)


def register_tractography(srcFile, targetFile, outDir):
    """
    Register srcFile to targetFile using wm_register_to_atlas_new.py
    """
    cmd = ['wm_register_to_atlas_new.py',
           srcFile,
           targetFile,
           outDir]

    __run_cmd(cmd)


def convert_vtp_to_vtk(path, outPath=None):
    """
    Converts a vtp file to vtk
    Expects full path to the file to convert. If outpath is supplied
    writes output there, otherwise creates in same location as input.
    """
    srcPath, fName = os.path.split(os.path.abspath(path))
    basename, ext = os.path.splitext(fName)

    if not outPath:
        outPath = srcPath

    outFile = basename + '.vtk'

    cmd = ['docker', 'run', '--rm',
           '-v', '{}:/srcDir'.format(srcPath),
           '-v', '{}:/dstDir'.format(outPath),
           'biomedia/mirtk',
           'convert-pointset',
           os.path.join('/srcDir', fName),
           os.path.join('/dstDir', outFile)]

    cmd = ['singularity', 'run',
           '-B', '{}:/input'.format(srcPath),
           '-B', '{}:/output'.format(outPath),
           CONTAINER_FILE,
           'convert-pointset',
           os.path.join('/input', fName),
           os.path.join('/output', outFile)]

    __run_cmd(cmd)
    return(os.path.join(outPath, outFile))


def convert_vtk_to_trk(path, anatFile, outPath=None):
    """
    Converts a vtk file to a trk file
    Expects full path to the file to convert. If outpath is supplied
    writes output there, otherwise creates a temp file.
    """
    srcPath, fName = os.path.split(os.path.abspath(path))
    basename, ext = os.path.splitext(fName)
    if not outPath:
        outPath = srcPath

    outFile = os.path.join(outPath, basename + '.trk')
    import pdb; pdb.set_trace()
    cmd = ['TractConverter.py',
           '-i', path,
           '-o', outFile,
           '-a', anatFile,
           '-f']
    __run_cmd(cmd)
    return(outFile)


def get_streamlines_from_trk(trkFile):
    """
    Extracts streamlines from a .trk file.
    Returns a list of streamlines.
    """
    streams, hdr = tv.read(trkFile)
    streamlines = [i[0] for i in streams]
    return(streamlines)


def convert_atlas_to_streams(atlasFile, anatFile=None, outDir=None):
    """
    Converts an atlas file to a list of streamlines
    """
    return(get_streams_from_file(atlasFile, anatFile=anatFile, outDir=outDir))


def get_most_advanced_file(filename):
    """
    Checks for the most advanced of .vtp, .vtk or .trk file.
    Returns the correct filename.
    """
    basename = os.path.splitext(filename)[0]
    for ext in ['.trk', '.vtk', '.vtp']:
        fname = basename + ext
        if os.path.exists(fname):
            return fname


def get_streams_from_file(fName, anatFile=None, outDir=None):
    """
    Process an input file to extract streamlines.
    Takes care of file type conversions if needed.
    Returns a list of streamlines
    """
    _, ext = os.path.splitext(fName)
    convert_to_vtk = False
    convert_to_trk = False
    tmpDir = None

    if ext == '.vtp':
        if not anatFile:
            msg = 'An anatomy file is required to convert vtp to vtk'
            logger.error(msg)
            sys.exit(msg)
        convert_to_vtk = True
        convert_to_trk = True
    elif ext == '.vtk':
        convert_to_trk = True
    elif ext != '.trk':
        logger.error('Unrecognised input file:{}'.format(fName))

    if convert_to_vtk or convert_to_trk:
        if not outDir:
            tmpDir = tempfile.mkdtemp()
            outDir = tmpDir
        #os.chmod(tmpdir, 0666)

    if convert_to_vtk:
        logger.info('Converting file to vtk')
        fName = convert_vtp_to_vtk(fName, outDir)

    if convert_to_trk:
        logger.info('Converting file to trk')
        fName = convert_vtk_to_trk(fName, anatFile, outDir)

    logger.info('Extracting streamlines from file')
    streams = get_streamlines_from_trk(fName)

    if tmpDir:
        logger.debug('Cleaning up:{}'.format(tmpDir))
        shutil.rmtree(tmpDir)

    return(streams)


def convert_clusters_to_streams(clusterDir,
                                pattern=None,
                                outDir=None,
                                anatFile=None):
    """
    Process a folder of cluster files, extracting the stream lines.

    Inputs:
        clusterDir - directory containing the cluster files (.vtp)
        pattern - regex pattern to limit which files are processed.
            default - '^.*cluster_\d{5}'
        outDir - directory to create working files
        anatFile - subject nifti file that was used for tractography
            This can be left out if processing has already been done

    Return:
        A dict {clustername: [streamlines]}
    """
    # set a default pattern, so other files in the folder don't get processed
    if not pattern:
        pattern = '^.*cluster_\d{5}'

    clusters = glob.glob(os.path.join(clusterDir, '*'))
    # just get the cluster filename
    clusters = [os.path.splitext(os.path.basename(f))[0]
                for f in clusters]
    clusters = set(clusters)

    p = re.compile(pattern)
    clusters = [f for f in clusters if p.match(f)]

    logger.info('Found {} cluster files.'.format(len(clusters)))

    cluster_streams = {}

    # loop through all the possible files in reverse order
    # of computational difficulty
    for cluster_id in clusters:
        target_f = os.path.join(outDir, cluster_id)
        target_f = get_most_advanced_file(target_f)
        if not target_f:
            # no processing done, start with the raw file
            target_f = os.path.join(clusterDir, cluster_id)
            target_f = get_most_advanced_file(target_f)

        logger.info('Converting file:{}'.format(target_f))
        streams = get_streams_from_file(target_f, anatFile, outDir=outDir)
        cluster_streams[cluster_id] = streams

    return(cluster_streams)


def match_fibers_to_clusters(fiber_streams, cluster_streams):
    """
    Matches fiber streamlines to cluster streamlines.
    Returns a vector of length fiber_streams, with keys from cluster_streams.
    """
    logger.info('Matching streams to clusters')
    matches = [None] * len(fiber_streams)
    cluster_streams_flat = []
    cluster_labels_flat = []
    for key, cluster in cluster_streams.iteritems():
        for cluster_stream in cluster:
            cluster_streams_flat.append(cluster_stream)
            cluster_labels_flat.append(key)

    logger.info('{} streams in {} clusters.'.format(len(cluster_streams_flat),
                                                    len(cluster_streams)))

    for i, stream in enumerate(fiber_streams):
        logger.info('Searching for fiber:{}'.format(i))
        idx = [np.array_equal(stream, cluster_stream)
               for cluster_stream in cluster_streams_flat].index(True)

        matches[i] = cluster_labels_flat[idx]

    # for i, stream in enumerate(fiber_streams):
    #     logger.info('Searching for fiber:{}'.format(i))
    #     for key, cluster in cluster_streams.iteritems():
    #         for cluster_stream in cluster:
    #             if stream == cluster:
    #                 matches[i] = key
    #                 continue
    return(matches)


def get_stream_ends(streamlines, tractMap):
    """
    Extracts start end endpoints of fibers and maps to clusters
    """
    count = len(tractMap)
    assert count == len(streamlines), ('All streamlines should'
                                       ' be defined in the tractMap')
    tract_ends = {}
    for i, stream in enumerate(streamlines):
        logger.info('Extracting ends for stream {} / {}'.format(i, count))
        start_coords = streamlines[i][0].tolist()
        end_coords = streamlines[i][-1].tolist()

        try:
            tract_ends[tractMap[i]]['starts'].append(start_coords)
            tract_ends[tractMap[i]]['ends'].append(end_coords)
        except KeyError:
            tract_ends[tractMap[i]] = {'starts': [],
                                       'ends': []}
            tract_ends[tractMap[i]]['starts'].append(start_coords)
            tract_ends[tractMap[i]]['ends'].append(end_coords)
    return(tract_ends)


def map_clusters_to_tracts(cluster_list, tract_map):
    """
    Takes a vector of clusters and a dict of tract membership.
    Returns a vector of same length as clusters with values
    replaced by tract membership
    """

    cluster_map = {}
    tracts = tract_map.keys()
    #  rearrange tract map so can index by cluster instead of tract
    for tract in tracts:
        clusters = [v[0] for v in tract_map[tract]]
        for cluster in clusters:
            cluster_map[cluster] = tract

    for i, cluster in enumerate(cluster_list):
        cluster_list[i] = cluster_map[cluster]

    return(cluster_list)


def convert_mm_to_voxels(coords, anat):
    """
    Uses measurements from the anatomy file to convert coordinates
    from mm to voxels.
    """
    img = nib.load(anat)
    affine = img.affine
    for key, val in coords.iteritems():
        coords_start = val['starts']
        coords_end = val['ends']
        voxels_start = nib.affines.apply_affine(npl.inv(affine), coords_start)
        voxels_end = nib.affines.apply_affine(npl.inv(affine), coords_end)
        # convert from numpy array back to list for json.dumps
        voxels_start = voxels_start.tolist()
        voxels_end = voxels_end.tolist()
        coords[key]['starts'] = voxels_start
        coords[key]['ends'] = voxels_end
    return coords


def process_atlas(atlas_file, subject_file, output_dir, anatFile=None):
    """
    Convert an atlas file to streamlines in subject space.

    Checks to see how much processing has been performed on the atlas.
    Inputs:
        .vtp (or .vtk) atlas and subject files
        output_dir - directory to create working files
        anatFile - subject nifti file that was used for tractography
            This can be left out if processing has already been done

    Return:
        Dict {'registered': A list of streamlines from the registered atlas,
              'raw': Streamlines from the unregistered atlas.}
    """
    atlas_name = os.path.splitext(os.path.basename(atlas_file))[0]

    # find the most processed version of the atlas file
    atlas_raw = get_most_advanced_file(os.path.join(output_dir, atlas_name))
    if not atlas_raw:
        atlas_raw = get_most_advanced_file(atlas_file)

    streams_raw = convert_atlas_to_streams(atlas_raw,
                                           anatFile=anatFile,
                                           outDir=output_dir)

    # check if registration has already been performed
    atlas_reg = os.path.join(output_dir,
                             '{}/output_tractography/{}_reg.vtk')

    atlas_reg = atlas_reg.format(atlas_name, atlas_name)

    if not os.path.isfile(atlas_reg):
        # need to register atlas to subject space
        register_tractography(atlas_file, subject_file, output_dir)

    # next check if the registered atlas has already been converted to .trk
    atlas_reg = get_most_advanced_file(atlas_reg)
    streams_reg = convert_atlas_to_streams(atlas_reg,
                                           anatFile=anatFile,
                                           outDir=os.path.dirname(atlas_reg))
    return {'registered': streams_reg,
            'raw': streams_raw}


def clean_working_dir(outputDir):
    shutil.rmtree(outputDir)


def main(atlas_fibers, atlas_clusters, cluster_pattern,
         subject_fibers, mrml_map, subject_anat, output_dir,
         cleanup):

    # create working directories
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    cluster_dir = os.path.join(output_dir, 'clusters')
    if not os.path.isdir(cluster_dir):
        os.mkdir(cluster_dir)

    # convert a tractography atlas to subject space and get the streamlines
    atlas_streams = process_atlas(atlas_fibers,
                                  subject_fibers,
                                  output_dir,
                                  subject_anat)

    # convert clustered fibers in atlas space to identified streamlines
    cluster_streams = convert_clusters_to_streams(atlas_clusters,
                                                  anatFile=subject_anat,
                                                  pattern=cluster_pattern,
                                                  outDir=cluster_dir)
    # get the mapping from cluster id to tract
    tract_map = MapTracts(mrml_map)

    # match the tracts identified in the unregistered atlas to
    # fibers in the registered atlas
    matches = match_fibers_to_clusters(atlas_streams['raw'], cluster_streams)
    matches = map_clusters_to_tracts(matches, tract_map.tract_map)

    # use tract -> fiber map to obtain fiber end points from registered atlas
    # check to see if this atlas has already been registered, create if not.
    tract_ends = get_stream_ends(atlas_streams['registered'], matches)
    if subject_anat:
        tract_ends_voxels = convert_mm_to_voxels(tract_ends, subject_anat)

    if cleanup:
        clean_working_dir(output_dir)

    return json.dumps(tract_ends_voxels)


if __name__ == "__main__":
    arguments = docopt(__doc__)
    atlasFile = arguments['--atlas_file']
    clusterDir = arguments['--cluster_dir']
    mrmlFile = arguments['--mrml_file']
    subjectFile = arguments['<subjectFile>']
    workingDir = arguments['--work_dir']
    anatFile = arguments['<anatFile>']
    outfile = arguments['--output']

    CONTAINER_FILE = arguments['--mirtk_file']

    pattern = arguments['--cluster-pattern']

    if arguments['--cleanup']:
        cleanup = True
    else:
        cleanup = False

    if arguments['--debug']:
        logger.setLevel(logging.DEBUG)
    elif arguments['--quiet']:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    script_dir = os.path.dirname(__file__)
    if not os.path.isabs(atlasFile):
        atlasFile = os.path.abspath(os.path.join(script_dir, atlasFile))
    if not os.path.isabs(clusterDir):
        clusterDir = os.path.abspath(os.path.join(script_dir, clusterDir))
    if not os.path.isabs(mrmlFile):
        mrmlFile = os.path.abspath(os.path.join(script_dir, mrmlFile))

    if workingDir:
        ends = main(atlasFile,
                    clusterDir,
                    pattern,
                    subjectFile,
                    mrmlFile,
                    anatFile,
                    workingDir,
                    cleanup)
    else:
        with tempdir.TempDir(prefix="tractmap_") as workingDir:
            ends = main(atlasFile,
                        clusterDir,
                        pattern,
                        subjectFile,
                        mrmlFile,
                        anatFile,
                        workingDir,
                        cleanup)

    if outfile:
        with open(outfile, 'w+') as outfile:
            outfile.writelines(ends)
    else:
        print(ends)
