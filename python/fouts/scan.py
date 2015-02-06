#!/usr/bin/env python

import sys
import subprocess
import os.path
import getopt
import threading
import shutil

global lock, cmds, stop

lock = threading.Lock()
cmds = []
stop = False

def usage():
    print 'scan.py:\n    arguments: samples_location, [other_options, echo_options, origin, axes, output]'

def buildScanCmds(samples, other_options='', echo_options='', origin='', axes='', output_dir='', verbose=False, svd=0):
    """Batch scans in multiple dimensions"""

    samples_ext = samples.split('.')[-1]

    if samples_ext == 'sst':

        state_ext = 'str'

    elif samples_ext == 'tst':

        state_ext = 'tct'

    else:

        raise Exception("Unrecognised samples extension '%s'" % samples_ext)


    path = samples.split('/')

    if path.count('mcmc'):

        ind = path.index('mcmc')

        sub_path = '/'.join(path[ind + 2:])
        sub_path = '.'.join(sub_path.split('.')[0:-1])

    else:
        sub_path = raw_input("Sub-path name: ")



    if not output_dir:
        output_dir = '/home/tclose/Data/Tractography/analysis/scan/%s' % sub_path




    if svd:

        #Round the value of svd up to the nearest multiple of 2.
        if svd % 2:
            svd += 1

        axes = '/home/tclose/Data/Tractography/analysis/svd/%s' % sub_path

        if os.path.isdir(axes):
            try:
                subprocess.call('svn delete --force %s' % axes)
            except OSError:
                shutil.rmtree(axes)

        os.makedirs(axes)

        svd_file = axes + '.' + samples_ext

        svd_cmd = '/home/tclose/Code/Tractography/bin/svd_fibres %s %s -max %d' % (samples, svd_file, svd)

        if verbose:
            print 'SVD Command: '
            print '\n\n%s\n' % svd_cmd

        subprocess.call(svd_cmd, shell=True)

        comp_i = 0

        for dir_i in range(0, svd, 2):

            subdir = 'principle-comp_%d-%d' % (dir_i, dir_i + 1)

            os.mkdir(os.path.join(axes, subdir))


        for comp_i in range(0, svd):

        #Round the value of svd down to the nearest multiple of 2.            
            if comp_i % 2:
                dir_i = comp_i - 1
            else:
                dir_i = comp_i

            subdir = 'principle-comp_%d-%d' % (dir_i, dir_i + 1)

            select_cmd = '/home/tclose/Code/Tractography/bin/select_fibres %s %s/%s/%d.%s -inc %d' % (svd_file, axes, subdir, comp_i, state_ext, comp_i)

            if verbose:
                print 'Select Command: '
                print '\n\n%s\n' % select_cmd

            subprocess.call(select_cmd, shell=True)


        os.remove(svd_file)
        try:
            os.remove(svd_file + 'x')
            os.remove(svd_file + 'xx')
        except OSError:
            None

    elif not axes:

        if state_ext == 'str':

            axes = '/home/tclose/Data/Tractography/analysis/scan/templates/strand/no_intens'

        elif state_ext == 'tct':

            axes = '/home/tclose/Data/Tractography/analysis/scan/templates/tract/no_intens'

        else:

            raise Exception("Unrecognised state extension '%s'" % state_ext)


    # If directory exists, clean all files withinness it first.
    if os.path.isdir(output_dir):
        try:
            subprocess.call('svn delete --force %s' % output_dir)
        except OSError:
            shutil.rmtree(output_dir)

    os.makedirs(output_dir)

    axes_subdir = os.listdir(axes)

    if '.svn' in axes_subdir:
        axes_subdir.remove('.svn')

    global cmds

    cmds = []

    for dir in axes_subdir:

        masks = os.listdir(os.path.join(axes, dir))

        if '.svn' in masks:
            masks.remove('.svn')

        for mask in masks[:]:
            if mask.split('.')[-1] <> state_ext:
                masks.remove(mask)

        mask1 = os.path.join(axes, dir, masks[0])
        mask2 = os.path.join(axes, dir, masks[1])

        if svd:
            strip_mask1 = 'comp_%s' % masks[0]
            strip_mask2 = 'comp_%s' % masks[1]
        else:
            strip_mask1 = '.'.join(masks[0].split('.')[0:-1])
            strip_mask2 = '.'.join(masks[1].split('.')[0:-1])

        echo_cmd = '/home/tclose/code/bin/echo_parameters %s -method scan -output %s/%s_%s.mif -axis1 %s -axis2 %s %s' % (samples, output_dir, strip_mask1, strip_mask2, mask1, mask2, ' '.join(echo_options))

        if verbose:
            print 'Echo Command: '
            print echo_cmd

        echo = subprocess.Popen(echo_cmd, stdout=subprocess.PIPE, shell=True).communicate()[0]

        if not echo:
            raise Exception ('Echo command failed: (%s)' % echo_cmd)

        cmd = '/home/tclose/code/bin/' + echo.strip() + ' ' + other_options

        cmds.append(cmd)

        if verbose:
            print 'Scan Command (echo): '
            print '\n%s\n\n' % cmd




    print "MATLAB plot command:\n\nplot_scan %s -state %s.iter.%s \n\n" % (output_dir, '.'.join(samples.split('.')[0:-1]), samples.split('.')[-1])


def runScanCmds():

    global cmds

    while True:
        lock.acquire()
        if not len(cmds):
            lock.release()
            return None
        cmd = cmds.pop()
        lock.release()

        subprocess.call(cmd, shell=True)


if __name__ == '__main__':

    try:
            opts, args = getopt.getopt(sys.argv[1:], "ho:a:p:e:vs:", ["help", "origin", "axes", "output", "echo_options", "verbose", "svd"])
    except getopt.GetoptError:
            usage()
            sys.exit(2)

    origin = ''
    axes = ''
    output_dir = ''
    echo_options = ''
    verbose = False
    stop = False
    svd = 0

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-o", "--origin"):
            origin = arg
        elif opt == ('-a', "--axes"):
            axes = arg
        elif opt in ("-p", "--output"):
            output_dir = arg
        elif opt in ("-e", "--echo_options"):
            output_dir = arg
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-s", "--svd"):
            svd = int(arg)

    if len(args) == 0:
        usage()
        sys.exit(2)

    samples = args[0]

    if len(args) > 1:
        other_options = ' '.join(args[1:])
    else:
        other_options = ''

    global cmds

    buildScanCmds(samples, other_options, echo_options, origin, axes, output_dir, verbose, svd)


    print 'Running %d scans...\n\n' % len(cmds)

    try: num_processors = os.sysconf('SC_NPROCESSORS_ONLN')
    except:
        try: num_processors = int(os.environ['NUMBER_OF_PROCESSORS'])
        except: num_processors = 1

    if verbose:
        print ''
        print 'launching ' + str(num_processors) + ' threads'
        print ''

    threads = []
    for i in range (1, num_processors):
        t = threading.Thread (target=runScanCmds);
        t.start()
        threads.append (t)

    exit_value = runScanCmds()

    for t in threads: t.join()

    exit (stop > 1)
