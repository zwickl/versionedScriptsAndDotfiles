#!/usr/bin/env python
import sys
#import re
from os import path
#import cPickle
from itertools import izip, cycle
import collections
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl

if __name__ == "__main__":
    import doctest
    doctest.testmod()


def read_files_and_split_columns_as_strings(filenames, skipRows=None, ignorePatts=None, allowMissing=False):
    if skipRows and ignorePatts:
        exit('pass either skipRows or ignorePatts to\nread_and_split_columns_as_strings, not both')
    data = []
    if ignorePatts:
        if isinstance(ignorePatts, str):
            ignorePatts = [ignorePatts]
        compiledPatts = []
        for patt in ignorePatts:
            compiledPatts.append(re.compile(patt)) 
    for f in filenames:
        if ignorePatts:
            rawData = [l for l in open(f, 'rb')]
            theseData = []
            for line in rawData:
                for cpat in compiledPatts:
                    if cpat.search(line):
                        break
                else:
                    theseData.append(line.strip().split())
        else:
            if not path.exists(f) and allowMissing:
                sys.stderr.write('WARNING: skipping missing file %s\n' % f)
            else:
                theseData = [l.strip().split() for l in open(f, 'rb')]
        if skipRows:
            theseData = theseData[skipRows:]
        data.append(theseData)
    return data


def translate_annotation_code(c):
    #code to indicate intron/exon regions
    '''
    G = full gap

    E, F = full exon (no gap, gap)
    I, J = full intron (no gap, gap)
    M, N = mixed intron/exon (no gap, gap)
    '''
    translate = {'E':5, 'F':5, 'I':-5, 'J':-5, 'M':0, 'N':0}
    return translate[c]


def translate_annotation_code_string(lines, column):
    translate = {'E':5, 'F':5, 'I':-5, 'J':-5, 'M':0, 'N':0}
    return [ translate[line[colnum]] for line in lines ]


def string_to_float_or_zero(val):
    try:
        return float(val)
    except ValueError:
        return 0.0


def path_to_plot_title(filename, sep='\n'):
    '''
    ../mafft.cds/
    ../mafft.cds.allBadAlignAA/
    ../mafft.cds.gblocks/
    ../mafft.genes/
    ../mafft.genes.gblocks/
    ../mafft.genes.intronsStripped/
    ../mafft.genes.noFullIntrons/
    ../mafft.genes.stripNs/
    
    all
    11T

    boot
    ml
    '''

    filename = filename.lower()

    plotTitle = ''
    if 'mafft' in filename:
        plotTitle += 'MAFFT'
    elif 'prank' in filename:
        plotTitle += 'PRANK'
    elif 'muscle' in filename:
        plotTitle += 'MUSCLE'

    plotTitle += sep

    if 'cds.allbadalignaa' in filename:
        plotTitle += 'CDS-BADCOLAA '
    elif 'cds.allbadalign' in filename:
        plotTitle += 'CDS-BADCOL '
    elif 'cds.maskedbadalign.onlyaa' in filename:
        plotTitle += 'CDS-MASKEDAA '
    elif 'cds.maskedbadalign.nobrach' in filename:
        plotTitle += 'CDS-MASKEDNOBRACH '
    elif 'cds.maskedbadalign' in filename:
        plotTitle += 'CDS-MASKED '
    elif 'cds.gblocks.maskedbadalign' in filename:
        plotTitle += 'CDS+GBLOCKS%s-MASKED' % sep
    elif 'cds.gblocks' in filename:
        plotTitle += 'CDS+GBLOCKS'
    elif 'cds' in filename:
        plotTitle += 'CDS'
    elif 'genes.gblocks.maskedbadalign' in filename:
        plotTitle += 'GENE+GBLOCKS%s-MASKED' % sep
    elif 'genes.gblocks' in filename:
        plotTitle += 'GENE+GBLOCKS'
    elif 'genes.intronsstripped' in filename:
        #plotTitle += 'GENE-INTRON'
        plotTitle += 'MASKED-INTRON'
    elif 'genes.nofullintrons' in filename:
        plotTitle += 'GENE-FULLI'
    elif 'genes.allbadalignAA' in filename:
        plotTitle += 'GENE-BADCOLAA '
    elif 'genes.allbadalign' in filename:
        plotTitle += 'GENE-BADCOL '
    elif 'genes.maskedbadalign.onlyaa' in filename:
        plotTitle += 'GENE-MASKEDAA '
    elif 'genes.maskedbadalign.nobrach' in filename:
        plotTitle += 'GENE-MASKEDNOBRACH '
    elif 'genes.maskedbadalign' in filename:
        plotTitle += 'GENE-MASKED '
    elif 'genes.stripns' in filename:
        plotTitle += 'GENE-STRIPNs '
    elif 'genes' in filename:
        plotTitle += 'GENE'
    else:
        plotTitle += 'UNKNOWN'

    if 'ssr' in filename:
        plotTitle += 'SSR'

    return plotTitle


def path_to_plot_title_random_supermatrix(filename, sep='\n'):
    plotTitle = path_to_plot_title(filename, sep)

    if 'correctaanoindica' in filename:
        plotTitle += '%scorrect AA-indica' % sep
    elif 'correctallnoindica' in filename:
        plotTitle += '%scorrect All-indica' % sep
    elif 'correctaa' in filename:
        plotTitle += '%scorrect AA' % sep
    elif 'correctall' in filename:
        plotTitle += '%scorrect All' % sep
    elif 'correctBBCC' in filename:
        plotTitle += '%scorrect AA-BB-CC' % sep

    return plotTitle


def prepare_plot_kwargs(kwargDict, optionList, fontscaler=1.0):
    '''
    parse and prepare kwargs in various formats, add new options,
    This is mainly for use with the PlottingArgumentParser defined below
    and would be called like this:
        
        superTitleKwargs = prepare_kwargs(options.additionalSuperTitleKwargs,
                [('fontsize', options.superTitleFontsize)],
                fontscaler=scaler)
        plt.suptitle(options.superTitle, **superTitleKwargs)

    where 
    -the first argument is any extra arbitrary kwargs that might have been
    massed on the command line
    -then a list of things to set specifically.  Defaults for those should 
    have been set when the parser was instantiated, but can be overridden
    
    kwargDict - can just be a list of strings, i.e. ['blah=2', 'poo=4']
        or a list of tuples [('blah', 2), ('poo', 4)]
        or an actual dictionary {'blah':2, 'poo':4}
    optionList - list of (key, value) tuples, added to returned dict if
        not already in it
    fontscaler - value to rescale fonts by, if desired
    returns outKwargs - dict of kwargs to pass to other functions
    '''

    outKwargs = {}
    if kwargDict:
        error = False
        if isinstance(kwargDict, dict):
            outKwargs = kwargDict
        elif isinstance(kwargDict, list) or isinstance(kwargDict, set) or isinstance(kwargDict, tuple):
            if isinstance(kwargDict[0], str):
                outKwargs = dict([ kw.split('=') for kw in kwargDict ])
            elif len(kwargDict[0]) == 2:
                outKwargs = dict(kwargDict)
            else:
                error = True
        else:
            error = True
        if error:
            exit('kwargDict must be a dictionary, or a list, set or tuple containing strings, or a list, set or tuple containing strings containing 2-tuples')
   
    #this will extract an actual dictionary out of the text list of kwargs specified on the command line
    #outKwargs = dict([ kw.split('=') for kw in kwargDict ])
    #this just ensures that if a kwarg is specified at the command line for something, it trumps any actual argparse option for it
    for arg, setting in optionList:
        outKwargs.setdefault(arg, setting)

    for key, val in outKwargs.iteritems():
        #try to convert values to numbers where possible
        try:
            outKwargs[key] = float(val)
        except ValueError:
            pass
        except TypeError:
            pass
        if val == 'False':
            outKwargs[key] = False
        elif val == 'True':
            outKwargs[key] = True
        
        if 'fontsize' in key or 'labelsize' in key:
            outKwargs[key] = outKwargs[key] * fontscaler

    #print 'reading kwargs'
    #print kwargDict
    #print optionList
    #print outKwargs

    return outKwargs

class PlottingStyleArgumentParser(ArgumentParser):
    '''abandoned this for the moment'''
    def __init__(self, **kwargs):
        #base constructor
        super(PlottingStyleArgumentParser, self).__init__(self, **kwargs)
        
        #argGroup = self.add_argument_group(title='style arguments')
        
        #argGroup.add_argument('-p', '--poo', dest='poop',type=str, default=None, 
        #                    help='test')


'''
class ArgparseAction_AppendToSpecifiedList(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
'''

class PlottingHistogramArgumentParser(ArgumentParser):
    '''A generic argparse parser that comes prepackaged with a bunch of common arguments
    for matplotlib/pyplot preloaded.
    Default values for the created parser can be passed in as keyword arguments.  Pass
    False to remove a built in option.'''

    def __init__(self, **kwargs):
        '''These are the defaults which can be overwridden with keyword arguments passed
        to the constructor.
        '''
        defaultSubplotKwargs = kwargs.pop('defaultSubplotmKwargs', [])
        defaultHistKwargs = kwargs.pop('defaultHistKwargs', [])
        defaultNumBins = kwargs.pop('defaultNumBins', 20)
        defaultBarFaceColor = kwargs.pop('defaultBarFaceColor', 'gray')
        
        if defaultNumBins is not False:
            self.add_argument('-n', '--num-bins', dest='numBins', type=int, default=defaultNumBins,
                                help='number of evenly spaced histogram bins (default %s)' % str(defaultNumBins))

        if defaultBarFaceColor is not False:
            self.add_argument('-c', '--bar-color', dest='barFaceColor', type=str, default=defaultBarFaceColor,
                                help='color of histogram bars (default %s)' % str(defaultBarFaceColor))

 
class PlottingArgumentParser(ArgumentParser):
    '''A generic argparse parser that comes prepackaged with a bunch of common arguments
    for matplotlib/pyplot preloaded.
    Default values for the created parser can be passed in as keyword arguments.  Pass
    False to remove a built in option.'''

    def __init__(self, **kwargs):
        '''These are the defaults which can be overwridden with keyword arguments passed
        to the constructor.
        '''
        defaultSubplotKwargs = kwargs.pop('defaultSuplotKwargs', [])
        
        #plot defaults
        defaultGreys = kwargs.pop('defaultGreys', ['0.0'])
        defaultColors = kwargs.pop('defaultColors', None)
        defaultMarkerSizes = kwargs.pop('defaultMarkerSizes', [12.0])
        defaultMarkers = kwargs.pop('defaultMarkers', 'osx*^h')
        defaultLineWidth = kwargs.pop('defaultLineWidth', [3.0])
        defaultCycleStylesWithinFiles = kwargs.pop('defaultCycleStylesWithinFiles', True)
        defaultPlotKwargs = kwargs.pop('defaultPlotKwargs', [])
        
        #axis defaults
        defaultYrange = kwargs.pop('defaultYrange', None)
        defaultXrange = kwargs.pop('defaultXrange', None)
        defaultAxisTickFontsize = kwargs.pop('defaultAxisTicFontsize', 16.0)
        defaultAxisLabelFontsize = kwargs.pop('defaultAxisLabelFontsize', 32.0)
        defaultXTickFunction = kwargs.pop('defaultXTickFunction', None)
        defaultTickKwargs = kwargs.pop('defaultTickKwargs', [])
        defaultXTickKwargs = kwargs.pop('defaultXTickKwargs', [])
        defaultYTickKwargs = kwargs.pop('defaultYTickKwargs', [])
        defaultXTickLabelLocation = kwargs.pop('defaultXTickLabelLocation', 'm')
        defaultYTickLabelLocation = kwargs.pop('defaultYTickLabelLocation', 'm')
        defaultDrawAxesAtZero = kwargs.pop('defaultDrawAxesAtZero', None)
        
        #axis labels
        defaultXlabel = kwargs.pop('defaultXlabel', None)
        defaultYlabel = kwargs.pop('defaultYlabel', None)
        defaultXlabelLocation = kwargs.pop('defaultXlabelLocation', 's')
        defaultYlabelLocation = kwargs.pop('defaultYlabelLocation', 's')
        defaultXlabelKwargs = kwargs.pop('defaultXlabelKwargs', [])
        defaultYlabelKwargs = kwargs.pop('defaultYlabelKwargs', [])
        
        #title defaults
        defaultTitleFunc = kwargs.pop('defaultTitleFunc', None)
        defaultTitles = kwargs.pop('defaultTitles', None)
        defaultTitleFontsize = kwargs.pop('defaultTitleFontsize', 18.0)
        defaultSuperTitle = kwargs.pop('defaultSuperTitle', None)
        defaultSuperTitleFontsize = kwargs.pop('defaultSuperTitleFontsize', 32.0)
        defaultTitleHalign = kwargs.pop('defaultTitleHAlign', 'center')
        defaultTitleValign = kwargs.pop('defaultTitleVAlign', 'top')
        defaultTitleKwargs = kwargs.pop('defaultTitleKwargs', [])
        defaultSuperTitleKwargs = kwargs.pop('defaultSuperTitleKwargs', [])
        
        defaultNcol = kwargs.pop('defaultNcol', 2)
        defaultOnePlot = kwargs.pop('defaultOnePlot', None)
        defaultSubplotPerFunction = kwargs.pop('defaultSubplotPerFunction', None)
        defaultInput = kwargs.pop('defaultInput', True)
        defaultOutput = kwargs.pop('defaultOutput', True)

        defaultMissingOk = kwargs.pop('defaultMissingOk', None)
        defaultDataColumnFunc = kwargs.pop('defaultDataColumnFunc', None)

        defaultHatches = kwargs.pop('defaultHatches', 'Xo+ ')

        defaultSkipRows = kwargs.pop('defaultSkipRows', 0)
        defaultRowIgnorePatterns = kwargs.pop('defaultRowIgnorePatterns', [])
        
        defaultHistogram = kwargs.pop('defaultHistogram', None)
        defaultHistogramKwargs = kwargs.pop('defaultHistogramKwargs', [])
        defaultHistogramBins = kwargs.pop('defaultHistogramBinsBins', 20)
       
        #base constructor
        super(PlottingArgumentParser, self).__init__(self, *kwargs)
      
        if defaultSubplotKwargs is not False:
            self.add_argument('--subplot-kwargs', dest='additionalSubplotKwargs', nargs='*', type=str, default=defaultSubplotKwargs,
                                help='kwargs to be passed on to matplotlib Figure.subplots_adjust function (default %s)' % ' '.join(defaultSubplotKwargs))
        
        ################
        #argument for plotting of lines/points
        styleArgs = self.add_argument_group('ARGUMENTS FOR POINT AND LINE STYLES')

        if defaultGreys is not False:
            styleArgs.add_argument('-gv', '--grey-values', dest='greyValues', nargs='*', type=str, default=defaultGreys,
                                help='floating point values in range 0.0 - 1.0 (black to white) indicating grey value cycle of lines (default %s)' % str(defaultGreys))

        if defaultColors is not False:
            styleArgs.add_argument('-cv', '--color-values', dest='colorValues', nargs='*', type=str, default=defaultColors,
                                help='valid matplotlib color specs indicating color cycle of lines (default %s)' % str(defaultColors))
        
        if defaultMarkers is not False:
            styleArgs.add_argument('-m', '--markers', dest='markers', type=str, default=defaultMarkers,
                                help='single string with marker designations for 1 2 3 (default %s)' % str(defaultMarkers))

        if defaultMarkerSizes is not False:
            styleArgs.add_argument('-ms', '--marker-sizes', dest='markerSizes', nargs='*', type=int, default=defaultMarkerSizes,
                                help='floating point values indicating size of markers in points for series 1 2 3, or a single value to apply to all (default %s)' % str(defaultMarkerSizes))

        if defaultLineWidth is not False:
            styleArgs.add_argument('-lw', '--line-width', dest='lineWidth', nargs='*', type=float, default=defaultLineWidth,
                                help='point size of lines to cycle through (default %s)' % str(defaultLineWidth))
 
        if defaultCycleStylesWithinFiles is not False:
            self.add_argument('--styles-per-file', dest='cycleStylesWithinFiles', default=defaultCycleStylesWithinFiles, action='store_false',
                                                    help='use the same line/point styles for all series within a given file, rather \
                                                    than cycling through styles for each series _within_ a file')

        if defaultPlotKwargs is not False:
            self.add_argument('--plot-kwargs', dest='additionalPlotKwargs', nargs='*', type=str, default=defaultPlotKwargs,
                                help='kwargs to be passed on to matplotlib axes.plot function (default %s)' % ' '.join(defaultPlotKwargs))
           
        ############
        #axis labels
        axisLabelArgs = self.add_argument_group('ARGUMENTS FOR AXIS LABELING')
        if defaultXlabel is not False:
            axisLabelArgs.add_argument('-xl', '--x-label', dest='xLabel', type=str, default=defaultXlabel,
                                help='label for x axis (default %s)' % str(defaultXlabel))

            if defaultXlabelLocation is not False:
                axisLabelArgs.add_argument('-xll', '--x-label-location', dest='xLabelLocation', type=str, default=defaultXlabelLocation, choices=['s', 'm', 'n'],
                                    help='where xlabels should appear in multiplot: s (single centered), m (one per plot on margins), n (none)')
        
            if defaultXlabelKwargs is not False:
                axisLabelArgs.add_argument('--x-label-kwargs', dest='additionalXlabelKwargs', nargs='*', type=str, default=defaultXlabelKwargs,
                                help='kwargs to be passed on to pyplot.xlabel function (default %s)' % ' '.join(defaultXlabelKwargs))
        
        if defaultYlabel is not False:
            axisLabelArgs.add_argument('-yl', '--y-label', dest='yLabel', type=str, default=defaultYlabel,
                                help='label for y axis (default %s)' % str(defaultYlabel))

            if defaultYlabelLocation is not False:
                axisLabelArgs.add_argument('-yll', '--y-label-location', dest='yLabelLocation', type=str, default=defaultYlabelLocation, choices=['s', 'm', 'n'],
                                    help='where ylabels should appear in multiplot: s (single centered), m (one per plot on margins), n (none)')
        
            if defaultYlabelKwargs is not False:
                axisLabelArgs.add_argument('--y-label-kwargs', dest='additionalYlabelKwargs', nargs='*', type=str, default=defaultYlabelKwargs,
                                help='kwargs to be passed on to pyplot.ylabel function (default %s)' % ' '.join(defaultYlabelKwargs))
        
        if defaultXlabel is not False or defaultYlabel is not False:
            if defaultAxisLabelFontsize is not False:
                axisLabelArgs.add_argument('-als', '--axis-label-fontsize', dest='axisLabelFontsize', type=int, default=defaultAxisLabelFontsize,
                                    help='font size of axis labels (default %d)' % defaultAxisLabelFontsize)
            
        ##########
        #axis range and tics
        axisArgs = self.add_argument_group('ARGUMENTS FOR AXIS RANGES AND TICKS')
        
        if defaultYrange is not False:
            axisArgs.add_argument('-yr', '--y-range', dest='yRange', nargs='*', type=float, default=defaultYrange,
                                help='min and max values on y axis pass one min max pair to have it shared, or a series of min max min max for each subplot (default is auto determined by pyplot from the data)')

        if defaultXrange is not False:
            axisArgs.add_argument('-xr', '--x-range', dest='xRange', nargs='*', type=float, default=defaultXrange,
                                help='min and max values on x axis pass one min max pair to have it shared, or a series of min max min max for each subplot (default is auto determined by pyplot from the data)')
        '''
        if defaultYrange is not False:
            axisArgs.add_argument('-yr', '--y-range', dest='yRange', nargs=2, type=float, default=defaultYrange,
                                help='min and max values on y axis (default is auto determined by pyplot from the data)')

        if defaultXrange is not False:
            axisArgs.add_argument('-xr', '--x-range', dest='xRange', nargs=2, type=float, default=defaultXrange,
                                help='min and max values on x axis (default is auto determined by pyplot from the data)')
        '''

        if defaultAxisTickFontsize is not False:
            axisArgs.add_argument('-ats', '--axis-tick-fontsize', dest='axisTickFontsize', type=int, default=defaultAxisTickFontsize,
                                help='font size of cumulative graph tic labels (default %d)' % defaultAxisTickFontsize)

        if defaultXTickFunction is not False:
            axisArgs.add_argument('--x-tick-function', dest='xTickFunction', type=str, default=defaultXTickFunction,
                                help='optional lambda or named function to get x-axis tick labels from input file')

        if defaultTickKwargs is not False:
            axisArgs.add_argument('--tick-kwargs', dest='additionalTickKwargs', nargs='*', type=str, default=defaultTickKwargs,
                                help='kwargs to be passed on to matplotlib axes.tick_params function, applied to both axes (default %s)' % ' '.join(defaultTickKwargs))
        
        if defaultXTickKwargs is not False:
            axisArgs.add_argument('--x-tick-kwargs', dest='additionalXTickKwargs', nargs='*', type=str, default=defaultXTickKwargs,
                                help='kwargs to be passed on to matplotlib axes.tick_params function, applied to only x axis (default %s)' % ' '.join(defaultXTickKwargs))
        
        if defaultYTickKwargs is not False:
            axisArgs.add_argument('--y-tick-kwargs', dest='additionalYTickKwargs', nargs='*', type=str, default=defaultYTickKwargs,
                                help='kwargs to be passed on to matplotlib axes.tick_params function, applied to only y axis (default %s)' % ' '.join(defaultYTickKwargs))
        
        if defaultDrawAxesAtZero is not False:
            axisArgs.add_argument('--draw-axes-at-zero', dest='drawAxesAtZero', default=defaultDrawAxesAtZero, action='store_true', 
                                help='whether to draw horizontal and vertical lines through 0, 0')
       
        if defaultXTickLabelLocation is not False:
            axisArgs.add_argument('--x-tick-label-location', dest='xTickLabelLocation', default=defaultXTickLabelLocation, choices=['m', 'e', 'n'],
                    help='whether to draw X tick labels only on (m)arginal plots, (n)one, or on (e)very plot. Default: %s' % defaultXTickLabelLocation)
       
        if defaultYTickLabelLocation is not False:
            axisArgs.add_argument('--y-tick-label-location', dest='yTickLabelLocation', default=defaultYTickLabelLocation, choices=['m', 'e', 'n'],
                    help='whether to draw Y tick labels only on (m)arginal plots, (n)one, or on (e)very plot. Default: %s' % defaultYTickLabelLocation)
       
        #########
        #titles
        titleArgs = self.add_argument_group('ARGUMENTS FOR PLOT AND SUBPLOT TITLES')
        if defaultTitleFunc is not False or defaultTitle is not False:
            titleType = titleArgs.add_mutually_exclusive_group()
            titleType.add_argument('--titles', dest='titles', nargs='*', type=str, default=defaultTitles,
                                help='plot titles, HACKY!, must supply one per SERIES, although if multiple series per plot then later will supercede earlier')
            
            titleType.add_argument('--title-func', dest='titleFunc', type=str, default=defaultTitleFunc,
                                help='function used to map arbitrary strings (datafile path names) to plot names')
            
            '''
            if isinstance(defaultTitle, str):
                titleArgs.add_argument('-t', '--title', dest='title', type=str, default=defaultTitle,
                                    help='plot title')
            else:
                titleArgs.add_argument('-t', '--title-func', dest='titleFunc', type=str, default=defaultTitle,
                                    help='function used to map arbitrary strings to plot names')
            '''

            if defaultTitleHalign is not False:
                titleArgs.add_argument('-tha', '--title-horiz-align', dest='titleHalign', type=str, default=defaultTitleHalign, choices='lcr',
                                    help='horizontal alignment of title (default c)')

            if defaultTitleValign is not False:
                titleArgs.add_argument('-tva', '--title-vert-align', dest='titleValign', type=str, default=defaultTitleValign, choices='tcb',
                                    help='vertical alignment of title (default t)')
        
            if defaultTitleKwargs is not False:
                titleArgs.add_argument('--title-kwargs', dest='additionalTitleKwargs', nargs='*', type=str, default=defaultTitleKwargs,
                                    help='kwargs to be passed on to matplotlib axes.setTitle function (default %s)' % ' '.join(defaultTitleKwargs))

            if defaultTitleFontsize is not False:
                titleArgs.add_argument('-ts', '--title-fontsize', dest='titleFontsize', type=int, default=defaultTitleFontsize,
                                    help='font size of title (default %d)' % defaultTitleFontsize)

        if defaultSuperTitle is not False:
            titleArgs.add_argument('-st', '--super-title', dest='superTitle', type=str, default=defaultSuperTitle,
                                help='plot super title on multiplot')

            if defaultSuperTitleFontsize is not False:
                titleArgs.add_argument('-stf', '--super-title-fontsize', dest='superTitleFontsize', type=int, default=defaultSuperTitleFontsize,
                                    help='font size of super title of multiplot (default %d)' % defaultSuperTitleFontsize)

            if defaultSuperTitleKwargs is not False:
                titleArgs.add_argument('--super-title-kwargs', dest='additionalSuperTitleKwargs', nargs='*', type=str, default=defaultSuperTitleKwargs,
                                    help='kwargs to be passed on to pyplot.suptitle function (default %s)' % ' '.join(defaultSuperTitleKwargs))
       
        ##############
        dataArgs = self.add_argument_group('ARGUMENTS FOR PARSING AND MANIPULATING DATA')
        
        if defaultInput is not False:
            dataArgs.add_argument('-i', '--input', dest='inFiles', nargs='*', type=str, default=None, 
                                help='Series of intput files with data to be plotted')

        if defaultDataColumnFunc is not False:
            dataArgs.add_argument('--data-column-func', dest='dataColumnFunc', nargs='*', type=str, default=defaultDataColumnFunc,
                                help='functions to select or convert column(s) in datafiles to a plotable series for pyplot.plot. \
                                To plot a single series of only x values, something like \'lambda rows:([float(r[1]) - float(r[3]) for r in rows])\'\
                                will work.  To plot x-y points, \'lambda rows: ([float(r[1]) for r in rows], [float(r[3])/float(r[2]) for r in rows])\'\
                                Multiple functions can be used, in which case each is applied to each datafile.  (default \'%s\')' % defaultDataColumnFunc )
        
        if defaultSkipRows is not False:
            dataArgs.add_argument('--skip-rows', dest='skipRows', type=int, default=defaultSkipRows,
                                help='number of rows to skip at start of data file (default %d)' % defaultSkipRows )
        
        if defaultRowIgnorePatterns is not False:
            dataArgs.add_argument('--row-ignore-patterns', dest='rowIgnorePatterns', nargs='*', type=str, default=defaultRowIgnorePatterns,
                                help='lines in the data file matching these regex patterns are ignored (default %s)' % defaultRowIgnorePatterns )
        
        if defaultMissingOk is not False:
            dataArgs.add_argument('--missing-ok', dest='missingOk', action='store_true', default=False, 
                                help='allow some input files to be missing')

        ###############
        plotType = self.add_mutually_exclusive_group()
        if defaultSubplotPerFunction is not False:
            plotType.add_argument('--subplot-per-function', dest='subplotPerFunction', default=False, action='store_true',
                                help='plot each data column function on a new subplot, rather than using one subplot per _file_')
        
        if defaultOnePlot is not False:
            plotType.add_argument('--one-plot', dest='onePlot', default=False, action='store_true',
                                help='plot all data from all files and data functions on one axis, rather than using one file per subplot')
        
        if defaultNcol is not False:        
            self.add_argument('-c', '--columns', dest='ncols', type=int, default=defaultNcol,
                                help='number of figures to place per row in a multiplot')

        if defaultOutput is not False:
            self.add_argument('-o', '--output', dest='outFile', type=str, default=None, 
                                help='file to write output to (default stdout or display plot to screen)')

        #this is specialized and not generally useful (but used for pies)
        if defaultHatches is not False:
            self.add_argument('-p', '--patterns', dest='hatchString', nargs='*', type=str, default=defaultHatches, 
                                help='string with four characters representing patterns (hatching) for the four wedges \
                                (e.g., default is %s) out of \"/ \\ | - + x X o O . *\" and blank' % defaultHatches)

        if defaultHistogram is not False:
            histogramArgs = self.add_argument_group('ARGUMENTS FOR PLOTTING OF HISTOGRAMS')

            histogramArgs.add_argument('--histogram', dest='histogram', default=False, action='store_true',
                                help='plot a histogram rather than line plot, disregarding or redefining some options')

            histogramArgs.add_argument('--histogram-bins', dest='histogramBins', default=defaultHistogramBins, type=int,
                                help='the number of bins to calculate, if this is a histogram plot (default %d)' % defaultHistogramBins)

            if defaultHistogramKwargs is not False:
                histogramArgs.add_argument('--histogram-kwargs', dest='additionalHistogramKwargs', nargs='*', type=str, default=defaultHistogramKwargs,
                                help='kwargs to be passed on to hist function (default %s), if --histogram is used' % ' '.join(defaultHistogramKwargs))

            histogramArgs.add_argument('--normalize-histogram', dest='normalizeHistogram', default=False, action='store_true',
                                help='normalilze histogram series, so that multiple series are comparable')


def make_figure_and_subplots(nrows, ncols):
    
    #set up the subplots
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True)

    #if only one subplot ax is a single axes object, otherwise a tuple of them
    #rather annoying
    axes = [axes] if not isinstance(axes, collections.Iterable) else np.ravel(axes)

    #if don't do this, aspect ratio of subplots won't be same as a single plot would be, 
    #but then need to scale fonts and points up as well (actually, not doing that anymore)
    oldSize = fig.get_size_inches()
    fig.set_size_inches((oldSize[0] * ncols), oldSize[1] * nrows)
    return fig, axes


def prepare_all_kwargs(options):
    #allKwargDict = collections.defaultdict(lambda: {})
    allKwargDict = {}
    
    allKwargDict['superTitleKwargs'] = prepare_plot_kwargs(options.additionalSuperTitleKwargs,
        [('fontsize', options.superTitleFontsize)])
    
    allKwargDict['subplotKwargs'] = prepare_plot_kwargs(options.additionalSubplotKwargs, [])
    
    allKwargDict['xLabelKwargs'], allKwargDict['yLabelKwargs'] = [ prepare_plot_kwargs(kw,  [('size', options.axisLabelFontsize)]) for kw in [options.additionalXlabelKwargs, options.additionalYlabelKwargs ] ]
    
    #various kwargs
    allKwargDict['titleKwargs'] = prepare_plot_kwargs(options.additionalTitleKwargs, 
            [('fontsize', options.titleFontsize), ('horizontalalignment', options.titleHalign), ('verticalalignment', options.titleValign)])

    allKwargDict['tickKwargs'] = prepare_plot_kwargs(options.additionalTickKwargs,
            [('labelsize', options.axisTickFontsize)])
    
    #these can only be combined like this because no specific values are being passed in
    allKwargDict['xTickKwargs'], allKwargDict['yTickKwargs'], allKwargDict['plotKwargs'], allKwargDict['histogramKwargs'] = [ prepare_plot_kwargs(kw, []) for kw in [ options.additionalXTickKwargs, options.additionalYTickKwargs, options.additionalPlotKwargs, options.additionalHistogramKwargs ] ]

    return allKwargDict
        

def full_plot_routine(options, fileData):
    '''This does a whole bunch of stuff related to making subplots, plotting data, cycling through styles,
    etc.  It is directed by the options returned from the parse_args call to a PlottingArgumentParser.
    There are a number of ways to line up plots/files/series, etc.
    A=axis, F=file, S=series(i.e., data function)
    1. Everything goes onto one axis. (options.onePlot = True)
        A   A   A   A
        F   F   F2  F2
        S   S2  S   S2
    2. Everything goes on its own axis. (options.subplotPerFunction = True)
        A1  A2  A3  A4
        F   F   F2  F2
        S   S2  S   S2
    3. Each file has its own axis, possibly with multiple series (default)
        A1  A1  A2  A2
        F   F   F2  F2
        S   S2  S   S2
    '''

    numFiles = len(fileData)
    numFuncs = len(options.dataColumnFunc)
    numSeries = numFiles * numFuncs

    if options.subplotPerFunction:
        numSubplots = numSeries
    elif options.onePlot:
        numSubplots = 1
    else:
        numSubplots = numFiles

    #################
    #now the plotting

    #make a multiplot (i.e., call pyplot.subplots) even if there is only one plot, to standardize below code
    multiplot = True
    nrows = int(np.ceil(numSubplots / float(options.ncols)))
    ncols = min(options.ncols, numSubplots)
    fig, axes = make_figure_and_subplots(nrows, ncols)

    allKwargDict = prepare_all_kwargs(options)

    fig.subplots_adjust(**allKwargDict['subplotKwargs'])
    
    if options.superTitle:
        plt.suptitle(options.superTitle, **allKwargDict['superTitleKwargs'])

    #if this isn't explicitly set it is figured out by pyplot and applied
    #equally to all subplots
    if options.yRange:
        plt.ylim(options.yRange)

    if options.xRange:
        plt.xlim(options.xRange)

    #think that setting this for axes[0] only works because set to share axes in multiplot
    #xtics = mpl.ticker.FixedLocator([0, 50, 100, 200, 400])
    #axes[0].get_xaxis().set_major_locator(xtics)

    #if location options is s(ingle), center label in overall figure
    #don't rescale fonts in that case
    xFontScaler = 1.0 if multiplot and options.xLabelLocation == 's' else ncols
    yFontScaler = 1.0 if multiplot and options.yLabelLocation == 's' else ncols

    #make a single centered label spanning individual plots
    if options.xLabel and options.xLabelLocation == 's':
        fig.text(allKwargDict['xLabelKwargs'].pop('x'), allKwargDict['xLabelKwargs'].pop('y'), options.xLabel, **allKwargDict['xLabelKwargs'])
    if options.yLabel and options.yLabelLocation == 's':
        fig.text(allKwargDict['yLabelKwargs'].pop('x'), allKwargDict['yLabelKwargs'].pop('y'), options.yLabel, **allKwargDict['yLabelKwargs'])

    #it will be easiest to just set up a bunch of lists of equal length, and then iterate through
    #these to get the proper association between data, functions, axes, etc.
    #########################
    def make_plot_iterator():
        if options.subplotPerFunction:
            #the # of axes should already be the # of plots, and will match up to d1f1 d1f2 d2f1 d2f2 etc
            axesToIterate = [ ax for num, ax in enumerate(axes) if num < numSeries ]
        elif options.onePlot:
            axesToIterate = [ ax for _, ax in zip(range(numSeries), cycle(axes))]
        else:
            #a single axis for each datafile, which will be equivalent to the above unless there are multiple functions per file 
            axesToIterate = np.ravel( [ [ ax for _ in options.dataColumnFunc ] for num, ax in enumerate(axes) if num * len(options.dataColumnFunc) < numSeries ] )

        #can't use ravel here, because lowest level list of lines will also be collapsed
        dataToIterate = []
        for dat in fileData:
            dataToIterate.extend([ dat for _ in xrange(numFuncs) ] )

        #this will be equivalent to the data
        inFilesToIterate = np.ravel( [ [ inf for _ in xrange(numFuncs) ] for inf in options.inFiles ] )

        #functions to iterate should just cycle through the functions for each file
        functionsToIterate = np.ravel( [ [ func for func in options.dataColumnFunc ] for _ in xrange(numFiles) ] )

        #style stuff
        if options.onePlot:
            #each series gets a sucessive style within the one plot
            markersToIterate = np.ravel( [ mark for mark, _ in zip(cycle(options.markers), axesToIterate) ] )
            markerSizesToIterate = np.ravel( [ size for size, _ in zip(cycle(options.markerSizes), axesToIterate) ] )
            lineWidthsToIterate = np.ravel( [ width for width, _ in zip(cycle(options.lineWidth), axesToIterate) ] )
            colorValuesToIterate = np.ravel( [ color for color, _ in zip(cycle(options.colorValues if options.colorValues else options.greyValues), axesToIterate) ] )
            hatchesToIterate = np.ravel( [ hatch for hatch, _ in zip(cycle(options.hatchString if options.hatchString else 'o'), axesToIterate) ] )
        elif options.cycleStylesWithinFiles:
            #this is the default really, sucessive styples within subplots
            markersToIterate = np.ravel( [ [ mark for mark, _ in zip(cycle(options.markers), options.dataColumnFunc) ] for _ in fileData ] )
            markerSizesToIterate = np.ravel( [ [ size for size, _ in zip(cycle(options.markerSizes), options.dataColumnFunc) ] for _ in fileData ] )
            lineWidthsToIterate = np.ravel( [ [ width for width, _ in zip(cycle(options.lineWidth), options.dataColumnFunc) ] for _ in fileData ] )
            colorValuesToIterate = np.ravel( [ [ color for color, _ in zip(cycle(options.colorValues if options.colorValues else options.greyValues), options.dataColumnFunc) ] for _ in fileData ] )
            hatchesToIterate = np.ravel( [ [ hatch for hatch, _ in zip(cycle(options.hatchString if options.hatchString else 'o'), options.dataColumnFunc) ] for _ in fileData ] )
        else:
            #use a single style for all series in a single file
            markersToIterate = np.ravel( [ [ mark for _ in options.dataColumnFunc ] for mark, f in zip( cycle(options.markers), fileData) ] )
            markerSizesToIterate = np.ravel( [ [ size for _ in options.dataColumnFunc ] for mark, f in zip( cycle(options.markerSizes), fileData) ] )
            lineWidthsToIterate = np.ravel( [ [ width for _ in options.dataColumnFunc ] for mark, f in zip( cycle(options.lineWidth), fileData) ] )
            colorValuesToIterate = np.ravel( [ [ color for _ in options.dataColumnFunc ] for mark, f in zip( cycle(options.colorValues if options.colorValues else options.greyValues), fileData) ] )
            hatchesToIterate = np.ravel( [ [ hatch for _ in options.dataColumnFunc ] for hatch, f in  zip(cycle(options.hatchString if options.hatchString else 'o'), fileData) ] )
            
        iterLen = numSeries

        assert(len(dataToIterate) == iterLen)
        assert(len(functionsToIterate) == iterLen)
        assert(len(markersToIterate) == iterLen)
        assert(len(markerSizesToIterate) == iterLen)
        assert(len(lineWidthsToIterate) == iterLen)
        assert(len(colorValuesToIterate) == iterLen)
        assert(len(axesToIterate) == iterLen)
        assert(len(hatchesToIterate) == iterLen)

        return izip(
                inFilesToIterate,
                dataToIterate,
                functionsToIterate,
                markersToIterate,
                markerSizesToIterate,
                lineWidthsToIterate,
                colorValuesToIterate,
                hatchesToIterate,
                axesToIterate)
    #####################

    #cycle through markers, marker sizes, colors, axes, etc., plotting one series per loop
    for num, (inFile, series, function, marker, markerSize, lineWidth, color, hatch, subplot) in enumerate(make_plot_iterator()):
        #turn off extra ticks
        subplot.get_xaxis().tick_bottom()
        subplot.get_yaxis().tick_left()
        
        #if location option is m(arginal), add a normal axis label for plots on left or bottom margin
        if options.yLabel:
            if options.yLabelLocation == 'm' or not multiplot:
                if num % ncols == 0:
                    subplot.set_ylabel(options.yLabel, **allKwargDict['yLabelKwargs'])
        if options.xLabel:
            if options.xLabelLocation == 'm' or not multiplot:
                if num >= numPlots - ncols:
                    subplot.set_xlabel(options.xLabel, **allKwargDict['xLabelKwargs'])

        #set the title - this could get set multiple times for a single file/plot with multiple funcs
        if options.titleFunc and options.titleFunc != 'None':
            subplot.set_title(options.titleFunc(inFile, sep=' '), **allKwargDict['titleKwargs'])
        elif options.titles:
            subplot.set_title(options.titles[num], **allKwargDict['titleKwargs'])
     
        #this is a little goofy, as it will be re-set for each subplot, despite being shared
        if options.xTickFunction:
            #xtickKwargs = {'rotation':90, 'size':'xx-small'}
            tickFunc = eval(options.xTickFunction)
            plt.xticks(range(len(series)), tickFunc(series))

        #change tick properties
        if allKwargDict['tickKwargs']:
            subplot.tick_params(**allKwargDict['tickKwargs'])
        if allKwargDict['xTickKwargs']:
            subplot.tick_params(axis='x', **allKwargDict['xTickKwargs'])
        if allKwargDict['yTickKwargs']:
            subplot.tick_params(axis='y', **allKwargDict['yTickKwargs'])
        
        #do the actual data evaluation and plotting
        dataFunc = eval(function)
        toPlot = dataFunc(series)
        '''to plot just x values, lambda looks like this:
        lambda rows:[float(r[2]) for r in rows]
        i.e., output is [x1, x2, ...]
        for x,y, like this:
        lambda rows:([float(r[2]) for r in rows], [float(r[3]) for r in rows])
        but, the plot function wants [x1, x2, ...], [y1, y2, ...] as the first
        two arguments, hence the * to remove the tuple containing the two lists,
        which I think is required for the lambda'''
        if isinstance(toPlot[0], list):
            if not options.histogram:
                subplot.plot(*dataFunc(series), 
                        marker=marker, 
                        markersize=markerSize, 
                        linewidth=lineWidth, 
                        color=color,
                        mfc=None,
                        mec=color,
                        **allKwargDict['plotKwargs'])
            else:
                #bins = np.linspace(0, 1.0, options.histogramBins)
                bins = np.linspace(options.xRange[0] if options.xRange else 0.0, options.xRange[1] if options.xRange else 1.0, options.histogramBins)
                outN, outBins, patches = subplot.hist(*dataFunc(series), 
                        bins=bins,
                        facecolor=color,
                        normed=options.normalizeHistogram,
                        **allKwargDict['histogramKwargs'])
                for pat in patches:
                    pat.set_hatch(hatch)

        else:
            if not options.histogram:
                subplot.plot(dataFunc(series), 
                        marker=marker, 
                        markersize=markerSize, 
                        linewidth=lineWidth, 
                        color=color,
                        **allKwargDict['plotKwargs'])
            else:
                #bins = np.linspace(0, 1.0, options.histogramBins)
                bins = np.linspace(options.xRange[0] if options.xRange else 0.0, options.xRange[1] if options.xRange else 1.0, options.histogramBins)
                outN, outBins, patches = subplot.hist(dataFunc(series), 
                        bins=bins,
                        facecolor=color,
                        normed=options.normalizeHistogram,
                        **allKwargDict['histogramKwargs'])
                for pat in patches:
                    pat.set_hatch(hatch)
                '''
                subplot.hist(dataFunc(series), 
                        bins=options.histogramBins,
                        facecolor=color,
                        normed=True,
                        **plotKwargs)
                '''
            
        #this plots x=0 and y=0 axis lines
        if options.drawAxesAtZero:
            subplot.axhline(xmin=-sys.maxint - 1, xmax=sys.maxint)
            subplot.axvline(ymin=-sys.maxint - 1, ymax=sys.maxint)

    #return the current figure, so more manipulations could potentially be done on it
    return plt.gcf()

    if options.outFile:
        plt.savefig(options.outFile, transparent=True, bbox_inches='tight')
    else:
        plt.show()


