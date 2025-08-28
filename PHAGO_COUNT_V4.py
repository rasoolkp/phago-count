from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate.features.edges import EdgeTargetAnalyzer

from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.util import LogRecorder
from fiji.plugin.trackmate.util import TMUtils

from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTracker
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking.jaqaman import LAPUtils

from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate.stardist import StarDistDetectorFactory
from fiji.plugin.trackmate.stardist import StarDistCustomDetectorFactory #for custom models
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.action import CaptureOverlayAction
from fiji.plugin.trackmate.action import IJRoiExporter

import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
import fiji.plugin.trackmate.features.track.TrackSpeedStatisticsAnalyzer as TrackSpeedStatisticsAnalyzer
import fiji.plugin.trackmate.features.track.TrackDurationAnalyzer as TrackDurationAnalyzer
import fiji.plugin.trackmate.features.edges.EdgeTargetAnalyzer as EdgeTargetAnalyzer
import fiji.plugin.trackmate.features.edges.EdgeSpeedAnalyzer as EdgeSpeedAnalyzer
import fiji.plugin.trackmate.features.edges.EdgeTimeLocationAnalyzer as EdgeTimeLocationAnalyzer
import fiji.plugin.trackmate.features.edges.DirectionalChangeAnalyzer as DirectionalChangeAnalyzer
#https://imagej.net/plugins/trackmate/scripting#Export_spot.2C_edge_and_track_numerical_features_after_tracking

from ij.plugin.frame import RoiManager
from ij.gui import PointRoi
#import components to work with results tables
from ij import WindowManager
from ij.text import TextWindow
from ij.measure import ResultsTable
from ij import IJ
from operator import itemgetter # allows easier sorting of list of lists
from datetime import datetime as dt
import sys
import time 
import os
import glob
import csv
import math

reload(sys)
sys.setdefaultencoding ('utf-8')

#image files to analyze
input_dir = IJ.getDirectory("Choose a directory containing image files")
if input_dir is None:
    raise Exception("No directory selected!")

file_paths = glob.glob(os.path.join(input_dir, "*.tif"))
if not file_paths:
    raise Exception("No .tif files found in the selected directory!")

# creating object for roiManager and opening it
RM = RoiManager()
rm = RM.getRoiManager()

show_output = True
closeToEdge = 10
closeSpots = True
outAppendix = '_vx'


# Get user input via dialog
from javax.swing import JPanel, JLabel, JTextField, JOptionPane, BoxLayout, JComboBox

# Unified input panel
panel = JPanel()
panel.setLayout(BoxLayout(panel, BoxLayout.Y_AXIS))

fields = {}

labels = [
    ("Target cell channel:", "b_cell_channel"),
    ("Macrophage channel:", "macrophage_channel"),
    ("Start frame for target cell analysis:", "bcell_start"),
    ("End frame for target cell analysis:", "bcell_end"),
    ("Start frame for macrophage analysis:", "macrophage_start"),
    ("End frame for macrophage analysis:", "macrophage_end"),
]

# Default values to prefill fields
defaults = {
    "b_cell_channel": "2",
    "macrophage_channel": "1",
    "bcell_start": "0",
    "bcell_end": "49",
    "macrophage_start": "0",
    "macrophage_end": "49"
}

# Add text fields with prefilled values
for text, key in labels:
    label = JLabel(text)
    field = JTextField(defaults.get(key, ""), 5)  # use default if defined
    panel.add(label)
    panel.add(field)
    fields[key] = field

# Add dropdown for model type
model_label = JLabel("Macrophage StarDist model type:")
model_dropdown = JComboBox(["default", "custom"])
panel.add(model_label)
panel.add(model_dropdown)

# Show dialog
result = JOptionPane.showConfirmDialog(None, panel, "Enter analysis parameters", JOptionPane.OK_CANCEL_OPTION)

if result != JOptionPane.OK_OPTION:
    raise Exception("User cancelled input.")

# Parse and validate input
try:
    b_cell_channel = int(fields["b_cell_channel"].getText())
    macrophage_channel = int(fields["macrophage_channel"].getText())
    bcell_start = int(fields["bcell_start"].getText())
    bcell_end = int(fields["bcell_end"].getText())
    macrophage_start = int(fields["macrophage_start"].getText())
    macrophage_end = int(fields["macrophage_end"].getText())
    macrophage_model = model_dropdown.getSelectedItem()  # "default" or "custom"
except Exception as e:
    raise Exception("Invalid input: " + str(e))

# Parse and validate input
try:
    b_cell_channel = int(fields["b_cell_channel"].getText())
    macrophage_channel = int(fields["macrophage_channel"].getText())
    bcell_start = int(fields["bcell_start"].getText())
    bcell_end = int(fields["bcell_end"].getText())
    macrophage_start = int(fields["macrophage_start"].getText())
    macrophage_end = int(fields["macrophage_end"].getText())
except Exception as e:
    raise Exception("Invalid input: " + str(e))

def run_trackmate_on_channel(file_path, channel, model_type, subthreshold_run=False):
    pixelSize= 0.16
    firstAreaNorm = 0.7 # track normalized area of first frame after division... maybe a little stringent
    lastAreaNorm = 0.8 # track normalized area of last frame before division
    firstAreaPx = 100
    lastAreaPx = 140
    firstNormaMean = 1.2
    lastNormMean = 1.0
    firstTotInt = 0.8 #prob better to ref rel to parent
    lastTotInt = 1.2
    minAvgArea = 0.9
    childTotIntFraction = 0.7
    maxDeltaNormTotIntCh1 = 0.5
    maxDeltaStartNormTotIntCh1 = 0.5
    maxBranchesNormArea = 0.8
    
    t0= time.time()*1000
    rm.reset() #reset ROI manager
    
    #open and process image
    print("\n=== Processing file: {} | Channel: {} | Model: {} ===".format(file_path, channel, model_type))
    imp  = IJ.openImage(file_path)
    cal = imp.getCalibration()
    channels = imp.getNChannels()
    slices = imp.getNSlices()
    frames = imp.getNFrames()
    height = imp.getHeight()
    width = imp.getWidth()
    t1 = time.time()*1000
    print("Channels: "+str(channels)+' slices:'+str(slices)+' frames:'+str(frames) + ' width: '+str(width)+'px height: '+str(height)+'px')
    dt1 = t1-t0
    print("Opened image in: "+str(round(dt1/1000,3))+" s")
    
    # logger -> saving content to xml file
    logger = LogRecorder( Logger.VOID_LOGGER )
    logger.log( 'TrackMate-StarDist analysis script\n' )
    dt_string = dt.now().strftime("%d/%m/%Y %H:%M:%S")
    logger.log( dt_string + '\n\n' )
    
    #instantiate roixporter built into trackmate gui
    roiExporter = IJRoiExporter(imp, logger)
    
    #preparing settings object
    settings = Settings(imp)
    setup = settings.toStringImageInfo
    
    #configure stardist detector
    if channel == b_cell_channel:
        settings.tstart = bcell_start
        settings.tend = bcell_end
        print("Using Target cell frame range: {}–{}".format(bcell_start, bcell_end))
    elif channel == macrophage_channel:
        settings.tstart = macrophage_start
        settings.tend = macrophage_end
        print("Using macrophage frame range: {}–{}".format(macrophage_start, macrophage_end))
    else:
        settings.tstart = 0
        settings.tend = imp.getNFrames() - 1
        print("Using default frame range: full video")
        
    if model_type == 'default':
        from fiji.plugin.trackmate.stardist import StarDistDetectorFactory
        detectorFactory = StarDistDetectorFactory()
        settings.detectorFactory = detectorFactory
        settings.detectorSettings = {'TARGET_CHANNEL': channel}
        if channel == b_cell_channel:
            settings.addSpotFilter(FeatureFilter('AREA', 40.0, True))
        elif channel == macrophage_channel:
            settings.addSpotFilter(FeatureFilter('AREA', 150.0, True))
            
    elif model_type == 'subthreshold':
        from fiji.plugin.trackmate.stardist import StarDistDetectorFactory
        detectorFactory = StarDistDetectorFactory()
        settings.detectorFactory = detectorFactory
        settings.detectorSettings = {'TARGET_CHANNEL': channel}
        if channel == b_cell_channel:
            settings.addSpotFilter(FeatureFilter('AREA', 40.0, False))
            
    elif model_type == 'custom':
        from fiji.plugin.trackmate.stardist import StarDistCustomDetectorFactory
        detectorFactory = StarDistCustomDetectorFactory()
        settings.detectorFactory = detectorFactory
        settings.detectorSettings = {
            'TARGET_CHANNEL': macrophage_channel,
            'MODEL_FILEPATH': "D:/Rasool/Stardist models/pCytes105_300ep_128rays_2grd_TF1.14/TF_SavedModel.zip",
            'SCORE_THRESHOLD': 0.1,
            'OVERLAP_THRESHOLD': 0.2
        }
        settings.addSpotFilter(FeatureFilter('AREA', 150.0, True))

    #tracker configuration
    settings.trackerFactory = SparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
    settings.trackerSettings['LINKING_MAX_DISTANCE'] = 15.0
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 12.0
    settings.trackerSettings['MAX_FRAME_GAP'] = 0
    settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
    settings.trackerSettings['ALLOW_TRACK_MERGING']  = False
    #settings.trackerSettings['SPLITTING_MAX_DISTANCE'] = 12.0
    settings.initialSpotFilterValue = -1.

    #tracker filters
    settings.addAllAnalyzers()
    settings.addTrackAnalyzer(TrackSpeedStatisticsAnalyzer())
    settings.addTrackAnalyzer(TrackDurationAnalyzer())
    settings.addEdgeAnalyzer(EdgeTargetAnalyzer())
    settings.addEdgeAnalyzer(EdgeSpeedAnalyzer())
    settings.addEdgeAnalyzer(EdgeTimeLocationAnalyzer())
    
    print "Spot filters added = ", settings.getSpotFilters()
    
    #instanstiate trackmate
    print ("Starting Trackmate")
    trackmate = TrackMate(settings)
    trackmate.computeSpotFeatures(True)
    trackmate.computeTrackFeatures(True)
    trackmate.getModel().setLogger(logger)
        
    #processing check
    if not trackmate.checkInput():
        print(str(trackmate.getErrorMessage()))
        return
    
    if not trackmate.process():
        print(str(trackmate.getErrorMessage()))
        return
        
    #Save results 
    from java.io import File
    
    # Build file paths
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    output_dir = os.path.dirname(file_path)
    xml_name = "{}_vx{}_TrackMate.xml".format(base_name, channel)
    xml_file = File(os.path.join(output_dir, xml_name))  # Java File object
    
    # Save XML
    writer = TmXmlWriter(xml_file, logger)
    writer.appendLog(logger.toString())
    writer.appendModel(trackmate.getModel())
    writer.appendSettings(trackmate.getSettings())
    writer.writeToFile()
    
    print("Results saved to: " + xml_file.toString() + '\n')
        
    #display results
    if show_output:
        model = trackmate.getModel()
        selectionModel = SelectionModel (model)
        ds = DisplaySettings()
        ds = DisplaySettingsIO.readUserDefault()
        ds.spotDisplayedAsRoi = True
        displayer =  HyperStackDisplayer( model, selectionModel, imp, ds )
        displayer.render()
        displayer.refresh()
            
    # capture overlay - RGB file
    image = trackmate.getSettings().imp
    
    t2 = time.time()*1000
    dt2 = (t2-t1)/1000
    spf = dt2/frames #seconds per frame
    spmp =1000000*spf/(width*height) #seconds/megapixel
    print("TrackMate took "+str(round(dt2/60,2))+" min at "+str(round(spf,3))+" s/frame "+  str(round(spmp,3)) + " s/megapixel")# capture overlay - RGB file
    
    #capture = CaptureOverlayAction.capture(image, -1, imp.getNFrames(), logger)
    #capture.setTitle("TracksOverlay")
    #capture.show()
    t2 = time.time()*1000
    dt2 = (t2-t1)/1000
    spf = dt2/frames #seconds per frame
    spmp =1000000*spf/(width*height) #seconds/megapixel
    print("TrackMate took "+str(round(dt2/60,2))+" min at "+str(round(spf,3))+" s/frame "+  str(round(spmp,3)) + " s/megapixel")
    
    #-----------------------------
    # Store Edge values 
    #-----------------------------
    
    #instantiate a table to store spots values
    MeasureTable = WindowManager.getWindow("Results")
    if MeasureTable == None:
        MeasureTable = ResultsTable()
    else:
        MeasureTable = WindowManager.getWindow("Results")
        MeasureTable = MeasureTable.getTextPanel().getOrCreateResultsTable()
    MeasureTable.show("Results")
            
    # The feature model, that stores edge and track features.
    fm = model.getFeatureModel()
    tm = model.getTrackModel()
    spaceUnit = model.getSpaceUnits() # extract the distance units
    timeUnit = model.getTimeUnits() # extract the time units
    
    edgeHolder = []
    
    edgesHeader = ["Track_ID","SPOT_SOURCE_ID", "SPOT_TARGET_ID","SOURCE_FRAME", "TARGET_FRAME",
    "EDGE_COST","DIRECTIONAL_CHANGE_RATE","SPEED","DISPLACEMENT","EDGE_TIME",
                    "EDGE_X_LOCATION","EDGE_Y_LOCATION"]
    
    #print(' '.join([str(elem) for elem in edgesHeader]))           
    for id in model.getTrackModel().trackIDs(True):     # Iterate over all the tracks that are visible.
    
        #print( "Exporting track " + str(id) )
        # Fetch the track feature from the feature model.
        v = fm.getTrackFeature(id, 'TRACK_MEAN_SPEED')
        model.getLogger().log('')
        model.getLogger().log('Track (ID), ' + str(id) + ': mean velocity = ' + str(v) + ' ' + model.getSpaceUnits() + '/' + model.getTimeUnits())
     
        # Get all the spots of the current track.
        edges = tm.trackEdges(id)
        for edge in edges:
            ed = fm.getEdgeFeature (edge, EdgeSpeedAnalyzer.DISPLACEMENT )
            et = fm.getEdgeFeature( edge, EdgeTimeLocationAnalyzer.TIME)
            ev = fm.getEdgeFeature( edge, EdgeSpeedAnalyzer.SPEED )
            startSpot = tm.getEdgeSource(edge)
            endSpot = tm.getEdgeTarget (edge)
            #startID = startSpot.ID()
            startID = tm.getEdgeSource(edge).ID()
            endID = endSpot.ID()
            startFrame = startSpot.getFeature ('FRAME')
            endFrame = endSpot.getFeature ('FRAME')
            
            #model.getLogger().log(str (id) +
            #print(str (id) + ',' + str(startID)+ ',' + str(endID) + ',' + str(startFrame) +
            #',' + str (endFrame) ',' + str(ev) + ',' + str(ed) +  ',' + str(et) )
            row = []
            row.append(id)
            row.append(tm.getEdgeSource(edge).ID())
            row.append(tm.getEdgeTarget(edge).ID())
            row.append(tm.getEdgeSource(edge).getFeature ('FRAME'))
            row.append(tm.getEdgeTarget(edge).getFeature ('FRAME'))
            row.append(fm.getEdgeFeature(edge, EdgeTargetAnalyzer.EDGE_COST))
            row.append(fm.getEdgeFeature(edge, DirectionalChangeAnalyzer.DIRECTIONAL_CHANGE_RATE))
            row.append(fm.getEdgeFeature(edge, EdgeSpeedAnalyzer.SPEED ))
            row.append(fm.getEdgeFeature(edge, EdgeSpeedAnalyzer.DISPLACEMENT))
            row.append(fm.getEdgeFeature(edge, EdgeTimeLocationAnalyzer.TIME))
            row.append(fm.getEdgeFeature(edge, EdgeTimeLocationAnalyzer.X_LOCATION))
            row.append(fm.getEdgeFeature(edge, EdgeTimeLocationAnalyzer.Y_LOCATION))
            #row.append(fm.getEdgeFeature( edge, EdgeTimeLocationAnalyzer.Z_LOCATION))
            edgeHolder.append(row)
            
    edgeHolder = sorted(edgeHolder, key=itemgetter(1))
    edgeHolder = sorted(edgeHolder, key=itemgetter(0))
    
    for i in range(len(edgeHolder)):
        MeasureTable.incrementCounter()
        for j in range(len(edgeHolder[i])):
            MeasureTable.addValue(edgesHeader[j], edgeHolder[i][j])
    
    SPOT_SOURCE_IDs = map(itemgetter(1), edgeHolder)
    SPOT_TARGET_IDs = map(itemgetter(2), edgeHolder)
    
    MeasureTable.show("Results")
    edges_path = file_path.replace(".tif","_trackMateEdges"+outAppendix+str(channel)+".csv")            #edges_path = file_path.replace(".nd2","_trackMateEdges"+outAppendix+".csv")
    IJ.saveAs("Results", edges_path); #Needs to be listed a "Result' to be saved as csv. Alternat table names can be saved as .txt.
    IJ.renameResults("Results","Edges") 
    IJ.selectWindow("Edges"); 
    IJ.run("Close");
    
    t3 = time.time()*1000
    dt3 = t3-t2
    print("Assessing and saving edges took "+str(round(dt3/1000,3))+" s")
    
    #-----------------------------
    # Store Spot values 
    #-----------------------------
    
    #instantiate a 2nd table to store spots values
    MeasureTable = WindowManager.getWindow("Results")
    if MeasureTable == None:
        MeasureTable = ResultsTable()
    else:
        MeasureTable = WindowManager.getWindow("Results")
        MeasureTable = MeasureTable.getTextPanel().getOrCreateResultsTable()
    MeasureTable.show("Results")                    
        
    spotFeatures = ['POSITION_X','POSITION_Y','POSITION_T','FRAME','QUALITY','RADIUS','SNR_CH1','MEAN_INTENSITY_CH1',
    'MIN_INTENSITY_CH1','MAX_INTENSITY_CH1','TOTAL_INTENSITY_CH1','AREA','ELLIPSE_MAJOR','ELLIPSE_MINOR','ELLIPSE_THETA','ELLIPSE_ASPECTRATIO',
    'PERIMETER','CIRCULARITY','SOLIDITY']
    ch2Features = ['SNR_CH2','MEAN_INTENSITY_CH2',  'MIN_INTENSITY_CH2','MAX_INTENSITY_CH2','TOTAL_INTENSITY_CH2']
    ch3Features = ['SNR_CH3','MEAN_INTENSITY_CH3',  'MIN_INTENSITY_CH3','MAX_INTENSITY_CH3','TOTAL_INTENSITY_CH3']
    #Add extra features to be logged/saved for multi-channel images
    if channels >1 :
        spotFeatures.extend(ch2Features)
    if channels >2 :
        spotFeatures.extend(ch3Features)
        
    
    #create a list to hold spot features used to assess branch names
    spotsHolder = []
    
    for id in model.getTrackModel().trackIDs(True):
        track = model.getTrackModel().trackSpots(id)
        #print( "Working on track " + str(id) + " with spots " +str(track)) 
        nSpots = 0 #keep track of number of rows in Results Table
        
        for spot in track:
            
            roiExporter.export(spot)  #export coordinates of spot to ROI manager
            nSpots +=1
    
    t3b = time.time()*1000
    dt3b = t3b-t2
    print("Time to export ROIs "+str(round(dt3b/1000,3))+" s")
            
    for id in model.getTrackModel().trackIDs(True):
        track = model.getTrackModel().trackSpots(id)
        #print( "Working on track " + str(id) + " with spots " +str(track)) 
        nSpots = 0 #keep track of number of rows in Results Table
        
        for spot in track:
            
    #           roiExporter.export(spot)  #export coordinates of spot to ROI manager
    #           nSpots +=1
            
            track = model.getTrackModel().trackSpots(id)
            #print( "Exporting spot " + str(id) + "."+str(spot))
            sid = spot.ID()
            # Fetch spot features directly from spot.
            # Note that for spots the feature values are not stored in the FeatureModel
            # object, but in the Spot object directly. This is an exception; for tracks
            # and edges, you have to query the feature model.
            #model.getLogger().log('\tspot ID = ' + str(sid) + ': x='+str(x)+', y='+str(y)+', t='+str(t)+', q='+str(q) + ', snr='+str(snr) + ', mean = ' + str(mean))
            
            # Append to custom Measurement table or create it if non existing
            MeasureTable.incrementCounter()
    #           MeasureTable.addValue('LABEL', 'ID'+str(sid))
    #           MeasureTable.addValue('ID', str(sid))
            MeasureTable.addValue('LABEL', 'ID'+'{:0>5}'.format(str(sid)))
            MeasureTable.addValue('ID', '{:0>5}'.format(str(sid)))
    
            MeasureTable.addValue('TRACK_ID', '{:0>5}'.format(str(id)))
            
            for i in spotFeatures:
                MeasureTable.addValue(i, spot.getFeature(i))
        
            #parse the edges list of source and target spots to find child(ren) and parent(s) of curent spot
            childIndex = [i for i in range(len(SPOT_SOURCE_IDs)) if SPOT_SOURCE_IDs[i]==sid] #short hand to check if edge source matches current sid
            parentIndex = [i for i in range(len(SPOT_TARGET_IDs)) if SPOT_TARGET_IDs[i]==sid]
            
            if len(childIndex) > 0:
                childID = [SPOT_TARGET_IDs[i] for i in childIndex ]
            else:
                childID = [-1]
            if len(parentIndex) > 0:
                parentID =  [SPOT_SOURCE_IDs[i] for i in parentIndex ]
            else:
                parentID = [-1] 
            #Add values to results table
            MeasureTable.addValue('nCHILDREN', len(childIndex))
            MeasureTable.addValue('nPARENTS', len(parentIndex))
            MeasureTable.addValue('CHILDREN', ",".join([str(i) for i in childID]))
            MeasureTable.addValue('PARENTS', ",".join([str(i) for i in parentID]))
    
            #alternate approach is to store results in a list of lists, then print all to results table
            row = []
            row.append(sid)
            row.append(id)
            for i in spotFeatures:
                row.append(spot.getFeature(i))
            row.append(len(childIndex)) 
            row.append(len(parentIndex))    
            row.append(",".join([str(i) for i in childID])) 
            row.append(",".join([str(i) for i in parentID]))    
            spotsHolder.append(row) 
            #print(row)
    
    #Sort the Results table by track and spotID
    MeasureTable.sort("ID") 
    MeasureTable.sort("TRACK_ID") 
    MeasureTable.show("Results") 
    
    #Sort spotsHolder by track and spotID to match the Results table
    spotsHolder = sorted(spotsHolder, key=itemgetter(0)) #sort by Spot ID or 'ID'
    spotsHolder = sorted(spotsHolder, key=itemgetter(1)) #sort by 'TRACK_ID'
    
    #Pull values from the results table
    idList = map(itemgetter(0), spotsHolder)
    trackList = map(itemgetter(1), spotsHolder)
    frameList = map(itemgetter(spotFeatures.index('FRAME')+2), spotsHolder)
    xList = map(itemgetter(spotFeatures.index('POSITION_X')+2), spotsHolder)
    yList = map(itemgetter(spotFeatures.index('POSITION_Y')+2), spotsHolder)
    ch1MeanList = map(itemgetter(spotFeatures.index('MEAN_INTENSITY_CH1')+2), spotsHolder)
    ch1IntDenList = map(itemgetter(spotFeatures.index('TOTAL_INTENSITY_CH1')+2), spotsHolder)
    areaList = map(itemgetter(spotFeatures.index('AREA')+2), spotsHolder)
    nChildList = map(itemgetter(len(spotFeatures)+2), spotsHolder)
    nParentList = map(itemgetter(len(spotFeatures)+3), spotsHolder)
    childrenList = map(itemgetter(len(spotFeatures)+4), spotsHolder)
    parentsList = map(itemgetter(len(spotFeatures)+5), spotsHolder)
    branchList = [0]*len(idList)
    branchLength = [0]*len(idList)
    nSiblings = [0]*len(idList)
    branchTypes = [0]*len(idList)
    
    if channels >1:
        ch2MeanList = map(itemgetter(spotFeatures.index('MEAN_INTENSITY_CH2')+2), spotsHolder)
        ch2IntDenList = map(itemgetter(spotFeatures.index('TOTAL_INTENSITY_CH2')+2), spotsHolder)
    if channels >2:
        ch3MeanList = map(itemgetter(spotFeatures.index('MEAN_INTENSITY_CH3')+2), spotsHolder)
        ch3IntDenList = map(itemgetter(spotFeatures.index('TOTAL_INTENSITY_CH3')+2), spotsHolder)
    
    MeasureTable.updateResults()
    
    
    t4 = time.time()*1000
    dt4 = (t4-t3b)/60000
    nROIs = rm.getCount()
    
    rate = nROIs / ((dt4 + 1e-6) * 60)
    print("Measuring and saving {} spots took {:.2f} min @ {:.2f} ROIs/s".format(nROIs, dt4, rate))
    
    
    t5 = time.time()*1000
    dt5 = (t5-t4)/60000
    print("Saving "+str(nROIs)+" ROIs took "+str(round(dt5,2))+" min")          
    
    t6 = time.time()*1000
    dt6 = t6-t5
    print("Summarizing tracks took "+str(round(dt6/1000,3))+" s")                                               
                                                                                            
    
    t7 = time.time()*1000
    dt7 = t7-t6
    print("Summarizing branches took "+str(round(dt7/1000,3))+" s")                                             
    
    
    MeasureTable.updateResults()
    t8 = time.time()*1000
    dt8 = t8-t7
    print("Analyzing special cases took "+str(round(dt8/1000,3))+" s")                                              
    
    if model_type == 'subthreshold':
        spots_path = file_path.replace(".tif","_subthreshold_trackMateSpots"+outAppendix+str(channel)+".csv")
        roi_path = file_path.replace(".tif","_subthreshold_trackMateROIs"+outAppendix+str(channel)+".zip")
    else:
        spots_path = file_path.replace(".tif","_trackMateSpots"+outAppendix+str(channel)+".csv")
        roi_path = file_path.replace(".tif","_trackMateROIs"+outAppendix+str(channel)+".zip")
    
    MeasureTable.updateResults()
    MeasureTable.show("Results")
    
    IJ.saveAs("Results", spots_path);
    IJ.renameResults("Spots");
    
    rm.runCommand("Save", roi_path);
   
    def get_mac_roi_bounds(mac_roi, padding = 20):
        bounds = mac_roi.getBounds()
        x_min = bounds.x - padding
        y_min = bounds.y - padding
        x_max = bounds.x + bounds.width + padding
        y_max = bounds.y + bounds.height + padding
        return x_min, x_max, y_min, y_max
    
    def circle_roi_overlap(bx, by, br, mac_roi, required_hits = 24):
        #checking if the centered bx,by b cell with radius is in the macrophage ROI
        import math
        hits = 0
        br = br * 0.5 
        for angle in range(0, 360, 15): #using 8 perimeter points by 45 degrees to check the B cell perimeter
            rad = math.radians(angle)
            px = bx + br * math.cos(rad)
            py = by + br * math.sin(rad)
            if mac_roi.contains(int(px), int(py)):
                hits += 1
            if hits >= required_hits:
                return True
        return False

    def circle_roi_inside_fraction(bx, by, br, mac_roi, step_deg=15):
        import math
        total = 0
        hits = 0
        br = br * 0.5 
        for angle in range(0, 360, step_deg): 
            rad = math.radians(angle)
            px = bx + br * math.cos(rad)
            py = by + br * math.sin(rad)
            total += 1
            if mac_roi.contains(int(px), int(py)):
                hits += 1
        return (hits / float(total)) if total else 0.0

    
    if channel == macrophage_channel: #only run this if we're now in macrophage channel to analyze
        phagocytic_frame_threshold = 10
        slip_off_threshold = 2
        slip_outside_fraction_threshold = 0.6
        mac_spots_path = spots_path
        bc_spots_path = spots_path.replace("_vx{}".format(macrophage_channel), "_vx{}".format(b_cell_channel))
        sub_bc_spots_path = file_path.replace(".tif", "_subthreshold_trackMateSpots" + outAppendix + str(b_cell_channel) + ".csv")
        sub_bc_roi_path = file_path.replace(".tif", "_subthreshold_trackMateROIs" + outAppendix + str(b_cell_channel) + ".zip")
        
        if not os.path.exists(mac_spots_path) or not os.path.exists(bc_spots_path):
            print("Missing input files for spatial overlap check.")
        else:
            import csv
            from collections import defaultdict
    
            def load_spots(path):
                with open(path, 'r') as f:
                    reader = csv.DictReader(f)
                    return [row for row in reader]
    
            mac_spots = load_spots(mac_spots_path)
            print("Loaded {} macrophage spot rows.".format(len(mac_spots)))
            if len(mac_spots) > 0:
                print("First mac_spot example:", mac_spots[0])
            bc_spots = load_spots(bc_spots_path)
            
            #grabbing the b cell ROIs here
            bc_roi_zip_path = bc_spots_path.replace("_trackMateSpots", "_trackMateROIs").replace(".csv", ".zip")
            print("Loading B cell ROIs from:", bc_roi_zip_path)
            
            bcell_track_rois = defaultdict(list)
            bc_id_to_track_frame = {}
            #mapping and making library of spots
            for row in bc_spots:
                try:
                    spot_id = int(row['ID'])
                    track_id = int(row['TRACK_ID'])
                    frame = int(row['FRAME'])
                    bc_id_to_track_frame[spot_id] = (track_id, frame)
                except:
                    continue
            #grab the spot IDs and take out the ID to keep the numbers only from the name        
            rm.runCommand("Open", bc_roi_zip_path)
            for i in range(rm.getCount()):
                name = rm.getName(i)
                if name.startswith("ID"):
                    try:
                        spot_id = int(name.replace("ID", ""))
                        roi = rm.getRoi(i)
                        if spot_id in bc_id_to_track_frame:
                            track_id, frame = bc_id_to_track_frame[spot_id]
                            bcell_track_rois[track_id].append((frame, roi))
                    except:
                        continue
            print("Loaded {} B cell ROIs.".format(len(bcell_track_rois)))
            
            #grabbing the mac ROIs here
            roi_zip_path = mac_spots_path.replace("_trackMateSpots", "_trackMateROIs")
            roi_zip_path = roi_zip_path.replace(".csv", ".zip")
            print("Loading ROIs from:", roi_zip_path)
            #doing same thing to save the spot IDs for the mac ROIs
            rm.reset()
            IJ.run("ROI Manager...", "")
            rm.runCommand("Open", roi_zip_path)
            roi_dict ={}
            for i in range(rm.getCount()):
                name = rm.getName(i)
                if name.startswith("ID"):
                    try:
                        spot_id = int(name.replace("ID", ""))
                        roi = rm.getRoi(i)
                        roi_dict[spot_id] = roi
                    except:
                        continue
            
            print("Loaded ROI count:", rm.getCount())
#            for i in range(min(5, rm.getCount())):
                #print("ROI name at index {}: {}".format(i, rm.getName(i)))
                
            def mean_std(values):
                if not values:
                    return (0.0, 0.0)
                n = len(values)
                mean = sum(values) / float (n)
                var = sum((x - mean) ** 2 for x in values) / float(n)
                std = var ** 0.5
                return mean, std
            
            target_values_f0 = []
            channel_key = "MEAN_INTENSITY_CH{}".format(b_cell_channel)
            for row in mac_spots:
                try:
                    frame = int(row['FRAME'])
                    target = float(row.get(channel_key, 0.0))
                    if frame == 0:
                        target_values_f0.append(target)
                except:
                    continue
            
            mean_f0, std_f0 = mean_std(target_values_f0)
            mac_threshold = mean_f0 + 2 * std_f0
            print("Empty macrophage (mean + 2*std): {:.2f}".format(mac_threshold))
            
            mac_tracks = defaultdict(list)
            #here we're grouping all the corresponding spots into the right macrophage tracks
            for row in mac_spots:
                try:
                    track_id = int(row['TRACK_ID'])
                    frame = int(row['FRAME'])
                    x = float(row['POSITION_X'])
                    y = float(row['POSITION_Y'])
                    mac_tracks[track_id].append((frame, x, y, row))
                except Exception as e:
                    print("Skipping row due to error:", e)
            
            print("Grouped {} tracks from mac_spots.".format(len(mac_tracks)))
            #now you group your b cells tracks too 
            bc_by_frame = defaultdict(list)
            for row in bc_spots:
                try:
                    frame = int(row['FRAME'])
                    x = float(row['POSITION_X'])
                    y = float(row['POSITION_Y'])
                    br = float(row['RADIUS'])
                    track_id = int(row['TRACK_ID'])
                    bc_by_frame[frame].append((x, y, br, track_id))
                except:
                    continue
            
            # === Grabbing the subthreshold B cell ROIs for trogocytosis ===
            sub_bc_spots_path = file_path.replace(".tif", "_subthreshold_trackMateSpots" + outAppendix + str(b_cell_channel) + ".csv")
            sub_bc_roi_path = file_path.replace(".tif", "_subthreshold_trackMateROIs" + outAppendix + str(b_cell_channel) + ".zip")
            
            print("Loading subthreshold target cell ROIs from:", sub_bc_roi_path)
            
            sub_bcell_track_rois = {}
            sub_bc_id_to_track_frame = {}
            
            trogocytic_tracks = []
            trogocytic_rows = []
            
            if not os.path.exists(sub_bc_spots_path) or not os.path.exists(sub_bc_roi_path):
                print("Missing subthreshold CSV or ROI file — skipping trogocytosis detection.")
            else:
                sub_bcell_track_rois = defaultdict(list)
                # Load mapping from spot ID → (track, frame)
                with open(sub_bc_spots_path, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        try:
                            spot_id = int(row['ID'])
                            track_id = int(row['TRACK_ID'])
                            frame = int(row['FRAME'])
                            sub_bc_id_to_track_frame[spot_id] = (track_id, frame)
                        except:
                            continue
            
                # Load ROIs and group them into sub_bcell_track_rois
                rm.reset()
                rm.runCommand("Open", sub_bc_roi_path)
            
                for i in range(rm.getCount()):
                    name = rm.getName(i)
                    if name.startswith("ID"):
                        try:
                            spot_id = int(name.replace("ID", ""))
                            roi = rm.getRoi(i)
                            if spot_id in sub_bc_id_to_track_frame:
                                track_id, frame = sub_bc_id_to_track_frame[spot_id]
                                bx = roi.getBounds().getCenterX()
                                by = roi.getBounds().getCenterY()
                                br = max(roi.getBounds().width, roi.getBounds().height) / 2.0
                                sub_bcell_track_rois[track_id].append((frame, bx, by, br, roi))
                        except:
                            continue
            
                print("Loaded {} subthreshold Target cell ROIs grouped into {} tracks.".format(rm.getCount(), len(sub_bcell_track_rois)))
            
            
    
            #make the output lists 
            phagocytic_tracks = set()
            phagocytic_rows = []
            multi_phagocytic_tracks = set()
            multi_phagocytic_rows = []
            edge_case_tracks = []
            edge_case_rows = []
            total_tracks = 0
            stalling_cup_tracks = []
            stalling_cup_rows = []
            trogocytic_tracks = []
            trogocytic_rows = []
            engulfed_bcell_rows_all = []            # NEW: accumulate engulfed B-cell rows across all macs
            engulfed_bcell_track_ids_all = set()    # NEW: accumulate the set of engulfed B-cell track IDs
            phago_slip_audit = []
            total_tracks = 0  
            
            for i, (track_id, mac_entries) in enumerate(mac_tracks.items()):
                if i >= 5:
                    break
                frame_set = set([entry[0] for entry in mac_entries])
                #print("Track {} has {} unique frames.".format(track_id, len(frame_set)))
                
            for track_id, mac_entries in mac_tracks.items():

                # overlaps_by_bcell now reset inside each mac track loop
                frame_set = set([entry[0] for entry in mac_entries])
                if len(frame_set) >= 15:
                    total_tracks += 1

                overlaps_by_bcell = defaultdict(list)
                mac_roi_by_frame = {}
                for frame, mx, my, row in mac_entries:
                    try:
                        mac_spot_id = int(row['ID'])
                    except:
                        continue
                    mac_roi = roi_dict.get(mac_spot_id, None)
                    if mac_roi is not None:
                        mac_roi_by_frame[int(frame)] = mac_roi
                overlap_frames = []
                for frame, mx, my, row in mac_entries: #look through all macrophage spots in the track
                    mac_spot_id = int(row['ID'])
                    mac_roi = roi_dict.get(mac_spot_id, None) #grab the ROI of this macrophage spot
                    if mac_roi is None:
                        print("Missing ROI for ID:", mac_spot_id)
                        continue
                    if frame not in bc_by_frame: #if no b cells in the frame, skip it 
                        continue
                    found = False
                    x_min, x_max, y_min, y_max = get_mac_roi_bounds(mac_roi)
                    
                    for b_track_id, frame_rois in bcell_track_rois.items():
                        for f, b_roi in frame_rois:
                            if f != frame: #check the b cell spots in the same frame
                                continue
                            #calculating b cell geometry here
                            bx = b_roi.getBounds().getCenterX()
                            by = b_roi.getBounds().getCenterY()
                            br = max(b_roi.getBounds().width, b_roi.getBounds().height) / 2.0
                            if x_min <= bx <= x_max and y_min <= by <= y_max:
                                if circle_roi_overlap(bx, by, br, mac_roi): #calls the perimeter point check
                                    #print("B cell track {} at frame {} overlaps mac track {}".format(b_track_id, frame, track_id))
                                    overlap_frames.append(frame) #if its overlapping, add this frame to the overlap_frames
                                    found = True
                                    break
                        if found:
                            break 
    
                overlap_frames = sorted(set(overlap_frames)) #sort everything in sequence
                threshold = phagocytic_frame_threshold #consecutive frame threshold
                is_phagocytic = False
                runs = []
                current_run = []
    
                for f in overlap_frames: #group everything in consecutive runs
                    if not current_run or f == current_run[-1] + 1:
                        current_run.append(f)
                    else:
                        runs.append(current_run)
                        current_run = [f]
                if current_run:
                    runs.append(current_run)
    
                for run in runs: #is the run long enough
                    if len(run) >= threshold:
                        is_phagocytic = True
                        break
    
                engulfed_bcell_rows = []
                for frame, mx, my, row in mac_entries:
                    mac_spot_id = int(row['ID'])
                    mac_roi = roi_dict.get(mac_spot_id, None)
                    if mac_roi is None:
                        continue
                    found = False 
                    if frame not in bc_by_frame:
                        continue
                    x_min, x_max, y_min, y_max = get_mac_roi_bounds(mac_roi)
                    for b_track_id, frame_rois in bcell_track_rois.items(): #checking every b cell track roi at a particular frame, this will let us track the b cells so we can check for more than 2 phagocytic events
                        for f, b_roi in frame_rois:
                            if f != frame:
                                continue
                            bx = b_roi.getBounds().getCenterX()
                            by = b_roi.getBounds().getCenterY()
                            br = max(b_roi.getBounds().width, b_roi.getBounds().height) / 2.0
                            if x_min <= bx <= x_max and y_min <= by <= y_max:
                                if circle_roi_overlap(bx, by, br, mac_roi):
                                    overlaps_by_bcell[b_track_id].append(f)
                                    found = True
                                    break
                        if found:
                            break
                        #if circle_ellipse_overlap(bx, by, br, mx, my, a, b, theta):
                        #roi_name = "ID{:05d}".format(int(row['ID']))
                        roi_name = "ID{}".format(int(row['ID']))  # No zero-padding
                      
            
                    threshold = phagocytic_frame_threshold
                    def has_consecutive_run(frames, min_length):
                        frames = sorted(set(frames))
                        current_run = []
                        for f in frames:
                            if not current_run or f == current_run[-1] + 1:
                                current_run.append(f)
                            else:
                                if len(current_run) >= min_length:
                                    return True
                                current_run = [f]
                        return len(current_run) >= min_length
                
                    qualified_bcell_tracks = [tid for tid, frames in overlaps_by_bcell.items()
                                              if has_consecutive_run(frames, threshold)]
                
                    
                    # --- SLIP-OFF FILTER ---
                    # Helpers for runs
                    def _runs(frames_sorted):
                        runs = []
                        cur = []
                        for f in frames_sorted:
                            if not cur or f == cur[-1] + 1:
                                cur.append(f)
                            else:
                                runs.append(cur)
                                cur = [f]
                        if cur:
                            runs.append(cur)
                        return runs
                
                    def _first_run_at_least(frames_sorted, L):
                        for r in _runs(frames_sorted):
                            if len(r) >= L:
                                return (r[0], r[-1])
                        return None
                
                    stable_bcell_tracks = []
                    for tid in qualified_bcell_tracks:
                        inside_frames_sorted = sorted(set(overlaps_by_bcell[tid]))
                        first_ok = _first_run_at_least(inside_frames_sorted, phagocytic_frame_threshold)
                        if first_ok is None:
                            continue  # safety; shouldn't happen given qualification
                        start_of_first_run, end_of_first_run = first_ok
                
                        # Collect frames where this B cell is OUTSIDE after the first qualifying run
                        outside_after = []
                        for f, b_roi in bcell_track_rois.get(tid, []):
                            f_int = int(f)
                            if f_int <= int(end_of_first_run):
                                continue
                            mac_roi = mac_roi_by_frame.get(f_int, None)  # same mac, same frame
                            if mac_roi is None:
                                continue
                            # perimeter-based containment check (same method you already use)
                            bx = b_roi.getBounds().getCenterX()
                            by = b_roi.getBounds().getCenterY()
                            br = max(b_roi.getBounds().width, b_roi.getBounds().height) / 2.0
                            
                            inside_frac = circle_roi_inside_fraction(bx, by, br, mac_roi, step_deg=15)
                            outside_frac = 1.0 - inside_frac
                            if outside_frac >= slip_outside_fraction_threshold:
                                outside_after.append(f_int)
                
                        # Did we get >= slip_off_threshold consecutive outside frames?
                        outside_after = sorted(set(outside_after))
                        slipped = False
                        if outside_after:
                            # check for a consecutive run >= slip_off_threshold
                            cur = []
                            for f in outside_after:
                                if not cur or f == cur[-1] + 1:
                                    cur.append(f)
                                else:
                                    if len(cur) >= slip_off_threshold:
                                        slipped = True
                                        break
                                    cur = [f]
                            if len(cur) >= slip_off_threshold:
                                slipped = True
                        
                        phago_slip_audit.append({
                        	'mac_track_id': int(track_id),
                        	'bcell_track_id': int(tid),
                        	'first_run_start': int(start_of_first_run),
                        	'first_run_end': int(end_of_first_run),
                        	'slip_check_from': int(end_of_first_run) + 1,
                        	'slipped': 1 if slipped else 0,
                        	'slip_frames': ','.join(str(f) for f in outside_after) if outside_after else ''
                        })
                
                        if not slipped:
                            stable_bcell_tracks.append(tid)
                
                    # Only now, if at least one B-cell track is stable, mark this macrophage track as phagocytic
                    if len(stable_bcell_tracks) > 0:
                        phagocytic_tracks.add(track_id)
                        for _, _, _, row in mac_entries:
                            phagocytic_rows.append(row)
                
                        # Use only the stable B-cell tracks for downstream export
                        engulfed_bcell_rows = []
                        for row in bc_spots:
                            try:
                                tid = int(row['TRACK_ID'])
                                if tid in stable_bcell_tracks:
                                    engulfed_bcell_rows.append(row)
                            except:
                                continue
                
                        if len(stable_bcell_tracks) >= 2:
                            multi_phagocytic_tracks.add(track_id)
                            for _, _, _, row in mac_entries:
                                multi_phagocytic_rows.append(row)
                                
                        engulfed_bcell_rows_all.extend(engulfed_bcell_rows)
                        engulfed_bcell_track_ids_all.update(stable_bcell_tracks)
    
            base_name = os.path.basename(spots_path).replace("_vx{}".format(macrophage_channel), "")
            summary_dir = os.path.dirname(spots_path)
    
            if phagocytic_rows:
                out_csv = os.path.join(summary_dir, base_name + "_phagocytic_spots.csv")
                with open(out_csv, 'w') as f:
                    writer = csv.DictWriter(f, fieldnames=phagocytic_rows[0].keys(), lineterminator='\n')
                    writer.writeheader()
                    writer.writerows(phagocytic_rows)
                print("Phagocytic spots saved to:", out_csv)
    
            excluded = set(phagocytic_tracks) | set(multi_phagocytic_tracks)
            if sub_bcell_track_rois and isinstance(sub_bcell_track_rois, dict):
                for track_id, mac_entries in mac_tracks.items():
                    if int(track_id) in excluded:
                        continue
                        
                    overlap_frames = []
                    
                    for frame, mx, my, row in mac_entries:
                        mac_spot_id = int(row['ID'])
                        mac_roi = roi_dict.get(mac_spot_id, None)
                        if mac_roi is None:
                            continue
                        
                        found = False
                        x_min, x_max, y_min, y_max = get_mac_roi_bounds(mac_roi)
                        for b_track_id, frame_rois in sub_bcell_track_rois.items():
                            for f, bx, by, br, b_roi in frame_rois:
                                if f != frame:
                                    continue
                                if x_min <= bx <= x_max and y_min <= by <= y_max:
                                    if circle_roi_overlap(bx, by, br, mac_roi):
                                        overlap_frames.append(frame)
                                        found = True
                                        break
                            if found:
                                break
                    
                    overlap_frames = sorted(set(overlap_frames))
                    threshold = phagocytic_frame_threshold
                    is_trogocytic = False
                    runs =[]
                    current_run = []
                    
                    for f in overlap_frames:
                        if not current_run or f == current_run[-1] + 1:
                            current_run.append(f)
                        else:
                            runs.append(current_run)
                            current_run = [f]
                    if current_run:
                        runs.append(current_run)
                        
                    for run in runs:
                        if len(run) >= threshold:
                            is_trogocytic = True
                            break
                            
                    if is_trogocytic:
                        trogocytic_tracks.append(track_id)
                        for _, _, _, row in mac_entries:
                            trogocytic_rows.append(row)
                print("Detected {} ROI-based trogocytic macrophage tracks.".format(len(trogocytic_tracks)))
                
            excluded_for_stalled = set([int(tid) for tid in phagocytic_tracks] + [int(tid) for tid in multi_phagocytic_tracks] + [int(tid) for tid in trogocytic_tracks])
            for track_id, mac_entries in mac_tracks.items():
                if int(track_id) in excluded_for_stalled:
                    continue
                
                frame_intensity = {}
                for frame, mx, my, row in mac_entries:
                    try:
                        target_ch_mean = float(row.get(channel_key, 0.0))
                        if target_ch_mean >= mac_threshold:
                            frame_intensity[int(frame)] = target_ch_mean
                    except:
                        continue

                above_thresh_frames = [f for f in frame_intensity if frame_intensity[f] >= mac_threshold]
                
                runs = []
                current_run = []
                for f in sorted(above_thresh_frames):
                    if not current_run or f == current_run[-1] +1:
                        current_run.append(f)
                    else:
                        runs.append(current_run)
                        current_run = [f]
                if current_run:
                    runs.append(current_run)
                
                stalling = False
                for run in runs:
                    if len(run) >= (phagocytic_frame_threshold):
                        stalling = True
                        break
                
                if stalling:
                    stalling_cup_tracks.append(track_id)
                    for _, _, _, row in mac_entries:
                        stalling_cup_rows.append(row)
                        
            if multi_phagocytic_rows:
                out_multi_csv = os.path.join(summary_dir, base_name + "_multi_phagocytic_spots.csv")
                with open(out_multi_csv, 'w') as f:
                    writer = csv.DictWriter(f, fieldnames=multi_phagocytic_rows[0].keys(), lineterminator='\n')
                    writer.writeheader()
                    writer.writerows(multi_phagocytic_rows)
                print("Multi-phagocytic spots saved to:", out_multi_csv)
    
            if stalling_cup_rows:
                out_stalling_csv = os.path.join (summary_dir, base_name + "_stalling_spots.csv")
                with open(out_stalling_csv, 'w') as f:
                    writer = csv.DictWriter(f, fieldnames=stalling_cup_rows[0].keys(), lineterminator='\n')
                    writer.writeheader()
                    writer.writerows(stalling_cup_rows)
                print("Stalling spots saved to:", out_stalling_csv)
            
            if trogocytic_rows:
                out_trogo_csv = os.path.join (summary_dir, base_name + "_trogocytic_spots.csv")
                with open(out_trogo_csv, 'w') as f:
                    writer = csv.DictWriter(f, fieldnames=trogocytic_rows[0].keys(), lineterminator='\n')
                    writer.writeheader()
                    writer.writerows(trogocytic_rows)
                print("Trogocytic spots saved to:", out_trogo_csv)
                
            if engulfed_bcell_rows_all:
            	engulfed_csv_path = os.path.join(summary_dir, base_name + "_engulfed_target_spots.csv")
            	with open(engulfed_csv_path, 'w') as f:
            		writer = csv.DictWriter(f, fieldnames=engulfed_bcell_rows_all[0].keys(), lineterminator='\n')
            		writer.writeheader()
            		writer.writerows(engulfed_bcell_rows_all)
            	print("Engulfed target spots saved to:", engulfed_csv_path) 
            
            if phago_slip_audit:
            	audit_csv_path = os.path.join(summary_dir, base_name + "_phago_slip_audit.csv")
            	audit_fields = ['mac_track_id', 'bcell_track_id', 'first_run_start', 'first_run_end', 'slip_check_from', 'slipped', 'slip_frames']
            	with open(audit_csv_path, 'w') as f:
            		writer = csv.DictWriter(f, fieldnames=audit_fields, lineterminator='\n')
            		writer.writeheader()
            		for row in phago_slip_audit:
            			writer.writerow(row)
            	print("Wrote slip audit:", audit_csv_path)
            
            t9 = time.time()*1000
            dt9= (t9-t0)/60000
            summary_csv = os.path.join(summary_dir, "phago_efficiency_summary.csv")
            write_header = not os.path.exists(summary_csv)
            with open(summary_csv, 'a') as f:
                writer = csv.writer(f, lineterminator="\n")
                if write_header:
                    writer.writerow([
                        "Image",
                        "Total_Macrophage_Tracks",
                        "Phagocytic_Tracks",
                        "Efficiency_Percent",
                        "Multi-Phagocytic_Tracks",
                        "Multi-Phagocytic_Percent",
                        "Trogocytic_Tracks",
                        "Trogocytic_Percent",
                        "Stalling_Tracks",
                        "Stalling_Percent",
                        "Total_Time_Min"
                    ])
                efficiency = 100.0 * len(phagocytic_tracks) / total_tracks if total_tracks > 0 else 0.0
                multi_percent = 100.0 * len(multi_phagocytic_tracks) / len(phagocytic_tracks) if len(phagocytic_tracks) > 0 else 0.0 if len(phagocytic_tracks) > 0 else 0.0
                stalling_percent = 100.0 * len(stalling_cup_tracks) / total_tracks if total_tracks > 0 else 0.0
                trogocytic_percent = 100.0 * len(trogocytic_tracks) / total_tracks if total_tracks > 0 else 0.0
                writer.writerow([
                    base_name,
                    total_tracks,
                    len(phagocytic_tracks),
                    "{:.2f}".format(efficiency),
                    len(multi_phagocytic_tracks),
                    "{:.2f}".format(multi_percent),
                    len(trogocytic_tracks),
                    "{:.2f}".format(trogocytic_percent),
                    len(stalling_cup_tracks),
                    "{:.2f}".format(stalling_percent),
                    "{:.2f}".format(dt9)
                ])
                
            summary_txt = os.path.join(summary_dir, base_name + "_phago_summary.txt")
            with open(summary_txt, 'w') as f:
                f.write("Phagocytic Efficiency Summary (Ellipse Overlap)\n")
                f.write("Image: {}\n".format(base_name))
                f.write("Total macrophage tracks: {}\n".format(total_tracks))
                f.write("Phagocytic tracks: {}\n".format(len(phagocytic_tracks)))
                f.write("Trogocytic tracks:{}\n".format(len(trogocytic_tracks)))
                f.write("Stalled tracks:{}\n".format(len(stalling_cup_tracks)))
                if total_tracks > 0:
                    efficiency = 100.0 * len(phagocytic_tracks) / total_tracks
                    f.write("Efficiency: {:.2f}%\n".format(efficiency))
                    f.write("Phagocytic macrophages with ≥2 B cells: {}\n".format(len(multi_phagocytic_tracks)))
                    if len(phagocytic_tracks) > 0:
                        f.write("Percent of phagocytic macrophages with ≥2 B cells: {:.2f}%\n".format(
                            100.0 * len(multi_phagocytic_tracks) / len(phagocytic_tracks)))
                    if len(trogocytic_tracks) > 0:
                    	f.write("Percent of trogocytic macrophages: {:.2f}%\n".format(
                            100.0 * len(trogocytic_tracks) / total_tracks))
                    if len(stalling_cup_tracks) > 0:
                    	f.write("Percent of stalling events: {:.2f}%\n".format(
                            100.0 * len(stalling_cup_tracks) / total_tracks))
    
            # (removed string conversion) phagocytic_tracks = set(phagocytic_tracks)
            print("Final list of phagocytic tracks to export:", phagocytic_tracks)
            # === Export phagocytic tracks to new TrackMate XML and ROIs ===
            from fiji.plugin.trackmate import Model as TrackMateModel
            from java.io import File
    
            if phagocytic_tracks:
                print("Creating filtered TrackMate model and ROI export...")
    
                original_model = trackmate.getModel()
                original_track_model = original_model.getTrackModel()
                original_spot_model = original_model.getSpots()
    
                filtered_model = TrackMateModel()
                filtered_model.beginUpdate()
    
                for track_id in original_track_model.trackIDs(True):
                    if track_id in phagocytic_tracks:
                        for spot in original_track_model.trackSpots(track_id):
                            filtered_model.addSpotTo(spot, int(spot.getFeature('FRAME')))
                        for edge in original_track_model.trackEdges(track_id):
                            source = original_track_model.getEdgeSource(edge)
                            target = original_track_model.getEdgeTarget(edge)
                            cost = original_model.getFeatureModel().getEdgeFeature(edge, EdgeTargetAnalyzer.EDGE_COST)
                            if cost is None:
                                cost = 1.0  # fallback weight
                            filtered_model.addEdge(source, target, cost)

    
                filtered_model.endUpdate()
    
                filtered_xml_path = os.path.join(summary_dir, base_name + "_phagocytic_TrackMate.xml")
                filtered_writer = TmXmlWriter(File(filtered_xml_path), logger)
                filtered_writer.appendModel(filtered_model)
                filtered_writer.appendSettings(settings)
                filtered_writer.writeToFile()
                print("Saved filtered TrackMate XML to:", filtered_xml_path)
    
                from java.awt import Color
                
                roiManager = RoiManager.getInstance()
                if roiManager is None:
                    roiManager = RoiManager()
    
                roiManager.reset()
                for track_id in original_track_model.trackIDs(True):
                    if track_id in phagocytic_tracks:
                        print("Exporting track ID:", track_id)
                        for spot in original_track_model.trackSpots(track_id):
                            roiExporter.export(spot)
                            if roi:
                            	roi.setName("ID{:05d}".format(spot_id))
                            	roi.setStrokeColor(Color.MAGENTA)
                            	roiManager.addRoi(roi)

                phago_roi_path = os.path.join(summary_dir, base_name + "_phagocytic_ROIs.zip")
                roiManager.runCommand("Save", phago_roi_path)
                print("Saved filtered ROIs to:", phago_roi_path)
                
            from java.awt import Color
            
            if stalling_cup_tracks:
                roiManager = RoiManager.getInstance()
                if roiManager is None:
                    roiManager = RoiManager()
                roiManager.reset()
                
                for track_id in stalling_cup_tracks:
                    for frame, mx, my, row in mac_tracks[track_id]:
                        spot_id = int(row['ID'])
                        roi = roi_dict.get(spot_id, None)
                        if roi:
                            roi.setName("ID{:05d}".format(spot_id))
                            roi.setStrokeColor(Color.BLUE)
                            roiManager.addRoi(roi)
                    
                    stalling_roi_path = os.path.join(summary_dir, base_name + "_stalling_cup_ROIs.zip")
                    roiManager.runCommand("Save", stalling_roi_path)
                    print("Saved stalling cup ROIs to:", stalling_roi_path)
                    
            from java.awt import Color
            
            if engulfed_bcell_track_ids_all:
                roiManager = RoiManager.getInstance()
                if roiManager is None:
                    roiManager = RoiManager()
                roiManager.reset()
                
                exported_count = 0
                
                for b_tid in sorted(list(engulfed_bcell_track_ids_all)):
                    for f, b_roi in bcell_track_rois.get(b_tid, []):
                    	try:
                    		b_roi.setStrokeColor(Color.MAGENTA)
                    	except:
                    		pass
                    	roiManager.addRoi(b_roi)
                    	exported_count += 1
                if exported_count >0:
                    engulfed_roi_path = os.path.join(summary_dir, base_name + "_engulfed_cell_ROIs.zip")
                    roiManager.runCommand("Save", engulfed_roi_path)
                    print("Saved engulfed cell ROIs to:", engulfed_roi_path)
            
    else:
        print("No CSV file in path:", spots_path)
    
    if closeSpots == True:
        IJ.selectWindow("Spots"); 
        IJ.run("Close");
    
    t9 = time.time()*1000
    dt9 = (t9-t0)/60000
    print("Image loop took "+str(round(dt9,1))+" min at "+str(round(spf,3))+" s/frame")
    
    imp.close()
    print("TrackMate finished for file {}, channel {}".format(file_path, channel))
    
fNumber = 0
for file_path in file_paths:
    fNumber = fNumber+1
    dt_string = dt.now().strftime("%d/%m/%Y %H:%M:%S")
    for channel_info in [ {'channel': b_cell_channel, 'model': 'subthreshold', 'start': bcell_start, 'end': bcell_end}, {'channel': b_cell_channel, 'model': 'default', 'start': bcell_start, 'end': bcell_end}, {'channel': macrophage_channel, 'model': macrophage_model, 'start': macrophage_start, 'end': macrophage_end}]:
        print("I'm in the loop with channel info: " + str(channel_info) )
        run_trackmate_on_channel(file_path, channel_info['channel'], channel_info['model'])