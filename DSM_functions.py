import numpy as np
import math
import pydicom
import pyvista as pv
from shapely.geometry import Polygon
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

####    -----------------------     ####
#   CODE TO GET DATA FROM DICOMS        #
####    -----------------------     ####
def get_ROI_Contour_Points(ROI_Name, RS_filepath):
    """
    Function: extracts the contour points of a specific ROI from an RS file

    Input:
    ROI_Name: name of the region of interest found in the specified RS file 
    RS_filepath: path from current directory to region where RS.dcm file is located

    Output:
    A 3D np array of the contour points of the ROI
    """
    rs = pydicom.read_file(RS_filepath)
    ROI_index,ROIlist = None,[]
    for index, item in enumerate(rs.StructureSetROISequence):
        ROIlist.append(item.ROIName)
        if item.ROIName == ROI_Name:
            ROI_index = index
            break
    if ROI_index == None:
        raise Exception('An ROI with the name you specified does not exist in the RT Structure file. The structures in this file are:\n',ROIlist)
    #Get all contour points from RS file (organized as: [[x0-N, y0-N, z0-N][x0-N, y0-N, z0-N]] )
    contour = []
    for item in rs.ROIContourSequence[ROI_index].ContourSequence:
        contour.append(item.ContourData)
    return np.array(contour)

def get_contourpoint_cloud(ROI_Name, RS_filepath):
    """
    Function: extracts the contour points of a specific ROI from an RS file and creates a 2D list of 
    vertex points that can be used to create a 3D surface mesh. It also generates a list of the 
    centroids of each slice of the contour that can serve as a centerline for a non-planar DSM

    Input:
    ROI_Name: name of the region of interest found in the specified RS file 
    RS_filepath: path from current directory to region where RS.dcm file is located

    Output:
    xdata,ydata,zdata: 1D arrays of the x,y,z locations of each vertex in the contour mesh
    CoMList:           The centers of mass (centroids) of each axial slice of the contour (can be used to create a centerline)
    """
    #get the raw mesh data
    mesh = get_ROI_Contour_Points(ROI_Name, RS_filepath)
    #prep output data arrays
    xdata = []
    ydata = []
    zdata = []
    Ztracking = [[],[]]     #used to check for duplicate Z (axial) indices
    CoMList = []
    for plane in mesh:
        #get points for point cloud
        xvals = (plane[0::3])
        yvals = (plane[1::3])
        zvals = (plane[2::3])
        #get details for centerline
        zval = plane[2]
        points = np.array([plane[0::3],plane[1::3]]).T
        poly = Polygon(points)
        #if a contour on same slice exists, use the bigger of the two (happens for eclipse contours if contourer accidentally missed or included some extra voxel)
        if zval in Ztracking[0] and poly.area > Ztracking[1][Ztracking[0].index(zval)]:
            print('taking the larger at',zval)
            #replace the old centerline values with new ones
            CoMList[Ztracking[0].index(zval)] = [poly.centroid.x,poly.centroid.y,float(zval)]
            Ztracking[1][Ztracking[0].index(zval)] = poly.area
            #replace the old pointcloud ones with new
            oldindex = zdata.index(zvals[0])
            xdata, ydata, zdata = xdata[:oldindex], ydata[:oldindex], zdata[:oldindex]
            xdata.extend(xvals)
            ydata.extend(yvals)
            zdata.extend(zvals)
        elif zval in Ztracking[0]:
            print('leaving the larger at',zval)
            continue    #dont add, just go to next iteration of the loop
        else:
            #add the data to the arrays
            Ztracking[0].append(zval)
            Ztracking[1].append(poly.area)
            CoMList.append([poly.centroid.x,poly.centroid.y,float(zval)])
            xdata.extend(xvals)
            ydata.extend(yvals)
            zdata.extend(zvals)
    #return three 1D arrays of all the coordinate info (which can be combined into a 2D numpy array later if desired), as well as the COM line
    return xdata,ydata,zdata, CoMList

def get_RT_Dose_info(DCM_folderpath,RD_file):
    """
    Function: extracts the dose grid and axes of that dose grid from a DICOM RD file
    Input:
    DCM_folderpath: path to a directory that contains the RD_dose file
    RD_file:        name of the RD_dose file to read

    Output:
    dosegrid:       A 3D np.array of RT dose values 
    xpos,ypos,zpos: 3 1D np.arrays that give the xyz coordinates of each grid voxel (eg. indicates nth voxel along y direction is at ypos[n] mm along y axis )
    """
    ds = pydicom.read_file(DCM_folderpath + '/' + RD_file)
    #STEP1: get the dose grid
        #ds.pixel_array is the dose grid, and ds.DoseGridScaling is the scaling factor that converts pixels to Gy
    dosegrid= ds.pixel_array*ds.DoseGridScaling
    #STEP2: make positional arrays that will be used to lookup dose values when creating DSMs
        #each array monotonically increases from the topcorner of the grid (imgpospat) with steps equal to the gird resolution (pixelspacing)
    xpos = np.arange(ds.Columns)*ds.PixelSpacing[0]+ds.ImagePositionPatient[0]
    ypos = np.arange(ds.Rows)*ds.PixelSpacing[1]+ds.ImagePositionPatient[1]
    zpos = np.asarray(ds.GridFrameOffsetVector) + ds.ImagePositionPatient[2]

    return dosegrid, xpos,ypos,zpos

####    -----------------------     ####
#   SMALL FUNCTIONS USED IN OTHER FUNCTIONS     #
####    -----------------------     ####
def find_nearest(array,value):
    """
    Function: returns the index of the element in the array closest in value to the value specified (provided array is sorted)
    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def plot_single_DSM(DSMarray,title,savename):
    """
    Function: Uses matplotlib to plot a single DSM
    """
    fig, ax = plt.subplots()
    f = ax.pcolormesh(DSMarray,cmap='inferno')
    ax.set_title(title)
    fig.colorbar(f, ax=ax)
    plt.savefig(savename+".png")
    return

def find_tangent_line(Pfix,P1,P2):
    """
    Function: estimates the tangent of a line in 3D space at the requested point, given its two nearest neighbours
    Input:
    Pfix:       the point at which we wish to find the tangent (x,y,z)
    P1,P2:      The points immediately before and after Pfix on the line (x,y,z coordinates)

    Output:
    direction_vector:   the tangent vector to the line at the point Pfix    
    """
    #Set an arbitrary start point for the optimization function
    X1 = Pfix + 2*(Pfix - P1)
    print(P1,Pfix,P2)
    #define the point-line distance equations for P1 and P2
    D1 = lambda X : np.linalg.norm(np.cross((P1 - Pfix),(P1 - X)) ) / np.linalg.norm((X - Pfix))  
    D2 = lambda X : np.linalg.norm(np.cross((P2 - Pfix),(P2 + X - 2*Pfix)) ) / np.linalg.norm((X - Pfix))
    Dtot = lambda X: D1(X)**2 + D2(X)**2
    #try to optimize
    res = minimize(Dtot,X1,method ='powell')
    print(res.message,'Success:',res.success)
    X1 = res.x
    #get the direction vector (AKA 3D slope surrogate) we need to define the line equation
    direction_vector = Pfix - X1
    return direction_vector
####    -----------------------     ####
#   CODE TO SAMPLE POINTS ALONG VARIOUS OBJECTS     #
####    -----------------------     ####
def sample_centerline(spline,step,steptype='dist'):
    """
    Function: given a spline, define equidistant points at which rings of the contour mesh will be sampled at
    Input:
    spline:         A pyvista.spline object that gives the centerline of the contour
    step:           Either the step size (in mm) or number of steps along the contour to be sampled at (float)
    steptype:       string that specifies if step represents mm between sampling points or total number of sampling points
                        Values: 'dist' or 'nsteps'

    Output:
    sampinds:       1D array of indices that give the points in spline.points at which rings of the contour mesh shoudl be defined at
    """
    #STEP1: get total length of the centerline spline
    totlen = spline['arc_length'][-1]
        #STEP1b: if using nsteps method, figure out how big step size should be
    if steptype=='nsteps':
        step = totlen/step
    #STEP2: set the sampling points along the spline and get the sampling points
    sampdists = np.arange(0,totlen,step)
    sampinds = [] 
    for samp in sampdists:
        sampinds.append(find_nearest(spline['arc_length'],samp))
    #STEP3: return the indicies of the sampling points (we give this instead of just the points themseves so can easily use to get points and path lengths at each sampling point)
    return sampinds
  
def sample_contour_surface(SurfaceObject,Centerline,SamplingPoints,NonPlanar=False):
    """
    Function: extracts the dose grid and axes of that dose grid from a DICOM RD file
    Input:
    SurfaceObject:      A pyvista volume or shell object we will be sampling in slices
    Centerline:         A pyvista spline object defining the centerline of the contour to be sampled
    SamplingPoints:     1D list of indicies of points along the centerline that define where sampling slices are anchored
    NonPlanar:          True or False flag that indicates if a planar or non-planar sampling style should be used (default=planar)

    Output:
    Slicedict:          Nested dictionary object containing the centroid, plane normal vector, and vertices (for visualization) of each sampling slice
                        Keys are integers from SamplingPoints (eg. first point = Slicedict[0])
    """
    #STEP1: Declare dictionary to store slice data
    Slicedict = {}
    #STEP2: Step through the sampling points list and begin taking slices
    for p in SamplingPoints:
        #determine the spline tangent at the sampling point to get the normal of the sampling slice
        if NonPlanar == False or p ==0 or p == len(Centerline):
            tan = [0,0,1]   #just take a slice on the axial plane (we do this if using planar style or for top and bottom slices)
        else:
            tan = find_tangent_line(Centerline[p],Centerline[p-1],Centerline[p+1])  #custom stepwise tangent function
            tan = tan/np.linalg.norm(tan)   #scale so that its a unit vector with length 1 (important for later)
        #use the pyvista slice function to slice the contour surface
        pslice = SurfaceObject.slice(normal=list(tan),origin=list(Centerline[p]))
        #add the slice and point data to the dictionary
        Slicedict[p] = {'point':Centerline[p], 'slice':pslice, 'planenormal':tan}
    #finish and return the slice dictionary
    return Slicedict
  
def interpolate_missed_points(samppoints,failvects,centroid):
    """
    Function: If ray-tracing failed to detect a collision at some points, uses a spline loop of contour slice to find 
                the closest point on the mesh where ray tracing should have found a collision
                NOTE: because the spline loop is created with the successful collision points, if ray-tracing failed too 
                many points (eg. more than half) the missed points may be interpolated incorrectly
    Input:
    samppoints:     2D np.array of sampling point coordinates. any missing points are given as np.nan
    failvects:      2D np.array of the directional vectors used in ray tracing that failed to detect a collision with the contour mesh
    centroid:       coordinate from which ray-tracing rays originate from

    Output:
    samppoints:     updated version of the original samppoints array that has filled in all the missing points with best interpolated values
    """
    #STEP1: Create a densely sampled spline with the successfully sampled points    
    successes = samppoints[~np.isnan(samppoints).any(axis=1)]
    failinds = np.argwhere(np.isnan(samppoints[:,0])).ravel()
    points = np.zeros((len(samppoints)-len(failvects)+1,3)) #make placeholder array 1 point longer than successful points
    points[:-1,:],points[-1,:]= successes, successes[0,:]   #fill with ordered list and repeat first point at end to make loop
    intspline = pv.Spline(points,720)   #using 720 points to try to make sure one point per degree at least
    ringpoints = intspline.points
    #STEP2: for each vector that didnt successfully ray trace, find the closest point on the spline
    print('Calculating points missed by ray tracing...')
    Dtoline = lambda X : np.linalg.norm(np.cross((X - centroid),(X - (centroid+vect))) ) / np.linalg.norm((vect)) #distance between point and line in 3D coordinate space
    ii = 0
    for vect in failvects:
        dists = np.apply_along_axis(Dtoline, 1, ringpoints) #calculate distances for all points on the spline
        #get the 3 closest dists to the vecorline, then check to make sure you grab one along correct direction of ray (not 180 of intended direction)
        idx = np.argpartition(dists, 3)
        closest = ringpoints[idx[:3],:]
        vdists = np.linalg.norm(closest - np.array([centroid+vect,centroid+vect,centroid+vect]), axis=1) #distances from point one unit step along ray in desired direction
        if (vdists[0]-vdists[1]) > 1 or (vdists[0]-vdists[2]) > 1:
            #if the first point more than a mm further from the vector end point than the other two, take the best of them
            if (vdists[1]-vdists[2]) > 1:
                closest = closest[2]     #if point 2 much further than point 3 (aka both points 1&2 were at the 180 position, take point3)
            else:
                closest = closest[1]    #otherwise take point 2, as the vectorline passes closer to it
        else:
            closest = closest[0]
        #overwrite the failed points (aka the [nan,nan,nan]s) in the original points array
        samppoints[failinds[ii],:] = closest
        ii += 1
    return samppoints

def sample_uniangular(plane,normal,centroid,nsamples,fixfailures=True):
    """
    Function: Finds N equiangularly spaced points on the contour for a given sampling plane.
                This done by using ray tracing from the centroid of the sampling slice
    Input:
    plane:
    normal:         a unit vector (x,y,z) perpendicular to the sampling plane
    centroid:       3D coordinate (x,y,z) at which rays originate from
    nsamples:       how many points to sample in the plane (note: for unwrapping purposes, points 1 and n will both be the point posterior to the centroid)
    fixfailures:    True or False flag that indicates if interpolation should be used to calculate points for which ray-tracing failed

    Output:
    samppoints:     2D np.array of sampling point coordinates. any missing points are given as np.nan
    failvects:      list of directional vectors for whichray-tracing failed to detect a collision with the contour (for troubleshooting)
    """
    #STEP1: create a list of the angles to sample at, then define the sampling vectors around a unit circle in an xy plane
    sampangs = np.linspace(0,2*math.pi,nsamples)
    sampvects = np.array([np.sin(sampangs),np.cos(sampangs),np.zeros(len(sampangs))]).T 
    #take the dotproduct of the vectors and the plane normain order to then calculate projection onto the plane
        #we will also multiply these vectors by 100 to increase their magnitude and ensure they reach the contour surface
    dotprod = np.sum(sampvects*np.tile(normal,(len(sampvects),1)),axis=1)   #I dont think this is doing anything
    rayvects = (sampvects - np.tile(normal,(len(sampvects),1))*dotprod[:, np.newaxis])*150
    #STEP2: begin iterating through the list and performing ray tracking
    failcount = 0
    failvects = []
    samppoints = np.full_like(rayvects,np.nan)
    for i in range(len(rayvects)):
        #project the sampling vector onto the slice plane
        inter,ind = plane.ray_trace(centroid,centroid+rayvects[i],first_point=True) #give start and end of ray, and only return first intersection
        if len(inter) == 0:
            failcount += 1
            failvects.append(rayvects[i])
            continue
        samppoints[i,:] = inter     #add to the list fo samp points (failed intersections will appear as nan)
    #STEP3: if ray tracing failed (tends to happen in less densely sapled regions of the mesh), run alternate method to get points
    if failcount > 0 and fixfailures==True:
        print('failed ',failcount,' of ',nsamples,'points for slice at',centroid,'. \nRunning interpolation fix now...')
        samppoints = interpolate_missed_points(samppoints,failvects,centroid)
    return samppoints, failvects

def sample_dose(samples,dosegrid,x,y,z):
    """
    Function: Calculates the dose at each point in a list of points to be sampled and returns a dose-surface map
    Input:
    samples:        3D array of coordinate points that define points on contour mesh to sample to create a DSM
    dosegrid:       3D dose grid retrived from RT Dose dicom file
    x,y,z:          Three 1D arrays that give the x,y,z coordinates of te voxels in the dose grid (output from get_RT_Dose_data)

    Output:
    DSM:            2D array of dose values (in Gy) AKA a dose-surface map
    """
    #STEP1: define the interpolation function we will use to sample the dose grid
        #please note dicome dose grids in np. array format have indexes (Z,Y,X) and not (X,Y,Z) order
    dose_interp = RegularGridInterpolator((z,y,x),dosegrid,fill_value=None) #if point outside the grid, dont extrapolate
    #STEP2: declare the array well store the dose data in
    DSM = []
    print('SAMPLING DSM NOW...')
    #STEP3: loop through the sampling list and take the dose at each point
    for plane in samples:
        row = []
        for point in plane:
            if np.isnan(point[0]):  #if coordinate missing for some reason, leave a hole in the map
                row.append(np.nan)
                continue
            dose = dose_interp(point[::-1]) #use interpolation function to find dose at the point
            row.append(dose)
        DSM.append(row)
    DSM = np.array(DSM,dtype=np.float)
    print('Done!')
    return DSM

def correct_slice_overlaps(Slices,samppoints,centerline):
    """
    Function: Checks if any slices of the contour mesh overlap with one another, and if they do correct them using the method
                of Witztum et al 2016 (doi: doi.org/10.1118/1.4964790) with a minor modification of ensuring the centroids of
                the corrected slices lie on the centerline (done to ensure origin used for ray tracing isnt too off-center)
    Input:
    Slices:         Slice dictionary created in sample_contour_surface
    samppoints:     Indices of points along the centerline spline where slices were defined at (also the keys of the slice dictionary)
    centerline:     A pyvista spline object defining the centerline passing through the contour of interest

    Output:
    Slices:         An updated version of the slice dictionary with adjusted normals and centroids for colliding slices
    """
    #STEP1: find all overlaps that exist
    collisions = []
    for ii in range(len(samppoints[:-4])):
        #check if collision with 4 nearest neighbours
        col1,col2  = Slices[samppoints[ii]]['slice'].boolean_difference(Slices[samppoints[ii+1]]['slice']), Slices[samppoints[ii]]['slice'].boolean_difference(Slices[samppoints[ii+2]]['slice'])
        col3,col4  = Slices[samppoints[ii]]['slice'].boolean_difference(Slices[samppoints[ii+3]]['slice']), Slices[samppoints[ii]]['slice'].boolean_difference(Slices[samppoints[ii+4]]['slice'])       
        #if there were collisions list the slice numbers where it happened
        if sum([len(col1.points),len(col2.points),len(col3.points),len(col4.points)])>0:
            cc = [ii]
            for col in [col1,col2,col3,col4]:
                if len(col.points)>0:
                    cc.append(ii+1+[col1,col2,col3,col4].index(col) )
            collisions.append(cc)

    print('collisions:',collisions)
    #STEP2: per collision, find four points at 0,90.180.270 on the nearest slices outside the collision (i call them crusts)
    for col in collisions:
        bcrust,tcrust = max(0,col[0]-1), min(len(samppoints),col[-1]+1)  #if the inclusion includes the first or last slice we only correct the other slices
        bLocks = sample_uniangular(shell,Slices[samppoints[bcrust]]['planenormal'],Slices[samppoints[bcrust]]['point'],4,fixfailures=False)
        tLocks = sample_uniangular(shell,Slices[samppoints[tcrust]]['planenormal'],Slices[samppoints[tcrust]]['point'],4,fixfailures=False)
        #STEP3: define vectors that connect corresponding angular points between the two crusts
        vects = tLocks[0] - bLocks[0]
        #STEP4: per slice in the collision that needs to be corrected, define new sampling plane that lies 1/(nfix+1)th along the vectors between the lock points on the crusts (n=number to be resampled)
        nfix = tcrust-bcrust-1
        fixus = np.arange(bcrust+1,tcrust)  #define slices to fix (different from col list in the event one of the colliding slices has to be used as a crust)
        lineseg = centerline.points[samppoints[bcrust]:samppoints[tcrust] ]
        for n in range(nfix):
            nLocks = bLocks[0] + vects*((n+1)/(nfix+1))     #eg. get points 1/3 of the way along the paths from bottom to top crusts
            norm = np.cross((nLocks[1]-nLocks[0]),(nLocks[2]-nLocks[0]))    #normal of plane is crossproduct of two vectors in the plane
            #find the point on the centerline that intersects or is closest to the plane (use dot product)
            Dplane = lambda X : abs(np.dot(norm,(X-nLocks[0])) )
            Dmin = np.argmin(np.apply_along_axis(Dplane, 1, lineseg) )
            #Update the plane normal vector and point for the sampling slice
            print('fixing colliding slice:',fixus[n],'from',col)
            Slices[samppoints[fixus[n]]]['point'] = lineseg[Dmin]
            Slices[samppoints[fixus[n]]]['planenormal'] = norm/np.linalg.norm(norm)
    #STEP5: return the fixed Slice array
    return Slices
