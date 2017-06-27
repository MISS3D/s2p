import gdal
import sys 
import math
import csv
from osgeo import osr

# def read3DPointCloudOfMOLAFromCSV(filename):
# 	"""read the 3D point cloud of the PEDR File MOLA"""
# 	listPoint=[]
# 	with open(filename,'r') as csvfile:
# 		spamreader=csv.reader(csvfile,delimiter='\n')
# 		for row in spamreader:
# 			listPoint.append(row)
# 		listPoint=[listPoint[i][0].split(",") for i in range(0,len(listPoint))]
# 	return listPoint

# def convert3DPointOFMolaForProcessing(listPoint):
# 	del listPoint[0]
# 	listPoint=[[float(listPoint[i][0]),float(listPoint[i][1]),float(listPoint[i][4])] for i in range(0,len(listPoint))]
# 	listPoint=[[listPoint[i][0]-360,listPoint[i][1],listPoint[i][2]-3396100] for i in range(0,len(listPoint))]
# 	#title=['long_East', 'lat_North', "heightFromMarsDatum"]
# 	#listPoint=[title,listPoint]
# 	return listPoint
# def grabMolaPoints(cloudPoint,long_min,long_max,lat_min,lat_max):
# 	long_max=float(long_max)
# 	long_min=float(long_min)
# 	lat_max=float(lat_max)
# 	lat_min=float(lat_min)
# 	localCloudPoint=[]
# 	for el in cloudPoint:
# 		if ((long_min <= el[0] <= long_max) and (lat_min <= el[1] <= lat_max)):
# 			localCloudPoint.append(el)
# 	return localCloudPoint

def heights_of_mola_from_DEM(DEM,long_min_DEM,lat_min_DEM,long_res_DEM,lat_res_DEM,long_min,long_max,lat_min,lat_max):
	"""compute the minimum and maximal height within a part of the DEM of MOLA, limited in the area: long_min,long_max,lat_min,lat_max 

	   args: DEM is a numpy.ndarray, all the other elements are of type float
	"""
	# long_max=float(long_max)
	# long_min=float(long_min)
	# lat_max=float(lat_max)
	# lat_min=float(lat_min)

	# in the image the lowest latitude is in its bottom, and the lowest longitude in the left
	# also here the index of the rows and collumns are
	max_row=len(DEM)-1-int(math.floor((lat_min-lat_min_DEM)/lat_res_DEM))
	min_row=len(DEM)-1-int(math.floor((lat_max-lat_min_DEM)/lat_res_DEM))

	min_col=int(math.floor((long_min-long_min_DEM)/long_res_DEM))
	max_col=int(math.floor((long_max-long_min_DEM)/long_res_DEM))

	# if the zone if outside the mola DEM, we shuld get the intersection
	try:
		assert (min_row >= 0 or min_col >= 0 or max_row < len(DEM) or max_col < len(DEM[0]))
	except AssertionError:
		print("The RPC functions don't cover the patch of the mola DEM ")
		# min_row=max(0,min_row)
		# max_row=min(len(DEM)-1,max_row)
		# min_col=max(0,min_col)
		# max_col=min(len(DEM[0])-1,max_col)
	# rows=len(array)
	# cols=len(array[0])
	heights=[]
	for col in range(min_col,max_col+1):
		for row in range(min_row,max_row+1):
			heights.append(DEM[row][col])
	maxHeight=max(heights)
	minHeight=min(heights)
	return [minHeight,maxHeight]
				# longAndlat=pixel2coord(col,row,conversionParams)/
				# cloudPoints.append([longAndlat[0],longAndlat[1],array[row][col]])
	# dx and dy, the resolution od the dem must be saved
	# dx=conversionParams[1];dy=conversionParams[5]

	# localCloudPoint=[]
	# for el in cloudPoint:
	# 	if ((ul[0] <= el[0] <= ul[1]) and (lat_min <= el[1] <= lat_max)):
	# 		localCloudPoint.append(el)
	# return localCloudPoint

# def closest_point(longitude_of_ray,latitude_of_ray,cloudPoint):
# 	"""search in a cloud point computed from a dem, the nearest point to a ray of longitude "longitude_of_ray" and latitude "latitude_of_ray
# 	   the cloudPoint must cover this range
# 	" 
# 	   args: cloudPoint is a list, and longitude_of_ray, latitude_of_ray are from type double
# 	"""
# 	closest_point=cloudPoint[0][0],cloudPoint[0][1]
# 	for el in cloudPoint:

def pixel2coord(x, y,conversionParams):
    """Returns global coordinates from pixel x, y coords, in the case of the DEM of s2p, it's the latitude and longitude coordinate"""
    xoff, a, b, yoff, d, e =conversionParams
    xp = a*(x+0.5)+b*(y+0.5)+xoff
    yp = d*(x+0.5)+e*(y+0.5)+yoff
    return[xp, yp]

def read_mola_DEM(filename):
	"""Tranform the DEM to a cloud points by considering the point at the center of the pixel and with the same elevation value as height and the value of the pixel od longitude and latitude
	    args: the full path to the dem
	    return: the dem expressed as an "numpy.ndarray", the upper left point of the dem expressed in latitude and longitude coordinates, and the resolution of the DEM

	"""
	ds=gdal.Open(filename)
	cols=ds.RasterXSize
	rows=ds.RasterYSize
	array=ds.ReadAsArray() #array[i] is the ith row 
	cloudPoints=list()
	ulx, xres, xskew, uly, yskew, yres=ds.GetGeoTransform()
	# for row in range(0,rows):
	# 	for col in range(0,cols):
	# 		if (not(math.isnan(array[row][col]))):
	# 			longAndlat=pixel2coord(col,row,conversionParams)
	# 			cloudPoints.append([longAndlat[0],longAndlat[1],array[row][col]])
	#dx and dy, the resolution od the dem must be saved
	#dx=conversionParams[1];dy=conversionParams[5]

	# the coordinates must be converted  to lat long
	src = osr.SpatialReference()
	tgt = osr.SpatialReference()
	src.ImportFromProj4("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=3396100 +b=3396100 +units=m +no_defs") # can be verifies with gdalsrsinfo, it's outpus should be put directly here 
	tgt.ImportFromUrl("http://spatialreference.org/ref/iau2000/49900/")
	transform = osr.CoordinateTransformation(src, tgt)
	# for i in range(0,len(cloudPoints)):
	# 	el=cloudPoints[i]
	# 	coords = transform.TransformPoint(el[0],el[1])[0:2]
	# 	alt=el[2]
	# 	el=[coords[0],coords[1],alt]
	# 	cloudPoints[i]=el
	ulx,uly = transform.TransformPoint(ulx,uly)[0:2]
	resx,resy = transform.TransformPoint(xres,yres)[0:2]
	[lrx,lry]=[ulx+cols*resx,uly+rows*resy]
	long_min=min(lrx,ulx)
	lat_min=min(lry,uly)
	long_res=math.fabs(resx)
	lat_res=math.fabs(resy)
	return array,[long_min,lat_min,long_res,lat_res]
# def computeNewElevation(x,y,w,h,min_elevation,max_elevation,cloudPoint,zoneError=0.05,elevationMolaError=0):
# 	"""x,y,w,h the zone of points near which we aim to find molaPoints,
# 	   min_elevation max_elevations are the maximal and minimal elevation for which the RPC function are computed,
# 	   cloudPoint is the cloud of MOLA Point 
# 	    """

# 	    a = np.array([x, x,   x,   x, x+w, x+w, x+w, x+w])
# 	    b = np.array([y, y, y+h, y+h,   y,   y, y+h, y+h])
# 	    c = np.array([m, M,   m,   M,   m,   M,   m,   M])

# 	    t1,t2,t3=rpc1.direct_estimate(a,b,c)
# 	    long_min=min(t1);long_max=max(t1);lat_min=min(t2);lat_max=max(t2)
# 	    long_max=float(long_max);long_min=float(long_min);lat_max=float(lat_max);lat_min=float(lat_min)
# 	    localCloudPoint=grabMolaPoints(cloudPoint,long_min-zoneError,long_max+zoneError,lat_min-zoneError,lat_max+zoneError)
# 	    return [min(localCloudPoint)[2]-elevationMolaError,max(localCloudPoint)[2]+elevationMolaError]