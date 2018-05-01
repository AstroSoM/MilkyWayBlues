import numpy as np
from PIL import Image

'''
Purpose:
--------
(1) Crop the original Milky Way roadmap
(2) Create a grayscale version of the cropped image
(3) Create color and grayscale images with transparent pixels outside the Milky Way radius

Resources:
----------
Milky Way Roadmap TIFF image by Robert Hurt (IPAC) available here:
http://www.spitzer.caltech.edu/images/1923-ssc2008-10a-A-Roadmap-to-the-Milky-Way
'''
#====================================================================================================
# Output filenames
fcrop  = "/Users/salvesen/outreach/asom/milkyway/data/MilkyWayRoadmapCrop.png"
fgray  = "/Users/salvesen/outreach/asom/milkyway/data/MilkyWayRoadmapGray.png"
fcolor = "/Users/salvesen/outreach/asom/milkyway/data/MilkyWayRoadmapColor.png"

# Input TIFF filename
fimg = "/Users/salvesen/outreach/asom/milkyway/data/MilkyWayRoadmap.tiff"

# NOTE: Python convention is for (0,0) to be the top left pixel of the image

# Read in the JPG file
im = Image.open(fimg)

#====================================================================================================
# Crop the Milky Way image to a square area centered on the Galactic center (by eye in SAOImage ds9)
ycroptop = 130   # [pixels]
ycropbot = 2530  # [pixels]
imcrop   = im.crop((0, ycroptop, im.width, ycropbot))
imcrop.save(fcrop)

#====================================================================================================
# Color version of the Milky Way image
imColor = imcrop.convert('RGBA')

# Make a Grayscale version of the Milky Way image
imGray = imcrop.convert('RGBA')
Nx = imGray.width   # [pixels]
Ny = imGray.height  # [pixels]
for i in np.arange(Nx):
    for j in np.arange(Ny):
        R,G,B,A = imGray.getpixel((i,j))
        Y = int(0.299*R + 0.587*G + 0.114*B)
        imGray.putpixel((i,j),(Y,Y,Y,255))

#====================================================================================================
# Make everything transparent outside of a radius spanning the smaller dimension of the images

# Galactic center (x,y) location (determined by eye in original image with SAOImage ds9)
xGC = 0.5 * Nx
yGC = 0.5 * Ny

# Radial extent of the Milky Way image
R = 0.5 * Nx

# Set pixels outside circle of radius "R" to be transparent
for i in np.arange(Nx):
    for j in np.arange(Ny):
        this_R = np.sqrt((i - xGC)**2 + (j - yGC)**2)
        if (this_R > R):
            # Grayscale image
            Rgray,Ggray,Bgray,Agray = imGray.getpixel((i,j))
            imGray.putpixel((i,j),(Rgray,Ggray,Bgray,0))
            # Color image
            Rcolor,Gcolor,Bcolor,Acolor = imColor.getpixel((i,j))
            imColor.putpixel((i,j),(Rcolor,Gcolor,Bcolor,0))

#====================================================================================================
# Save the results
imGray.save(fgray)
imColor.save(fcolor)

#====================================================================================================
# I DON'T THINK THE STUFF BELOW IS NEEDED
'''
# (x,y) coordinates of the Galactic center and the Sun
rcirc = imcrop.width         # "Radius" of the Milky Way [pixels]
xGC   = 0.5 * imcrop.width   # Galactic center x-location [pixels]
yGC   = 0.5 * imcrop.height  # Galactic center y-location [pixels]
xSun  = xGC                  # Sun x-location (in cropped image) [pixels]
ySun  = ySun - ycroptop      # Sun y-location (in cropped image) [pixels]

# Circle projection
xmin = -1.0  # Min x-value
xmax =  1.0  # Max x-value
ymin = -1.0  # Min y-value
ymax =  1.0  # Max y-value
Nx = imcrop.size[0]
Ny = imcrop.size[1]
x  = np.linspace(xmin, xmax, Nx)
y  = np.linspace(ymin, ymax, Ny)
lon_grid = np.zeros([Nx, Ny])  # [lonmin, lonmax] = [-180, +180] degrees
rad_grid = np.zeros([Nx, Ny])  # [radmin, radmax] = [-sqrt(2), +sqrt(2)]
rgb_grid = np.zeros([Nx, Ny, 3])

# Additive factor to transform from yGC to ySun
yGC2Sun = (ySun - yGC) / float(Ny)

# Loop through the pixel in the JPG image
for i in np.arange(Nx):
    for j in np.arange(Ny):
    
        # Check that we are inside the circle
        if ((x[i]**2 + y[j]**2) < 1.0):

            # Transform from (x,y) in the image to (xlon,ylon) with the Sun as the origin
            xlon = x[i]
            ylon = y[j] #+ yGC2Sun
            
            # The function arctan2(y,x) chooses the correct quadrant (see docs)
            lon = np.arctan2(ylon, xlon)
            
            # Convert from the arctan2() convention to Galactic longitude convention
            if (lon < 0): lon = np.abs(lon)
            if (lon > 0): lon = 360.0 - lon
            
            # Populate the grid of Galactic longitudes for the Milky Way image
            lon_grid[i,j] = lon
            
            # Populate the grid of radial coordinates for the Milky Way image
            rad = np.sqrt(xlon**2 + ylon**2)
            rad_grid[i,j] = rad

            # Populate the grid of (R,G,B) colors for the Milky Way image
            r, g, b = rgb_imcrop.getpixel((i,j))
            rgb_grid[i,j,:] = [r/255.0, g/255.0, b/255.0]
        
        # Set values to NaN if outside the circle
        else:
            lon_grid[i,j]   = np.nan
            rad_grid[i,j]   = np.nan
            rgb_grid[i,j,:] = np.nan

# Flatten the RGB grid for plotting with pcolormesh
colorTuple = tuple(np.array([rgb_grid[:,:,0].flatten(), rgb_grid[:,:,1].flatten(), rgb_grid[:,:,2].flatten()]).transpose().tolist())

# Output colorTuple to an HDF5 file
f = h5py.File(fh5, 'w')
f.create_dataset('lon_grid',   data=lon_grid)
f.create_dataset('rad_grid',   data=rad_grid)
f.create_dataset('rgb_grid',   data=rgb_grid)
f.create_dataset('colorTuple', data=colorTuple)
f.close()
'''
