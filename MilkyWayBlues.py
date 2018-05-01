import numpy as np
import pylab as plt
import os
import subprocess
import matplotlib.font_manager as fm
from astropy.io import ascii
from PIL import Image
from matplotlib import patches

'''
Purpose:
--------
Generate the visualization to accompany the Milky Way Blues sonification by Mark Heyer (UMass).

Notes:
------
PIL Image tips: http://matthiaseisen.com/pp/patterns/p0202/
'''

#====================================================================================================
# LOAD IN THE IMAGES AND DO SOME INITIALIZATIONS

# Milky Way Roadmap filenames (grayscale and color)
fgray = "/Users/salvesen/outreach/asom/milkyway/data/MilkyWayRoadmapGray.png"
fmway = "/Users/salvesen/outreach/asom/milkyway/data/MilkyWayRoadmapColor.png"

# Read in the Grayscale image
imGray = Image.open(fgray).convert('RGBA')

# Read in the Color image
imColor = Image.open(fmway).convert('RGBA')

# Get image width and height
w = imGray.width   # [pixels]
h = imGray.height  # [pixels]

# Galactic center (x,y) location (determined by eye in original image with SAOImage ds9)
xGC = 0.5 * w
yGC = 0.5 * h

# Sun (x,y) location (determined by eye in cropped image with SAOImage ds9)
xSun = 0.5 * w  # Sun x-location [pixels]
ySun = 716      # Sun y-location [pixels]

# Radial extent of the Milky Way image
R = 0.5 * w

# Unit conversion between pixels and kpc
Dkpc = 8.4         # Distance from Sun to Galactic center [kpc]
Dpix = yGC - ySun  # Distance from Sun to Galactic center [pixels]
kpc2pix = Dpix / Dkpc  # [pixels kpc^-1]

#====================================================================================================
# PLOTTING HELPER FUNCTIONS

#----------------------------------------------------------------------------------------------------
# For YouTube production, determine the appropriate (xsize,ysize) for the figures (16:9 aspect)
def XYsize(pixres, dpi):
    '''
    Inputs:
    -------
    pixres - Resolution (240, 360, 720, 1080, etc.)
    dpi    - Dots per inch
    '''
    aspect = 16.0 / 9.0
    Nxpix  = float(pixres) * aspect
    Nypix  = float(pixres)
    xsize  = Nxpix / dpi
    ysize  = Nypix / dpi
    return xsize, ysize

#----------------------------------------------------------------------------------------------------
# For an input Galactic longitude, get the corresponding (x,y) on the edge of the Milky Way circle
def GalLon(lon=0.0):

    # Use the law of cosines to figure out the distance to the edge of a circle
    acoeff = 1.0
    bcoeff = -2.0 * (yGC - ySun) * np.cos(lon * np.pi/180.0)
    #bcoeff = -2.0 * (yGC - ySun) * np.cos((360.0 - lon) * np.pi/180.0)
    ccoeff = (yGC - ySun)**2 - R**2
    coeffs = [acoeff, bcoeff, ccoeff]
    roots  = np.roots(coeffs)

    # Choose the first positive root
    if (roots[0].real >= 0): c = roots[0].real
    elif (roots[1].real >= 0): c = roots[1].real

    # Get the unit vector direction
    dx = np.cos((lon + 90.0) * np.pi/180.0)  # Unit vector x-pointing
    dy = np.sin((lon + 90.0) * np.pi/180.0)  # Unit vector y-pointing
    
    # Calculate the (x,y) location on the circle
    xlon = c * dx + R
    ylon = c * dy + ySun

    return xlon, ylon

#----------------------------------------------------------------------------------------------------
# Convert an input Galactic longitude [degrees] and distance from the Sun [kpc] into (x,y) [pixels]
def ld2xy(lon, d):

    # Convert distance from [kpc] to [pixels]
    dpix = d * kpc2pix  # [pixels]

    # Get the unit vector direction
    dx = np.cos((lon + 90.0) * np.pi/180.0)  # Unit vector x-pointing
    dy = np.sin((lon + 90.0) * np.pi/180.0)  # Unit vector y-pointing

    # Calculate the (x,y) location on the circle
    xpix = dpix * dx + R
    ypix = dpix * dy + ySun

    return xpix, ypix

#----------------------------------------------------------------------------------------------------
# Convert an input (x,y) location [pixels] into an angle for plotting the arcs
def arcangle(x, y, x0=0, y0=0):
    X = x - x0
    Y = y - y0
    theta = np.arctan2(Y, X) * 180.0 / np.pi
    return theta

#----------------------------------------------------------------------------------------------------
# From the color Milky Way full image, excise a square of side length "r" [pixels] centered on the
# specified location (lon, d) and gray-mask outside of a circular region of radius "r".
def imgMask(img, lon, d, r=50.0, circle=True, xyPoly=None):

    # Determine the (x,y) pixel location corresponding to the input (lon,d) coordinates
    xpix, ypix = ld2xy(lon=lon, d=d)

    # Crop out the square region centred on (xpix,ypix)
    imgCrop = img.crop((xpix-r,h-ypix-r,xpix+r,h-ypix+r))  # Careful: Need to convert y!

    # Apply a grayscale mask outside of a circular region in the cropped image
    Nx   = imgCrop.width
    Ny   = imgCrop.height
    xmin = -1.0  # Min x-value
    xmax =  1.0  # Max x-value
    ymin = -1.0  # Min y-value
    ymax =  1.0  # Max y-value
    xarr = np.linspace(xmin, xmax, Nx)
    yarr = np.linspace(ymin, ymax, Ny)
    # Loop over all pixels in the cropped image
    for i in np.arange(Nx):
        for j in np.arange(Ny):
            # Get the (R,G,B,A) values for this pixel
            R,G,B,A = imgCrop.getpixel((i,j))
            # Create a circular mask (color inside, invisible outside)
            if ((xarr[i]**2 + yarr[j]**2) > 1.0): imgCrop.putpixel((i,j),(R,G,B,0))
            # Convert (R,G,B) to grayscale
            #Y = int(0.299*R + 0.587*G + 0.114*B)
            # Create a circular mask (color inside, grayscale outside)
            #if ((xarr[i]**2 + yarr[j]**2) > 1.0): imgCrop.putpixel((i,j),(Y,Y,Y,A))

    return imgCrop, xpix, ypix

#----------------------------------------------------------------------------------------------------
# Extract a regular polygon from an image, given the (x,y) locations of its vertices
# https://stackoverflow.com/questions/22588074/polygon-crop-clip-using-python-pil
'''
# convert to numpy (for convenience)
imArray = numpy.asarray(im)
# create mask
polygon = [(444,203),(623,243),(691,177),(581,26),(482,42)]
maskIm = Image.new('L', (imArray.shape[1], imArray.shape[0]), 0)
ImageDraw.Draw(maskIm).polygon(polygon, outline=1, fill=1)
mask = numpy.array(maskIm)
# assemble new image (uint8: 0-255)
newImArray = numpy.empty(imArray.shape,dtype='uint8')
# colors (three first columns, RGB)
newImArray[:,:,:3] = imArray[:,:,:3]
# transparency (4th column)
newImArray[:,:,3] = mask*255
'''

#----------------------------------------------------------------------------------------------------
# Return the area "A" of an n-sided regular polygon of circumradius "r"
def PolygonArea(r, n):
    A = 0.5 * n * r**2 * np.sin(2.0 * np.pi / n)
    return A

#----------------------------------------------------------------------------------------------------
# Return the circumradius "r" of an n-sided regular polygon of area "A"
def PolygonRadius(A, n):
    r = np.sqrt(2.0 * A / (n * np.sin(2.0 * np.pi / n)))
    return r

#----------------------------------------------------------------------------------------------------
# Return the (x,y) vertices of an n-sided regular polygon of circumradius "r" centered on the point (x0, y0)
def PolygonVertices(r, n, x0=0.0, y0=0.0):
    xyPoly = []
    for i in np.arange(n):
        x = x0 + r * np.cos(2.0 * np.pi * i / n)
        y = y0 + r * np.sin(2.0 * np.pi * i / n)
        xyPoly.append((x, y))
    return xyPoly

#----------------------------------------------------------------------------------------------------
# Map the duration (x) of a note to the circumradius (r) of the n-sided regular polygon (or circle)
# NOTE: Rest notes are mapped to negative radii
# NOTE: Only works for 1D input "x" array or list --> returns an array
def dur2rad(x, rmin=25.0, circle=True, n=None):
    x    = np.array(x)
    ipos = np.where(x >= 0)[0]
    ineg = np.where(x < 0)[0]
    Nx   = np.size(x)
    r    = np.zeros(Nx)
    
    # Calculate the radius of the circle with duration x
    if (circle is True):
        # Minimum area "Amin" corresponding to a circle with radius "rmin" and duration x=1
        Amin = np.pi * rmin**2  # [pixels^2]
        if (len(ipos) > 0): r[ipos] = np.sqrt(Amin * x[ipos] / np.pi)  # [pixels]
        if (len(ineg) > 0): r[ineg] = -1

    # Calculate the circumradius of the regular polygon with duration x
    if (circle is False):
        # Minimum area "Amin" corresponding to a regular polygon with circumradius "rmin" and duration x=1
        Amin = PolygonArea(r=rmin, n=n) # [pixels^2]
        if (len(ipos) > 0): r[ipos] = PolygonRadius(A=Amin*x[ipos], n=n)  # [pixels]
        if (len(ineg) > 0): r[ineg] = -1

    return r

#====================================================================================================
# PLOT A FRAME

# Plotting preferences
dpi    = 200   # Resolution [dots per inch]
pixres = 1080  # Number of y-pixels
fres   = pixres / 1080.0
lw     = 0.5 * fres
fsL    = 48 * fres
fsS    = 14 * fres
fsXS   = 8 * fres
left, right, bottom, top = 0.10, 0.90, 0.10, 0.90
xsize, ysize = XYsize(pixres=pixres, dpi=dpi)
plt.rcParams['axes.linewidth'] = lw
plt.rcParams['axes.edgecolor'] = 'black'
deg2rad = np.pi / 180.0

# Fonts
#jsansR  = fm.FontProperties(fname='../fonts/JosefinSans-Regular.ttf')
#jsansB  = fm.FontProperties(fname='../fonts/JosefinSans-Bold.ttf')
#jsansSB = fm.FontProperties(fname='../fonts/JosefinSans-SemiBold.ttf')
#jsansBI = fm.FontProperties(fname='../fonts/JosefinSans-BoldItalic.ttf')
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.bf'] = 'Arial:bold'
plt.rcParams.update({'font.family': 'Arial'})
plt.rc('font', family='Arial')

# Colors corresponding to the instrument parts [1,2,3,4]
#colors = ['red', 'yellow', 'blue', 'green']
rgb     = [(240,0,255,255), (255,231,0,255), (77,238,234,255), (116,238,21,255)]  # http://www.color-hex.com/color-palette/8618
rgbA    = [(240,0,255,127), (255,231,0,127), (77,238,234,127), (116,238,21,127)]  # http://www.color-hex.com/color-palette/8618
#rgb     = [(255,255,255,255), (255,255,255,255), (255,255,255,255), (255,255,255,255)]
#rgbA    = [(255,255,255,255), (255,255,255,255), (255,255,255,255), (255,255,255,255)]
#rgb     = [(0,0,0,255), (255,231,0,255), (255,255,255,255), (116,238,21,255)]  # http://www.color-hex.com/color-palette/8618
#rgbA    = [(0,0,0,127), (255,231,0,127), (255,255,255,127), (116,238,21,127)]  # http://www.color-hex.com/color-palette/8618
colors  = []
colorsA = []
Ncolors = len(rgb)
for i in np.arange(Ncolors):
    colors.append(tuple(ti/255.0 for ti in rgb[i]))
    colorsA.append(tuple(ti/255.0 for ti in rgbA[i]))

SunYellow = tuple(ti/255.0 for ti in (255,231,0,127))
GrayTrans = (0.5,0.5,0.5,0.5)

# Linestyles corresponding to the instrument parts []
# NOTE: Parts 1/2 always overlap in longitude, so choose one solid, one dashed
linestyles = ['solid','dashed','solid','solid']
zorders    = [1, 2, 1, 1]

# Shapes corresponding to the instrument parts
#nPoly = [0, 4, 3, 5]  # [circle, square, triangle, pentagon]
nPoly = [0, 0, 0, 0]  # [circle, circle, circle, circle]

# "Clock" hour (x,y) locations
Nhours = 12
hours  = np.linspace(1, Nhours, Nhours)
lonhr  = np.linspace(360, 0, Nhours+1)[1:Nhours+1]
xpixhr = np.zeros(Nhours)
ypixhr = np.zeros(Nhours)
texthr = []
for i in np.arange(Nhours):
    xpixhr[i], ypixhr[i] = GalLon(lon=lonhr[i])
    texthr.append(str(int(hours[i])))

# Plot the frame, include all parts at a given moment in time ("granular" duration)
def plot_frame(fout, lon=None, d=None, r=None, v=None, lon_meas=None, d_meas=None, r_meas=None, v_meas=None, lon_all=None, d_all=None, r_all=None, v_min=-150.0, v_max=150.0, MWcolor=False):
    '''
    Inputs:
    -------
    fout - Output filename
    lon  - 1D array of Galactic longitudes for each instrumental part [degrees]
    d    - 1D array of distances for each instrumental part [kpc]
    r    - 1D array of radii (circle glyph size) for each instrumental part [pixels]
    v    - 1D array of velocities for each instrumental part [km s^-1]
    lon_meas - 2D list of lists of Galactic longitudes for each part <-- Longitudes for a given part should all be the same
    d_meas   - 2D list of lists of distances for each part
    r_meas   - 2D list of lists of radii for each part
    v_meas   - 2D list of lists of velocities for each part
    v_min    - Minimum velocity (used for color range) [km s^-1]
    v_max    - Maximum velocity (used for color range) [km s^-1]

    Notes:
    ------
    When (lon_meas, d_meas, r_meas) are supplied, all of the notes in the measure will be displayed.
    Also, (lon_meas, d_meas, r_meas) are generated directly from the original data files...
    ...while (lon, d, r) are generated by matching the original data files to a common "granularity" in duration
    '''
    # Initialize the figure and axes objects
    fig = plt.figure(figsize=(xsize, ysize), dpi=dpi)
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
    fig.patch.set_facecolor('black')
    ax = fig.add_subplot(111)

    # Suppress everything having to do with the x- and y-axes
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.axis('off')
    
    # Set the axes limits
    axpad = 5  # Allow a little extra room so that the line forming the Milky Way circle is not cutoff
    ax.set_xlim([-axpad,w+axpad])
    ax.set_ylim([-axpad,h+axpad])

    # Plot the title and credits
    MWBtxt = "Milky" + "\n" + "Way" + "\n" + "Blues"
    Ctxt = "Sonification: " + u"\u00A9" + " Mark Heyer (UMass)" + "\n" + \
           "Visualization: Greg Salvesen (UCSB)" + "\n" + \
           "Image: Robert Hurt (JPL/NASA)" + "\n" + \
           "Data: Anderson et al. (2011); Kalberla et al. (2005); Dame et al. (2001)"
    fig.text(0.133, 0.5, MWBtxt, fontsize=fsL, color='white', ha='center', va='center')#, fontproperties=jsansSB)
    fig.text(0.05, 0.1, Ctxt, fontsize=fsXS, color=GrayTrans, ha='left', va='center')#, fontproperties=jsansSB)
    
    # Plot the legend
    emdash  = u"\u2014"
    bigskip = "\n\n\n\n\n"
    L1txt = "Acoustic Bass / Atomic"
    L2txt = bigskip + "Woodblocks / Molecular"
    L3txt = bigskip + bigskip + "Saxophone / Ionized"
    L4txt = bigskip + bigskip + bigskip + "Piano / Molecular"
    fig.text(0.85, 0.65, L1txt, fontsize=fsS, color=colors[0], ha='center', va='center')#, fontproperties=jsansSB)
    fig.text(0.85, 0.65, L2txt, fontsize=fsS, color=colors[1], ha='center', va='center')#, fontproperties=jsansSB)
    fig.text(0.85, 0.65, L3txt, fontsize=fsS, color=colors[2], ha='center', va='center')#, fontproperties=jsansSB)
    fig.text(0.85, 0.65, L4txt, fontsize=fsS, color=colors[3], ha='center', va='center')#, fontproperties=jsansSB)
    #ax.annotate(Ctxt, (0.05,0.1), xycoords='figure fraction', fontsize=fsXS, color=GrayTrans, ha='left', va='center', fontproperties=jsansSB)

    # Plot the Milky Way image
    if (MWcolor is False): ax.imshow(imGray, origin='lower', extent=[0,w,h,0], alpha=0.5)  # Notice y-extent is flipped
    if (MWcolor is True): ax.imshow(imColor, origin='lower', extent=[0,w,h,0], alpha=0.5)  # Notice y-extent is flipped

    # Plot a circle around the Milky Way image
    circle = plt.Circle((xGC, yGC), R, color=GrayTrans, linewidth=lw*2, fill=False, zorder=3)
    ax.add_artist(circle)

    # Plot the hours of the "clock"
    #for i in np.arange(Nhours):
    #    ax.annotate(texthr[i], (xpixhr[i],ypixhr[i]), xycoords='data', fontsize=fsXS, color='w', ha='center', va='center', fontproperties=jsansSB, zorder=4)
    
    # Plot the Sun
    markSize  = 100
    markColor = SunYellow #'white'
    markShape = '*'
    ax.scatter(xSun, ySun, s=markSize, c=markColor, marker=markShape, alpha=1, linewidths=0, zorder=4)
    ax.scatter(xSun, ySun, s=markSize, c=markColor, marker=markShape, alpha=1, linewidths=lw*2, zorder=4, facecolors='none', edgecolors='black')

    # Plot a circle of the Sun's orbit around the Galactic Center
    rSun    = 8.5 * kpc2pix
    circSun = plt.Circle((xGC, yGC), rSun, color=GrayTrans, linewidth=lw, fill=False, zorder=0, alpha=1)
    ax.add_artist(circSun)

    # Save the frame and return
    if (lon is None):
        fig.savefig(fout, facecolor=fig.get_facecolor(), edgecolor='none', dpi=dpi)
        return
    
    # Loop through instrument parts
    Nparts = np.size(lon)
    for i in np.arange(Nparts):

        # This longitude, distance, radius, velocity
        this_lon = lon[i]
        this_d   = d[i]
        this_r   = r[i]
        this_v   = v[i]

        # This shape (i.e., instrument)
        this_nPoly = nPoly[i]

        # Plot the Galactic longitude line
        xlon, ylon = GalLon(lon=this_lon)
        if (linestyles[i] == 'solid'):
            ax.plot([xSun,xlon], [ySun,ylon], color=colors[i], linestyle=linestyles[i], linewidth=lw, zorder=zorders[i])
        if (linestyles[i] == 'dashed'):
            ax.plot([xSun,xlon], [ySun,ylon], color=colors[i], linestyle=linestyles[i], linewidth=lw, zorder=zorders[i], dashes=(6, 6))


        # Plot the circle, unless it is a rest note
        if (this_r > 0):
        
            # Plot the masked/highlighted region
            imgCrop, xpix, ypix = imgMask(img=imGray, lon=this_lon, d=this_d, r=this_r)
            extent = [xpix-this_r, xpix+this_r, ypix+this_r, ypix-this_r]  # Notice the y-extent is flipped
            ax.imshow(imgCrop, origin='lower', extent=extent, zorder=5)
            
            # Plot a shape in color around the masked/highlighted region to indicate the instrument
            Nlw = 2.5
            if (this_nPoly < 3): shape = plt.Circle((xpix, ypix), this_r+(5*lw*Nlw), color=colors[i], linewidth=lw*Nlw, fill=False, zorder=7)
            if (this_nPoly >= 3): shape = patches.RegularPolygon((xpix, ypix), numVertices=this_nPoly, radius=this_r, orientation=this_lon*deg2rad, color=colors[i], linewidth=lw*2, fill=False, zorder=7)
            #shape = plt.Circle((xpix, ypix), this_r, color=colors[i], linewidth=lw*2, fill=False, zorder=7)
            ax.add_artist(shape)
            
            rim1 = plt.Circle((xpix, ypix), this_r, color='black', linewidth=lw, fill=False, zorder=7)
            rim2 = plt.Circle((xpix, ypix), this_r+2*(5*lw*Nlw), color='black', linewidth=lw, fill=False, zorder=7)
            ax.add_artist(rim1)
            ax.add_artist(rim2)

            # Fill the shape with color to indicate the velocity (w.r.t. local standard of rest)
            this_ic   = int((this_v - v_min) / (v_max - v_min) * 255)
            this_c    = list(plt.cm.RdBu_r(this_ic))
            this_c[3] = 1.0
            if (this_nPoly < 3): shapeFill = plt.Circle((xpix, ypix), this_r, color=tuple(this_c), linewidth=0, fill=True, zorder=6)
            if (this_nPoly >= 3): shapeFill = patches.RegularPolygon((xpix, ypix), numVertices=this_nPoly, radius=this_r, orientation=this_lon*np.pi/180.0, color=tuple(this_c), linewidth=0, fill=True, zorder=6)
            #shapeFill = plt.Circle((xpix, ypix), this_r, color=colorsA[i], linewidth=0, fill=True, zorder=6)
            ax.add_artist(shapeFill)

            '''
            # Plot the orbit arc, unless the note lies at the Galactic center
            this_dpix  = this_d * kpc2pix  # [pixels]
            this_rGal  = np.sqrt(this_dpix**2 + rSun**2 - 2.0 * this_dpix * rSun * np.cos(this_lon * np.pi/180.0))  # [pixels]
            this_angle = arcangle(x=xpix, y=ypix, x0=xGC, y0=yGC)  # [degrees]
            if (this_rGal > 0):
                arc = patches.Arc(xy=(xGC, yGC), width=2*this_rGal, height=2*this_rGal, angle=this_angle, theta1=-22.5, theta2=22.5, color=colors[i], linewidth=lw*2, fill=False, zorder=4, alpha=1)
                ax.add_patch(arc)
            '''

        # For the current part "i"...
        # Loop through each note and plot the masked/highlighted region in grayscale
        if (lon_meas is not None):
            # Number of notes in the measure for the current part "i"
            Nnotes = np.size(lon_meas[i])
            # This shape (i.e., instrument)
            nPoly_now = nPoly[i]
            for j in np.arange(Nnotes):
                # This longitude, distance, radius, velocity
                lon_now = lon_meas[i][j]
                d_now   = d_meas[i][j]
                r_now   = r_meas[i][j]
                v_now   = v_meas[i][j]
                # Plot the circle, unless it is a rest note
                if (r_now > 0):
                    # Plot the masked/highlighted region in grayscale
                    imgCrop, xpix, ypix = imgMask(img=imGray, lon=lon_now, d=d_now, r=r_now)
                    extent = [xpix-r_now, xpix+r_now, ypix+r_now, ypix-r_now]
                    ax.imshow(imgCrop, origin='lower', extent=extent, zorder=2)
                    # Fill the shape with color to indicate the velocity (w.r.t. local standard of rest)
                    ic_now   = int((v_now - v_min) / (v_max - v_min) * 255)
                    c_now    = list(plt.cm.RdBu_r(ic_now))
                    c_now[3] = 0.5
                    if (nPoly_now < 3): shapeFill = plt.Circle((xpix, ypix), r_now, color=tuple(c_now), linewidth=0, fill=True, zorder=3)
                    if (nPoly_now >= 3): shapeFill = patches.RegularPolygon((xpix, ypix), numVertices=nPoly_now, radius=r_now, orientation=lon_now*np.pi/180.0, color=tuple(c_now), linewidth=0, fill=True, zorder=3)
                    ax.add_artist(shapeFill)
                    # Plot a shape around the masked/highlighted region
                    if (nPoly_now < 3): shape = plt.Circle((xpix, ypix), r_now, color=colors[i], linewidth=lw, fill=False, zorder=3)
                    if (nPoly_now >= 3): shape = patches.RegularPolygon((xpix, ypix), numVertices=nPoly_now, radius=r_now, orientation=lon_now*np.pi/180.0, color=colors[i], linewidth=lw, fill=False, zorder=3)
                    ax.add_artist(shape)


        # For the current part "i"...
        # Plot all masked/highlighted circle regions so far in color
        # !!! NOTE: This approach replots everything every time.
        # !!! NOTE: Better approach is to clear only certain plotting elements and not close the plot environment.
        if (lon_all is not None):
            Nall = len(lon_all)  # <-- "lon_all" is a list of lists: [[lon1part1,lon1part2,lon1part3,lon1part4], [lon2part1,lon2part2,lon2part3,lon2part4], ...]
            for j in np.arange(Nall):
                # This longitude, distance, radius
                lon_now = lon_all[j][i]
                d_now   = d_all[j][i]
                r_now   = r_all[j][i]
                if (r_now > 0):
                    # Plot the region residue in color
                    imgCrop, xpix, ypix = imgMask(img=imColor, lon=lon_now, d=d_now, r=r_now)
                    extent = [xpix-r_now, xpix+r_now, ypix+r_now, ypix-r_now]
                    ax.imshow(imgCrop, origin='lower', extent=extent, zorder=2)
    

    # Save the figure
    fig.savefig(fout, facecolor=fig.get_facecolor(), edgecolor='none', dpi=dpi)

    # Close the plot object
    plt.close()

#====================================================================================================
# HELPER FUNCTIONS FOR WORKING WITH MULTIPLE INSTRUMENTAL PARTS

#----------------------------------------------------------------------------------------------------
# Calculate the number of notes per measure for a given instrumental part
def notes_per_measure(fpart):
    divs  = int(os.popen("grep Divisions " + fpart).read().split()[-1])
    beats = int(os.popen("grep Beats " + fpart).read().split()[-1])
    npm   = divs * beats
    return npm

#----------------------------------------------------------------------------------------------------
# Determine the duration scaling factor for an instrumental part
def scale_dur(fpart, fpart_list):

    # For the current part, determine the number of notes per measure
    npm_part = notes_per_measure(fpart)

    # Considering all parts, determine the maximum number of notes per measure
    Nparts  = np.size(fpart_list)
    npm_max = 0
    for f in fpart_list:
        npm = notes_per_measure(fpart=f)
        if (npm > npm_max): npm_max = npm

    # Calculate the duration scaling factor
    fscl = int(npm_max / npm_part)
    return fscl

#----------------------------------------------------------------------------------------------------
# Output the original data for all parts in 2D lists [[part1],[part2],...]
# NOTE: We do a 2D list of lists because each part is of different length (arrays will not work)
def original_granularity(fpart_list):

    # Number of instrumental parts
    Nparts = np.size(fpart_list)

    # Initialize the relevant data lists to be returned
    measure_list   = []
    longitude_list = []
    distance_list  = []
    duration_list  = []
    velocity_list  = []

    # Loop through each instrumental part
    for i in np.arange(Nparts):

        # Filename for the current part
        fpart = fpart_list[i]

        # Collect the relevant data for the current part
        data      = ascii.read(fpart)
        measure   = data['Measure']
        longitude = data['X']
        distance  = data['Distance']
        duration  = data['Duration']
        velocity  = data['Vlsr']

        # Collect the results for the current part
        measure_list.append(measure.tolist())
        longitude_list.append(longitude.tolist())
        distance_list.append(distance.tolist())
        duration_list.append(duration.tolist())
        velocity_list.append(velocity.tolist())

    # Return the results
    return measure_list, longitude_list, distance_list, duration_list, velocity_list

#----------------------------------------------------------------------------------------------------
# Match all instrumental parts according to a common duration or "granularity"
def common_granularity(fpart_list):

    # Sanity Check: All parts must have the same number of measures "Nmeas"
    data  = ascii.read(fpart_list[0])
    Nmeas = data['Measure'][-1]
    for f in fpart_list:
        data   = ascii.read(f)
        Nmeasf = data['Measure'][-1]
        if (Nmeasf != Nmeas):
            print "\nERROR: The data files must all have the same number of measures.\n"
            quit()

    # Get the maximum number of notes per measure
    npm_max = 0
    for f in fpart_list:
        npm = notes_per_measure(fpart=f)
        if (npm > npm_max): npm_max = npm

    # Calculate the total number of notes at the finest "granularity"
    Nnotes = Nmeas * npm_max

    # Number of instrumental parts
    Nparts = np.size(fpart_list)

    # Initialize the relevant data arrays to be returned
    measure_arr   = np.zeros([Nparts, Nnotes])
    longitude_arr = np.zeros([Nparts, Nnotes])
    distance_arr  = np.zeros([Nparts, Nnotes])
    duration_arr  = np.zeros([Nparts, Nnotes])
    velocity_arr  = np.zeros([Nparts, Nnotes])

    # Loop through each instrumental part
    for i in np.arange(Nparts):

        # Filename for the current part
        fpart = fpart_list[i]

        # Collect the relevant data for the current part
        data      = ascii.read(fpart)
        measure   = data['Measure']
        longitude = data['X']
        distance  = data['Distance']
        duration  = data['Duration']
        velocity  = data['Vlsr']

        # Find duration scale factor for the current part
        fscl = scale_dur(fpart=fpart, fpart_list=fpart_list)
        
        # Remap the relevant data to the finest granularity
        Nnotes_part     = len(duration)       # Number of notes (lines) in the data file for the current part
        measure_rescl   = []
        longitude_rescl = []
        distance_rescl  = []
        duration_rescl  = []
        velocity_rescl  = []
        for j in np.arange(Nnotes_part):        # Loop through each note (line) for the current part
            abs_duration = np.abs(duration[j])  # <-- This allows us to handle durations of rest notes (i.e., negative durations)
            for k in np.arange(abs_duration):   # Repeat each note a number of times according to its coarse duration
                for l in np.arange(fscl):       # Repeat each note an additional number of times according to the fine granularity scaling
                    measure_rescl.append(measure[j])
                    longitude_rescl.append(longitude[j])
                    distance_rescl.append(distance[j])
                    duration_rescl.append(duration[j] * fscl)
                    velocity_rescl.append(velocity[j])

        # Collect the results for the current part
        measure_arr[i,:]   = measure_rescl
        longitude_arr[i,:] = longitude_rescl
        distance_arr[i,:]  = distance_rescl
        duration_arr[i,:]  = duration_rescl
        velocity_arr[i,:]  = velocity_rescl

    # Return the results
    return measure_arr, longitude_arr, distance_arr, duration_arr, velocity_arr

#====================================================================================================
# MERGE THE SONIFICATION AND VISUALIZATION (A.K.A. ffmpeg WITCHCRAFT)

#----------------------------------------------------------------------------------------------------
# Get the duration of an input audio, video, movie file in [seconds]
def duration(fin):
    # Here is the most disgusting Unix command I have ever seen in my life
    cmd = "ffmpeg -i " + fin + \
          " 2>&1 | grep 'Duration'| cut -d ' ' -f 4 | sed s/,// |" + \
          " awk '{ split($1, A, " + '":"' + ");" + \
          " print 3600*A[1] + 60*A[2] + A[3] }'"
    tsec = float(os.popen(cmd).read())
    return tsec

#----------------------------------------------------------------------------------------------------
# Cut an audio file from time t0 to tf [seconds]
def cut_audio(fout, faud, t0, tf):
    '''
    Inputs:
    -------
    fout - Output filename for the cut audio file (AIFF, MP3)
    faud - Audio file (AIFF, MP3)
    t0   - Start time [seconds]
    tf   - End time [seconds]
    '''
    # Remove any pre-existing cut audio file
    cmd_rm = "rm " + fout
    subprocess.call(cmd_rm, shell=True)

    # Cut the audio file
    cmd = "ffmpeg -ss " + str(t0) + " -t " + str(tf) + " -i " + faud + " " + fout
    subprocess.call(cmd, shell=True)

#----------------------------------------------------------------------------------------------------
# Combine all static frames into a video
def make_video(fout, froot, fps):
    '''
    Inputs:
    -------
    fout  - Output filename for the video (MP4)
    froot - Root name (with path) for the static frames (PNG) to stitch together into a video
            *Do not include the assumed suffix of froot, which must follow the format: 0000.png
    fps   - Frames per second
    '''
    # Collect the width X height for the video
    img = Image.open(froot + "0000.png")
    Nx  = str(int(img.size[0]))  # [pixels]
    Ny  = str(int(img.size[1]))  # [pixels]

    # Remove any pre-existing video file
    cmd_rm = "rm " + fout
    subprocess.call(cmd_rm, shell=True)

    # Make the video
    cmd = "ffmpeg -f image2 -framerate " + str(fps) + " -i '" + froot + "%04d.png' -s " + Nx + "X" + Ny + " -pix_fmt yuv420p " + fout
    subprocess.call(cmd, shell=True)

#----------------------------------------------------------------------------------------------------
# Merge the audio and video
def make_movie(fout, faud, fvid):
    '''
    Inputs:
    -------
    fout - Output filename for the movie (MP4)
    faud - Audio file (AIFF, MP3)
    fvid - Video file (MP4)
    '''
    # Remove any pre-existing movie file
    cmd_rm = "rm " + fout
    subprocess.call(cmd_rm, shell=True)

    # Make the movie
    cmd = "ffmpeg -i " + fvid + " -i " + faud + " -c:v copy -c:a aac -strict experimental " + fout
    subprocess.call(cmd, shell=True)

#====================================================================================================
# CREATE THE VISUALIZATION

# Set input/output filenames
faud  = "/Users/salvesen/outreach/asom/milkyway/data/MWB_V3.2.aiff"
fcut  = "/Users/salvesen/outreach/asom/milkyway/data/MWB_V3.2_cut.aiff"
dout  = "/Users/salvesen/outreach/asom/milkyway/results/"
froot = dout + "MWB_V3.2_"
fvid  = dout + "MWB_V3.2_Video.mp4"
fmov  = dout + "MWB_V3.2_circles_colored.mp4"

# Minimum circle radis [pixels]
rmin = 25.0

# Min/max velocities for color scaling
v_min = -100.0  # [km s^-1]
v_max = 100.0   # [km s^-1]

# Data files for each instrumental part
ddir   = "/Users/salvesen/outreach/asom/milkyway/data/"
fpart1 = ddir + "MWB_V3.2.1.txt"  # HI 21 cm (Acoustic Bass)
fpart2 = ddir + "MWB_V3.2.2.txt"  # ^12CO (Wood Blocks)
fpart3 = ddir + "MWB_V3.2.3.txt"  # Radio Recombination Lines (Baritone Sax)
fpart4 = ddir + "MWB_V3.2.4.txt"  # ^12CO (Piano)
fpart_list = [fpart1, fpart2, fpart3, fpart4]

# Find duration scale factor for each part
Nparts    = len(fpart_list)
fscl_orig = []
for i in np.arange(Nparts):
    fscl = scale_dur(fpart=fpart_list[i], fpart_list=fpart_list)
    fscl_orig.append(fscl)

# Collect the data in its original granularity
measure_orig, longitude_orig, distance_orig, duration_orig, velocity_orig = original_granularity(fpart_list=fpart_list)

# Remap all parts to the finest granularity and collect the relevant data
# --> 2D arrays with dimensions [Nparts, Nnotes]
measure_arr, longitude_arr, distance_arr, duration_arr, velocity_arr = common_granularity(fpart_list=fpart_list)

# Map duration to radius of the masked/highlighted circles
# (Loop through all parts)
Nparts      = len(measure_arr[:,0])
Nnotes      = len(measure_arr[0,:])
radius_orig = []
radius_arr  = np.zeros([Nparts, Nnotes])
for i in np.arange(Nparts):
    x_orig = np.array(duration_orig[i]) * np.array(fscl_orig[i]) # <-- NOTE: Had to multiply duration by the approprite scale factor!
    x_arr  = duration_arr[i,:]
    radius_orig.append(dur2rad(x=x_orig, rmin=rmin).tolist())
    radius_arr[i,:] = dur2rad(x=x_arr, rmin=rmin)

# Total number of frames for the visualization
Nframes = len(measure_arr[0,:])

# Populate an alternating True/False list for each changing measure
MWcolor = []
for j in np.arange(Nframes):
    if (measure_arr[0,j] % 2 == 1): MWcolor.append(False)  # Odd measures --> Milky Way in black and white
    if (measure_arr[0,j] % 2 == 0): MWcolor.append(True)   # Even measures --> Milky Way in color

# Keep track of all longitudes, distances, radii played so far
# (These are used to plot the "residue" left behind after each note is played)
all_longitudes = []
all_distances  = []
all_radii      = []

# Loop through all frames (i.e., the finest granularity level for duration)
jcnt = 0
for j in np.arange(Nframes):
    
    # Current longitude, distance, circle radius, velocity (for each part)
    this_longitude = longitude_arr[:,j]  # [degrees]
    this_distance  = distance_arr[:,j]   # [kpc]
    this_radius    = radius_arr[:,j]     # [pixels]
    this_velocity  = velocity_arr[:,j]   # [km s^-1]

    # Up until and including now, all longitudes, distances, and circle radii (for each part)
    all_longitudes.append(this_longitude.tolist())
    all_distances.append(this_distance.tolist())
    all_radii.append(this_radius.tolist())

    # For the current measure, get the sets of Longitudes, Distances, Radii, Velocities for all parts
    this_measure = measure_arr[0,j]  # Measure number
    # OLD WAY: Caused the same circle to be overplotted many times
    # NEW WAY: Uses the original data files, rather than the data re-mapped to a common duration granularity
    # CAREFUL: np.where() does not work on lists
    longitude_set = []
    distance_set  = []
    radius_set    = []
    velocity_set  = []
    for i in np.arange(Nparts):
        # Find the indices corresponding to the current measure for the current part
        kset = np.where(np.array(measure_orig[i]) == this_measure)[0]
        # Pick out the longitudes, distances, radii contained in the current measure for the current part
        longitude_list = np.array(longitude_orig[i])[kset].tolist()
        distance_list  = np.array(distance_orig[i])[kset].tolist()
        radius_list    = np.array(radius_orig[i])[kset].tolist()
        velocity_list  = np.array(velocity_orig[i])[kset].tolist()
        # Collect longitudes, distances, radii, velocities for the current measure for ALL parts
        longitude_set.append(longitude_list)
        distance_set.append(distance_list)
        radius_set.append(radius_list)
        velocity_set.append(velocity_list)

    # Create the current frame
    ftag = str(jcnt).zfill(4)
    fout = froot + ftag + ".png"
    plot_frame(fout=fout, lon=this_longitude, d=this_distance, r=this_radius, v=this_velocity, lon_meas=longitude_set, d_meas=distance_set, r_meas=radius_set, v_meas=velocity_set, v_min=v_min, v_max=v_max, MWcolor=True)#MWcolor[j])
    #plot_frame(fout=fout, lon=this_longitude, d=this_distance, r=this_radius, v=this_velocity, lon_meas=longitude_set, d_meas=distance_set, r_meas=radius_set, v_meas=velocity_set, lon_all=all_longitudes, d_all=all_distances, r_all=all_radii, v_min=v_min, v_max=v_max, MWcolor=MWcolor[j])

    # Progress report
    pctComplete = float(jcnt+1) / float(Nframes) * 100
    print "    % Complete: ", '{:.1f}'.format(pctComplete)
    jcnt = jcnt + 1
    

# Tempo and duration of the original audio file
bpm   = 120                 # Tempo of the audio file [beats per minute] <-- 120 bpm is the default in Finale
bps   = bpm / 60.0          # Tempo of the audio file [beats per second]
tsec  = duration(fin=faud)  # Duration of the audio file [seconds]

# Number of measures in the original audio file, the visualization, and their difference
Nmeas_tot = tsec / bps             # Total number of measures in the audio file
Nmeas_vis = np.max(measure_arr)    # Number of measures in the visualization
Nmeas_end = Nmeas_tot - Nmeas_vis  # Number of extra rest measures tacked on to the end of the audio file

# Unseful information about beat and measure durations
bpmeas    = 4.0                                # Beats per measure (need to tweak things if != 4)
meas_gran = float(Nmeas_vis) / float(Nframes)  # Granularity of a measure (e.g., 0.0625 = 16th note)
beat_gran = meas_gran * bpmeas                 # Granularity of a beat (e.g., 0.25 --> 1 beat = 1 quarter note)

# Number of beats...
Nbeats_vis = beat_gran * Nframes      # ...in the entire visualization
Nbeats_tot = bps * tsec               # ...in the entire audio file
Nbeats_end = Nbeats_tot - Nbeats_vis  # ...tacked on to the end of the audio file as rests
Nbeats_add = np.floor(Nbeats_end)     # ...to add to the visuzalization
Nbeats_cut = Nbeats_end - Nbeats_add  # ...to cut off the end of the audio file

# Cut the audio file to exactly match an integer number of beats (finest granularity)
tcut  = Nbeats_cut / bps  # Time to cut from the end of the audio file [seconds]
tgood = tsec - tcut       # Time to keep from the beginning of the audio file [seconds]
cut_audio(fout=fcut, faud=faud, t0=0, tf=tgood)

# Make the additional frames
fpb   = Nframes / Nbeats_vis   # Frames per beat in the visualization
Nadd  = int(Nbeats_add * fpb)  # Number of frames to add to the end of the visualization
icnt  = int(Nframes)           # Initialize the counter
fcopy = froot + str(icnt-1).zfill(4) + ".png"  # <-- Copy this frame "Nadd" times
for i in np.arange(Nadd):
    ftag = str(icnt).zfill(4)
    fout = froot + ftag + ".png"
    cmd  = "cp " + fcopy + " " + fout
    subprocess.call(cmd, shell=True)
    icnt = icnt + 1

# Frames per second for the visualization
fps = float(Nframes + Nadd) / tgood  # [fps]

# Make the video
make_video(fout=fvid, froot=froot, fps=fps)

# Make the movie
make_movie(fout=fmov, faud=fcut, fvid=fvid)

# Create the first "empty" frame
# CAREFUL!!! <-- Adding extra frames at this stage screws up things later with matching notes to frames. Do this at the very end.
#jcnt = 0
#ftag = str(jcnt).zfill(4)
#fout = froot + ftag + ".png"
#plot_frame(fout=fout)
#jcnt = jcnt + 1

