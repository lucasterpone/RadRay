import os
import numpy as np 
import numpy.f2py as myf2py 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys 
import numpy
import random
import math
from matplotlib import pyplot
from matplotlib.colors import rgb2hex
from matplotlib.collections import PolyCollection
import matplotlib.gridspec as gridspec
from matplotlib import cm

def randrange(n, vmin, vmax):
    '''
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    '''
    return (vmax - vmin)*np.random.rand(n) + vmin

fig_1 = plt.figure('GDS 3D view')
ax_1 = fig_1.gca(projection='3d')
#ax_1.set_aspect("equal") #3D plot : geometrical representation

#
# Definition of points: to be defined before any elaboration
#
MAX_POINTS = 1000
MAX_CUBE = 300
MAX_LAYER = 300
RES_STEP = 10

if len (sys.argv) != 9:
	print('')
	print('Usage: ')
	print('RadRay_Tool [GDS txt filename] [Metals Values] [RadRay Energy] [T/C] [#] [V] [F] [HS]')
	print(' ')
	print('[GDS txt filename]: GDS-II Textual description of the Cell Under Analysis ')
	print('[Metal Values]: Data Values of the Metals with a Uniform Time Range ')
	print('[RadRay Energy]: Energy Profile for each Cell Layer ')
	print('[T/C]: Transient or Cumulative mode')
	print('[#] : Number of iterations ')
	print('[V] : Y for printing detailed PDF reports')
	print('[F] : V for 0 degree tilting (mimic radiation facility)')
	print('     R for random degree tilting (mimic realistic radiation environment)')
	print('     S for screen mode ')
	print('     Default Resolution 10 nm ')
	print('[HS]: Y enable the generation of the HSPICE simulation file')
	print('      Layer Number - Hspice test reference - Hspice reference port')
	
	
	
            
	
	
	exit()

#Initialization
cnt = 0
#Definition of layers cube LAYERS, CUBES, VERTICES, [0.x0,1.y0,2.z0]
layers_cube = np.zeros ((250,MAX_CUBE,100,6)) #250,MAX_CUBE,20,6 or normally #250,MAX_CUBE,50,6
cube_layer_id = np.zeros ((MAX_CUBE))
layer_energ_id = np.zeros ((MAX_LAYER))

#Reference for Hspice analysis
hspice_ref = []
hspice_sig = []
hspice_port = []
hspice_full_cmd = []


GEO_tlinex = []
GEO_tliney = []
GEO_tlinez = []
GEO_xlinex1 = []
GEO_yliney1 = []
GEO_zlinez1 = []
GEO_xlinex2 = []
GEO_yliney2 = []
GEO_zlinez2 = []



#DEFINITION OF THE 3D MESH FOR ELABORATION
max_cube_size = np.zeros ((MAX_CUBE,6))
#Maximal Cube Size 0.xmin 1.xmax 2.ymin 3.ymax 4.zmin 5.zmax
energy_values = np.zeros ((MAX_CUBE,MAX_POINTS))
data_points = np.zeros ((MAX_CUBE,MAX_POINTS,3))
example_values = np.zeros ((MAX_CUBE,MAX_POINTS))





#DEFINITION OF THE LAYERS TECHNOLOGICAL HEIGHT and THICKNESS - 45nm FreePDK
#[height, thickness]
layers_tech = np.zeros((300,2)) 
#layers_tech[255][0]= -900
#layers_tech[255][1]= 1200
layers_tech[1][0]= 0 #N-Well -300
layers_tech[1][1]= 520 #520
layers_tech[2][0]= 220 #P-Plus
layers_tech[2][1]= 100
layers_tech[3][0]= 220 #N-Plus
layers_tech[3][1]= 110
layers_tech[4][0]= 320 #Active
layers_tech[4][1]= 20
layers_tech[5][0]= 320 #Poly 
layers_tech[5][1]= 100
layers_tech[9][0]= 620 
layers_tech[9][1]= 150
layers_tech[10][0]= 770 #Via-1
layers_tech[10][1]= 200
layers_tech[11][0]= 970 #Metal-1
layers_tech[11][1]= 250
layers_tech[12][0]= 1220 #Via-2
layers_tech[12][1]= 200
layers_tech[13][0]= 1420 #Metal-2
layers_tech[13][1]= 250



layer_list = []
vertex_cube_list = []

en_xy = 0
index_cube = 0
index_vertex = 0
cu_layer = 0


print('--------------------------------------')
print('- Rad Ray Tool 3.0 - Politecnico di Torino  ')
print('- Author: Luca Sterpone ')
print('                                        []         ')
print('                                       []          ')
print('     ######              ######       []                     ')
print('     ##  ###             ##  ###     []                   ')
print('     ##   ##             ##   ##    []                  ')
print('     ##  ###             ##  ###   []                    ')
print('     ######              ######   []                      ')
print('     ####     ###  ####  ####    []## ##  ##                            ')
print('     ## ##   ## ## ## ## ## ##  []# ## ####                               ')
print('     ##  ##  ##### ## ## ##  ##[]#####  ##                             ')
print('     ##   ## ## ## ## ## ##   [] ## ##  ##                             ')
print('     ##    #### ## ####  ##* []#### ##  ##                           ')
print('                          **[]                    ')
print('   ======================*****================       ')
print('                        *******                     ')
print('   ======================*****================                 ')
print('                        []***  ')
print('                       []  *  ')
print('  FreePDK45nm Version []           ')
print('- ')
print('- v0.1: loading all GDS vertices')#24.08.2018
print('- v0.2: generation of cube values')
print('- v0.3: visualizaion of data values')
print('- v0.4: cumulative analysis ')#21.09.2018
print('- v0.5: transient analysis ')
print('- v0.7: cumulative energy computation')
print('- v1.0: maximal voltage peak computation')
print('- v2.0: generation of the transient pulse on GDS output layer')
print('- v2.1: generation of the amplitude pulse distribution ')
print('- v3.0: generation of the hspice simulation commands ')
print('- v3.1: inclusion of the FreePDK 45nm z-section ')
print('- v3.2: sensitivity heatmap XY view ')
print('- v3.3: sensitivity heatmap update ')#11.11.2020


print('--------------------------------------')

print('-- Elaboration: ',sys.argv[1])
print('-- Cube points: ', MAX_POINTS)
print('-- Reading GDS boundaries...')

#Update status 
#-------------
#24.08.2018 - Generation of graphical area




#Reading GDS TXT file
f = open(sys.argv[1],"r")
for line in f:
	line = line.strip()
	if line == 'BOUNDARY':
		in_b = 1
		cnt = cnt + 1
		while (in_b == 1):
			#print ('In boundary - ', cnt)
			data = f.readline() #Readling LAYER ID
			check = data.split(' ')
			data = data.strip()
			#print(data)
			#print(check)
			#Layer Storage
			if check[0] == 'LAYER':
				cu_layer = int(check[1])
				#add_lay = 1
				#for item in layer_list:
				#	if item == cu_layer:
				#		add_lay = 0
				#if add_lay == 1:
				layer_list.append(cu_layer)

			if data == 'ENDEL':
			    #A cube is closed
			    vertex_cube_list.append(index_vertex)
			    #print(layers_cube)
			    index_cube = index_cube + 1
			    index_vertex = 0
			    #print(index_vertex)
			    en_xy = 0
			    in_b = 0
			    #d = input()
			    #break

			if en_xy == 1:
			    #Inserting coordinates once already recognized XY
			    #print('-- Loading layer ', cu_layer)
			    	
			    layers_cube[cu_layer][index_cube][index_vertex][0] = check[0].replace(":","")
			    layers_cube[cu_layer][index_cube][index_vertex][1] = check[1].replace("\n","")
			    layers_cube[cu_layer][index_cube][index_vertex][2] = layers_tech[cu_layer][0]
			    layers_cube[cu_layer][index_cube][index_vertex][3] = check[0].replace(":","")
			    layers_cube[cu_layer][index_cube][index_vertex][4] = check[1].replace("\n","")
			    layers_cube[cu_layer][index_cube][index_vertex][5] = layers_tech[cu_layer][0] + layers_tech[cu_layer][1]
			    index_vertex = index_vertex + 1


			if check[0] == 'XY':
			    layers_cube[cu_layer][index_cube][index_vertex][0] = check[1].replace(":","")
			    layers_cube[cu_layer][index_cube][index_vertex][1] = check[2].replace("\n","")
			    layers_cube[cu_layer][index_cube][index_vertex][2] = layers_tech[cu_layer][0]
			    layers_cube[cu_layer][index_cube][index_vertex][3] = check[1].replace(":","")
			    layers_cube[cu_layer][index_cube][index_vertex][4] = check[2].replace("\n","")
			    layers_cube[cu_layer][index_cube][index_vertex][5] = layers_tech[cu_layer][0] + layers_tech[cu_layer][1]
			    index_vertex = index_vertex + 1
			    en_xy = 1

total_cube = index_cube		

print('-- Total Number of Cubes : ', total_cube)
print('-- Done')
print('-- Generating Layers Cube ')

index_cube = 0
for item in layer_list:
    #print ('The Current Layer is:',item)
    #print ('-- Cube ', index_cube, ' Vertices: ', vertex_cube_list[index_cube])
    cube_layer_id[index_cube] = item
    n = 0
    xline_1 = []
    yline_1 = []
    zline_1 = []
    xline_2 = []
    yline_2 = []
    zline_2 = []
    tlinex = []
    while n < vertex_cube_list[index_cube]:
    	#print (n, ' - vertex: ',layers_cube[item][index_cube][n])
    	tlinex = []
    	tliney = []
    	tlinez = []
    	xline_1.append(layers_cube[item][index_cube][n][0])
    	yline_1.append(layers_cube[item][index_cube][n][1])
    	zline_1.append(layers_cube[item][index_cube][n][2])
    	xline_2.append(layers_cube[item][index_cube][n][3])
    	yline_2.append(layers_cube[item][index_cube][n][4])
    	zline_2.append(layers_cube[item][index_cube][n][5])
    	tlinex.append(layers_cube[item][index_cube][n][0])
    	tliney.append(layers_cube[item][index_cube][n][1])
    	tlinez.append(layers_cube[item][index_cube][n][2])
    	tlinex.append(layers_cube[item][index_cube][n][3])
    	tliney.append(layers_cube[item][index_cube][n][4])
    	tlinez.append(layers_cube[item][index_cube][n][5])
    	ax_1.plot3D(tlinex, tliney, tlinez, 'gray')  #'gray'
    	n = n + 1
    ax_1.plot3D(xline_1, yline_1, zline_1, 'blue')
    ax_1.plot3D(xline_2, yline_2, zline_2, 'blue')  
    n = 0
    index_cube = index_cube + 1
    #d = input()

plt.show()

fig_2 = plt.figure('Static Signal 3D view')
fig_3 = plt.figure('Voltage Pulse Plot')
fig_t = plt.figure('Test Figure')
ax_2 = fig_2.gca(projection='3d')
ax_t = fig_t.gca(projection='3d')
#ax_2.set_aspect("equal") #3D plot : static signal view
#ax_t.set_aspect("equal") #3D plot : developing testing plot


print('-- Generation XYZ Longitudinal Distribution')

np.random.seed(564500)

# example data
mu = 0 # mean of distribution
sigma = 130  # standard deviation of distribution
#signa =  80 distribution -+220 Ang (for Ion C)
#sigma = 100 distribution -+300 Ang (for Ion Ar)
#signa = 130 distribution -+400 Ang (for Ion Kr)
#signa = 250 distribution -+800 Ang (for Ion Xe)
x = mu + sigma * np.random.randn(1000) #437

num_bins = 60

fig, ax = plt.subplots()
#fig, bx = plt.subplots()

# the histogram of the data
#n, bins, patches = ax.hist(x, num_bins, density=1)
n, bins, patches = ax.hist(x, num_bins, density=1)

# Energy Fit Lie
y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
     np.exp(-0.5 * (1 / sigma * (bins - mu))**2))

ener_scale_pos = bins / 10
ener_scale_val = (1 / max(y))*y

#print(ener_scale_pos)
#print(ener_scale_val)

ax.plot(bins, y, '--')
#bx.plot(bins, ener_scale_val, '-')
ax.set_xlabel('Ang')
ax.set_ylabel('Probability density')
ax.set_title('XY Ion Range Longitudinal Distribution')

min_x_ray = min(abs(ener_scale_pos))
max_x_ray = max(abs(ener_scale_pos))
min_y_ray = min(abs(ener_scale_pos))
max_y_ray = max(abs(ener_scale_pos))

print('-- Longitudinal Minimal and Maximal Distance (nm) ', '%.2f' % min_x_ray, '%.2f' % max_x_ray)



# Tweak spacing to prevent clipping of ylabel
fig.tight_layout()
fig.savefig('beam_ernergy_profile.pdf', format='pdf')
plt.close(fig)

print('-- Loading Energy Profile')
#Reading GDS TXT file
f = open(sys.argv[3],"r")
cnt = 0
for line in f:
	if cnt > 0:
		#print(line)
		check = line.split(' ')
		#print(int(check[0]))
		layer_energ_id[int(check[0])] = int(check[1])
		#print('here 1', check[0])
		#print('here 2', check[1])
	cnt = cnt + 1


print('-- Generating Data Points')

#Loading Hspice Profile
if sys.argv[8] == 'Y':
    print('-- Loading Hspice reference ports')
    h_spice_log = open("log_hspice_test.txt","a+")
    h = open(sys.argv[2],"r")
    cnt = 0
    for line in h:
        if cnt > 0:
            #print(line)
            check = line.split('-')
            #print(check)
            #print(len(check))
            total = len(check)-2
            for i in range(total):
                #print(i)
                hspice_ref.append(int(check[i]))
                hspice_sig.append(check[total])
            hspice_port.append(check[total])
        cnt = cnt + 1
    print('-- Hspice reference loaded ')
    #print(hspice_port)
    #print(hspice_ref)
    #print(hspice_sig)



i_ray = 0
index_cube = 0
for item in layer_list:
    #print ('The Current Layer is:',item)
    #print ('Cube ', index_cube, ' Vertices: ', vertex_cube_list[index_cube])
    n = 0
    xline_1 = []
    yline_1 = []
    zline_1 = []
    xline_2 = []
    yline_2 = []
    zline_2 = []
    tlinex = []
    while n < vertex_cube_list[index_cube]:
    	#print (n, ' - vertex: ',layers_cube[item][index_cube][n])
    	tlinex = []
    	tliney = []
    	tlinez = []
    	xline_1.append(layers_cube[item][index_cube][n][0])
    	yline_1.append(layers_cube[item][index_cube][n][1])
    	zline_1.append(layers_cube[item][index_cube][n][2])
    	xline_2.append(layers_cube[item][index_cube][n][3])
    	yline_2.append(layers_cube[item][index_cube][n][4])
    	zline_2.append(layers_cube[item][index_cube][n][5])
    	tlinex.append(layers_cube[item][index_cube][n][0])
    	tliney.append(layers_cube[item][index_cube][n][1])
    	tlinez.append(layers_cube[item][index_cube][n][2])
    	tlinex.append(layers_cube[item][index_cube][n][3])
    	tliney.append(layers_cube[item][index_cube][n][4])
    	tlinez.append(layers_cube[item][index_cube][n][5])
    	n = n + 1
    #individuating min and max (x,y,z) ranges
    min_x = xline_1[0]
    max_x = xline_1[0]
    for i in xline_1:
        if i < min_x:
            min_x = i
        if i > max_x:
            max_x = i
    dif_x = max_x - min_x
    min_y = yline_1[0]
    max_y = yline_1[0]
    for i in yline_1:
        if i < min_y:
            min_y = i
        if i > max_y:
            max_y = i
    dif_y = max_y - min_y
    min_z = zline_1[0]
    max_z = zline_2[0]
    dif_z = max_z - min_z
    if dif_z < 0:
        t = max_z
        max_z = min_z
        min_z = t

    if dif_x != 0 and dif_y != 0 and dif_z != 0:
        #print('Elaboration of the cube')
        if dif_x > dif_y and dif_x > dif_z:
            #print(' values on x')
            i = 0
            while i < MAX_POINTS:
                xs = random.randint(min_x, max_x)
                ys = random.randint(min_y, max_y)
                zs = random.randint(min_z, max_z)
                #print('i: ',i,' ',xs,ys,zs)
                data_points[index_cube][i][0] = xs
                data_points[index_cube][i][1] = ys
                data_points[index_cube][i][2] = zs
                i = i + 1


        if dif_y > dif_x and dif_y > dif_z:
            #print(' values on y')
            i = 0
            while i < MAX_POINTS:
                xs = random.randint(min_x, max_x)
                ys = random.randint(min_y, max_y)
                zs = random.randint(min_z, max_z)
                #print('i: ',i,' ',xs,ys,zs)
                data_points[index_cube][i][0] = xs
                data_points[index_cube][i][1] = ys
                data_points[index_cube][i][2] = zs
                i = i + 1

        if dif_z > dif_x and dif_z > dif_y:
            #print(' values on z')
            i = 0
            while i < MAX_POINTS:
                xs = random.randint(min_x, max_x)
                ys = random.randint(min_y, max_y)
                zs = random.randint(min_z, max_z)
                #print('i: ',i,' ',xs,ys,zs)
                data_points[index_cube][i][0] = xs
                data_points[index_cube][i][1] = ys
                data_points[index_cube][i][2] = zs
                i = i + 1
        
        #print('maxx ', max_x, "minx ", min_x)
        #print('maxy ', max_y, "minx ", min_y)
        #print('maxz ', max_z, "minx ", min_z)
        #print('difx ', dif_x, "dify ", dif_y, "difz ", dif_z)
        max_cube_size[index_cube][0] = min_x - max_x_ray
        max_cube_size[index_cube][1] = max_x - max_x_ray
        max_cube_size[index_cube][2] = min_y - max_y_ray
        max_cube_size[index_cube][3] = max_y - max_y_ray
        max_cube_size[index_cube][4] = min_z
        max_cube_size[index_cube][5] = max_z

        #Defintion of Ray Elaboration Area
        if i_ray == 0:
            A_MAX_XR = max_x
            A_MAX_YR = max_y
            A_MAX_ZR = max_z
            A_MIN_XR = min_x
            A_MIN_YR = min_y
            A_MIN_ZR = min_z
        else:
            if max_x > A_MAX_XR:
                A_MAX_XR = max_x
            if max_y > A_MAX_YR:
                A_MAX_YR = max_y
            if max_z > A_MAX_ZR:
                A_MAX_ZR = max_z
            if min_x < A_MIN_XR:
                A_MIN_XR = min_x
            if min_y < A_MIN_YR:
                A_MIN_YR = min_y
            if min_z < A_MIN_ZR:
                A_MIN_ZR = min_z
        i_ray = i_ray + 1

        
    n = 0
    index_cube = index_cube + 1
















    


if sys.argv[4] == 'C':
    print('-- Cumulative Mode for TID analysis ... ')
    #Ray Area
    if A_MIN_ZR > 0:
        A_MAX_ZR = A_MAX_ZR + A_MIN_ZR
        A_MIN_ZR = 0
        
    print('-- Ray Elaboration area X(', A_MAX_XR,',', A_MIN_XR,') Y(', A_MAX_YR,',', A_MIN_YR,') Z(', A_MAX_ZR,',', A_MIN_ZR,')')
    #Ray Elaboration on Z axis from top to bottom
    xs = []
    ys = []
    zs = []
    xs.append(A_MAX_XR)
    xs.append(A_MIN_XR)
    ys.append(A_MAX_YR)
    ys.append(A_MIN_YR)
    zs.append(A_MAX_ZR)
    zs.append(A_MIN_ZR)
    #ax_2.plot(xs,ys,zs, "none")

    
    
    
    print('-- Propagation Test ')

    itera = 0
    while itera < int(sys.argv[5]):    
        #Propagation Test
        rad_x_start = random.randint(A_MIN_XR, A_MAX_XR)
        rad_x_end = random.randint(A_MIN_XR, A_MAX_XR)
        rad_y_start = random.randint(A_MIN_YR, A_MAX_YR)
        rad_y_end = random.randint(A_MIN_YR, A_MAX_YR)
        rad_z_start = A_MAX_ZR
        rad_z_end = A_MIN_ZR
        print('-- ',itera,' - Ray Track X(', rad_x_start,',', rad_x_end,') Y(', rad_y_start,',', rad_y_end,') Z(', rad_z_start,',', rad_z_end,')')
        xs = []
        ys = []
        zs = []
        xs.append(rad_x_start)
        xs.append(rad_x_end)
        ys.append(rad_y_start)
        ys.append(rad_y_end)
        zs.append(rad_z_start)
        zs.append(rad_z_end)
        #Visualization of the RAY
        ax_1.plot(xs,ys,zs, "red")
        #ax_2.plot(xs,ys,zs, "red")
        
        
        d_x = (rad_x_start - rad_x_end)/(rad_z_start-rad_z_end)
        d_y = (rad_y_start - rad_y_end)/(rad_z_start-rad_z_end)
        steps = rad_z_start - rad_z_end
        #print('-- Delta XY (', '%.2f' % d_x, ',','%.2f' % d_y,') - Steps ', steps,'')
        n = 0
        px = rad_x_start
        py = rad_y_start
        pz = rad_z_start
        #
        # EXECUTION OF THE RADIATION RAY PROPAGATION 
        #
        while n < steps:
            #print('-- P(x,y,z) = ','%.2f' % px,'%.2f' % py,'%.2f' % pz)
            px = px - d_x
            xs = int(px)
            py = py - d_y
            ys = int(py)
            pz = pz - 1
            zs = int(pz)
            #POINT DEFINED HERE
            #print('-- P(x,y,z) = ', xs, ys, zs)
            #
            index_cube = 0
            while index_cube < total_cube:
                i = 0

                if ((xs >= max_cube_size[index_cube][0] and xs <= max_cube_size[index_cube][1]) and (ys >= max_cube_size[index_cube][2] and ys <= max_cube_size[index_cube][3]) and (zs >= max_cube_size[index_cube][4] and zs <= max_cube_size[index_cube][5])):
                    while i < MAX_POINTS:
                        vxs = data_points[index_cube][i][0]
                        vys = data_points[index_cube][i][1]
                        vzs = data_points[index_cube][i][2]
                        #Maximal Cube Size 0.xmin 1.xmax 2.ymin 3.ymax 4.zmin 5.zmax
                        #print('-- P(x,y,z) = ', xs, ys, zs)
                        #print('-- PC(x,y,z) = ','%.2f' % vxs,'%.2f' % vys,'%.2f' % vzs)
                        dx = abs(vxs - px)
                        dy = abs(vys - py)
                        #dz = abs(vzs - pz)
                        #print('-- D(x,y,z) = ','%.2f' % dx,'%.2f' % dy,'%.2f' % dz)
                        #print('-- D(x) = ','%.2f' % dx)
                        #energy_values[index_cube][i] = 0
                        if (dx < max_x_ray) and (dy < max_y_ray):
                            #print('-- D(x) = ','%.2f' % dx)
                            #print('-- D(y) = ','%.2f' % dy)
                            ind = 0
                            dist_x = max_x_ray
                            dist_y = max_y_ray
                            when = 0
                            when_x = 0
                            when_y = 0
                            for c in ener_scale_pos:
                                if (abs(dx - c)< dist_x):
                                    dist_x = dx - c
                                    when_x = ind                                    
                                ind = ind + 1
                            ind = 0
                            for c in ener_scale_pos:
                                if (abs(dy - c)< dist_y):
                                    dist_y = dy - c
                                    when_y = ind                                    
                                ind = ind + 1
                            if when_x <= when_y:
                            	when = when_x
                            else:
                            	when = when_y
                            print(when_x,when_y,when)
                            energy_values[index_cube][i] = energy_values[index_cube][i] + ener_scale_val[when]
                        i = i + 1
                index_cube = index_cube + 1
            #cs = [1,0,0]
            #m = '^'
            #ax_1.scatter(xs, ys, zs, c=cs, depthshade=True, marker=m)
            n = n + 1
        itera = itera + 1




    print('-- Generating Material Mesh (TBT) ')
    i = 0
    mas = np.amax(energy_values)
    print ('-- Maximal Normalized Energy Value : %.2f' % mas)
    


    
    i = 0
    m = 'o'
    while i < total_cube:
        j = 0
        while j < MAX_POINTS:
            xs = []
            ys = []
            zs = []
            cs = [1,0,0]
            val = 1
            if mas != 0:
                val = 1 - energy_values[i][j] / mas
                #print(val)
            cs = [1, val, val]
            #print (energy_values[i][j])
            if val != 1:
                xs.append(data_points[i][j][0])
                ys.append(data_points[i][j][1])
                zs.append(data_points[i][j][2])
                #ax_2.scatter(xs, ys, zs, c=cs, depthshade=False, marker=m)
            j = j + 1


        #ax_2.scatter(xs, ys, zs, c=cs, depthshade=False, marker=m)
        i = i + 1



    name = 'rad_exposure.pdf'
    #fig_2.savefig(name, format='pdf')
    #plt.show()











    exit()


    
    
    index_cube = 0
    n = 0

    ave_means, ave_std = (20, 35, 30, 35, 27), (2, 3, 4, 1, 2)
    peak_means, peak_std = (25, 38, 34, 40, 45), (3, 5, 2, 3, 3)

    ind = np.arange(len(ave_means))  # the x locations for the groups
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind - width/2, ave_means, width, yerr=ave_std,
                    color='SkyBlue', label='Average')
    rects2 = ax.bar(ind + width/2, peak_means, width, yerr=peak_std,
                    color='IndianRed', label='Peak')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('eV')
    ax.set_title('Energy Loss Distribution')
    ax.set_xticks(ind)
    ax.set_xticklabels(('L1', 'L2', 'L3', 'L4', 'L5'))
    ax.legend()


    def autolabel(rects, xpos='center'):
        """
        Attach a text label above each bar in *rects*, displaying its height.

        *xpos* indicates which side to place the text w.r.t. the center of
        the bar. It can be one of the following {'center', 'right', 'left'}.
        """

        xpos = xpos.lower()  # normalize the case of the parameter
        ha = {'center': 'center', 'right': 'left', 'left': 'right'}
        offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                    '{}'.format(height), ha=ha[xpos], va='bottom')


    autolabel(rects1, "left")
    autolabel(rects2, "right")
 

    plt.show()











    exit()
    
if sys.argv[4] == 'T':
    print('-- Transient Mode analysis ... ')

    
    file_track_log = open("log_radray_track.txt","a+")
    file_track_log.write("XS XE YS YE ZS ZE\n")
    file_log_all = open("log_radray_execution.txt","a+")
    file_cube_data = open("log_radray_data.txt","a+")


    file_v_point = open("log_voltage_points.txt","a+")
    file_data = open("log_max_voltages.txt","a+")
    file_data.write("Particle[#] Volts[V] Duration[ps]\n")
    file_v_point.write("")

    file_h_log = open("log_heatmap.txt","a+")
    

    max_dist = 0
    #Ray Area
    if A_MIN_ZR > 0:
        A_MAX_ZR = A_MAX_ZR + A_MIN_ZR
        A_MIN_ZR = 0
        
    print('-- Ray Elaboration area X(', A_MAX_XR,',', A_MIN_XR,') Y(', A_MAX_YR,',', A_MIN_YR,') Z(', A_MAX_ZR,',', A_MIN_ZR,')')

    #Generation of the XY map
    print('-- Init of XY sensitive map')
    hm_y, hm_x = np.meshgrid(np.linspace(0, int(A_MAX_YR/RES_STEP)+1, 400), np.linspace(0,int(A_MAX_XR/RES_STEP)+1,400))
    hm_z = 0

    #sens_top_XY_map = np.zeros((int(A_MAX_XR/RES_STEP)+1,int(A_MAX_YR/RES_STEP)+1))
    print('-- XY(',int(A_MAX_XR/RES_STEP)+1,',',int(A_MAX_YR/RES_STEP)+1,')')

    
    #Ray Elaboration on Z axis from top to bottom
    xs = []
    ys = []
    zs = []
    xs.append(A_MAX_XR)
    xs.append(A_MIN_XR)
    ys.append(A_MAX_YR)
    ys.append(A_MIN_YR)
    zs.append(A_MAX_ZR)
    zs.append(A_MIN_ZR)
    ax_2.plot(xs,ys,zs, "none")

    
    
    
    print('-- Propagation Test ')

    

    itera = 0
    if sys.argv[7] != 'S':
        max_itera = int(sys.argv[5])
    if sys.argv[7] == 'S':
        max_itera = (int(A_MAX_XR/RES_STEP)+1)*(int(A_MAX_YR/RES_STEP)+1) + 1

    print('-- Maximal Iterations %d' % max_itera)

        
    cube_data_all_en = np.zeros ((max_itera,total_cube))
    cube_data_per_en = np.zeros ((max_itera,total_cube))
    ave_means = np.zeros ((max_itera))
    tot_dist_ener = np.zeros ((max_itera))
    volt_pulse = np.zeros((max_itera))
    monte_energy = np.zeros ((max_itera))
    #Histogram Report Plotting Information
    narrow_amp_SET = np.zeros (10)  #max 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 Volts
    middle_amp_SET = np.zeros (10)  #max 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0 Volts
    large_amp_SET = np.zeros (10)   #max 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0 Volts
    h_gen = 1
    monte_carlo_energy = 0
    monte_carlo_delta = 0
    rad_x_start = random.randint(A_MIN_XR, A_MIN_XR)
    i_rad_x = -1
    i_rad_y = 0
    rad_x_end = random.randint(A_MIN_XR, A_MIN_XR)
    rad_y_start = random.randint(A_MIN_YR, A_MIN_YR)
    rad_y_end = random.randint(A_MIN_YR, A_MIN_YR)

    en_comp = 1
    
    while en_comp == 1:
        #Initialization
        #Maximal Cube Size 0.xmin 1.xmax 2.ymin 3.ymax 4.zmin 5.zmax
        #Maximal Cube Size 0.xmin 1.xmax 2.ymin 3.ymax 4.zmin 5.zmax
        energy_values.fill(0)
        #energy_values = np.zeros ((MAX_CUBE,MAX_POINTS))
        #example_values = np.zeros ((MAX_CUBE,MAX_POINTS))
        hspice_full_cmd = []


        #Propagation Test
        if sys.argv[7] == 'V':
            print('-- Vertical Mode')
            rad_x_start = random.randint(A_MIN_XR, A_MAX_XR)
            rad_x_end = rad_x_start
            rad_y_start = random.randint(A_MIN_YR, A_MAX_YR)
            rad_y_end = rad_y_start
        if sys.argv[7] == 'R':
            print('-- Random Mode')
            rad_x_start = random.randint(A_MIN_XR, A_MAX_XR)
            rad_x_end = random.randint(A_MIN_XR, A_MAX_XR)
            rad_y_start = random.randint(A_MIN_YR, A_MAX_YR)
            rad_y_end = random.randint(A_MIN_YR, A_MAX_YR)
        if sys.argv[7] == 'S':
            print('-- Screen Mode')
            rad_x_start = rad_x_start + RES_STEP
            i_rad_x = i_rad_x + 1
            rad_x_end = rad_x_start
            if (rad_x_start > A_MAX_XR):
                rad_x_start = A_MIN_XR
                rad_x_end = rad_x_start
                i_rad_x = 0
                rad_y_start = rad_y_start + RES_STEP
                rad_y_end = rad_y_start
                i_rad_y = i_rad_y + 1
                if (rad_y_start > (A_MAX_YR-RES_STEP)):
                    en_comp = 0
               
        rad_z_start = A_MAX_ZR
        rad_z_end = A_MIN_ZR
        print('-- ',itera,' - Ray Track X(', rad_x_start,',', rad_x_end,') Y(', rad_y_start,',', rad_y_end,') Z(', rad_z_start,',', rad_z_end,')')
        file_cube_data.write("-- %d START \n" % itera)
        file_log_all.write("-- %d - Ray Track X(%d,%d) Y(%d,%d) Z(%d,%d)\n" % (itera,rad_x_start, rad_x_end, rad_y_start, rad_y_end, rad_z_start, rad_z_end))
        file_track_log.write("%d %d %d %d %d %d\n" % (rad_x_start, rad_x_end, rad_y_start, rad_y_end, rad_z_start, rad_z_end))     
        file_v_point.write("\n%d "% itera)

        xs = []
        ys = []
        zs = []
        xs.append(rad_x_start)
        xs.append(rad_x_end)
        ys.append(rad_y_start)
        ys.append(rad_y_end)
        zs.append(rad_z_start)
        zs.append(rad_z_end)
        #Visualization of the RAY
        ax_1.plot(xs,ys,zs, "red")
        ax_2.plot(xs,ys,zs, "red")
        
        
        d_x = (rad_x_start - rad_x_end)/(rad_z_start-rad_z_end)
        d_y = (rad_y_start - rad_y_end)/(rad_z_start-rad_z_end)
        steps = rad_z_start - rad_z_end
        #print('-- Delta XY (', '%.2f' % d_x, ',','%.2f' % d_y,') - Steps ', steps,'')
        n = 0
        px = rad_x_start
        py = rad_y_start
        pz = rad_z_start
        #
        # EXECUTION OF THE RADIATION RAY PROPAGATION 
        #
        while n < steps:
            #print('-- P(x,y,z) = ','%.2f' % px,'%.2f' % py,'%.2f' % pz)
            px = px - d_x
            xs = int(px)
            py = py - d_y
            ys = int(py)
            pz = pz - 1
            zs = int(pz)
            #POINT DEFINED HERE
            #print('-- P(x,y,z) = ', xs, ys, zs)
            #
            index_cube = 0
            while index_cube < total_cube:
                i = 0

                if ((xs >= max_cube_size[index_cube][0] and xs <= max_cube_size[index_cube][1]) and (ys >= max_cube_size[index_cube][2] and ys <= max_cube_size[index_cube][3]) and (zs >= max_cube_size[index_cube][4] and zs <= max_cube_size[index_cube][5])):
                    while i < MAX_POINTS:
                        vxs = data_points[index_cube][i][0]
                        vys = data_points[index_cube][i][1]
                        vzs = data_points[index_cube][i][2]
                        #Maximal Cube Size 0.xmin 1.xmax 2.ymin 3.ymax 4.zmin 5.zmax
                        #print('-- P(x,y,z) = ', xs, ys, zs)
                        #print('-- PC(x,y,z) = ','%.2f' % vxs,'%.2f' % vys,'%.2f' % vzs)
                        dx = abs(vxs - px)
                        dy = abs(vys - py)
                        eu_dis = math.sqrt((vxs-px)**2+(vys-py)**2)
                        if (eu_dis < max_x_ray):
                            #print('-- D(x) = ','%.2f' % dx)
                            #print('-- D(y) = ','%.2f' % dy)
                            ind = 0
                            dist_x = max_x_ray
                            when_x = 0
                            for c in ener_scale_pos:
                                if (abs(eu_dis - c)< dist_x):
                                    dist_x = dx - c
                                    when_x = ind                                    
                                ind = ind + 1
                            ind = 0
                            #print(eu_dis, when_x, ener_scale_val[when_x])
                            #print(ener_scale_val)
                            energy_values[index_cube][i] = energy_values[index_cube][i] + ener_scale_val[when_x]
                            #sens_top_XY_map[i_rad_x][i_rad_y] = sens_top_XY_map[i_rad_x][i_rad_y] + ener_scale_val[when]
                            #Test
                        i = i + 1

                        #dz = abs(vzs - pz)
                        #print('-- D(x,y,z) = ','%.2f' % dx,'%.2f' % dy,'%.2f' % dz)
                        #print('-- D(x) = ','%.2f' % dx)
                        #energy_values[index_cube][i] = 0
                        #if (dx < max_x_ray) and (dy < max_x_ray):
                            #print('-- D(x) = ','%.2f' % dx)
                            #print('-- D(y) = ','%.2f' % dy)
                            #ind = 0
                            #dist = max_x_ray
                            #when = 0
                            #for c in ener_scale_pos:
                            #    if (abs(dx - c)< dist):
                            #        dist = dx - c
                            #        when = ind                                    
                            #    ind = ind + 1
                        #if (dx < max_x_ray) and (dy < max_y_ray):
                            #print('-- D(x) = ','%.2f' % dx)
                            #print('-- D(y) = ','%.2f' % dy)
                            #ind = 0
                            #dist_x = max_x_ray
                            #dist_y = max_y_ray
                            #when = 0
                            #when_x = 0
                            #when_y = 0
                            #for c in ener_scale_pos:
                            #    if (abs(dx - c)< dist_x):
                            #        dist_x = dx - c
                            #        when_x = ind                                    
                            #    ind = ind + 1
                            #ind = 0
                            #for c in ener_scale_pos:
                            #    if (abs(dy - c)< dist_y):
                            #        dist_y = dy - c
                            #        when_y = ind                                    
                            #    ind = ind + 1
                            #if when_x >= when_y:
                            #	when = when_x
                            #else:
                            #	when = when_y
                            #print(when_x,when_y,when)
                            #energy_values[index_cube][i] = energy_values[index_cube][i] + ener_scale_val[when]
                            #sens_top_XY_map[i_rad_x][i_rad_y] = sens_top_XY_map[i_rad_x][i_rad_y] + ener_scale_val[when]
                            #Test
                        #i = i + 1
                index_cube = index_cube + 1
            #cs = [1,0,0]
            #m = '^'
            #ax_1.scatter(xs, ys, zs, c=cs, depthshade=True, marker=m)
            n = n + 1

            
        #print('-- Generating Material Dynamic Mesh')
        i = 0

        
        mas = np.amax(energy_values)
        print ('-- Maximal Energy Value : %.2f ev/Ang' % mas)
        file_log_all.write("-- Maximal Energy Value : %.2f ev/Ang \n" % mas)
        i = 0
        m = 'o'
        #cube_data_all_en = np.zeros (total_cube)
        total_dis_en = 0
        total_v_cub = 0
        pulse_xs = []
        pulse_ys = []
        time_dist = 0
        max_v_pulse = 0
        pulse_xs.append(time_dist)
        pulse_ys.append(max_v_pulse)
        time_dist = 100
        max_v_pulse = 0
        dur_max_v_pulse = 0
        pulse_xs.append(time_dist)
        pulse_ys.append(max_v_pulse)




        while i < total_cube:
            j = 0
            cube_ener = 0
            va_cu = 0 #Number of "energetic" point
            while j < MAX_POINTS:
                xs = []
                ys = []
                zs = []
                cs = [1,0,0]
                val = 1
                if mas != 0:
                    val = 1 - energy_values[i][j] / mas
                    cube_ener = cube_ener + energy_values[i][j]
                    #print(val)
                cs = [1, val, val]
                if val != 1:
                	va_cu = va_cu + 1
                #print (energy_values[i][j])
                if val != 1:
                    xs.append(data_points[i][j][0])
                    ys.append(data_points[i][j][1])
                    zs.append(data_points[i][j][2])
                    #
                    # Heatmap coefficient storage
                    #
                    file_h_log.write("%d %d %d %f\n" % (data_points[i][j][0], data_points[i][j][1], data_points[i][j][2], val))
                                     

                    
                    #ax_2.scatter(xs, ys, zs, c=cs, depthshade=False, marker=m)
                    ax_2.scatter(xs, ys, zs, c=rgb2hex(cs), depthshade=False, marker=m)
                j = j + 1
            if cube_ener > 0:
                en = layer_energ_id[int(cube_layer_id[i])]
                if (en == 0):
                    print('-- Error Energy on Layer Undefined')
                    en = 1
                normal_en = (((va_cu/MAX_POINTS)*(cube_ener/10))* en) / 1000
                total_dis_en = total_dis_en + normal_en
                #print('-- Cube ',i,'Layer ', int(cube_layer_id[i]), 'DB Energy', int(en),' KeV/ang Energy %.2f eV/PartialCube' % cube_ener, ' Distributed Energy Percentage %.2f ' % ((va_cu/MAX_POINTS)*100), ' Distributed Energy on Cube %.2f KeV/Cube ' % normal_en)
                file_log_all.write("-- Cube %d Layer %d DB Energy %d KeV/ang Energy %.2f eV/PartialCube Distributed Energy Percentage %.2f Distributed Energy on Cube %.2f KeV/Cube\n" % (i,int(cube_layer_id[i]), int(en), cube_ener, ((va_cu/MAX_POINTS)*100), normal_en))
                file_cube_data.write("%d. %d %d %d %.2f %.2f %.2f" % (itera,i,int(cube_layer_id[i]),int(en),cube_ener, ((va_cu/MAX_POINTS)*100), normal_en))

                vol = (max_cube_size[i][1]-max_cube_size[i][0])*(max_cube_size[i][3]-max_cube_size[i][2])*(max_cube_size[i][5]-max_cube_size[i][4])
                el = 1.602188E-19
                # Related to the normalized value - Scale Factor 0.001 um - UNITS 0.001 1e-9
                vol = vol / 10E10
                # vol = vol / 10E11
                sil_v = 210E-12
                n_el_sil = 14
                cube_col = (vol/sil_v)*el*n_el_sil
                max_v_pulse = normal_en*1000*el/cube_col
                #print('-- Cube Volume %.2f' % vol, 'Ang ', el*14)
                #print('-- Cube Charge: ', cube_col, ' Coulomb - Maximal Pulse %.5f' % max_v_pulse, ' Volts - Duration %.2f' % (max_cube_size[i][5]-max_cube_size[i][4]),' ps ' )
                file_log_all.write("-- Cube Maximal Pulse %.5f Volts - Duration %.2f ps\n" % (max_v_pulse, (max_cube_size[i][5]-max_cube_size[i][4])))
                #Hspice Full Command Generation
                #Reference layer: int(cube_layer_id[i])
                #Duration (max_cube_size[i][5]-max_cube_size[i][4])
                #print(hspice_ref)
                #print(hspice_sig)
                hspice_full_cmd.append([int(cube_layer_id[i]),max_v_pulse,(max_cube_size[i][5]-max_cube_size[i][4])])

                file_cube_data.write(" %.5f %.2f \n" % (max_v_pulse, (max_cube_size[i][5]-max_cube_size[i][4])))
                time_dist = time_dist + (max_cube_size[i][5]-max_cube_size[i][4])
                file_v_point.write("%.5f %.2f " % (max_v_pulse, time_dist))
                pulse_xs.append(time_dist)
                pulse_ys.append(max_v_pulse)
                if ((cube_ener/MAX_POINTS) > max_dist):
                    max_dist = cube_ener/MAX_POINTS
                if (max_v_pulse > total_v_cub):
                    total_v_cub = max_v_pulse
                    dur_max_v_pulse = (max_cube_size[i][5]-max_cube_size[i][4])

                cube_data_all_en[itera][i] = cube_ener
                    
                if (va_cu/MAX_POINTS)*100 > ave_means[itera]:
                    ave_means[itera] = int((va_cu/MAX_POINTS)*100)
                cube_data_per_en[itera][i] = (va_cu/MAX_POINTS)*100


                #    dur_max_v_pulse = (max_cube_size[i][5]-max_cube_size[i][4])
            #ax_2.scatter(xs, ys, zs, c=cs, depthshade=False, marker=m)
            i = i + 1

        #print (cube_data_all_en)
        print ('-- Total Distributed Energy on Cubes : %.2f KeV/Cubes' % total_dis_en, ' Maximal Peak %.5f' % total_v_cub,' [V]')
        file_cube_data.write("-- %d END -- %.2f %.5f \n" % (itera, total_dis_en, total_v_cub))
        file_log_all.write("-- Total Distributed Energy on Cubes : %.2f KeV/Cubes Maximal Peak %.5f [V] \n" % (total_dis_en, total_v_cub))
        #update sensitive XY map

        if sys.argv[7] == 'S':
            print('-- Sensitive Map XY(',i_rad_x,',',i_rad_y,')')

        hm_z_new = total_v_cub * np.exp(-((hm_x-i_rad_x)*0.1)** 2 - ((hm_y-i_rad_y)*0.1)**2)
        hm_zp = np.maximum(hm_z_new, hm_z)
        hm_z = hm_zp


        
        #sens_top_XY_map[i_rad_x][i_rad_y] = total_v_cub
        
        
        #HSpice Effect Generation
        #print(hspice_port)
        #print(hspice_ref)
        #print(hspice_sig)
        t_ref = 50
        if sys.argv[8] == 'Y':
            if total_v_cub >= 0.1:
                print('-- Generating HSpice Command')
                #print(hspice_port)
                #print(hspice_ref)
                #print(hspice_sig)
                h_spice_log.write("*Testing Pattern from RadRay Tool\n")
                h_spice_log.write(".alter case %d: ref_radray_%d \n" % (h_gen, itera))
                h_gen = h_gen + 1
                for i in hspice_port:
                    #Generation of the command for the port i
                    #print('%s GND PWL (0n,0m 50n,0m ' % i)
                    h_spice_log.write("%s GND PWL (0n,0m 50n,0m " % i)
                    for j in range(len(hspice_full_cmd)):
                        #print (hspice_full_cmd[j][0])
                        #print (hspice_port[(hspice_ref.index(hspice_full_cmd[j][0]))])
                        ref = hspice_ref.index(hspice_full_cmd[j][0])
                        if hspice_sig[ref]==i:
                            #inside the port
                            #print('inside %s' % (hspice_sig[ref]))                        
                            t_ref = t_ref + (hspice_full_cmd[j][2]/1000) + 0.1
                            #print ('%.1f' % (hspice_full_cmd[j][2]/1000))
                            #print ('%.1fn,%.0fm ' % (t_ref,hspice_full_cmd[j][1]*1000))
                            h_spice_log.write("%.1fn,%.0fm " % (t_ref,hspice_full_cmd[j][1]*1000))
                        #print (hspice_sig[ref])
                    t_ref = t_ref + 0.1
                    h_spice_log.write("%.1fn,0m" % (t_ref))
                    h_spice_log.write(") \n")
                    #print('\n')
                    #for j
                    #print(hspice_ref.index[int(hspice_full_cmd[i][0])])
                print('-- Hspice test command done ')
                
                
                





        
        #print(hspice_full_cmd)
        #print(len(hspice_full_cmd))
                           
        #exit()
        



        #Maximal Voltage Peak Classification by Amplitude
        #Narrow Amplitude
        if (total_v_cub >= 0 and total_v_cub < 0.1):
            narrow_amp_SET[0] = narrow_amp_SET[0] + 1
        if (total_v_cub >= 0.1 and total_v_cub < 0.2):
            narrow_amp_SET[1] = narrow_amp_SET[1] + 1
        if (total_v_cub >= 0.2 and total_v_cub < 0.3):
            narrow_amp_SET[2] = narrow_amp_SET[2] + 1
        if (total_v_cub >= 0.3 and total_v_cub < 0.4):
            narrow_amp_SET[3] = narrow_amp_SET[3] + 1
        if (total_v_cub >= 0.4 and total_v_cub < 0.5):
            narrow_amp_SET[4] = narrow_amp_SET[4] + 1
        if (total_v_cub >= 0.5 and total_v_cub < 0.6):
            narrow_amp_SET[5] = narrow_amp_SET[5] + 1
        if (total_v_cub >= 0.6 and total_v_cub < 0.7):
            narrow_amp_SET[6] = narrow_amp_SET[6] + 1
        if (total_v_cub >= 0.7 and total_v_cub < 0.8):
            narrow_amp_SET[7] = narrow_amp_SET[7] + 1
        if (total_v_cub >= 0.8 and total_v_cub < 0.9):
            narrow_amp_SET[8] = narrow_amp_SET[8] + 1
        if (total_v_cub >= 0.9 and total_v_cub < 1.0):
            narrow_amp_SET[9] = narrow_amp_SET[9] + 1

        #Middle Amplitude
        if (total_v_cub >= 1.0 and total_v_cub < 1.4):
            middle_amp_SET[0] = middle_amp_SET[0] + 1
        if (total_v_cub >= 1.4 and total_v_cub < 1.8):
            middle_amp_SET[1] = middle_amp_SET[1] + 1
        if (total_v_cub >= 1.8 and total_v_cub < 2.2):
            middle_amp_SET[2] = middle_amp_SET[2] + 1
        if (total_v_cub >= 2.2 and total_v_cub < 2.6):
            middle_amp_SET[3] = middle_amp_SET[3] + 1
        if (total_v_cub >= 2.6 and total_v_cub < 3.0):
            middle_amp_SET[4] = middle_amp_SET[4] + 1
        if (total_v_cub >= 3.0 and total_v_cub < 3.4):
            middle_amp_SET[5] = middle_amp_SET[5] + 1
        if (total_v_cub >= 3.4 and total_v_cub < 3.8):
            middle_amp_SET[6] = middle_amp_SET[6] + 1
        if (total_v_cub >= 3.8 and total_v_cub < 4.2):
            middle_amp_SET[7] = middle_amp_SET[7] + 1
        if (total_v_cub >= 4.2 and total_v_cub < 4.6):
            middle_amp_SET[8] = middle_amp_SET[8] + 1
        if (total_v_cub >= 4.6 and total_v_cub < 5.0):
            middle_amp_SET[9] = middle_amp_SET[9] + 1


        #Large Amplitude
        if (total_v_cub >= 5.0 and total_v_cub < 6.0):
            large_amp_SET[0] = large_amp_SET[0] + 1
        if (total_v_cub >= 6.0 and total_v_cub < 8.0):
            large_amp_SET[1] = large_amp_SET[1] + 1
        if (total_v_cub >= 8.0 and total_v_cub < 10.0):
            large_amp_SET[2] = large_amp_SET[2] + 1
        if (total_v_cub >= 10.0 and total_v_cub < 12.0):
            large_amp_SET[3] = large_amp_SET[3] + 1
        if (total_v_cub >= 12.0 and total_v_cub < 14.0):
            large_amp_SET[4] = large_amp_SET[4] + 1
        if (total_v_cub >= 14.0 and total_v_cub < 16.0):
            large_amp_SET[5] = large_amp_SET[5] + 1
        if (total_v_cub >= 16.0 and total_v_cub < 18.0):
            large_amp_SET[6] = large_amp_SET[6] + 1
        if (total_v_cub >= 18.0 and total_v_cub < 20.0):
            large_amp_SET[7] = large_amp_SET[7] + 1
        if (total_v_cub >= 20.0 and total_v_cub < 22.0):
            large_amp_SET[8] = large_amp_SET[8] + 1
        if (total_v_cub >= 22.0 and total_v_cub < 24.0):
            large_amp_SET[9] = large_amp_SET[9] + 1


            







     # narrow_amp_SET = np.zeros (10)  #max 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 Volts
     #middle_amp_SET = np.zeros (10)  #max 1.4, 1.8, 2.2, 2.6, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0 Volts
      #large_amp_SET = np.zeros (10)   #max 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0 Volts








        
        monte_carlo_delta = abs((monte_carlo_energy/(itera+1) - ((monte_carlo_energy + total_dis_en)/(itera+1))))
        monte_carlo_energy = monte_carlo_energy + total_dis_en
        file_data.write("%d %f %f\n" % (itera,total_v_cub, dur_max_v_pulse))

        print ('-- Monte Carlo Delta %.2f' % monte_carlo_delta, ' Monte Carlo Energy Value %.2f' % (monte_carlo_energy/(itera+1)))
        file_log_all.write("-- Monte Carlo Delta %.2f Monte Carlo Energy Value %.2f \n" % (monte_carlo_delta,(monte_carlo_energy/(itera+1))))
        monte_energy[itera] = monte_carlo_energy/(itera+1)
        tot_dist_ener[itera] = int(total_dis_en)
        volt_pulse[itera] = total_v_cub

        index_cube = 0
        for item in layer_list:
        	n = 0
        	xline_1 = []
        	yline_1 = []
        	zline_1 = []
        	xline_2 = []
        	yline_2 = []
        	zline_2 = []
        	tlinex = []
        	while n < vertex_cube_list[index_cube]:
        		tlinex = []
        		tliney = []
        		tlinez = []
        		xline_1.append(layers_cube[item][index_cube][n][0])
        		yline_1.append(layers_cube[item][index_cube][n][1])
        		zline_1.append(layers_cube[item][index_cube][n][2])
        		xline_2.append(layers_cube[item][index_cube][n][3])
        		yline_2.append(layers_cube[item][index_cube][n][4])
        		zline_2.append(layers_cube[item][index_cube][n][5])
        		tlinex.append(layers_cube[item][index_cube][n][0])
        		tliney.append(layers_cube[item][index_cube][n][1])
        		tlinez.append(layers_cube[item][index_cube][n][2])
        		tlinex.append(layers_cube[item][index_cube][n][3])
        		tliney.append(layers_cube[item][index_cube][n][4])
        		tlinez.append(layers_cube[item][index_cube][n][5])
        		if (cube_data_all_en[itera][index_cube] > 0):
        			ax_2.plot3D(tlinex, tliney, tlinez, 'gray')  #'gray'
        		n = n + 1
        	if (cube_data_all_en[itera][index_cube] > 0):
        		ax_2.plot3D(xline_1, yline_1, zline_1, 'blue')
        		ax_2.plot3D(xline_2, yline_2, zline_2, 'blue')  
        	n = 0
        	index_cube = index_cube + 1






        #print(pulse_xs)
        #print(pulse_ys)
        time_dist = time_dist + 200
        max_v_pulse = 0
        pulse_xs.append(time_dist)
        pulse_ys.append(max_v_pulse)
        time_dist = time_dist + 200
        max_v_pulse = 0
        pulse_xs.append(time_dist)
        pulse_ys.append(max_v_pulse)

        fig_3, ax_3 = plt.subplots()
        
        ax_3.set_ylabel('[V]')
        ax_3.set_title('Voltage Pulse Plot')
        ax_3.set_xlabel('[ps]')
        ax_3.plot(pulse_xs,pulse_ys,"b-")
        name = 'pulse_plot_%d.pdf' % itera
        if sys.argv[6] == 'Y':
            fig_3.savefig(name, format='pdf')
        #plt.show(fig_2)
        fig_3.clear()
        plt.close(fig_3)


        name = 'rad_exposure_%d.pdf' % itera
        if sys.argv[6] == 'Y':
            fig_2.savefig(name, format='pdf')
        #plt.show(fig_2)
        fig_2.clear()
        fig_2 = plt.figure('Static Signal 3D view')
        ax_2 = fig_2.gca(projection='3d')
        #ax_2.set_aspect("equal") #3D plot : static signal view
        xs = []
        ys = []
        zs = []
        xs.append(A_MAX_XR)
        xs.append(A_MIN_XR)
        ys.append(A_MAX_YR)
        ys.append(A_MIN_YR)
        zs.append(A_MAX_ZR)
        zs.append(A_MIN_ZR)
        ax_2.plot(xs,ys,zs, "none")
            
        
        itera = itera + 1
        
        if itera == int(sys.argv[5]) and sys.argv[7] != 'S':
            en_comp = 0


    #print(cube_data_all_en)

    fig, ax = plt.subplots()
    im = ax.imshow(cube_data_all_en, cmap=plt.get_cmap('hot'), interpolation='nearest',vmin=0, vmax=max_dist)
    fig.colorbar(im)
    fig.savefig('heatmap_all_en.pdf', format='pdf')
    #plt.show(im)
    #plt.imshow(cube_data_all_en)
    #plt.show()

    #ave_means, ave_std = (20, 35, 30, 35, 27), (2, 3, 4, 1, 2)

    ind = np.arange(len(ave_means))  # the x locations for the groups
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind - width/2, ave_means, width, color='Blue', label='Max')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('%')
    ax.set_title('Maximal Distributed Energy Percentage / Particle')
    ax.set_xticks(ind)
    ax.set_xlabel('Particle ID')
    ax.legend()


    def autolabel(rects, xpos='center'):
        """
        Attach a text label above each bar in *rects*, displaying its height.

        *xpos* indicates which side to place the text w.r.t. the center of
        the bar. It can be one of the following {'center', 'right', 'left'}.
        """

        xpos = xpos.lower()  # normalize the case of the parameter
        ha = {'center': 'center', 'right': 'left', 'left': 'right'}
        offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

        for rect in rects:
            height = rect.get_height()
            #ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
            #        '{}'.format(height), ha=ha[xpos], va='bottom')


    autolabel(rects1, "left")

    fig.savefig('energy_percentage.pdf', format='pdf')



    ind = np.arange(len(tot_dist_ener))  # the x locations for the groups
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind - width/2, tot_dist_ener, width, color='Red', label='Max')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('KeV/Cubes')
    ax.set_title('Distributed Energy on Cubes')
    ax.set_xticks(ind)
    ax.set_xlabel('Particle ID')
    ax.legend()


    def autolabel(rects, xpos='center'):
        """
        Attach a text label above each bar in *rects*, displaying its height.

        *xpos* indicates which side to place the text w.r.t. the center of
        the bar. It can be one of the following {'center', 'right', 'left'}.
        """

        xpos = xpos.lower()  # normalize the case of the parameter
        ha = {'center': 'center', 'right': 'left', 'left': 'right'}
        offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

        for rect in rects:
            height = rect.get_height()
            #ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
            #        '{}'.format(height), ha=ha[xpos], va='bottom')


    autolabel(rects1, "left")

    fig.savefig('distributed_energy_profile.pdf', format='pdf')



    ind = np.arange(len(monte_energy))  # the x locations for the groups
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind - width/2, monte_energy, width, color='Red', label='Max')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('KeV/Cubes')
    ax.set_title('MonteCarlo Distributed Energy')
    ax.set_xticks(ind)
    ax.set_xlabel('Particle ID')
    ax.legend()


    def autolabel(rects, xpos='center'):
        """
        Attach a text label above each bar in *rects*, displaying its height.

        *xpos* indicates which side to place the text w.r.t. the center of
        the bar. It can be one of the following {'center', 'right', 'left'}.
        """

        xpos = xpos.lower()  # normalize the case of the parameter
        ha = {'center': 'center', 'right': 'left', 'left': 'right'}
        offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

        for rect in rects:
            height = rect.get_height()
            #ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                   # '{}'.format(height), ha=ha[xpos], va='bottom')


    autolabel(rects1, "left")

    fig.savefig('monte_carlo_profile.pdf', format='pdf')



    ind = np.arange(len(volt_pulse))  # the x locations for the groups
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind - width/2, volt_pulse, width, color='Red', label='Max')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('[V]')
    ax.set_title('Maximal Voltage Peak')
    ax.set_xticks(ind)
    ax.set_xlabel('Particle ID')
    ax.legend()


    def autolabel(rects, xpos='center'):
        """
        Attach a text label above each bar in *rects*, displaying its height.

        *xpos* indicates which side to place the text w.r.t. the center of
        the bar. It can be one of the following {'center', 'right', 'left'}.
        """

        xpos = xpos.lower()  # normalize the case of the parameter
        ha = {'center': 'center', 'right': 'left', 'left': 'right'}
        offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

        for rect in rects:
            height = rect.get_height()
            #ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                   # '{}'.format(height), ha=ha[xpos], va='bottom')


    autolabel(rects1, "left")

    fig.savefig('voltage_pulse_profile.pdf', format='pdf')


    fig, ax = plt.subplots()
    print('-- Narrow Amp Pulses ', narrow_amp_SET)
    file_log_all.write("-- Narrow Amp Pulses \n")
    np.savetxt(file_log_all, narrow_amp_SET, fmt="%d")
    items = ('0 < 0.1','< 0.2','< 0.3','< 0.4','< 0.5','< 0.6','< 0.7','< 0.8','< 0.9', '< 1.0')
    # An "interface" to matplotlib.axes.Axes.hist() method
    ax.bar(items,narrow_amp_SET,align='center') # A bar chart
    ax.set_title('Narrow Amplitude Pulse Distribution')
    ax.set_xlabel('Max Voltage [V]')
    ax.set_ylabel('Event [#]')
    
    fig.savefig('narrow_amp_SET.pdf', format='pdf')


    fig, ax = plt.subplots()
    print('-- Middle Amp Pulses ', middle_amp_SET)
    file_log_all.write("-- Middle Amp Pulses \n")
    np.savetxt(file_log_all, middle_amp_SET, fmt="%d")
    items = ('1 < 1.4','< 1.8','< 2.2','< 2.6','< 3.0','< 3.4','< 3.8','< 4.2','< 4.6', '< 5.0')
    # An "interface" to matplotlib.axes.Axes.hist() method
    ax.bar(items,middle_amp_SET,align='center') # A bar chart
    ax.set_title('Middle Amplitude Pulse Distribution')
    ax.set_xlabel('Max Voltage [V]')
    ax.set_ylabel('Event [#]')
    
    fig.savefig('middle_amp_SET.pdf', format='pdf')


    fig, ax = plt.subplots()
    print('-- Large Amp Pulses ', large_amp_SET)
    file_log_all.write("-- Large Amp Pulses \n")
    np.savetxt(file_log_all, large_amp_SET, fmt="%d")
    items = ('5 < 6','< 8','< 10','< 12','< 14','< 16','< 18','< 20','< 22', '< 24')
    # An "interface" to matplotlib.axes.Axes.hist() method
    ax.bar(items,large_amp_SET,align='center') # A bar chart
    ax.set_title('Large Amplitude Pulse Distribution')
    ax.set_xlabel('Max Voltage [V]')
    ax.set_ylabel('Event [#]')

    fig.savefig('large_amp_SET.pdf', format='pdf')








 



    file_h_log.close()
    file_track_log.close()

    file_log_all.close()
    file_cube_data.close()

    file_data.close()
    file_v_point.close()
    print('-- Transient Analysis Report done for ', max_itera,' particles ')


    #
    #    Heatmap analysis
    #   sens_top_XY_map[rad_x_start][rad_y_start] = total_v_cub
    max_dist = 1    
    if sys.argv[7] == 'S':
            hm_z = hm_z[:-1, :-1]
            z_max = np.abs(hm_z).max()
            fig, asx = plt.subplots()
            #figu, asx = plt.subplots()
            c = asx.pcolormesh(hm_x, hm_y, hm_z, cmap='magma_r', vmin = 0, vmax=z_max)
            asx.set_title('Cell Sensitivity Heatmap')
            asx.axis([hm_x.min(), hm_x.max(), hm_y.min(), hm_y.max()])
            fig.colorbar(c, ax=asx)
            #plt.show()

            


            
            #im = ax.imshow(sens_top_XY_map, cmap=plt.get_cmap('hot'), interpolation='nearest',vmin=0, vmax=max_dist)
            #fig.colorbar(im)
            fig.savefig('heatmap_XY_vertical.pdf', format='pdf')

        














    




exit()




    





    

    















#n = 1000
#for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
#    xs = randrange(n, 0, 100)
#    ys = randrange(n, 0, 100)
#    zs = randrange(n, zlow, zhigh)
#    ax_t.scatter(xs, ys, zs, c=c, depthshade=True, marker=m)

cs = [[1,1,1],[1,0.9,0.9],[1,0.8,0.8],[1,0.7,0.7],[1,0.6,0.6],[1,0.5,0.5],[1,0.4,0.4],[1,0.3,0.3],[1,0.2,0.2],[1,0.1,0.1],[1,0,0]]
m = 'o'
#cs = [0.175059,0.175059,0.175059,0.175059,0.175059,0.175059,0.175059,0.175059,0.175059,0.175059,0.175059]
#cs = [4.375059,4.375059,4.375059,4.375059,4.375059,4.375059,4.375059,4.375059,0.375059,0.375059,0.375059]
print('-- Generated colors : ', cs)

xs = [1,2,3,4,5,6,7,8,9,10,0]
ys = [2,2,2,2,2,2,2,2,2,2,0]
zs = [5,5,5,5,5,5,5,5,5,5,0]
ax_t.scatter(xs, ys, zs, c=cs, depthshade=False, marker=m)


print('-- Importing Energy Released (TBC)')




print('-- Importing Metals Values (TBC)')


print('-- Importing Cell Dynamic Condition (TBC)')

print('-- Generating Cell Static Condition (TBC)')

print('-- Generating Cell Dynamic Condition (TBC)')



plt.show()


# Data for a three-dimensional line
#zline = np.linspace(0, 50, 1000)
#xline = np.sin(zline)
#yline = np.cos(zline)
#print(' xline :', xline)
#print(' yline :', yline)
#print(' zline :', zline)
