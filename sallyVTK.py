import matplotlib.pyplot as plt
import numpy as np

#from evtk.hl import gridToVTK

# identify variables
resolution = 128
res_arr=np.linspace(0,resolution)
xs=resolution; ys=resolution; zs=resolution
n_times=1
final_time=1
initial_time = 1; SMALL = 0.00000000001

# define xyz values
xdelta = 1.0 / (xs-1.0)
ydelta = 1.0 / (ys-1.0)
zdelta = 1.0 / (zs-1.0)

# SALLY
SALLY = np.zeros((n_times,resolution,resolution,resolution,3))
times=np.linspace(initial_time,final_time,n_times)
xarr=np.arange(0,resolution,1)
yarr=np.arange(0,resolution,1)
zarr=np.arange(0,resolution,1)
time_counter=0
# define the ix, iy, iz
ix=0; iy=0; iz=0

for time in (times):
    # for iz in res_arr:
    for iz in (zarr):
        print(iz)
    # while iz < zs:
        z = iz * zdelta
        xc = 0.5 + 0.1*np.sin(0.04*time+10.0*z)   # For each z-slice, determine the spiral circle.
        yc = 0.5 + 0.1*np.cos(0.03*time+3.0*z)    #   (xc,yc) determine the center of the circle.
        r = 0.1 + 0.4 * z*z + 0.1 * z * np.sin(8.0*z) #  The radius also changes at each z-slice.
        r2 = 0.2 + 0.1*z #   r is the center radius, r2 is for damping
        #lets move to iy
        # for iy in res_arr:
        for iy in (yarr):
        # while iy < ys:
            y = iy * ydelta
            #lets move to ix
            # for ix in res_arr:
            for ix in (xarr):
            # while ix < xs:
                x = ix * xdelta
                temp = np.sqrt((y-yc)*(y-yc) + (x-xc)*(x-xc))
                #lets go to scale
                scale = np.abs(r - temp)
                if scale > r2:
                    scale = 0.8 - scale
                else:
                    scale = 1.0
                z0 = 0.1 * (0.1 - temp*z)
                if z0 < 0.0:
                    z0 = 0.0
                temp = np.sqrt( temp*temp + z0*z0 )
                scale = (r + r2 - temp) * scale / (temp + SMALL)
                scale = scale / (1+z)
                # create sally vector fields 
                sallyX = scale * (y-yc) + 0.1*(x-xc)
                sallyY = scale * -(x-xc) + 0.1*(y-yc)
                sallyZ = scale * z0
                # append into final data array
                SALLY[time_counter][iz] [iy][ix][0]=sallyX
                SALLY[time_counter][iz] [iy][ix][1]=sallyY
                SALLY[time_counter][iz] [iy][ix][2]=sallyZ
                # print('ix')
                # print(ix)
                ix +=1
            # print('iy')
            # print(iy)
            iy += 1
        # print('iz')
        # print(iz)
        iz += 1
    time_counter=time_counter+1
            
    
#print (SALLY)

print(np.shape(SALLY))

#only makes VTK file form first time iteration at the moment
#gridToVTK("./structured", x, y, z, cellData = {"pressure" : pressure}, pointData = {"temp" : temp})
#gridToVTK("./structured",)
#gridToVTK("./image",xdata={"xdata":SALLY[][][][]},,)
vtkreplacement=open("sally128.txt","w")
#write the data tofile
#work out no of lines of data
#loop over this and write each line into file
resolutioniterate=np.arange(0,resolution,1)
#label columns
vtkreplacement.write("t\t x\t y\t z\t xmag\t ymag\t zmag\n")
time_counter=0
for time in times:
    for zdata in resolutioniterate:
        for ydata in resolutioniterate:
            for xdata in resolutioniterate:
                #vtkreplacement.write(xdata, ydata, zdata, SALLY[time_counter][zdata][ydata][xdata][0], SALLY[time_counter][zdata][ydata][xdata][1], SALLY[time_counter][zdata][ydata][xdata][2])
                print("line:",((time_counter*resolution*resolution*resolution)+(zdata*resolution*resolution)+(ydata*resolution)+xdata))#for 1 time only
                print("{0[0]}\t {0[1]}\t {0[2]}\t {0[3]}\t {0[4]}\t {0[5]}\t {0[6]}\n".format((time,xdata*xdelta, ydata*ydelta, zdata*zdelta, SALLY[time_counter][zdata][ydata][xdata][0], SALLY[time_counter][zdata][ydata][xdata][1], SALLY[time_counter][zdata][ydata][xdata][2])))
                #placeholder="{}\t {}\t {}\t {}\t {}\t {}\n",format((xdata, ydata, zdata, SALLY[time_counter][zdata][ydata][xdata][0], SALLY[time_counter][zdata][ydata][xdata][1], SALLY[time_counter][zdata][ydata][xdata][2]))
                vtkreplacement.write("{0[0]}\t {0[1]}\t {0[2]}\t {0[3]}\t {0[4]}\t {0[5]}\t {0[6]}\n".format((time,xdata*xdelta, ydata*ydelta, zdata*zdelta, SALLY[time_counter][zdata][ydata][xdata][0], SALLY[time_counter][zdata][ydata][xdata][1], SALLY[time_counter][zdata][ydata][xdata][2])))
#add time data also
    time_counter=time_counter+1
vtkreplacement.close()

