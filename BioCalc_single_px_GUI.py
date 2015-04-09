################################################################################
### Script to calculate the cavity thicknes of a Biosensor fora single pixel ###
################################################################################

#############
### INPUT ###
#############


# enter average deviation of experiment to simulation in nanometer, "1" is a good value to start

tolerance = 1

# change version of the release here which will be included in the results files
version = 'BioCalc single px GUI 1.0'


#############################
#### start of the program ###
#############################

# load all the python modules needed
#import cython_all_fit as Fit # import self-written cython code
import matplotlib
matplotlib.use('TkAgg')
import cython_GUI as Fit
import numpy as np
import os 
from PIL import Image as im
from scipy import ndimage
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.filters import gaussian_filter
from skimage import transform
from matplotlib.image import AxesImage

matplotlib.use('TkAgg')

#####################
## GUI for BioCalc ##
#####################

from Tkinter import *
import ttk
import tkFileDialog




# define parameters for minima detection  

lookahead_min = 5 # something like peak width for the minima, 5 is a good value
lookahead_max = lookahead_min-1 # for the maxima --> should not be larger than lookahead_min
def f1(*args):
    for i in range(10):
        progress.step(5)

# define function to convert the image-string (read from file) to an array
def image2array(Img):
    
    newArr= np.fromstring(Img.tostring(), np.uint16)
    newArr= np.reshape(newArr, (Image_height,Image_width))

def onpick(event):
    x_pos.set(int(event.xdata))
    y_pos.set(int(event.ydata))
    print('onpick image', event.xdata,event.ydata,event.x,event.y)


def read_images(*args):

    progress = ttk.Progressbar(mainframe, orient=HORIZONTAL, length=500, mode='determinate')
    progress.grid(column=1,row=9, columnspan=4)    

    #progress = ttk.Progressbar(subframe_bottom, orient=HORIZONTAL, length=500, mode='determinate')
    #progress.grid(column=1,row=2)

    wave_s = int(wave_start_var.get())
    wave_e = int(wave_end_var.get())

    binning = binning_var.get()
    print binning_var.get()

    ###################
    # read image data #
    ###################

    data = tkFileDialog.askdirectory()

    # make wavelength list

    wave_step = 1       # [nm], this can be adjusted if you don't have an image evrey nanometer
    
    # create an empty list which will later contain all wavelengths, e.g. 550,551,552,...,750
    global waves
    waves=[]

    # write wavelengths into the list
    waves=[wave_s + i*wave_step for i in xrange((wave_e-wave_s)/wave_step + 1)]
    # make and sort list of file in the folder
    files = os.listdir(data)
    files.sort()

    # get size, bit-depth of the images
    for i in range(len(files)):
        # only consider files with the ending tiff, tif
        if files[i][-5:]=='.tiff' or files[i][-4:]=='.tif': 
            Img=im.open(data + '/' + files[i])
            break # stop the loop if an image was found


    Image_width = Img.size[0]/binning
    Image_height = Img.size[1]/binning

    # get colour and bit depth of image, this is important to know the range of the values (e.g. 8-bit is 0-255, 16-bit is 0-65535)
    Image_mode = Img.mode 
    if Image_mode == 'RGB' or Image_mode == 'P':
        Image_bit = '8'
    elif Image_mode == 'I;16' or Image_mode == 'I;12' or Image_mode=='I;16B':
        Image_bit = '16'
    else:
        print 'Image mode is:', Img.mode
        Image_bit = '16B'

    # create an empty array which will contain all image data 
    global all_images
    all_images = np.zeros(((wave_e-wave_s)/wave_step + 1,Image_height,Image_width),np.uint16)
    


    # read every image in folder and check if it is in the wavelength range --> write grey values into array
    
    # set a counter to check how many images have been processed
    counter=0
    # print some information for the user what is done
    print 'reading images from folder: ', data
    
    # start iterating over the files
    for i in xrange(len(files)):
        # update progress bar
        if i%20 == 0 and i!=0:
            progress.step(20.0/len(files)*100)
            root.update()
        # only consider files with the ending tiff, tif
        if files[i][-5:]=='.tiff' or files[i][-4:]=='.tif':
            # check if the current file is in the wavelength range the user specified at the beginning
            if int(files[i][:3]) >= wave_s and int(files[i][:3]) <= wave_e:
                #print files[i]
                #print counter

                # assign the current image name to a variable
                Img = data + '/' + files[i]
                # check if its 8-bit, convert, load
                if Image_bit == '8':
                    Img = im.open(data  + '/' + files[i])
                    Img = Img.convert('L')
                    all_images[counter]=transform.rescale(np.asarray(Img),1.0/binning,preserve_range=True).round().astype(np.uint16)

                # not sure, but I think I used "imread" because it is more platform independent
                else:
                    all_images[counter]=transform.rescale(imread(Img, as_grey=True),1.0/binning,preserve_range=True).round().astype(np.uint16)

                
                counter+= 1



####################################################
## Find new delta to deal with different dynamics ##
####################################################
               
    # use different delta scaling for 8,12,16 bit, delta something like the peak height of the minima, it differs of course significantly for 8-bit images and 16-bit images, just because of the different range

    # for 16 or 12 bit images do the following
    global delta
    if Image_bit == '16':
        # the proper delta shall be a 10th of the mean value of all images
        new_delta = int(all_images.mean()/10)
        # in case delta is larger than 7, take it
        if new_delta > 7:
            delta = new_delta

    # for 8-bit images just take a delta of 7 as this is sufficient
    if Image_bit == '8':
        delta = 7    

    # calculate how much the delta should be varied in case no value is found
    global delta_vary
    delta_vary = int(delta/5-delta/20) # this has be found empirically to lead to a good fit of the minima
    print 'delta = ', delta
    print "All images read"
    Image_h.set(Image_height)
    Image_w.set(Image_width)
    fig0 = plt.figure('wavelength image '+str(waves[len(waves)/2]) + ' nm')
    ax = fig0.add_subplot(111)
    ax.set_title('click on points')
    ax.imshow(all_images[len(waves)/2],cmap='Greys_r')

    fig0.canvas.mpl_connect('button_press_event', onpick)
    root.update()
    plt.show()

def callback():
    print "This is a callback"

def read_simulation_file(*args):

    # read simulation file 
    print 'read simulated data'

    # chose elastomer thickness range , the smaller the range the faster the program. If you are not sure, just take d_min = 1000, d_max = 19000

    d_min = 1000 # [nm]
    d_max = 20000 # [nm]

    # get the values of the thickness

    d_min = int(d_min_var.get())
    d_max = int(d_max_var.get())
    print d_min
    print d_max

    wave_s = wave_start_var.get()
    wave_e = wave_end_var.get()

    print wave_s
    print wave_e
    # open file
    p = tkFileDialog.askopenfile()
    head,tail = os.path.split(p.name)
    Simulation_file.set(tail)

    # read the whole content of the file into a string
    string=p.read()

    # close the file
    p.close()

    # list which will contain, thickness, length of a set of wavelengths, starting position of the block --> this is all useful to have to improve speed
    global thickness_len_pos
    thickness_len_pos=[]
    

    # current set of minima for one thickness
    minima_block=[]

    # list which contains all minima_blocks arrays one after another
    global list_minima_blocks
    list_minima_blocks = []


    # variable which defines at which position a set of minima starts
    position = 0

    # read every line of the string
    for thisline in string.split('\n'):
        # check if the current line contains a thickness
        if ('\t' in thisline) == False and len(thisline) != 0:
            thickness = thisline
        # check if this is a line with no thickness, but one is still in the thickness range
        if ('\t' in thisline) == True and int(thickness) >= d_min and int(thickness) <= d_max:
            # split line into words
            for word in thisline.split('\t'):
                # check word if it is a wavelengt, only append it if its in the wavelength range +- lookahead
                if len(word)<6 and float(word) >= wave_s + lookahead_min and float(word)<= wave_e - lookahead_min: # use only minima which are in the wave-range + lookahead_min
                    # append the wavelength to the current block of wavelengths
                    minima_block.append(float(word))

        # check if the line is empty and inside the thickness range
        if len(thisline) == 0 and int(thickness) >= d_min and int(thickness) <= d_max:

            # append thickness, length of waveblock, position of block to a list
            thickness_len_pos.append([np.uint16(thickness),np.uint16(len(minima_block)),np.uint32(position)]) # calculate length of the waveblock since it will be needed later

            # append waveblock to an array
            list_minima_blocks.append(np.array(minima_block,dtype=np.float))

            # update the current starting position of the next waveblock
            position += len(minima_block)
            # clear the waveblock to write new wavelengths into it
            minima_block=[]
    print "done"


def calculate(*args):
    #plt.close()
    # create sub array to fit only pixel in the middle
    y = int(y_pos.get())
    x = int(x_pos.get())
    global sim_thickness_error
    global single_px_3d_array 
    global result
    single_px_3d_array = all_images[:,y-0:y+1,x-0:x+1]

    tolerance = int(tolerance_var.get())
    ##########################
    # Single-Core Processing #
    ##########################

    # start row
    start = 0
    # last row
    ende = 1

    use_thickness_limits = False
    thickness_limit = 50
    area_avrg = 3
    # call the external cython/c++ function with all the parameters
    result, sim_thickness_error = Fit.c_Fit_Pixel(start,ende,single_px_3d_array, thickness_len_pos, waves, tolerance, lookahead_min, lookahead_max, delta,delta_vary,list_minima_blocks, use_thickness_limits, thickness_limit,area_avrg)

    print result
    thickness_result.set(result)
    ##############
    # plot datta #
    ##############

    # make a new figure
    #fig = plt.figure("row = " + str(y) + " column = " + str(x) )
    fig = plt.figure('intensity (lambda)')
    #ax1 = fig.add_subplot(1,1,1)

    plt.plot(waves,all_images[:,y,x],'.',label='x='+str(x)+', y='+str(y)+', thickness =' + str(int(result))+' nm')
    plt.legend(loc=4)
    if swp_var.get() == True:
        plt.vlines(list_minima_blocks[zip(*thickness_len_pos)[0].index(int(result))],plt.ylim()[0],plt.ylim()[1], color='r')
    # create plot of the results
    #plt.imshow(all_images[100])
    # set the color scale to the limits provided
    #plt.clim(all_images[100].mean()-color_min,all_images[100].mean()+color_max)
    # plot a color bar
    
    #plt.colorbar()

    # remove "#" to show the plot after the calculation
    if error_var.get() == True:
        print "Error plot should be plotted!"
        plt.figure('Error Plot')
        plt.plot(sim_thickness_error[0][1:],sim_thickness_error[1][1:],'.')

    

    plt.show()

def calculate_all(*args):

    global result

    # start row
    start = 0
    # last row
    ende = Image_h.get()


    # get values of variables

    use_thickness_limits = use_thickness_limits_var.get()

    thickness_limit = int(thickness_limit_var.get())

    tolerance = int(tolerance_var.get())

    area_avrg = 3

    print "use thickness limits: ", use_thickness_limits
    print "thickness_limit: ", thickness_limit
    print "tolerance: ", tolerance
    print "start: ", start
    print "ende: ", ende
    print "area_avrg: ", area_avrg
    print all_images
    # call external python function 

    result, sim_thickness_error = Fit.c_Fit_Pixel(start,ende,all_images, thickness_len_pos, waves, tolerance, lookahead_min, lookahead_max, delta,delta_vary,list_minima_blocks, use_thickness_limits, thickness_limit,area_avrg)

    # print
    plt.figure("Thickness Map")
    plt.imshow(result)
    plt.show()
# function to save calculated data into files
def save_results(*args):

    y = int(y_pos.get())
    x = int(x_pos.get())
    p = open("Intensity_Profile_"+str(x)+"_"+str(y)+".txt",'w')
    p.write('Wavelength [nm]\tIntensity\n')
    for i in range(len(waves)):
        p.write(str(waves[i])+"\t"+str(all_images[:,y,x][i])+'\n')
    p.close()

    if swp_var.get() == True:
        np.savetxt("Sim_Minima_Pos_"+str(int(result))+"nm_"+str(x)+"_"+str(y)+".txt",list_minima_blocks[int(result)-1],fmt='%1.1f')

    if error_var.get() == True:
        # array to write thickness vs error as columns 
        a = np.array((sim_thickness_error[0][1:],sim_thickness_error[1][1:]))
        np.savetxt("Thickness_vs_Error_"+str(x)+"_"+str(y)+".txt",a.transpose(),fmt=['%i','%1.4f'],delimiter="\t")

    if len(result)>2:
        np.savetxt("result.txt",result,fmt='%i')

def plot_wave_image(*args):

    wave_image = int(plot_wave_entry.get())

    if (wave_image in waves) == True: # check if entered wavelength is in the list
        fig0 = plt.figure('wavelength image '+str(wave_image) + ' nm')
        ax = fig0.add_subplot(111)
        ax.set_title('click on points')
        ax.imshow(all_images[waves.index(wave_image)],cmap='Greys_r')

        fig0.canvas.mpl_connect('button_press_event', onpick)
        root.update()
        plt.show()
    else:
        print wave_image


root = Tk()
root.title('BioCalc')

# make the main frame
mainframe = ttk.Frame(root,padding="3 3 12 12")
mainframe.pack(side=LEFT)
mainframe.columnconfigure(0,weight=1)
mainframe.rowconfigure(0,weight=1)

# make a subframe at the right side
subframe_right = ttk.Frame(root,padding=" 3 3 12 12")
subframe_right.pack(anchor=E)

# make the subframe at the bottom
#subframe_bottom = ttk.Frame(root,padding=" 3 3 12 12")
#subframe_bottom.pack(side=BOTTOM)



# make variables which can be changed by entering a string
wave_start_var = IntVar()
wave_end_var = IntVar()    
wave_step_var = IntVar()
# set default values
wave_start_var.set(550)
wave_end_var.set(750)
wave_step_var.set(1)

x_pos = IntVar()
y_pos = IntVar()    

x_pos.set(0)
y_pos.set(0)

thickness_result = IntVar()

# build entry fields and labels
ttk.Entry(mainframe, width=3, textvariable=wave_start_var).grid(column=2,row=1,sticky=(W,E))
ttk.Label(mainframe, text="wave start [nm]").grid(column=1, row=1, sticky=W)


ttk.Entry(mainframe, width=3, textvariable=wave_end_var).grid(column=2,row=2,sticky=(W,E))
ttk.Label(mainframe, text="wave start [nm]").grid(column=1, row=2, sticky=W)

ttk.Entry(mainframe, width=3, textvariable=wave_step_var).grid(column=2,row=3,sticky=(W,E))
ttk.Label(mainframe, text="wave step [nm]").grid(column=1, row=3, sticky=W)

ttk.Entry(mainframe, width=7, textvariable=x_pos).grid(column=2,row=4,sticky=(W,E))
ttk.Label(mainframe, text="x").grid(column=3, row=4, sticky=W)

ttk.Entry(mainframe, width=7, textvariable=y_pos).grid(column=2,row=5,sticky=(W,E))
ttk.Label(mainframe, text="y").grid(column=3, row=5, sticky=W)


# build buttons
ttk.Button(mainframe, text="2. Load Image Data", command=read_images).grid(column=1, row=5, sticky=W)

ttk.Button(mainframe, text="1. Load Simulation File", command=read_simulation_file).grid(column=1, row=4, sticky=W)


Simulation_file = StringVar()

ttk.Label(mainframe, textvariable=Simulation_file).grid(column=1,row=8, columnspan = 3)
# ttk.Label(subframe_bottom, textvariable=Simulation_file).grid(column=1,row=1)

# check box for plotting fitted thickness wavelengths, error plot

swp_var = BooleanVar() # simulated wavelength position variable
Checkbutton(mainframe, text="Wavelength positions", variable=swp_var).grid(column=4, row=5, sticky=W)
error_var = BooleanVar()
Checkbutton(mainframe, text="Error plot", variable=error_var).grid(column=4, row=4, sticky=W)

ttk.Label(mainframe, textvariable=thickness_result).grid(column=2, row=6, sticky=(E))
ttk.Label(mainframe, text='[nm]').grid(column=3, row=6, sticky=(W, E))
ttk.Button(mainframe, text="3. Calculate single", command=calculate).grid(column=1, row=6, sticky=W)

ttk.Button(mainframe, text="3. Calculate all", command=calculate_all).grid(column=1, row=7, sticky=W)

Image_h = IntVar()
Image_w = IntVar()

ttk.Label(mainframe, textvariable= Image_h).grid(column=3,row=1)
ttk.Label(mainframe, text= '[px] Image height').grid(column=4,row=1,sticky=W)
ttk.Label(mainframe, textvariable= Image_w).grid(column=3,row=2)
ttk.Label(mainframe, text= '[px] Image width').grid(column=4,row=2, sticky=W)


# build the progressbar in the subframe_bottom
progress = ttk.Progressbar(mainframe, orient=HORIZONTAL, length=500, mode='determinate').grid(column=1,row=9, columnspan=4)
# progress = ttk.Progressbar(subframe_bottom, orient=HORIZONTAL, length=500, mode='determinate').grid(column=1,row=2)
#progress.grid(column=1,row=10)



ttk.Button(mainframe, text="4. Save results", command=save_results).grid(column=4, row=6, sticky=W)

plot_wave_var = IntVar()
plot_wave_entry = ttk.Entry(mainframe, width=4, textvariable=plot_wave_var)
plot_wave_entry.grid(column=3,row=3,sticky=W)
ttk.Button(mainframe, text="Plot Wave", command=plot_wave_image).grid(column=4, row=3, sticky=W)


# subframe right

d_min_var = IntVar()
d_max_var = IntVar()

d_min_var.set(1000)
d_max_var.set(15000)

use_thickness_limits_var = BooleanVar()

binning_var = IntVar()
ttk.Label(mainframe, text= 'Binning').grid(column=5,row=1,sticky=W)
binning_box = Spinbox(mainframe, width = 3,values = (1,2,4,8,16),textvariable=binning_var).grid(column=6,row=1)

ttk.Label(mainframe, text= 'Min thickness [nm]').grid(column=5,row=2,sticky=W)
d_min_entry = ttk.Entry(mainframe, width=5, textvariable=d_min_var).grid(column=6,row=2)

ttk.Label(mainframe, text= 'Max thickness [nm]').grid(column=5,row=3,sticky=W)
d_max_entry = ttk.Entry(mainframe, width=5, textvariable=d_max_var).grid(column=6,row=3)

Checkbutton(mainframe, text="Use thickness limits", variable=use_thickness_limits_var).grid(column=5, row=4, sticky=W)

thickness_limit_var = IntVar()   
thickness_limit_var.set(50)
ttk.Entry(mainframe, width=5, textvariable=thickness_limit_var).grid(column=6,row=5,sticky=(W,E))
ttk.Label(mainframe, text="thickness limit [nm]").grid(column=5, row=5, sticky=W)

tolerance_var = IntVar()
tolerance_var.set(1)
ttk.Entry(mainframe, width=5, textvariable=tolerance_var).grid(column=6,row=6,sticky=(W,E))
ttk.Label(mainframe, text="tolerance [nm]").grid(column=5, row=6, sticky=W)


# make boardeers around fields in the main frame
for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)


#wave_start_var_entry.focus()
root.bind('<Return>', f1)

root.mainloop()
