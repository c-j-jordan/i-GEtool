import numpy as np
import math
import matplotlib.pyplot as plt

def read(tstep):
    cell_no = 100
    cell_size = 25
    max_iso = 9
    


    H = np.zeros((cell_no,cell_no),dtype=float)
    He = np.zeros((cell_no,cell_no),dtype=float)
    metallicity = np.zeros((cell_no,cell_no),dtype=float)

    fname = 'gas_t' + str(tstep) + '.dat'
    data = None
    file = open(fname, 'r')
    x = np.zeros(cell_no,dtype=float)
    y = np.zeros(cell_no,dtype=float)
    z = np.zeros(cell_no,dtype=float)
    mass = np.zeros((cell_no,cell_no),dtype=float)
    logmass = np.zeros((cell_no,cell_no),dtype=float)
    masst = np.zeros((cell_no,cell_no),dtype=float)
    logmasst = np.zeros((cell_no,cell_no),dtype=float)
    count = 0


    # Reads in z y x
    data  = np.loadtxt(fname ,usecols=(0,1,2,3,4,5,6,7,8), dtype=float)
    total_mass=0
    total_masst=0
    for i in range(0,cell_no):
        for j in range(0,cell_no):
            H[i][j] = data[count][0]
            He[i][j] = data[count][1]
            metallicity[i][j] = data[count][8]
            x[i] = cell_size*i
            y[j] = cell_size*j
            count+=1

    for i in range(0,cell_no):
        for j in range(0,cell_no):
            mass[i][j] = H[i][j] + He[i][j] + metallicity[i][j]
            total_mass += mass[i][j]
            logmass[i][j] = math.log10(mass[i][j])

    fname = 'gas_t' + str(tstep+1) + '.dat'
    datat  = np.loadtxt(fname ,usecols=(0,1,2,3,4,5,6,7,8), dtype=float)
    count = 0
    for i in range(0,cell_no):
        for j in range(0,cell_no):
            H[i][j] = datat[count][0]
            He[i][j] = datat[count][1]
            metallicity[i][j] = datat[count][8]
            count+=1

    for i in range(0,cell_no):
        for j in range(0,cell_no):
            masst[i][j] = H[i][j] + He[i][j] + metallicity[i][j]
            total_masst += masst[i][j]
            logmasst[i][j] = math.log10(masst[i][j])
    asdf = 0
    for i in range(79,99):
        for j in range(59,79): 
            asdf += masst[i][j] - mass[i][j]
    qwer = masst[50][50] - mass[50][50]
    print('Mass difference at centre %d' % qwer)
    qwer = masst[99][99] - mass[99][99]
    print('Mass difference at edge %d' % qwer)

    import pylab as py
    print ('Mass difference of %d' % asdf)
    py.title("Mass Distribution")
    py.xlabel("$x$ (pc)")
    py.ylabel("$y$ (pc)")
    py.pcolor(x,y,(masst[:][:]-mass[:][:]))
    py.colorbar(label="Mass Difference (M$_\odot$)")
    py.set_cmap('Blues')
    py.xlim(0, 2485)
    py.ylim(0,2485)
    title = 'Mass Distribution Between Times ' + str(tstep) + ' and ' + str(tstep+1) + ' Myr'
    py.title(title)
    py.show()
    py.clf()
