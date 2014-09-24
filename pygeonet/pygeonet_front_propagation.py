import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from osgeo import gdal
from time import clock
import scipy.signal as conv2
import numpy.ma as npma

'''
%   [D,S] = perform_front_propagation_2d_slow(W,start_points,\
%            end_points,nb_iter_max,H);
%
%   [The mex function is perform_front_propagation_2d]
%
%   'D' is a 2D array containing the value of the distance function to seed.
%	'S' is a 2D array containing the state of each point : 
%		-1 : dead, distance have been computed.
%		 0 : open, distance is being computed but not set.
%		 1 : far, distance not already computed.
%	'W' is the weight matrix (inverse of the speed).
%	'start_points' is a 2 x num_start_points matrix where k is
%                       the number of starting points.
%	'H' is an heuristic (distance that remains to goal). This is a 2D matrix.
%   
%   Copyright (c) 2004 Gabriel Peyr
The code is modified to suit syntax as per Python by
Harish Sangireddy
'''

def perform_fast_marching_step(D,O,S,W):
    print 'fast_matching step'
    kClose = -1
    kOpen = 0
    kFar = 1
    n,p = D.shape
    # step size
    h = 1/n
    if O.size ==0:
        D1 = D
        #break
    # selected vertex
    i = O[0]
    j = O[1]
    print O
    S[i,j] = kClose # now its close
    # its neighbor
    nei = np.array([[1,0],[0,1], [-1,0],[0,-1]])
    print nei
    for k in range(0,3):
        ii = i+ nei[k,0]
        jj = j + nei[k,1]
        print ii,jj
        if ii > 0 and jj > 0 and ii <= n and jj <=p:
            f = np.array([0,0]) # current father
            # update using upwind resolution
            P = h/W[ii,jj]
            # neighbors values
            a1 = np.Inf
            if ii < n:
                a1 = D[ii+1,jj]
                f[0] = const1
            if ii > 1 and D[ii-1,jj] < a1:
                a1 = D[ii-1,jj]
                f[1] = const2
            a2 = np.Inf
            if jj < n:
                a2 = D[ii,jj+1]
                f[1] = const3
            if jj > 1 and D[ii,jj-1]<a2:
                a2 = D[ii,jj-1]
                f[2] = const4
            # swap to reorder
            if a1 > a2:
                tmp = a1
                a1 = a2
                a2 = tmp
                f = f[1,0]
            # Now the equation is
            # (a-a1)^2+(a-a2)^2 = P, with a >= a2 >= a1.
            if P**2 > (a2-a1)**2:
                delta = 2 * P**2 - (a2-a1)**2
                A1 = (a1+a2+np.sqrt(delta))/2
            else:
                A1 = a1 + P
                f[2] = 0
            # Implementing the switch case as a if statement
            if S[ii,jj]==kClose:
                # check if action has change.
                #Should not append for FM
                if A1 < D[ii,jj]:
                    if False:
                        O[:,]
            elif S[ii,jj]==kOpen:
                # check if action has change
            elif S[ii,jj]== kFar:
                # add to open
            else:
                print 'Unknown state'
            
                
                
        
    return D1

def perform_front_propagation(W,start_points,nb_iter_max):
    print 'Performing fast marching'
    print 'W is the weight matrix (inverse of the speed)'
    print ' '
    D = W*0 + np.Inf
    print type(W),W
    D[start_points[0],start_points[1]]=0
    O = np.array(start_points)

    # S=1 : far,  S=0 : open,  S=-1 : close
    S = np.ones((W.shape))
    S[start_points[0],start_points[1]] = 0

    i = 0
    while i< nb_iter_max and O.size!=0:
        i = i +1
        data = perform_fast_marching_step(D,O,S,W)
    # End of while
    
    
    plt.figure(1)
    plt.imshow(D,cmap=cm.coolwarm)
    plt.plot(start_points[1],start_points[0],'or')
    plt.title('W')
    plt.colorbar()
    plt.show()
    
    

def main():
    outlettif = 'C:\\Mystuff\\grassgisdatabase\\ikawa_roi1_nutm54_clipped_costfunction.tif'
    dsout = gdal.Open(outlettif, gdal.GA_ReadOnly)
    aryfdrout = dsout.GetRasterBand(1).ReadAsArray()
    nanDemArray=np.array(aryfdrout.T)
    xS,yS = nanDemArray.shape
    W = nanDemArray
    start_points = [127,1006]
    #end_points = []# we basically don't use this
    nb_iter_max = 1.2* np.max(np.array([xS, yS]))**2
    # H = hueristic# we are also not going to use this
    perform_front_propagation(W,start_points,nb_iter_max)
    



if __name__ == '__main__':
    t0 = clock()
    main()
    t1 = clock()
    print "time taken to complete the script is::",t1-t0," seconds"
    print "script complete"
