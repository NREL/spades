import numpy as np
import matplotlib.pyplot as plt
import copy
from sys import argv
from IPython import display
import time

nulleventval=6666666.0

def get_total_energy(lattice,B,N):
    
    totalenergy=0.0
    #do x sweep
    for i in range(N):
        for j in range(N):
            totalenergy+=B*(lattice[i][j-1] & lattice[i][j])
    #do y sweep
    for j in range(N):
        for i in range(N):
            totalenergy+=B*(lattice[i-1][j] & lattice[i][j])
    return(totalenergy)

def getrates(loc,lattice,B,N):
    i=int(loc/N)
    j=loc%N
    rates=np.zeros(4)+nulleventval;
    rateind=0
    #left,right and bottom,top
    for dim in range(2):
        for m in (-1,1):
            index1=np.array([i,j])
            index1[dim]=(index1[dim]+m)%N
            index1=index1.astype(int)
            if(lattice[tuple(index1)]==1):
                energy_at_faces1=0.0
                ii=index1[0]
                jj=index1[1]
                energy_at_faces1+=B*(lattice[ii-1][jj] & lattice[ii][jj])
                energy_at_faces1+=B*(lattice[(ii+1)%N][jj] & lattice[ii][jj])
                energy_at_faces1+=B*(lattice[ii][jj-1] & lattice[ii][jj])
                energy_at_faces1+=B*(lattice[ii][(jj+1)%N] & lattice[ii][jj])
        
                energy_at_faces2=0.0
                for dim1 in range(2):
                    for m1 in (-1,1):
                        index2=np.array([i,j])
                        index2[dim1]=(index2[dim1]+m1)%N
                        index2=index2.astype(int)
                        iii=index2[0]
                        jjj=index2[1]
                        if(iii!=ii or jjj!=jj):
                            energy_at_faces2+=B*(lattice[iii%N][jjj%N] & 1)

                rates[rateind]=energy_at_faces2-energy_at_faces1

            rateind+=1

    return(rates)

N=int(argv[1])
f=float(argv[2])
niters=int(argv[3])
B=-1.0
lattice=np.random.rand(N,N)
lattice[lattice<f]=0.0
lattice[lattice>=f]=1.0

lattice=lattice.astype(int)
#find number of zero locations
pause_time=0.1
for it in range(niters):
    plt.imshow(lattice)
    plt.pause(pause_time)
    #plt.savefig("fig%3.3d.png"%(it))
    zerolocs=np.array([])
    for i in range(N):
        for j in range(N):
            if(lattice[i][j]==0):
                zerolocs=np.append(zerolocs,N*i+j)

    nzeros=len(zerolocs)
    zerolocs=zerolocs.astype(int)
    eventlist=np.zeros((nzeros,4))
    for i in range(nzeros):
        rates=getrates(zerolocs[i],lattice,B,N)
        eventlist[i,:]=rates

    #find probable events
    eventprobs=np.array([])
    zlocs=np.array([])
    zlocmove=np.array([])
    for i in range(nzeros):
        for j in range(4):
            if(eventlist[i][j]!=nulleventval):
                eventprobs=np.append(eventprobs,eventlist[i][j])
                zlocs=np.append(zlocs,zerolocs[i])
                zlocmove=np.append(zlocmove,j)

    minrate_loc=np.argmin(eventprobs.flatten())
    minloc_id=zlocs[minrate_loc]
    minloc_move=int(zlocmove[minrate_loc])
    print("it,nzeros,nevents,min:",it,nzeros,len(eventprobs),eventprobs[minrate_loc])
    if(eventprobs[minrate_loc]==0):
        break

    swap_loc1=np.array([int(minloc_id/N),int(minloc_id%N)])
    swap_loc2=np.array([int(minloc_id/N),int(minloc_id%N)])

    minloc_move_dim=int(minloc_move/2)
    minloc_m=2*(minloc_move%2)-1
    swap_loc2[minloc_move_dim]=(swap_loc2[minloc_move_dim]+minloc_m)%N

    if(lattice[tuple(swap_loc1)]!=0 or lattice[tuple(swap_loc2)]!=1):
        print("Problems:",lattice[tuple(swap_loc1)],lattice[tuple(swap_loc2)])
    tempval=lattice[tuple(swap_loc1)]
    lattice[tuple(swap_loc1)]=lattice[tuple(swap_loc2)]
    lattice[tuple(swap_loc2)]=tempval
    

plt.show()
