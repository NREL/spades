import argparse
from ase.build import bulk
from ase.neighborlist import NeighborList
import numpy as np
import random
import matplotlib.pyplot as plt
from copy import copy as copy
from itertools import accumulate
from bisect import bisect
import pathlib
import inspect


# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from timeit import default_timer as timer
import kmc_diffusion as kmcd

module_dir = pathlib.Path(inspect.getfile(kmcd)).parents[0]
base_dir = module_dir.parents[0]

def make_lattice(Nx,Ny,Nz,f_vacant):
    # f_vacant = fraction of vacant sites
    # Use ASE to create a simple cubic lattice and make N on the bottom
    # for what ultimately is going to be an alternating N-B-N-B surface
    spacing = 1
    fullcell = bulk('O','sc',a=spacing,cubic=True) # Here 'O' means "occupied"
    lattice = fullcell.repeat((Nx,Ny,Nz))
    lattice.set_pbc((True,True,True))

    xvec = np.arange(0,Nx,spacing)
    yvec = np.arange(0,Ny,spacing)
    zvec = np.arange(0,Nz,spacing)

    #print(xvec)
    #print(yvec)
    # Generate a neighborlist for each site
    cutoff_list = [0.4*spacing/np.sqrt(2.) for i in range(len(lattice))]
    #print(cutoff_list)
    nl = NeighborList(cutoff_list,self_interaction=False,bothways=True)
    nl.update(lattice)
    print("Neighborlist:")
    #for i in range(Nx*Ny*Nz):
    #  indices, offsets = nl.get_neighbors(i)
    #  print(i, indices)

    # Initialize the defect-containing cell by randomly setting sites to 'Ar' with probability f_vacant
    # Initialize random number generator
    for i in range(Nx*Ny*Nz):
      if (random.random() < f_vacant):
        lattice.symbols[i] = 'Ar'

    # Uncomment the stuff below to all sites in the left fraction
    # to be vacant, the rest filled
    '''
    ndiv = int(f_vacant*float(Nx))
    for nx in range(Nx):
      if (nx <= ndiv):
        site = 'Ar'
      else:
        site = 'O'
      for ny in range(Ny):
        for nz in range(Nz):
          i = Ny*Nz*nx +  Nz*ny + nz
          lattice.symbols[i] = site
    '''

    return lattice,nl,xvec,yvec,zvec

def build_site_events(i,lattice,nl,ratelist,beta,nn_delta_e,ndim):
    
    hop_rate_not = ratelist[0]
    psum = 0.0
    events_i = []
    itype = lattice.symbols[i]
    # The block below just looks for diffusion events - that's all there is right now
    event_type = 1
    if (itype == 'Ar'):
        # See if there is any unoccupied site next to i to which to stick
        indices,offsets = nl.get_neighbors(i)
        for j in indices:
            if (lattice.symbols[j] == 'O'): # position j is 'O' so we can try to move the vacancy there
                nneighbors = 0
                # See if site j has more than 1 non-Ar neighbor to keep us from "diffusing" into thin air
                jindices,joffsets = nl.get_neighbors(j)
                for k in jindices:
                    if (lattice.symbols[k] != 'Ar'): #Could speed this up with a pythonic count-y thing
                        nneighbors = nneighbors + 1
                if (nneighbors > 1): # Site j is a candidate for the vacancy to go here
                    ineighbors = 0
                    for l in indices:
                        if (lattice.symbols[l] != 'Ar'):
                          ineighbors = ineighbors + 1
                    # Define nearest-neighbor energy for energy model:
                    # e(i) = nn(i)*nn_delta_e 
                    # ndim = 6 in 3d, 4 in 2d, 2 in 1d
                    # With this definition, a site with fewer than ndim neighbors will have a positive
                    # energy if nn_delta_e is positive.
                    # So, positive nn_delta_e will cause clusters of vacancies to be more stable 
                    # Note that hop_rate_not = attempt_freq*np.exp(-beta*barrier_energy)
                    # so that the hop_rate = attempt_freq*exp(-beta*(barrier_energy - e_i))
                    # this means that the ratio of hopping rates i->j to that for j->i obeys
                    # detailed balance
                    site_i_energy = nn_delta_e*ineighbors 
                    boltzmann = np.exp(beta*(site_i_energy)) #Minus sign is exp(-beta*-site_i_energy)
                    #print("Site {} hopping to site {}".format(i,j))
                    rate = hop_rate_not*boltzmann
                    psum += rate
                    events_i.append([event_type,rate,i,j])
    
    return psum, events_i


def atom_counts(lattice,type1,type2):
    one_count = 0
    two_count = 0
    for i in range(len(lattice)):
        atomtype = lattice.symbols[i]
        if atomtype == type1:
            one_count = one_count + 1
        elif atomtype == type2:
            two_count = two_count + 1

    print("{} occupied and {} unoccupied sites".format(one_count,two_count))
    return one_count,two_count


def KMC_find_site(propensity,kmc_algorithm='rfKMC'):
    '''
    Three KMC methods are available: 
    =================================================================================
    First reaction MC method (used by Shirazi and Elliott for their simulations)
    kmc_algorithm = 'FirstReaction'
     
                Form a list with each event rate (all_rates, above) and a rate called k_i

                Choose a random number, rho_i, [0,1) for each event and form deltat_i = ln(rho_i)/k_i

                Perform the event with the smallest delta_t and repeat the loop

    =================================================================================
    Standard rejection-free KMC method (this already exists in SPPARKS as just KMC)
    kmc_algorithm = 'rfKMC'
      
    =================================================================================
    Rejection KMC method
    kmc_algorithm = 'R_KMC' (looks different to avoid rfKMC confusion/typos)
        
        Form a list of all events, i, with total of N events possible
        
        Choose a random number uniformly on [0,N) and select that event to perform.
        Call this event i.

        Accept or reject this event with probability r_i/rmax where rmax is some maximum possible
        rate. Here we will take rmax to be the maximum of all the rates in our list. This appears
        to be an arbitrary choice...so I'll do this.
    '''

    if (kmc_algorithm != 'FirstReaction' and kmc_algorithm != 'rfKMC' and kmc_algorithm != 'R_KMC'):
        print("Error - calling apply_KMC with an incorrect kmc_algorithm specification")
        exit(1)

    if (kmc_algorithm=='FirstReaction'):
        #Do a single-shot formation of dt_array and find the event with the minimum time
        dt_array = np.asarray([-np.log(random.random())/prop for prop in propensity])
        enumber = np.argmin(dt_array)
        dt = dt_array[enumber]
        #print("doing event {} with timestep {}: {}".format(enumber,dt,all_rates[enumber]))
        #Update the surface by enacting "do_this_event"
        return enumber,dt
    if (kmc_algorithm=='rfKMC'):
        cumulative_rates = np.cumsum(propensity)
        total_rate = cumulative_rates[-1]
        #Now find the event for whice cumulative_rates(event) < cum_choice < cumulative_rates(event+1)
        u = random.random()
        cum_choice = u*total_rate
        #return position to insert cum_choice into cumulative_rates in order to maintain sort
        #this can return a value greater than the largest element in the array if u=1.0 so we need to check for this #with the if statement below and correct by just taking the final event
        enumber = bisect(cumulative_rates,cum_choice)
        #Bisect may return len() if we get unlucky and u==1. Fix for this is to pick the final event
        if (enumber==len(cumulative_rates)):
            enumber = len(cumulative_rates)-1
        #print("Choosing event {} ({})".format(enumber,do_this_event))
        #Update time: delta_t = (1/total_rate)*ln(1/u') = -ln(u')/total_rate
        uprime = random.random()
        delta_t = (-1.0/total_rate)*np.log(uprime)
 
        return enumber,delta_t
    if (kmc_algorithm=='R_KMC'):
        print("In R_KMC")
        rmax = np.max(propensity)
        enumber = random.randint(0,len(propensity)-1)
        #rate for this event = propensity(enumber). Do this event with probability rate/rmax
        prob = propensity[enumber]/rmax # This is a made up probability - there are other ways to set the probability (e.g., an energy criterion) so this is a placeholder to test the R_KMC functionality
        if (random.random() < prob):
            #print("R_KMC - event at {} accepted".format(enumber))
            uprime = random.random()
            delta_t = - np.log(uprime)/(len(propensity)*rmax)
            return enumber, delta_t
        else:
            #print("R_KMC - event at {} rejected".format(enumber))
            return -1,0.0

def do_event_at_site(i,lattice,nl,ratelist,beta,nn_delta_e,ndim,site_events,propensity):
    '''
    Pick from the possible events at site i, update the lattice with the event
    that was chosen, then update the site propensities for site i, the other site (j)
    and the unique neighbors of i and j 

    This function assumes diffusion of vacancies ('Ar' site -> 'O' and vice versa)
    To extend this, use the event_type info in event[0]
    '''
    site_props = np.asarray([event[1] for event in site_events[i]])
    cumulative_props = np.cumsum(site_props)
    #print(cum_props)
    u = random.random()
    cumulative_choice = u*cumulative_props[-1]
    enumber = bisect(cumulative_props,cumulative_choice)
    if (enumber == len(cumulative_props)):
        enumber = len(cumulative_props)-1
    event = site_events[i][enumber]
    print("...doing event {} at site {}: {}".format(enumber,i,event))

    # Update lattice sites i and j from site_events[i][enumber
    site_i = event[2]
    site_j = event[3]
    lattice.symbols[site_i] = 'O'
    lattice.symbols[site_j] = 'Ar'

    # Update propensity and events at site_i, site_j, and the unique list of neighbors of
    # both site_i and site_j
    site_update_list = [site_i,site_j]
    ineighbors, offsets = nl.get_neighbors(site_i)
    for k in ineighbors:
      site_update_list.append(k)
    jneighbors, offsets = nl.get_neighbors(site_j)
    for l in jneighbors:
      site_update_list.append(l)

    unique_update_list = list(set(site_update_list))
    #print(unique_update_list)
    for site in unique_update_list:
        psum,events = build_site_events(site,lattice,nl,ratelist,beta,nn_delta_e,ndim)
        propensity[site] = psum
        site_events[site] = events
    return propensity,site_events

def latcuts_xyz_planes(lattice,Nx,Ny,Nz):
    latcut_xy = np.zeros((Nx,Ny))
    latcut_xz = np.zeros((Nx,Nz))
    latcut_yz = np.zeros((Ny,Nz))
    for nx in range(Nx):
      for ny in range(Ny):
        index = Ny*Nz*nx +  Nz*ny
        if (lattice.symbols[index] == 'Ar'):
          ltype = 0
        else:
          ltype = 1
        latcut_xy[nx][ny] = ltype
      for nz in range(Nz):
        index = Ny*Nz*nx + nz
        if (lattice.symbols[index] == 'Ar'):
          ltype = 0
        else:
          ltype = 1
        latcut_xz[nx][nz] = ltype
    for ny in range(Ny):
      for nz in range(Nz):
        index =  Nz*ny + nz 
        if (lattice.symbols[index] == 'Ar'):
          ltype = 0
        else:
          ltype = 1
        latcut_yz[ny][nz] = ltype

    return latcut_xy,latcut_xz,latcut_yz

def latcuts_xy_planes(lattice,nz,Nx,Ny,Nz):
    latcut_xy = np.zeros((Nx,Ny))
    for nx in range(Nx):
      for ny in range(Ny):
        index = Ny*Nz*nx +  Nz*ny + nz
        if (lattice.symbols[index] == 'Ar'):
          ltype = 0
        else:
          ltype = 1
        latcut_xy[nx][ny] = ltype
    return latcut_xy


def main():
    '''
    Simple prototype to make a lattice with vacancies and vacancy diffusion
    Occupied = 'O', Unoccupied or vacant = 'Ar' - lets us use ase to make the neighbor list
    '''
    start = timer()
    random.seed(34950435)

    parser = argparse.ArgumentParser(description="A simple KMC diffusion application")
    parser.add_argument(
        "-p", "--do_plots", help="Make plots", action="store_true"
    )
    parser.add_argument(
        "-s", "--save_plots", help="Save plots", action="store_true"
    )
    args = parser.parse_args()


    #Set microscopic parameters that give various rates and set sticking probabilities for deposition

    # Set externally applied parameter(s): T
    # Use SI units so deposition parameters are easy to compare to experiment: meters, kg, seconds
    # enter m_O and m_C in molar units then convert to kg/atom
    # and enter P_O and P_C in bar and convert to SI unit.
    # convert g/mol to kg/atom
    avogadro = 6.02214076e+23
    gpermol_to_kgperatom = 1.0e-3/avogadro
    bar_to_pascals = 1.0e5 # 1 bar = 10^5 Pa
    kB = 1.380649e-23 # Boltzmann's constant in J/K
    ev_to_J = 1.602176634e-19 # Convert eV to J

    data_dir = base_dir / "data"
    if args.save_plots:
        data_dir.mkdir(parents=True, exist_ok=True)

    # User-chosen parameters are below
    T = 300.0 # temperature in Kelvin
    kT = kB*T/ev_to_J # Convert kT to eV (means we should use energies in eV everywhere)
    beta = 1.0/kT
    print("T = 300 K, so kT = {} eV (beta = {} / eV)".format(kT,beta))

    # Set some lattice parameters
    Nx = 9
    Ny = 9
    # Nz = 9 or a multiple is nice for later because this code will display a 3x3 grid of XY slice
    Nz = 9
    Nsites = Nx*Ny*Nz
    ndim = 6 # Number of nearest neighbors depending on the number of active dimensions (2, 4, 6)
    f_vacant = 0.2 # Fraction of sites wtih vacancies

    # Define nearest-neighbor energy for energy model:
    # e(i) = nn(i)*nn_delta_e 
    # ndim = 6 in 3d, 4 in 2d, 2 in 1d
    # With this definition, a site with fewer than ndim neighbors will have a negative
    # energy if nn_delta_e is positive. So positive nn_delta_e will cause clusters of vacancies
    # to be more stable 
    nn_delta_e = 0.03  # absolute value of 1/2 the bond energy in eV (bond en = negative)
    print("T = 300 K so nn_delta_e ({} eV)//kT = {}".format(nn_delta_e,nn_delta_e/kT))

    # attempt frequency (say 10^12 / s) ~vibrational frequencies
    # and a hopping energy barrier, E_hop_not
    attempt_freq = 1.0e12

    E_hop_not = 0.2 # interbasin hopping barrier in eV. This could be 
    hop_rate_not = attempt_freq*np.exp(-beta*E_hop_not)
    print("Rate of Ar hops with no site energy changes would be {:e}".format(hop_rate_not))

    kmc_algorithm = 'rfKMC'
    #kmc_algorithm = 'FirstReaction'
    #kmc_algorithm = 'R_KMC'

    # Set either a maximum time or a maximum number of steps to take
    end_criterion = 'TIME'
    maxval = 2.0e-12 # Time in whatever units 1/rate gives.

    #end_criterion = 'STEPS'
    #maxval = 20

    print("End criterion is {}, maxval = {}".format(end_criterion,maxval))

    # Randomly place vacancies (see function for a commented out example of making a bilayer)
    lattice,nl,X,Y,Z = make_lattice(Nx,Ny,Nz,f_vacant)
    lattice_info = [X,Y,Nx,Ny,Nz] # Convenient list of necessary lattice parameters for some functions
    
    Xmg,Ymg = np.meshgrid(X,Y)

    print("Initial lattice:")
    nO,nAr = atom_counts(lattice,'O','Ar')
    print(lattice.symbols)

    '''
    Set things up to diffuse vacancies with a TST model: aattempt frequency and barrier.
    The rates will be a list of lists with each entry being a possible event and each event being a list with a format [event_type,rate,*], where * will be a any additional information
    such as reaction site or sites involved in the event. Note that in this diffusion model,
    "event_type" will always be "diffusion"

    --------------------------------------------------------------------

    The rates we define are for the following event types (can always add new event types)

    1. Rate for Ar to move to a neighboring empty site (must be an empty 'Ar' site
	that is itself a neighbor of an occupied site 'O').
    '''

    # Pick the kmc algorithm used in apply_KMC
    # default is kmc_algorithm = 'FirstReaction'
    #kmc_algorithm = 'R_KMC' # R_KMC is pretty inefficient for the parameter sets I've tested in this application
    
    spacing = 4.0e-10 # site spacing between sites in m. This is arbitrary, only used if we compute diffusion coefficient

    print("The site spacing is {} cm, so this corresponds to D = rate*d^2 m^2/s = {} for vacancies".format(spacing,spacing*spacing*hop_rate_not))

    ratelist = [hop_rate_not]
    maxrate = max(ratelist)
    print("The largest rate coefficient from ratelist is is {:.2e} / s".format(maxrate))

    print("Rate parameters: {}".format(ratelist))

    #Set up a run for some length of time (define it in terms of number of diffusers and in terms of
    # 1/maxrate?

    # Here do stuff that will be in a loop 
    # Set the max time to be the expected number of vacancies divided by the maximum rate, time
    # the number of hopping events we would expect
    #tmax = (Nx*Ny*Nz)*f_vacant/maxrate
    print("1/maxrate = {}".format(1.0/maxrate))
    time = 0.0

    output_time = 0.0;
    otime = 0.0
    print("Put out some information every {} seconds".format(output_time))	


    # Plot a 3x3 grid of XY lattice planes spread over 9 z values
    skip = int(float(Nz)/9.)
    print("skip = ",skip)
    zvals = range(0,Nz,skip)
    print(*zvals)
    if args.do_plots:
        fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(nrows=3,ncols=3,sharex=True,sharey=True)#,constrained_layout=True)
        plt.subplots_adjust(hspace=0.4,wspace=0.)
        ax_set = np.asarray(((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)))
        ax_imgs = []
        
        for nz, ax in enumerate(np.ndarray.flatten(ax_set)):
          z = zvals[nz]
          latcut_xy = latcuts_xy_planes(lattice,nz,Nx,Ny,Nz)
          if (ax == ax2):
            ax.set_title("time = {:.2e}\nnz = {}".format(0.0,z))
          else:
            ax.set_title("nz = {}".format(z))
          img = ax.imshow(latcut_xy)
          ax_imgs.append(img)
        if args.save_plots:
            plt.savefig(data_dir / f"{0:05d}.png")
        plt.pause(1)

    # Define a list with Nsites elements. Each element is itself a list of the possible events at
    # the site: [type, rate, i, j, k,...]. So site_events[i] = [[event_1],[event_2],...,[final_event] ]
    # If there are no events at site i, then site_events[i] = [] - a list with no elements
    site_events = [[] for i in range(Nsites)]
    #for i in range(Nsites):
    #  site_events.append([])

    # Define an array that contains the propensity for each site
    propensity = np.zeros(Nsites)

    #print("propensity: ",propensity)
    #print("site_events: ",site_events)

    # Now fill in propensity and site_events
    ptot = np.sum(propensity)
    for i in range(Nsites):
      psum,events = build_site_events(i,lattice,nl,ratelist,beta,nn_delta_e,ndim)
      propensity[i] = psum
      site_events[i] = events

    #print("Filled propensity: ",propensity)
    #print("Filled site_events: ",site_events)

    #print("Events by site:")
    #for i in range(len(site_events)):
    #  if site_events[i] != []:
    #    print(i,propensity[i],site_events[i],"\n")

    timestep_array = []
    nstep = 0
    #while (time < tmax):
    #while (nstep < nmax):
    run_state = True
    while (run_state):
      site, dt = KMC_find_site(propensity,'rfKMC')
      if (site != -1):
        print("We will do an event at site {} and increment dt by {}".format(site,dt))
        propensity,site_events = do_event_at_site(site,lattice,nl,ratelist,beta,nn_delta_e,ndim,site_events,propensity)
        time += dt
        nstep += 1
        timestep_array.append(dt)
        print(time)
        if (end_criterion == 'TIME' and time >= maxval):
            run_state = False
        if (end_criterion == 'STEPS' and nstep >= maxval):
            run_state = False

      otime += dt
      if (otime >= output_time):
        print("Time = {}".format(time))
        otime = 0.0

        if args.do_plots:
            ax_set = np.asarray(((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)))
            for nz, ax in enumerate(np.ndarray.flatten(ax_set)):
              z = zvals[nz]
              latcut_xy = latcuts_xy_planes(lattice,nz,Nx,Ny,Nz)
              if (ax == ax2):
                ax.set_title("nn_delta_e = {}, time = {:.3e}, nstep = {}\nnz = {}".format(nn_delta_e,time,nstep,z))
              else:
                ax.set_title("nz = {}".format(z))
              ax_imgs[nz].set_data(latcut_xy)
            if args.save_plots:
                plt.savefig(data_dir / f"{nstep:05d}.png")
            plt.pause(0.01)

    if args.do_plots:
        ax_set = np.asarray(((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)))
        for nz, ax in enumerate(np.ndarray.flatten(ax_set)):
          z = zvals[nz]
          latcut_xy = latcuts_xy_planes(lattice,nz,Nx,Ny,Nz)
          if (ax == ax2):
            ax.set_title("nn_delta_e = {}, time = {:.3e}, nstep = {}\nnz = {}".format(nn_delta_e,time,nstep,z))
          else:
            ax.set_title("nz = {}".format(z))
            ax_imgs[nz].set_data(latcut_xy)

    av_dt = np.mean(timestep_array)
    print("{} time steps".format(len(timestep_array)))
    print("Average time step = {:e}".format(av_dt))

    end = timer()
    print(f"Elapsed time {end - start:.2f} s")
        
if __name__=="__main__":
    main()

