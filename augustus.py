import yt
import caesar as cs
import numpy as np
from utility import AugustusUtility

# This class holds an array of augustus objects indexable by redshift or snapshot number
class Augustus:
    
    def __init__(self, snap_list, sim_namebase, sim_dirpath,
                 caesar_namebase, caesar_dirpath, snapsAsRedshift=False,
                 boxspaceFile=None):
        # Inputs
        #   + snap_list --> list of snapshot numbers
        #   + sim_namebase --> naming convention of simulation snapshots (without numbers)
        #   + sim_dirpath  --> path to directory containing all simulation snapshots
        #   + caesar_namebase --> naming convention of caesar snapshots (without numbers)
        #   + caesar_dirpath  --> path to directory containing all caesar snapshots
        #   + snapsAsRedshift --> if True then snaplist is interpreted as a list of redshifts
        #                         and the closest snapshots will be loaded, must also provide
        #                         boxspaceFile
        #   + boxspaceFile    --> path to boxspace file output from simulation. This contains
        #                         two columns, the first being snapshot number, the second
        #                         being the cosmological scale factor at this snapshot
        
        # Input variables
        if snapsAsRedshift:
            self.userRedshifts = snap_list
        else:
            self.snap_list = snap_list
        self.sim_namebase = sim_namebase
        self.sim_dirpath = sim_dirpath
        self.caesar_namebase = caesar_namebase
        self.caesar_dirpath = caesar_dirpath
        if snapsAsRedshift:
            if boxspaceFile == None:
                print("A boxspace file must be supplied when passing redshifts instead of snapshot numbers!")
                exit(0)
            self.snap_list = AugustusUtility.find_snapshots_from_redshift(self.userRedshifts, boxspaceFile)
            
        # Ensure dirpaths are correctly formatted
        if self.sim_dirpath[-1] != '/': self.sim_dirpath += '/'
        if self.caesar_dirpath[-1] != '/': self.caesar_dirpath += '/'
            
        # Load and store augustus objects
        self.data  = {}
        self.redshifts = []
        self._load_augustus()

    # Method to create and store augustus objects for each snapshot number in snap_list
    def _load_augustus(self):
        
        for snapnum in self.snap_list:
            snapnum = str(snapnum).zfill(3)
            simfile = f"{self.sim_dirpath}{self.sim_namebase}{snapnum}.hdf5"
            caesarfile = f"{self.caesar_dirpath}{self.caesar_namebase}{snapnum}.hdf5"
            
            augObj = AugustusObject()
            augObj.load_sim_snapshot(simfile)
            augObj.load_caesar_snapshot(caesarfile)
            redshift = augObj.cosmo_params['redshift']
            self.data[redshift] = augObj
            self.redshifts.append(redshift)
    
    # Method to find galaxies identified by Caesar at each redshift and store their particle data
    def find_galaxies(self, redshift):
        for z in self.redshifts:
            # Call object method for each redshift
            self.data[z].identify_galaxies()
            
    # Method to find halos identified by Caesar at each redshift and store their particle data
    def find_halos(self, redshift):
        for z in self.redshifts:
            # Call object method for each redshift
            self.data[z].identify_halos()
            
    

    
# This class is each individual augustus object
class AugustusObject:
    
    def __init__(self, yt_obj=None, caesar_obj=None):
        # Inputs
        #   + Can supply pre-loaded simulation yt and caesar objects with yt_obj and caesar_obj
        
        self.yt_obj = yt_obj
        self.caesar_obj = caesar_obj
        self.cosmo_params = {}
        
        # >>> These are populated by the identify_galaxies method >>>
        # List of galaxy IDs to index the dictionaries at
        self.galaxyIDs = []
        # List of all caesar galaxy objects
        self.galaxies  = None
        # Dictionaries are indexed by galaxy ID
        self.gasInGalaxies        = {}
        self.starsInGalaxies      = {}
        self.darkMatterInGalaxies = {}
        self.blackHolesInGalaxies = {}
        # <<< These are populated by the identify_galaxies method <<<
        
        # >>> These are populated by the identify_halos method >>>
        # List of galaxy IDs to index the dictionaries at
        self.haloIDs = []
        # List of all caesar galaxy objects
        self.halos  = None
        # Dictionaries are indexed by galaxy ID
        self.gasInHalos        = {}
        self.starsInHalos      = {}
        self.darkMatterInHalos = {}
        self.blackHolesInHalos = {}
        # <<< These are populated by the identify_halos method <<<
        
    # Method to load snapshot file using yt
    def load_sim_snapshot(self,sim_path):
        # Inputs
        #   + sim_path --> path to simulation snapshot to load
        # Returns
        #   + None. Saves yt object internally.
        
        self.yt_obj = yt.load(sim_path)
        self.get_cosmological_parameters()
        
    # Method to load caesar file using yt
    def load_caesar_snapshot(self,caesar_path):
        # Inputs
        #   + caesar_path --> path to caesar file to load
        # Returns
        #   + None. Saves caesar object internally.
        
        self.caesar_obj = cs.load(caesar_path)
        
    # Method to get various cosmological parameters from the yt object
    def get_cosmological_parameters(self):
        self.cosmo_params['redshift'] = self.yt_obj.current_redshift
        self.cosmo_params['omega_lambda'] = self.yt_obj.omega_lambda
        self.cosmo_params['omega_matter'] = self.yt_obj.omega_matter
        self.cosmo_params['hubble_constant'] = self.yt_obj.hubble_constant
        
    # Method to find all galaxies and store the corresponding particle data for each
    def identify_galaxies(self):
        self.galaxies = self.caesar_obj.galaxies
        # Loop over each galaxy
        for gal in self.galaxies:
            # Galaxy ID to index dictionaries at
            galID = gal.GroupID
            self.galaxyIDs.append(galID)
            # Get particle lists for each galaxy
            self.gasInGalaxies[galID]        = gal.glist
            self.starsInGalaxies[galID]      = gal.slist
            self.darkMatterInGalaxies[galID] = gal.dmlist
            self.blackHolesInGalaxies[galID] = gal.bhlist
            
    # Method to find all halos and store the corresponding particle data for each
    def identify_halos(self):
        self.halos = self.caesar_obj.halos
        # Loop over each halo
        for halo in self.halos:
            # Halo ID to index dictionaries at
            haloID = halo.GroupID
            self.haloIDs.append(haloID)
            # Get particle lists for each halo
            self.gasInHalos[haloID]        = halo.glist
            self.starsInHalos[haloID]      = halo.slist
            self.darkMatterInHalos[haloID] = halo.dmlist
            self.blackHolesInHalos[haloID] = halo.bhlist