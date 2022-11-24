import yt
import caesar as cs
import numpy as np
import os, sys

# Class containing utility methods to be used internally and externally with augustus
class AugustusUtility:
    
    # Method to find snapshot numbers at given redshifts
    @staticmethod
    def convert_redshift_to_snapnum(redshifts, boxspaceFile):
        # INPUTS
        #   + redshits --> list of redshifts to find closest snapshot to
        #   + boxspaceFile    --> path to boxspace file output from simulation. This contains
        #                         two columns, the first being snapshot number, the second
        #                         being the cosmological scale factor at this snapshot
        # OUTPUTS
        #   + list of snapshot numbers corresponding to closest snapshots to requested redshifts
            # Load each row of the file
            with open(boxspaceFile, 'r') as fptr:
                lines = fptr.readlines()
                # List of (snapshot,redshift) calculated from boxspace file
                boxspaceData = np.zeros((len(lines), 2))
                for snap, line in enumerate(lines):
                    # Get snapshot number and scale factor
                    a = float(line.strip())
                    z = 1/float(a) - 1
                    # Calculate the redshift
                    boxspaceData[snap,:] = snap, z
            # Find closest redshift in boxspaceData to supplied redshift
            snap_list=[]
            for z_user in redshifts:
                ind = np.argmin(np.abs(boxspaceData[:,1] - z_user))
                snapnum = boxspaceData[ind,0]
                snap_list.append(int(snapnum))
            return snap_list
        
    # Method to convert simulation snapshots into caesar snapshots by running a standard
    # Caesar script to identify galaxies and halos. These caesar snapshots are necessary
    # for running Augustus.
    @staticmethod
    def createCaesarSnapshots(snap_range, snap_namebase, snap_dir, new_namebase=None,
                              caesar_dir=None, fof6D_dir=None, caesar_search_options=None):
        # INPUTS
        #   + snap_range --> list of snapshot numbers
        #   + snap_namebase --> string of the snapshot name excluding snapshot number
        #   + snap_dir --> directory where snapshots are located
        #   + new_namebase --> namebase to use when writing Caesar and fof6D files
        #   + caesar_dir --> directory where caesar files will be saved
        #   + fof6D_dir --> directory where fof6D files will be saved
        #   + caesar_search_options --> dictionary containing keyword arguments to be passed to
        #                               caesar's member_search function
        # OUTPUTS
        #   + None. Will save files to directories specified.
        
        # Directories
        if caesar_dir == None:
            caesar_dir = str(os.path.join(snap_dir, 'caesar_files'))
        if fof6D_dir == None:
            fof6D_dir  = str(os.path.join(caesar_dir, 'fof6D'))
        # Files
        if new_namebase == None:
            new_namebase = snap_namebase
        caesar_namebase = f'caesar_{new_namebase}_'
        fof6D_namebase  = f'fof6D_{new_namebase}_'
        # caesar.member_search options
        if caesar_search_options == None:
            caesar_search_options = dict(haloid='snap', fsps_bands='uvoir', ssp_model='FSPS')

        # Loop over snapshots
        for snap_num in snap_range:
            # Check for directories and create them if needed
            if not os.path.isdir(caesar_dir):
                os.makedirs(caesar_dir)
            if not os.path.isdir(fof6D_dir):
                os.makedirs(fof6D_dir)

            # Create names of all files
            snum = str(snap_num).zfill(3)
            snap_name   = f'{snap_dir}/{snap_namebase}{snum}.hdf5'
            caesar_name = f'{caesar_dir}/{caesar_namebase}{snum}.hdf5'
            fof6D_name  = f'{fof6D_dir}/{fof6D_namebase}{snum}.hdf5'
            
            # Add fof6D name to member search options
            caesar_search_options['fof6d_file'] = fof6D_name

            # Load snapshot into caesar object
            snap_yt = yt.load(snap_name)
            snap_caesar = cs.CAESAR(snap_yt)

            # Execute Caesar's member_search method to identify astrophysical objects
            snap_caesar.member_search(**caesar_search_options)

            # Save Caesar catalogue
            snap_caesar.save(caesar_name)