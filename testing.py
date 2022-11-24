from utility import AugustusUtility

args = dict(haloid='snap', fsps_bands='uvoir', ssp_model='FSPS',
            ssp_table_file='/home/ejones/local/caesar/SSP_Chab_EL.hdf5',
            nproc=16)

snaprange = list(range(67))
snapnamebase = 'snap_m12-5n128_'
snapdir = '/disk01/ejones/sim_runs/gizgrain_production'

AugustusUtility.createCaesarSnapshots(snaprange, snapnamebase, snapdir, new_namebase='m12-5n128_',
                                        caesar_search_options=args)
