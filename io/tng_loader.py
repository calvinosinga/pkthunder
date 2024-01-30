import illustris_python as il
import numpy as np
import pickle as pkl
from typing import Union

class TNGLoader(object):

    _def_simnames = [
        'tng100',
        'tng300',
        'Mtng'
    ]

    pkpath = "/Users/cosinga/code/hcolor/fig_md_analysis/1-23_pkdatasort.pkl"
    def __init__(self, 
        path : str,
        simname : str,
        snap : int
    ) -> None:
        self.path = path
        if simname not in self._def_simnames:
            msg = "not an implemented simulation name\n"
            msg += "accepted simnames: " + str(self._def_simnames)
            raise NotImplementedError(msg)
        
        self.simname = simname
        self.snap = snap

        self.loadGal()
        self.loadPk()
        return
    

    def loadGal(self) -> None:
        path = self.path
        snap = self.snap

        sfields = [
            'SubhaloMassType',
            'SubhaloMass',
            'SubhaloGrNr',
            'SubhaloStellarPhotometrics',
            'SubhaloPos'
        ]

        gfields = [
            'GroupFirstSub',
            'GroupPos',
            'Group_M_Mean200',
            'Group_R_Mean200'
        ]
        self.head = il.groupcat.loadHeader(path, snap)
        h = self.head['HubbleParam']
        sdata = il.groupcat.loadSubhalos(path, snap, sfields)
        hdata = il.groupcat.loadHalos(path, snap, gfields)
        ngals = sdata[sfields[0]].shape[0]
        
        # label which galaxies are centrals
        is_cent = np.zeros(ngals, dtype = bool)
        is_cent[hdata['GroupFirstSub']] = True

        # save the host data for each galaxy
        host_mass = np.zeros_like(sdata['SubhaloMass'])
        host_pos = np.zeros((host_mass.shape[0], 3))
        host_rad = np.zeros_like(host_mass)
        host_mass = hdata['Group_M_Mean200'][sdata['SubhaloGrNr']]
        host_pos = hdata['GroupPos'][sdata['SubhaloGrNr'], :]
        host_rad = hdata['Group_R_Mean200'][sdata['SubhaloGrNr']]
        # calculate g - r
        photo = sdata['SubhaloStellarPhotometrics']
        gr = photo[:, 4] - photo[:, 5]
        
        # save arrays
        self.gdata = {
            'stmass' : sdata['SubhaloMassType'][:, 4] * 1e10/h,
            'total_mass' : sdata['SubhaloMass'] * 1e10/h,
            'is_central' : is_cent,
            'host_mass' : host_mass * 1e10/h,
            'gr' : gr,
            'host_x' : host_pos / 1e3,
            'x' : sdata['SubhaloPos'] / 1e3,
            'host_r' : host_rad / 1e3 * self.head['Time'],
            'host_idx' : sdata['SubhaloGrNr']
        }

        self.hdata = {
            'mass' : hdata['Group_M_Mean200'] * 1e10 / h
        }
        return
    
    def loadPk(self) -> None:
        ds = pkl.load(open(self.pkpath, 'rb'))

        ip = {
            'is_auto':True, 'fieldname':'galaxy',
            'axis':0, 'grid_resolution':800, 'simname': self.simname, 
            'gal_res':'diemer', 'snapshot':[self.snap],
            'gal_species':'stmass', 'sim_resolution':'high', 
            'path':'fiducial', 'mas':'CICW', 'post_process':'no key found',
            'space':['real']
        }
        self.pks = {}
        ip['color'] = 'resolved'
        matches = ds.getMatching(ip)
        if len(matches) > 1:
            print("found too many matches...")  
        self.pks['all'] = matches[0].getData()
        
        return
    
    def loadHI() -> None:
        return
    
    # data getters

    def getBox(self):
        # convert to Mpc / h
        return self.head['BoxSize'] / 1e3
    
    def getHaloMass(self):
        return self.hdata['mass'].copy()
        
    def getGalMass(self):
        return self.gdata['stmass'].copy()

    def getZ(self):
        return self.head['Redshift']
    
    def getResolvedGalMask(self):
        if self.simname == 'tng100':

            return self.gdata['stmass'] > 2e8
        else:
            raise NotImplementedError(
                'resolution definitions for other sims not defined'
            )
        
    def _getGRCut(self):
        if self.snap == 99:
            return 0.60
        elif self.snap == 67:
            return 0.55
        elif self.snap == 50:
            return 0.5
        else:
            raise NotImplementedError(
                'No GR Cut defined for snapshot {:d}'.format(self.snap)
            )
    
    def getBlueMask(self):
        grcut = self._getGRCut()
        return self.gdata['gr'] <= grcut
    
    def getCentralMask(self):
        return self.gdata['is_central'].copy()