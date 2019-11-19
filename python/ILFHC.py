import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline as RBS
import matplotlib.pyplot as plt

class InertialLiftForceHelper(object):
    """This class helps with parsing/utilising computed inertial lift  
    force and migration data for a neutrally buoyant spherical  
    particle suspended in flow through curved ducts having a  
    square or rectangular (with aspect ratio 2:1 or 4:1) cross-section.
    
    The class will interpolate the raw data for a desired bend radius of
    the duct. Interpolation of particle size is not currently supported.
    As such one must pick a particle size matching the provided data.
        
    Note: use of this helper class requires the files 
    'curved_duct_lift_data/square_lift_data.npz'  
    'curved_duct_lift_data/rect_2x1_lift_data.npz' 
    'curved_duct_lift_data/rect_4x1_lift_data.npz' 
    to be located relative to the working directory.
    
    Please ensure you cite our JFM paper (https://doi.org/10.1017/jfm.2019.323) 
    if you use this code/data. This code is provided under an MIT license 
    (see https://opensource.org/licenses/MIT). However, I would appreciate it
    if you contact me and to let me know if you use this code/data.
    Please also don't hesitate to contact me if you have any questions/queries.
    
    Brendan Harding, 2019."""
    def __init__(self,a,R,cs='square'):
        """Initialise the helper class for a given particle radius (a) and bend radius (R).
        Both a and R should be read as being relative to half of the duct height, 
        i.e. read a as 2a/H and R as 2R/H (H being the duct height).
        The cs argument specifies cross-section shape: 'square' for a square 
        cross-section, or 'rect_2x1' for a rectangular duct (having aspect ratio 2),
        or 'rect_4x1' for a rectangular duct (having aspect ratio 4).
        (Data for additional cross-sections may be added in the future)
        """
        self._a = a
        self._R = R
        self._cs = cs
        self._load_cs_data()
        self._update_a_data()
        self._update_R_data()
        self._setup_interpolants()
        return
    def _load_cs_data(self):
        """Load data associate with the desired cross-section"""
        if self._cs=='square':
            # If the  data file is not found, a FileNotFoundError exception will be raised
            self._raw_data = np.load('curved_duct_lift_data/square_lift_data.npz')
            self._ar = 1.0
        elif self._cs=='rect_2x1':
            # If the data file is not found, a FileNotFoundError exception will be raised
            self._raw_data = np.load('curved_duct_lift_data/rect_2x1_lift_data.npz')
            self._ar = 2.0
        elif self._cs=='rect_4x1':
            # If the data file is not found, a FileNotFoundError exception will be raised
            self._raw_data = np.load('curved_duct_lift_data/rect_4x1_lift_data.npz')
            self._ar = 4.0
        else:
            raise ValueError("The requestied cross-section is not one of 'square' or 'rect_2x1' or 'rect_4x1'")
        return
    def _update_a_data(self):
        """Pre-process the raw data associated with the desired a"""
        aR_lookup = self._raw_data['aR']
        ais = np.where(aR_lookup[:,0]==self._a)[0]
        if len(ais)==0:
            raise ValueError('The desired particle size is not present in the raw data. '+\
                             'Check what is available with the method get_available_a().')
        self._R_values = aR_lookup[ais,1]
        a_data = np.array([self._raw_data['data'][index] for index in ais])
        # strip NaN data around the edges (there's possibly a cleaner way)
        r = a_data[0,0,:,20]
        z = a_data[0,1,20,:]
        rf = np.where(np.isfinite(r))[0]
        zf = np.where(np.isfinite(z))[0]
        self._rs = r[rf[0]:rf[-1]+1]
        self._zs = z[zf[0]:zf[-1]+1]
        self._a_data = a_data[:,:,rf[0]:rf[-1]+1,zf[0]:zf[-1]+1]
        return
    def _update_R_data(self):
        """Pre-process the raw data associated with the desired R"""
        if self._R<self._ar:
            raise ValueError("The requested bend radius R is non-physical. Choose R>{:d}".format(self.ar_)+ 
                             "(ideally R>>{:d}) where R is interpreted to be relative".format(self.ar_)+
                             "to half the duct height, i.e. read R as 2R/H")
        if self._R<20.0 or self._R>1280.0:
            print("Warning: Results may be innacurate for the given bend radius, \n"+
                  "\t(20<=R<=1280 is best, where R should be read as 2R/H)")
        if self._R in self._R_values:
            # use the respective raw data as is
            index = np.where((self._R_values==self._R))[0][0]
            self._aR_data = self._a_data[index]
        else: 
            # construct an interpolant
            Rs = np.array(self._R_values)
            XS = self._a_data[0][0]
            ZS = self._a_data[0][1]
            aR_data = [XS,ZS]
            for fi in range(2,11):
                F_array = np.array([self._a_data[ai][fi] for ai in range(len(self._a_data))])
                F_interp = interp1d(1.0/Rs[::-1],F_array[::-1,:,:],axis=0,kind='cubic',
                                    bounds_error=False,fill_value="extrapolate",assume_sorted=True)
                aR_data.append(F_interp(1.0/self._R))
            self._aR_data = np.array(aR_data)
        return
    def _setup_interpolants(self):
        """Constructs interpolants of the appropriate fields"""
        nhH,nhW = 1.0/self._a-1.0,self._ar/self._a-1.0
        bbox = [-nhW,nhW,-nhH,nhH]
        self._Cr_RBS = RBS(self._rs,self._zs,self._aR_data[2,:,:], bbox=bbox)
        self._Cz_RBS = RBS(self._rs,self._zs,self._aR_data[3,:,:], bbox=bbox)
        self._Lr_RBS = RBS(self._rs,self._zs,self._aR_data[4,:,:], bbox=bbox)
        self._Lz_RBS = RBS(self._rs,self._zs,self._aR_data[5,:,:], bbox=bbox)
        self._Sr_RBS = RBS(self._rs,self._zs,self._aR_data[6,:,:], bbox=bbox)
        self._Sz_RBS = RBS(self._rs,self._zs,self._aR_data[7,:,:], bbox=bbox)
        self._Up_RBS = RBS(self._rs,self._zs,self._aR_data[8,:,:], bbox=bbox)
        self._Wr_RBS = RBS(self._rs,self._zs,self._aR_data[9,:,:], bbox=bbox)
        self._Wz_RBS = RBS(self._rs,self._zs,self._aR_data[10,:,:],bbox=bbox)
        self._kappa = 4.0/(self._R*self._a**3)
        self._Fr_RBS = RBS(self._rs,self._zs,self._aR_data[4,:,:]+self._kappa*self._aR_data[6,:,:],bbox=bbox)
        self._Fz_RBS = RBS(self._rs,self._zs,self._aR_data[5,:,:]+self._kappa*self._aR_data[7,:,:],bbox=bbox)
        return
    def _plot_aR_component(self,ind):
        """Useful for testing/debugging purposes"""
        RS,ZS = np.meshgrid(self._rs,self._zs,indexing='ij')
        plt.figure(figsize=(1+3*self._ar,3))
        plt.contourf(RS,ZS,self._aR_data[ind,:,:],17)
        plt.colorbar()
        nhH,nhW = 1.0/self._a,self._ar/self._a
        plt.plot([-nhW+1,nhW-1,nhW-1,-nhW+1,-nhW+1],[-nhH+1,-nhH+1,nhH-1,nhH-1,-nhH+1],'r--')
        plt.xlim(-nhW,nhW)
        plt.ylim(-nhH,nhH)
        plt.gca().set_aspect(1.0)
        plt.show()
    def get_cross_section(self):
        """Get the current cross-section"""
        return self._cs
    def get_bounds(self):
        """Get the bounds of the current cross-section 
        (given in the form [r_min,r_max,z_min,z_max])"""
        return [-self._ar/self._a,self._ar/self._a,-1.0/self._a,1.0/self._a]
    def get_particle_bounds(self):
        """Get the bounds for the particle centre in the current 
        cross-section (given in the form [r_min,r_max,z_min,z_max])"""
        return [-self._ar/self._a+1.0,self._ar/self._a-1.0,-1.0/self._a+1.0,1.0/self._a-1.0]
    def set_cross_section(self,cs):
        """Change the cross-section (and update the interpolants)"""
        self._cs = cs
        self._load_cs_data()
        self._update_a_data()
        self._update_R_data()
        self._setup_interpolants()
        return
    def get_available_a(self):
        """Get a list of available particle radii"""
        return np.unique(self._raw_data['aR'][:,0])
    def get_a(self):
        """Get the current particle radius"""
        return self._a
    def set_a(self,a):
        """Change the particle radius (and update the interpolants)"""
        self._a = a
        self._update_a_data()
        self._update_R_data()
        self._setup_interpolants()
        return
    def get_R(self):
        """Get the current bend radius"""
        return self._R
    def set_R(self,R):
        """Change the bend radius of the duct (and update the interpolants)"""
        self._R = R
        self._update_R_data()
        self._setup_interpolants()
    def migration_force(self,r,z):
        """Get the (net) migration force for a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via rho U_m^2 a^4 / H^2 )"""
        return np.squeeze([self._Fr_RBS(r,z),self._Fz_RBS(r,z)])
    def migration_velocity(self,r,z):
        """Get the migration velocity for a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via U_m a / H )"""
        return np.squeeze([self._Fr_RBS(r,z)/self._Cr_RBS(r,z),
                           self._Fz_RBS(r,z)/self._Cz_RBS(r,z)])
    def drag_coefficient(self,r,z):
        """Get the drag coefficients in the r,z directions of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via mu a )"""
        return np.squeeze([self._Cr_RBS(r,z),self._Cz_RBS(r,z)])
    def secondary_flow_drag(self,r,z):
        """Get the drag coefficients in the r,z directions of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via rho U_m^2 a^4 / H^2 )"""
        return np.squeeze([self._Sr_RBS(r,z),self._Sz_RBS(r,z)])
    def axial_velocity(self,r,z):
        """Get the terminal/steady axial velocity of a neutrally buoyant 
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via U_m a / H )"""
        return np.squeeze(self._Up_RBS(r,z))
    def spin_components(self,r,z):
        """Get the terminal/steady r,z spin components of a neutrally buoyant 
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via U_m / H )"""
        return np.squeeze([self._Wr_RBS(r,z),self._Wz_RBS(r,z)])
    def plot_migration_force(self):
        """Produces a rough sketch of the magnitude of the migration force including
        the zero contours of the horizontal and vertical components."""
        Fr = self._aR_data[4,:,:]+self._kappa*self._aR_data[6,:,:]
        Fz = self._aR_data[5,:,:]+self._kappa*self._aR_data[7,:,:]
        RS,ZS = np.meshgrid(self._rs,self._zs,indexing='ij')
        plt.figure(figsize=(1+3*self._ar,3))
        plt.contourf(RS,ZS,(Fr**2+Fz**2)**0.5,17)
        plt.colorbar()
        plt.contour(RS,ZS,Fr,[0.0],colors=['k'])
        plt.contour(RS,ZS,Fz,[0.0],colors=['w'])
        nhH,nhW = 1.0/self._a,self._ar/self._a
        plt.plot([-nhW+1,nhW-1,nhW-1,-nhW+1,-nhW+1],[-nhH+1,-nhH+1,nhH-1,nhH-1,-nhH+1],'r--')
        plt.xlim(-nhW,nhW)
        plt.ylim(-nhH,nhH)
        plt.gca().set_aspect(1.0)
        plt.show()
        return
    # end of class
