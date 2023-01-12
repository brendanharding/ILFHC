import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline as RBS
import matplotlib.pyplot as plt

class MDNInertialMigrationHelper(object):
    """This class helps with parsing/utilising computed inertial
    migration data for a neutrally buoyant spherical particle
    suspended in flow through curved ducts having a
    rectangular cross-section with aspect ratio 2:1 or 4:1.
    
    The class will interpolate the raw data for a desired bend radius of
    the duct. Interpolation of particle size is not currently supported.
    As such one must pick a particle size matching the provided data.
        
    Note: use of this helper class requires the files
    'MDN_data/rect_2x1_MDN_data.npz'
    'MDN_data/rect_4x1_MDN_data.npz'
    to be located relative to the working directory.
    
    Please ensure you cite our JFM paper (TODO: add DOI once provided)
    if you use this code/data. This code is provided under an MIT license 
    (see https://opensource.org/licenses/MIT). However, I would appreciate it
    if you contact me and to let me know if you use this code/data.
    Please also don't hesitate to contact me if you have any questions/queries.
    
    Brendan Harding, 2023."""
    def __init__(self,a,R,K,aspect=2):
        """Initialise the helper class for a given particle radius (a) and bend radius (R).
        Both a and R should be read as being relative to half of the duct height, 
        i.e. read a as 2a/H=alpha and R as 2R/H=1/epsilon (H being the duct height).
        The aspect argument specifies cross-section aspect ratio,
        which can be 2 or 4 and indicates the width as a multiple of the height.
        (Data for additional aspect ratios may be added in the future)
        """
        self._alpha = a # equivalent to alpha
        self._R = R # equivalent to 1/eps
        self._eps = 1./R
        self._K = K
        self._aspect = aspect
        self._load_aspect_data()
        self._update_alpha_data()
        self._update_eps_data()
        self._update_K_data()
        self._setup_interpolants()
        return
    def _load_aspect_data(self):
        """Load data associate with the desired aspect ratio"""
        if self._aspect==2:
            # If the data file is not found, a FileNotFoundError exception will be raised
            self._raw_data = np.load('MDN_data/rect_2x1_MDN_data.npz')
        elif self._aspect==4:
            # If the data file is not found, a FileNotFoundError exception will be raised
            self._raw_data = np.load('MDN_data/rect_4x1_MDN_data.npz')
        else:
            raise ValueError("The requestied aspect is not one of 2 or 4")
        return
    def _update_alpha_data(self):
        """Pre-process the raw data associated with the desired a"""
        alpha_lookup = self._raw_data['alphas']
        ai = np.where(alpha_lookup==self._alpha)[0]
        if len(ai)==0:
            raise ValueError('The desired particle size is not present in the raw data. '+\
                             'Check what is available with the method get_available_a().')
        
        self._eps_values = self._raw_data['epss'][ai]
        self._alpha_data = self._raw_data['alpha_data'+str(ai)]
        Rs = self._alpha_data[0,:,:,0] # assume samples are same for all eps values
        Zs = self._alpha_data[0,:,:,1]
        self._rs = Rs[:,0]
        self._zs = Zs[0,:]
        return
    def _update_eps_data(self):
        """Pre-process the raw data associated with the desired R"""
        if self._R<self._aspect:
            raise ValueError("The requested bend radius R is non-physical. Choose R>{:d}".format(self._aspect)+
                             "(ideally R>>{:d}) where R is interpreted to be relative".format(self._aspect)+
                             "to half the duct height, i.e. read R as 2R/H")
        if self._R<20.0 or self._R>1280.0:
            print("Warning: Results may be innacurate for the given bend radius, \n"+
                  "\t(20<=R<=1280 is best, where R should be read as 2R/H)")
        if self._eps in self._eps_values:
            # use the respective raw data as is
            ei = np.where((self._eps_values==self._eps))[0][0]
            self._eps_data = self._alpha_data[ei]
        else: 
            # interpolate the raw data
            interpolant = interp1d(self._eps_values,self._alpha_data,axis=0,kind='cubic',
                                   bounds_error=False,fill_value="extrapolate",assume_sorted=True)
            self._eps_data = interpolant(self._eps)
        return
    def _update_K_data(self):
        """Updates the forces based on the current value of K"""
        cc-ignore = False
        pre_inverted = True
        K = self._K
        eps = self._eps
        alpha = self._alpha
        sfs = (K*eps)**0.5 # secondary flow scale
        Lda = 2./alpha
        nt = 6
        # Unpack and reshape the data
        Rs = self._eps_data[:,:,0]
        Zs = self._eps_data[:,:,1]
        UyWxWz0_Kpert = self._eps_data[:,:, 2:11].reshape(Rs.shape+(3,3))
        UxUzWy0_Kpert = self._eps_data[:,:,11:20].reshape(Rs.shape+(3,3))
        self._Aa = self._eps_data[:,:,20:29].reshape(Rs.shape+(3,3))
        self._As = self._eps_data[:,:,29:38].reshape(Rs.shape+(3,3))
        bg_iner = self._eps_data[:,:, 38:101].reshape(Rs.shape+(3,21))
        il_coef = self._eps_data[:,:,101:335].reshape(Rs.shape+(3,78))
        # First construct the leading order motion
        UyWxWz0 = UyWxWz0_Kpert[:,:,:,0]\
                  + K*UyWxWz0_Kpert[:,:,:,1]\
                  + K**2*UyWxWz0_Kpert[:,:,:,2]
        UxUzWy0 = UxUzWy0_Kpert[:,:,:,0]\
                  + K*UxUzWy0_Kpert[:,:,:,1]\
                  + K**2*UxUzWy0_Kpert[:,:,:,2]
        UxUzWy0 *= sfs # appropriate scale for the secondary flow
        Theta0 = UyWxWz0[:,:,0]*(eps*alpha/(1.0+eps*alpha*Rs))
        # Account for centrifugal/centripetal forces
        if cc_ignore:
            UW1 = np.zeros(Rs.shape+(6,))
        else:
            UW1 = np.empty((SX.shape[0],SX.shape[1],6))
            UW1[:,:,1] = +2*(rho_p*4.0/ 3.0*np.pi)*UxUzWy0[:,:,0]*Theta0
            UW1[:,:,3] = -  (rho_p*8.0/15.0*np.pi)*UxUzWy0[:,:,2]*Theta0
            UW1[:,:,0] = -(rho_p*4.0/ 3.0*np.pi)*Theta0*UyWxWz0[:,:,0]
            UW1[:,:,4] = +(rho_p*8.0/15.0*np.pi)*UyWxWz0[:,:,1]*Theta0
            UW1[:,:,[2,5]] = 0
            if pre_inverted:
                UW1[:,:,1::2] = np.linalg.solve(-Aa,UW1[:,:,1::2])
                UW1[:,:,0::2] = np.linalg.solve(-As,UW1[:,:,0::2])
            ind = 0
            for i in range(0,nt):
                for j in range(i,nt):
                    if (i%2)==(j%2):
                        factor = -K**((i//2)+(j//2))
                        if i%2==1: # i,j are either both even, or both odd, in the subsequent loop
                            factor *= K*eps # note K*eps=sfs**2
                        UW1[:,:,0::2] += factor*bg_iner[:,:,:,ind]
                    else:
                        factor = -sfs*K**((i//2)+(j//2))
                        UW1[:,:,1::2] += factor*bg_iner[:,:,:,ind]
                    ind += 1
        # Now re-construct the inertial lift contribution
        dat = [UxUzWy0[:,:,0],UyWxWz0[:,:,0],UxUzWy0[:,:,1],\
               UyWxWz0[:,:,1],UxUzWy0[:,:,2],UyWxWz0[:,:,2]]
        for i in range(nt):
            dat += [Lda*K**(i//2) if i%2==0 else Lda*sfs*K**(i//2)]
        ind = 0
        for i in range(len(dat)):
            for j in range(i,len(dat)):
                if i<6:
                    #UW1 += (dat[i]*dat[j])[:,:,np.newaxis]*il_coef[:,:,:,ind]
                    if (j<6 and (i%2)==(j%2)) or (j>=6 and (i%2)!=(j%2)):
                        UW1[:,:,0::2] += (dat[i]*dat[j])[:,:,np.newaxis]*il_coef[:,:,:,ind]
                    else:
                        UW1[:,:,1::2] += (dat[i]*dat[j])[:,:,np.newaxis]*il_coef[:,:,:,ind]
                else:
                    #UW1 += (dat[i]*dat[j])*il_coef[:,:,:,ind]
                    if (j<6 and (i%2)==(j%2)) or (j>=6 and (i%2)!=(j%2)):
                        UW1[:,:,0::2] += (dat[i]*dat[j])*il_coef[:,:,:,ind]
                    else:
                        UW1[:,:,1::2] += (dat[i]*dat[j])*il_coef[:,:,:,ind]
                ind +=1
        if pre_inverted:
            UyWxWz1 = UW1[:,:,1::2]
            UxUzWy1 = UW1[:,:,0::2]
        else:
            UyWxWz1 = np.linalg.solve(-Aa,UW1[:,:,1::2])
            UxUzWy1 = np.linalg.solve(-As,UW1[:,:,0::2])
        # Lastly, put the data into arrays to be called upon...
        self._UxUzWy0 = UxUzWy0
        self._UyWxWz0 = UyWxWz0
        #self._Theta0 = Theta0
        self._UyWxWz1 = UyWxWz1
        self._UxUzWy1 = UxUzWy1
        return
    def _setup_interpolants(self,include_spin=False):
        """Constructs interpolants of the appropriate fields"""
        nhH,nhW = 1.0/self._a-1.0,self._aspect/self._a-1.0
        bbox = [-nhW,nhW,-nhH,nhH]
        rs,zs = self._rs,self._zs
        self._Ux0_RBS = RBS(rs,zs,self._UxUzWy0[:,:,0], bbox=bbox)
        self._Uy0_RBS = RBS(rs,zs,self._UyWxWz0[:,:,0], bbox=bbox)
        self._Uz0_RBS = RBS(rs,zs,self._UxUzWy0[:,:,1], bbox=bbox)
        if include_spin:
            self._Wx0_RBS = RBS(rs,zs,self._UyWxWz0[:,:,1], bbox=bbox)
            self._Wy0_RBS = RBS(rs,zs,self._UxUzWy0[:,:,2], bbox=bbox)
            self._Wz0_RBS = RBS(rs,zs,self._UyWxWz0[:,:,2], bbox=bbox)
        self._Ux1_RBS = RBS(rs,zs,self._UxUzWy1[:,:,0], bbox=bbox)
        self._Uy1_RBS = RBS(rs,zs,self._UyWxWz1[:,:,0], bbox=bbox)
        self._Uz1_RBS = RBS(rs,zs,self._UxUzWy1[:,:,1], bbox=bbox)
        if include_spin:
            self._Wx1_RBS = RBS(rs,zs,self._UyWxWz1[:,:,1], bbox=bbox)
            self._Wy1_RBS = RBS(rs,zs,self._UxUzWy1[:,:,2], bbox=bbox)
            self._Wz1_RBS = RBS(rs,zs,self._UyWxWz1[:,:,2], bbox=bbox)
        Re = (self._K/self._eps)**0.5
        Rep = 0.5*Re*self._alpha**2
        self._Ux_RBS = RBS(rs,zs,self._UxUzWy0[:,:,0]+Rep*self._UxUzWy1[:,:,0], bbox=bbox)
        self._Uy_RBS = RBS(rs,zs,self._UyWxWz0[:,:,0]+Rep*self._UyWxWz1[:,:,0], bbox=bbox)
        self._Uz_RBS = RBS(rs,zs,self._UxUzWy0[:,:,1]+Rep*self._UxUzWy1[:,:,1], bbox=bbox)
        if include_spin:
            self._Wx_RBS = RBS(rs,zs,self._UyWxWz0[:,:,1]+Rep*self._UyWxWz1[:,:,1], bbox=bbox)
            self._Wy_RBS = RBS(rs,zs,self._UxUzWy0[:,:,2]+Rep*self._UxUzWy1[:,:,2], bbox=bbox)
            self._Wz_RBS = RBS(rs,zs,self._UyWxWz0[:,:,2]+Rep*self._UyWxWz1[:,:,2], bbox=bbox)
        #self._Cr_RBS = RBS(rs,zs,self._aR_data[2,:,:], bbox=bbox)
        #self._Cz_RBS = RBS(rs,zs,self._aR_data[3,:,:], bbox=bbox)
        #self._Lr_RBS = RBS(rs,zs,self._aR_data[4,:,:], bbox=bbox)
        #self._Lz_RBS = RBS(rs,zs,self._aR_data[5,:,:], bbox=bbox)
        #self._Sr_RBS = RBS(rs,zs,self._aR_data[6,:,:], bbox=bbox)
        #self._Sz_RBS = RBS(rs,zs,self._aR_data[7,:,:], bbox=bbox)
        #self._Up_RBS = RBS(rs,zs,self._aR_data[8,:,:], bbox=bbox)
        #self._Wr_RBS = RBS(rs,zs,self._aR_data[9,:,:], bbox=bbox)
        #self._Wz_RBS = RBS(rs,zs,self._aR_data[10,:,:],bbox=bbox)
        #self._kappa = 4.0/(self._R*self._a**3)
        #self._Fr_RBS = RBS(rs,zs,self._aR_data[4,:,:]+self._kappa*self._aR_data[6,:,:],bbox=bbox)
        #self._Fz_RBS = RBS(rs,zs,self._aR_data[5,:,:]+self._kappa*self._aR_data[7,:,:],bbox=bbox)
        return
    def _plot_aR_component(self,ind):
        """Useful for testing/debugging purposes"""
        # TODO: requires update
        RS,ZS = np.meshgrid(self._rs,self._zs,indexing='ij')
        plt.figure(figsize=(1+3*self._aspect,3))
        plt.contourf(RS,ZS,self._aR_data[ind,:,:],17)
        plt.colorbar()
        nhH,nhW = 1.0/self._a,self._aspect/self._a
        plt.plot([-nhW+1,nhW-1,nhW-1,-nhW+1,-nhW+1],[-nhH+1,-nhH+1,nhH-1,nhH-1,-nhH+1],'r--')
        plt.xlim(-nhW,nhW)
        plt.ylim(-nhH,nhH)
        plt.gca().set_aspect(1.0)
        plt.show()
    def get_asepct_ratio(self):
        """Get the current cross-section aspect ratio"""
        return self._aspect
    def get_bounds(self):
        """Get the bounds of the current (dimensionless) cross-section
        (given in the form [r_min,r_max,z_min,z_max])"""
        return [-self._aspect/self._a,self._aspect/self._a,-1.0/self._a,1.0/self._a]
    def get_particle_bounds(self):
        """Get the bounds for the particle centre in the current (dimensionless)
        cross-section (given in the form [r_min,r_max,z_min,z_max])"""
        return [-self._aspect/self._alpha+1.0,self._aspect/self._alpha-1.0,\
                -1.0/self._alpha+1.0,1.0/self._alpha-1.0]
    def set_aspect_ratio(self,aspect):
        """Change the cross-section (and update the interpolants)"""
        self._aspect = aspect
        self._load_aspect_data()
        self._update_alpha_data()
        self._update_eps_data()
        self._update_K_data()
        self._setup_interpolants()
        return
    def get_available_alpha(self):
        """Get a list of available particle radii"""
        return self._raw_data['alphas']
    def get_alpha(self):
        """Get the current particle radius"""
        return self._alpha
    def set_alpha(self,a):
        """Change the particle radius (and update the interpolants)"""
        self._alpha = a
        self._update_alpha_data()
        self._update_eps_data()
        self._update_K_data()
        self._setup_interpolants()
        return
    def get_R(self):
        """Get the current bend radius"""
        return self._R
    def set_R(self,R):
        """Change the bend radius of the duct (and update the interpolants)"""
        self._R = R
        self._eps = 1./self._R
        self._update_eps_data()
        self._update_K_data()
        self._setup_interpolants()
    def get_K(self):
        """Get the current Dean number K=eps*Re^2"""
        return self._K
    def get_Re(self):
        """Get the current (channel) Reynolds number Re=(rho/mu)*U*(H/2)
        (with U the maximum axial velocity and H the duct height)"""
        return (self._K/self._eps)**0.5
    def get_Rep(self):
        """Get the current particle Reynolds number Re_p=Re*alpha^2
        (with U the maximum axial velocity and H the duct height)"""
        return (0.5*self._alpha)**2*(self._K/self._eps)**0.5
    def set_K(self,K):
        """Change the dean number"""
        self._K = K
        self._update_K_data()
        self._setup_interpolants()
    def migration_velocity(self,r,z):
        """Get the migration velocity for a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via U_m a / H, not by
        ( rho / mu ) U_m^2 a^3 / H^2 as was the case previously)"""
        return np.squeeze([self._Ur_RBS(r,z),
                           self._Uz_RBS(r,z)])
    def axial_velocity(self,r,z):
        """Get the terminal/steady axial velocity of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via U_m a / H )"""
        return np.squeeze(self._Uy_RBS(r,z))
    def spin_components(self,r,z):
        """Get the terminal/steady r,z spin components of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via U_m / H )"""
        return np.squeeze([self._Wr_RBS(r,z),self._Wz_RBS(r,z)])
    # TODO: Everything below needs checking/updating
    def migration_force(self,r,z):
        raise NotImplementedError # need to multiply Aa,As by the velocity fields
        """Get the (net) migration force for a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via rho U_m^2 a^4 / H^2 )"""
        return np.squeeze([self._Fr_RBS(r,z),self._Fz_RBS(r,z)])
    def migration_force_jacobian(self,r,z):
        raise NotImplementedError # need to multiply Aa,As by the velocity fields
        """Get the jacobian of the (net) migration force for a neutrally 
        buoyant spherical particle centred at (r,z) within the 
        cross-section (non-dimensionalised via rho U_m^2 a^3 / H^2 )"""
        return np.squeeze([[self._Fr_RBS(r,z,dx=1),self._Fr_RBS(r,z,dy=1)],
                           [self._Fz_RBS(r,z,dx=1),self._Fz_RBS(r,z,dy=1)]])
    def drag_coefficient(self,r,z):
        raise NotImplementedError # need to multiply Aa,As by the velocity fields
        """Get the drag coefficients in the r,z directions of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via mu a )"""
        return np.squeeze([self._Cr_RBS(r,z),self._Cz_RBS(r,z)])
    def secondary_flow_drag(self,r,z):
        raise NotImplementedError # need to multiply Aa,As by the velocity fields
        """Get the drag coefficients in the r,z directions of a neutrally buoyant
        spherical particle centred at (r,z) within the cross-section
        (non-dimensionalised via rho U_m^2 a H^2 / ( 4 R ) )"""
        return np.squeeze([self._Sr_RBS(r,z),self._Sz_RBS(r,z)])
    def plot_migration_force(self):
        raise NotImplementedError # need to multiply Aa,As by the velocity fields
        """Produces a rough sketch of the magnitude of the migration force including
        the zero contours of the horizontal and vertical components."""
        Fr = self._aR_data[4,:,:]+self._kappa*self._aR_data[6,:,:]
        Fz = self._aR_data[5,:,:]+self._kappa*self._aR_data[7,:,:]
        RS,ZS = np.meshgrid(self._rs,self._zs,indexing='ij')
        plt.figure(figsize=(1+3*self._aspect,3))
        plt.contourf(RS,ZS,(Fr**2+Fz**2)**0.5,17)
        plt.colorbar()
        plt.contour(RS,ZS,Fr,[0.0],colors=['k'])
        plt.contour(RS,ZS,Fz,[0.0],colors=['w'])
        nhH,nhW = 1.0/self._a,self._aspect/self._a
        plt.plot([-nhW+1,nhW-1,nhW-1,-nhW+1,-nhW+1],[-nhH+1,-nhH+1,nhH-1,nhH-1,-nhH+1],'r--')
        plt.xlim(-nhW,nhW)
        plt.ylim(-nhH,nhH)
        plt.gca().set_aspect(1.0)
        plt.show()
        return
    # end of class
