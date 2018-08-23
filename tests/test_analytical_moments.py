import galsim
import unittest
import numpy as np
import os, sys
realizer_path = os.path.join(os.environ['SLREALIZERDIR'], 'slrealizer')
sys.path.insert(0, realizer_path)
from utils.utils import *

class AnalyticalTest(unittest.TestCase):  
    """Tests the analytical equations used for calculating moments."""
    
    @classmethod
    def setUpClass(cls):
        truth_params = {
            'gal_sigma': 4.1,
            'gal_flux' : 1.5,
            'gal_e1' : 0.3,
            'gal_e2' : 0.4,
            'gal_x' : -2.0,
            'gal_y' : 2.4,
            'qso_sigma' : 5.3, 
            'qso_flux' : 0.5,
            'qso_x' : -1.0,
            'qso_y' : 1.0,
            'psf_sigma' : 0.4,
            'nx' : 257,
            'pixel_scale' : 0.2,
        }

        for k, v in truth_params.items():
            setattr(cls, k, v)

    def test_gaussian_galaxy(self):
        """ Compares sigma and flux of galaxy """
        gal = galsim.Gaussian(sigma=self.gal_sigma, flux=self.gal_flux)
        galsim_img = gal.drawImage(scale=self.pixel_scale, nx=self.nx, ny=self.nx, method='no_pixel')

        hsm_output = galsim_img.FindAdaptiveMom(guess_sig=self.gal_sigma/self.pixel_scale)
        num_sigma_from_det_hsm = hsm_output.moments_sigma*self.pixel_scale
        num_sigma_from_trace = galsim_img.calculateMomentRadius(center=hsm_output.moments_centroid, rtype='trace')
        num_sigma_from_det = galsim_img.calculateMomentRadius(center=hsm_output.moments_centroid, rtype='det')
        
        assert np.isclose(num_sigma_from_det_hsm, self.gal_sigma)
        assert np.isclose(num_sigma_from_trace, self.gal_sigma)
        assert np.isclose(num_sigma_from_det, self.gal_sigma)
        assert np.isclose(hsm_output.moments_amp, self.gal_flux)
        
    def test_ellip_gaussian_galaxy(self):
        """Compares trace and det of second moment matrix and total flux"""
        gal = galsim.Gaussian(sigma=self.gal_sigma, flux=self.gal_flux).shear(e1=self.gal_e1, e2=self.gal_e2)
        galsim_img = gal.drawImage(scale=self.pixel_scale, nx=self.nx, ny=self.nx, method='no_pixel')

        # Analytical
        e, phi = e1e2_to_ephi(e1=self.gal_e1, e2=self.gal_e2)
        q = np.sqrt((1.0-e)/(1.0+e))
        lam1 = self.gal_sigma**2.0/q
        lam2 = self.gal_sigma**2.0*q
        Ixx = lam1*np.cos(phi)**2.0 + lam2*np.sin(phi)**2.0
        Ixy = (lam1 - lam2)*np.cos(phi)*np.sin(phi)
        Iyy = lam1*np.sin(phi)**2.0 + lam2*np.cos(phi)**2.0
        trace = Ixx + Iyy
        det = Ixx*Iyy - Ixy*Ixy

        # Numerical
        hsm_output = galsim_img.FindAdaptiveMom(guess_sig=self.gal_sigma/self.pixel_scale)
        num_det_hsm = (hsm_output.moments_sigma*self.pixel_scale)**4.0
        num_trace = 2.0*(galsim_img.calculateMomentRadius(center=hsm_output.moments_centroid, rtype='trace')**2.0)
        num_det = galsim_img.calculateMomentRadius(center=hsm_output.moments_centroid, rtype='det')**4.0

        assert np.isclose((Ixx - Iyy)/(Ixx + Iyy), self.gal_e1)
        assert np.isclose((2*Ixy)/(Ixx + Iyy), self.gal_e2)
        assert np.isclose(num_det_hsm, det)
        assert np.isclose(num_det, det)
        assert np.isclose(num_trace, trace)
        assert np.isclose(hsm_output.moments_amp, self.gal_flux)
        
    def test_ellip_gaussian_galaxy_with_psf(self):
        """Compares trace and det of second moment matrix and total flux"""
        gal = galsim.Gaussian(sigma=self.gal_sigma, flux=self.gal_flux).shear(e1=self.gal_e1, e2=self.gal_e2)
        psf = galsim.Gaussian(sigma=self.psf_sigma, flux=1.0)
        total = galsim.Convolve([gal, psf])
        galsim_img = total.drawImage(scale=self.pixel_scale, nx=self.nx, ny=self.nx, method='no_pixel')
        
        # Analytical
        e, phi = e1e2_to_ephi(e1=self.gal_e1, e2=self.gal_e2)
        q = np.sqrt((1.0-e)/(1.0+e))
        lam1 = self.gal_sigma**2.0/q
        lam2 = self.gal_sigma**2.0*q
        Ixx = lam1*np.cos(phi)**2.0 + lam2*np.sin(phi)**2.0 + self.psf_sigma**2.0
        Ixy = (lam1 - lam2)*np.cos(phi)*np.sin(phi)
        Iyy = lam1*np.sin(phi)**2.0 + lam2*np.cos(phi)**2.0 + self.psf_sigma**2.0
        trace = Ixx + Iyy
        det = Ixx*Iyy - Ixy*Ixy
        
        # Numerical
        hsm_output = galsim_img.FindAdaptiveMom(guess_sig=self.gal_sigma/self.pixel_scale)
        num_det_hsm = (hsm_output.moments_sigma*self.pixel_scale)**4.0
        num_trace = 2.0*(galsim_img.calculateMomentRadius(center=hsm_output.moments_centroid, rtype='trace')**2.0)
        num_det = galsim_img.calculateMomentRadius(center=hsm_output.moments_centroid, rtype='det')**4.0

        num_x = pixel_to_physical(hsm_output.moments_centroid.x,
                                  canvas_size=self.nx,
                                  pixel_scale=self.pixel_scale)
        num_y = pixel_to_physical(hsm_output.moments_centroid.y,
                                         canvas_size=self.nx,
                                         pixel_scale=self.pixel_scale)
        
        assert np.allclose([num_x, num_y], [0.0, 0.0])
        assert np.isclose(hsm_output.observed_shape.e1, (Ixx - Iyy)/(Ixx + Iyy))
        assert np.isclose(hsm_output.observed_shape.e2, (2*Ixy)/(Ixx + Iyy))
        assert np.isclose(num_det_hsm, det, rtol=1.e-3)
        assert np.isclose(num_det, det, rtol=1.e-3)
        assert np.isclose(num_trace, trace, rtol=1.e-3)
        assert np.isclose(hsm_output.moments_amp, self.gal_flux)
        
    def test_ellip_shifted_gaussian_galaxy_with_psf(self):
        """Compares trace and det of second moment matrix and total flux"""
        gal = galsim.Gaussian(sigma=self.gal_sigma, flux=self.gal_flux)
        gal = gal.shear(e1=self.gal_e1, e2=self.gal_e2)
        gal = gal.shift(dx=self.gal_x, dy=self.gal_y)
        psf = galsim.Gaussian(sigma=self.psf_sigma, flux=1.0)
        total = galsim.Convolve([gal, psf])
        galsim_img = total.drawImage(scale=self.pixel_scale, nx=self.nx, ny=self.nx, method='no_pixel')
        
        # Analytical
        e, phi = e1e2_to_ephi(e1=self.gal_e1, e2=self.gal_e2)
        q = np.sqrt((1.0-e)/(1.0+e))
        # 1. Centered moments (no PSF)
        lam1 = self.gal_sigma**2.0/q
        lam2 = self.gal_sigma**2.0*q
        Ixx = lam1*np.cos(phi)**2.0 + lam2*np.sin(phi)**2.0 + self.psf_sigma**2.0
        Ixy = (lam1 - lam2)*np.cos(phi)*np.sin(phi)
        Iyy = lam1*np.sin(phi)**2.0 + lam2*np.cos(phi)**2.0 + self.psf_sigma**2.0
        # 2. First moments
        Ix = self.gal_x
        Iy = self.gal_y
        det = Ixx*Iyy - Ixy**2.0
        trace = Ixx + Iyy
        
        # Numerical
        hsm_output = galsim_img.FindAdaptiveMom(guess_sig=self.gal_sigma/self.pixel_scale)
        num_det_hsm = (hsm_output.moments_sigma*self.pixel_scale)**4.0
        num_trace = 2.0*(galsim_img.calculateMomentRadius(center=hsm_output.moments_centroid, rtype='trace')**2.0)
        num_det = galsim_img.calculateMomentRadius(center=hsm_output.moments_centroid, rtype='det')**4.0

        num_x = pixel_to_physical(hsm_output.moments_centroid.x,
                                  canvas_size=self.nx,
                                  pixel_scale=self.pixel_scale)
        num_y = pixel_to_physical(hsm_output.moments_centroid.y,
                                         canvas_size=self.nx,
                                         pixel_scale=self.pixel_scale)
        
        assert np.allclose([num_x, num_y], [Ix, Iy])
        assert np.isclose(hsm_output.observed_shape.e1, (Ixx - Iyy)/(Ixx + Iyy))
        assert np.isclose(hsm_output.observed_shape.e2, (2*Ixy)/(Ixx + Iyy))
        assert np.isclose(num_det_hsm, det, rtol=1.e-3)
        assert np.isclose(num_det, det, rtol=1.e-3)
        assert np.isclose(num_trace, trace, rtol=1.e-3)
        assert np.isclose(hsm_output.moments_amp, self.gal_flux)
        
    def test_ellip_gaussian_galaxy_with_psf_and_qso(self):
        """Compares first moments, trace and det of second moment matrix, and total flux"""
        gal = galsim.Gaussian(sigma=self.gal_sigma, flux=self.gal_flux).shear(e1=self.gal_e1, e2=self.gal_e2)
        qso = galsim.Gaussian(sigma=self.qso_sigma, flux=self.qso_flux).shift(dx=self.qso_x, dy=self.qso_y)
        psf = galsim.Gaussian(sigma=self.psf_sigma)
        total = galsim.Convolve([gal + qso, psf])
        #total = gal + qso
        galsim_img = total.drawImage(scale=self.pixel_scale, nx=self.nx, ny=self.nx, method='no_pixel')
        image_array = galsim_img.array
        
        gal_ana = {}
        qso_ana = {}
        # Analytical
        # 1. Galaxy: centered moments (no PSF)
        gal_ana['e'], gal_ana['phi'] = e1e2_to_ephi(e1=self.gal_e1, e2=self.gal_e2)
        gal_ana['q'] = np.sqrt((1.0 - gal_ana['e'])/(1.0 + gal_ana['e']))
        gal_ana['lam1'] = self.gal_sigma**2.0/gal_ana['q']
        gal_ana['lam2'] = self.gal_sigma**2.0*gal_ana['q']
        gal_ana['Ixx_centered'] = gal_ana['lam1']*np.cos(gal_ana['phi'])**2.0 + gal_ana['lam2']*np.sin(gal_ana['phi'])**2.0
        gal_ana['Ixy_centered'] = (gal_ana['lam1'] - gal_ana['lam2'])*np.cos(gal_ana['phi'])*np.sin(gal_ana['phi'])
        gal_ana['Iyy_centered'] = gal_ana['lam1']*np.sin(gal_ana['phi'])**2.0 + gal_ana['lam2']*np.cos(gal_ana['phi'])**2.0
        # 2. QSO: centered moments (no PSF)
        qso_ana['Ixx_centered'] = self.qso_sigma**2.0
        qso_ana['Ixy_centered'] = 0.0
        qso_ana['Iyy_centered'] = self.qso_sigma**2.0
        # 3. Define flux ratios for weighting
        total_flux = self.gal_flux + self.qso_flux
        gal_flux_ratio = self.gal_flux/total_flux
        qso_flux_ratio = self.qso_flux/total_flux
        # 4. First moments
        Ix = qso_flux_ratio*self.qso_x
        Iy = qso_flux_ratio*self.qso_y
        # 6. Off-center moments of galaxy (no PSF)
        gal_ana['Ixx_uncentered'] = gal_ana['Ixx_centered'] + Ix**2.0 # TODO: check sign
        gal_ana['Ixy_uncentered'] = gal_ana['Ixy_centered'] + Ix*Iy # TODO: check sign
        gal_ana['Iyy_uncentered'] = gal_ana['Iyy_centered'] + Iy**2.0 # TODO: check sign
        # 7. Off-center mometns of QSO (no PSF)
        qso_x_to_Ix = self.qso_x - Ix
        qso_y_to_Iy = self.qso_y - Iy
        qso_ana['Ixx_uncentered'] = qso_ana['Ixx_centered'] + (qso_x_to_Ix)**2.0
        qso_ana['Ixy_uncentered'] = qso_ana['Ixy_centered'] + (qso_x_to_Ix)*(qso_y_to_Iy)
        qso_ana['Iyy_uncentered'] = qso_ana['Iyy_centered'] + (qso_y_to_Iy)**2.0
        # 8. Weighted moments of galaxy and QSO (no PSF)
        Ixx = gal_flux_ratio*gal_ana['Ixx_uncentered'] + qso_flux_ratio*qso_ana['Ixx_uncentered']
        Ixy = gal_flux_ratio*gal_ana['Ixy_uncentered'] + qso_flux_ratio*qso_ana['Ixy_uncentered']
        Iyy = gal_flux_ratio*gal_ana['Iyy_uncentered'] + qso_flux_ratio*qso_ana['Iyy_uncentered']
        # 9. Add PSF
        Ixx += self.psf_sigma**2.0
        Iyy += self.psf_sigma**2.0
        
        # Numerical
        num_Ix, num_Iy = get_first_moments_from_image(image_array, pixel_scale=self.pixel_scale)
        num_Ixx, num_Ixy, num_Iyy = get_second_moments_from_image(image_array, pixel_scale=self.pixel_scale)
        
        assert np.allclose([num_Ix, num_Iy], [Ix, Iy], rtol=1.e-3)
        assert np.isclose(num_Ixx, Ixx, rtol=1.e-3)
        assert np.isclose(num_Ixy, Ixy, rtol=1.e-3)
        assert np.isclose(num_Iyy, Iyy, rtol=1.e-3)
        assert np.isclose(np.sum(galsim_img.array), total_flux)

if __name__ == '__main__':
    unittest.main()