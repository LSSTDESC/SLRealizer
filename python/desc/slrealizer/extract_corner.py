import numpy as np
import desc.slrealizer
import math

"""
This file contains helper methods to draw cornerplot.
Each method extracts the information user requested (ex. size, position, ellipticity, etd)
"""

def extract_features(df, names):
    """
    Parameters
    ----------
    df: csv toy catalog
    names : csv toy catalog column names

    Returns
    ---------
    features: float, ndarray
    Features requested for the cornerplot
    labels: strings, list
    Corresponding labels
    """
    features = np.array([])
    labels = []

    p = len(names)
    n = len(df)

    for name in names:
        features = np.append(features, df[name])
        labels.append(axis_labels[name])

    return features.reshape(p,n).transpose(), labels

def calculate_size(df):

    """
    Parameters
    ----------
    df : csv toy catalog

    Returns
    ----------
    features: float, ndarray
    Sizes for each filter

    labels: string, list
    Corresponding axis labels
    """

    features = np.array([])
    labels = []
    
    for filter in ['u', 'g', 'r', 'i', 'z']:
        features = np.append(features, df[filter+'_size'])
        labels.append(axis_labels[filter+'size'])

    return features, labels


def calculate_phi(df):

    """                                                                         
    Parameters                                                                  
    ----------                                                                  
    df : csv toy catalog                                                        
                                                                                
    Returns                                                                     
    ----------                                                                  
    features: float, ndarray                                                    
    rotation angle for each filter                                                       
                                                                                
    labels: string, list                                                        
    Corresponding axis labels                                                   
    """

    features = np.array([])
    labels = []

    for filter in ['u', 'g', 'r', 'i', 'z']:
        features = np.append(features, df[filter+'_phi'])
        labels.append(axis_labels[filter+'phi'])

    return features, labels
    #return features.reshape(5, len(df)).transpose(), labels  

def calculate_ellipticity(df):

    """
    Parameters
    ----------
    df : csv toy catalog
    
    Returns
    ----------
    features: float, ndarray
    ellipticities for each filter
    
    labels: string, list
    Corresponding axis labels
    """
    features = np.array([])
    labels = []

    for filter in ['u', 'g', 'r', 'i', 'z']:
        e = df[filter+'_e']
        features = np.append(features,e)
        labels.append(axis_labels[filter+'e'])
    return features, labels

def calculate_magnitude(df):

    """
    Parameters
    ----------
    df : csv toy catalog
    
    Returns
    ----------
    features: float, ndarray
    Magnitudes for each filter
    
    labels: string, list
    Corresponding axis labels
    """
        
    features = np.array([])
    labels = []

    for filter in ['u', 'g', 'r', 'i', 'z']:
        filter_flux = df[filter+'_flux']
        filter_mag = desc.slrealizer.return_zeropoint()-2.5*np.log10(filter_flux)
        features = np.append(features, filter_mag)
        labels.append(axis_labels[filter+'mag'])

    return features, labels

def calculate_color(df):

    labels = []
    features = np.array([])

    for filters in [['u', 'g'], ['g', 'r'], ['r', 'i'], ['i', 'z']]:
        filter_flux_1 = df[filters[0]+'_flux']
        filter_mag_1 = desc.slrealizer.return_zeropoint()-2.5*np.log10(filter_flux_1)
        filter_flux_2 = df[filters[1]+'_flux']
        filter_mag_2 = desc.slrealizer.return_zeropoint()-2.5*np.log10(filter_flux_2)
        features = np.append(features, (filter_mag_1 - filter_mag_2))
        labels.append(axis_labels[filters[0]+filters[1]])

    return features, labels

def calculate_position(df):
    # using pythagorean theorem calculate the distance between centroids
    x_distance, _ = desc.slrealizer.calculate_x_position(df)
    y_distance, _ = desc.slrealizer.calculate_y_position(df)
    x_distance_squared = np.multiply(x_distance, x_distance)
    y_distance_squared = np.multiply(y_distance, y_distance)
    distance = np.sqrt(x_distance_squared+y_distance_squared)
    labels = []
    for filter in ['u', 'g', 'i', 'z']:
        labels.append(axis_labels[filter+'pos'])
    return distance, labels

def calculate_x_position(df):
    # reference filter : r

    """
    Parameters
    ----------
    df : csv toy catalog
    
    Returns
    ----------
    features: float, ndarray
    Center of moments for x axis for each filter
    
    labels: string, list
    Corresponding axis labels
    """
    features = np.array([])
    labels = []

    for filter in ['u', 'g', 'i', 'z']:
        name = filter+'_x'
        filter_pos = (df[name] - df['r_x'])
        features = np.append(features, filter_pos*desc.slrealizer.get_pixel_arcsec_conversion())
        labels.append(axis_labels[filter+'xpos'])
    
    return features, labels

def calculate_y_position(df):
    # reference filter : i
    """
        Parameters
        ----------
        df : csv toy catalog
        
        Returns
        ----------
        features: float, ndarray
        Center of moments for y axis for each filter
        
        labels: string, list
        Corresponding axis labels
        """

    features = np.array([])
    labels = []

    for filter in ['u', 'g', 'i', 'z']:
        name = filter+'_y'
        filter_pos = (df[name] - df['r_y'])
        features = np.append(features, filter_pos*desc.slrealizer.get_pixel_arcsec_conversion())
        labels.append(axis_labels[filter+'ypos'])

    return features, labels

#============================================================================================

axis_labels = {}
axis_labels['g_flux'] = '$g$'
axis_labels['z_flux'] = '$z$'
axis_labels['r_flux'] = '$r$'
axis_labels['u_flux'] = '$u$'
axis_labels['i_flux'] = '$i$'
axis_labels['g_x'] = '$g_x / arcsec$'
axis_labels['z_x'] = '$z_x / arcsec$'
axis_labels['r_x'] = '$r_x / arcsec$'
axis_labels['u_x'] = '$u_x / arcsec$'
axis_labels['i_x'] = '$i_x / arcsec$'
axis_labels['g_y'] = '$g_y / arcsec$'
axis_labels['z_y'] = '$z_y / arcsec$'
axis_labels['r_y'] = '$r_y / arcsec$'
axis_labels['u_y'] = '$u_y / arcsec$'
axis_labels['i_y'] = '$i_y / arcsec$'
axis_labels['gsize'] = '$size_g / arcsec$'
axis_labels['zsize'] = '$size_z / arcsec$'
axis_labels['rsize'] = '$size_r / arcsec$'
axis_labels['usize'] = '$size_u / arcsec$'
axis_labels['isize'] = '$size_i / arcsec$'
axis_labels['ge'] = '$e_g$'
axis_labels['ze'] = '$e_z$'
axis_labels['re'] = '$e_r$'
axis_labels['ue'] = '$e_u$'
axis_labels['ie'] = '$e_i$'
axis_labels['gmag'] = '$mag_g$'
axis_labels['zmag'] = '$mag_z$'
axis_labels['rmag'] = '$mag_r$'
axis_labels['umag'] = '$mag_u$'
axis_labels['imag'] = '$mag_i$'
axis_labels['gr'] = '$g-r$'
axis_labels['ri'] = '$r-i$'
axis_labels['iz'] = '$i-z$'
axis_labels['ug'] = '$u-g$'
axis_labels['gxpos'] = '$x_g / arcsec$'
axis_labels['zxpos'] = '$x_z / arcsec$'
axis_labels['rxpos'] = '$x_r / arcsec$'
axis_labels['uxpos'] = '$x_u / arcsec$'
axis_labels['ixpos'] = '$x_i / arcsec$'
axis_labels['gypos'] = '$y_g / arcsec$'
axis_labels['zypos'] = '$y_z / arcsec$'
axis_labels['rypos'] = '$y_r / arcsec$'
axis_labels['uypos'] = '$y_u / arcsec$'
axis_labels['iypos'] = '$y_i / arcsec$'
axis_labels['zpos'] = '$pos_z / arcsec$'
axis_labels['gpos'] = '$pos_g / arcsec$'
axis_labels['upos'] = '$pos_u / arcsec$'
axis_labels['ipos'] = '$pos_i / arcsec$'
axis_labels['gphi'] = '$phi_g$'
axis_labels['zphi'] = '$phi_z$'
axis_labels['rphi'] = '$phi_r$'
axis_labels['uphi'] = '$phi_u$'
axis_labels['iphi'] = '$phi_i$'
