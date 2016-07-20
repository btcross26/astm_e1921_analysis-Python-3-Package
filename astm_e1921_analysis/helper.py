def ksi_to_MPa(stress_in_ksi):
    '''
    Convert a measurement in ksi to units of MPa.
    
    Arguments
    ---------
    stress_in_ksi: float
        Measurement in ksi.
        
    Returns
    -------
    Measurement in MPa.
    '''
    return 6.89476 * stress_in_ksi

    
def in_to_mm(length_in_inches):
    '''
    Convert a measurement in inches to millimeters.
    
    Arguments
    ---------
    length_in_inches: float
        Measurement in inches.
        
    Returns
    -------
    Measurement in millimeters (mm).
    '''
    return 25.4 * length_in_inches

    
def lbf_to_newtons(force_in_lbf):
    '''
    Convert a measurement in lbf to newtons.
    
    Arguments
    ---------
    force_in_lbf: float
        Measurement in lbf.
        
    Returns
    -------
    Measurement in newtons (N).
    '''
    return 4.44822 * force_in_lbf

    
def stress_intensity_in_SI_to_english(K_in_MPa_root_meter):
    '''
    Convert a stress_intensity measurement in MPa root meters to ksi root
    inches.
    
    Arguments
    ---------
    K_in_MPa_root_meter: float
        Measurement in MPa root meters.
        
    Returns
    -------
    Measurement in ksi root inches.
    '''
    return 0.9100477 * K_in_MPa_root_meter

    
def adjust_KJc(KJc, B0, BX):
    '''
    Adjust a stress_intensity measurement in MPa root meters to account for
    size.
    
    Arguments
    ---------
    KJc: float
        KJc in MPa root meters.    
    B0: float
        Actual specimen size (length).      
    BX: float
        Size for which stress intensity should be calculated (length).
        
    Returns
    -------
    Adjusted KJc in MPa root meters.
    '''
    return 20.0 + (KJc - 20.0) * (B0 / BX)**(0.25)

    
def celsius_to_fahrenheit(temp_in_celsius):
    '''
    Convert a temperature in degrees celsius to degrees fahrenheit.
    
    Arguments
    ---------
    temp_in_celsius: float
        Measurement in degrees celsius.
        
    Returns
    -------
    Measurement in degrees fahrenheit.
    '''
    return 1.8 * temp_in_celsius + 32.0

def fahrenheit_to_celsius(temp_in_fahrenheit):
    '''
    Convert a temperature in degrees fahrenheit to degrees celsius.
    
    Arguments
    ---------
    temp_in_fahrenheit: float
        Measurement in degrees fahrenheit.
        
    Returns
    -------
    Measurement in degrees celsius.
    '''
    return (temp_in_fahrenheit - 32.0) / 1.8