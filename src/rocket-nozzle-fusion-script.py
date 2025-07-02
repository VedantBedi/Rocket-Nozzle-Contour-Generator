import adsk.core, adsk.fusion, adsk.cam, traceback
import math

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
        design = app.activeProduct
        rootComp = design.rootComponent
        
        gamma = 1.4
        R = 287.0
        
        sketches = rootComp.sketches
        xzPlane = rootComp.xZConstructionPlane
        sketch = sketches.add(xzPlane)
        
        (returnValue, cancelled) = ui.inputBox(
            'Enter parameters as: P0,T0,mdot,Me,Scale,Thickness',
            'Nozzle Parameters',
            '2000000,3000,100,3.5,0.5,2')
        
        if cancelled:
            return
        
        try:
            params = returnValue.split(',')
            P0 = float(params[0])
            T0 = float(params[1])
            mdot = float(params[2])
            Me = float(params[3])
            scale_factor = float(params[4])
            thickness = float(params[5])
        except:
            ui.messageBox('Invalid input format. Using defaults.')
            P0, T0, mdot, Me = 2e6, 3000, 100, 3.5
            scale_factor = 0.0000005
            thickness = 2.0
        
        mdot_star = P0 * math.sqrt(gamma/(R*T0)) * ((gamma+1)/2) ** (-(gamma+1)/(2*(gamma-1)))
        A_throat = mdot / mdot_star
        R_throat = math.sqrt(A_throat / math.pi)
        
        x_nozzle, r_nozzle = create_nozzle_profile(Me, R_throat, gamma)
        
        x_mm = [x * 1000 * scale_factor for x in x_nozzle]
        r_mm = [r * 1000 * scale_factor for r in r_nozzle]
        
        small_offset = 0.1
        points_outer = adsk.core.ObjectCollection.create()
        for i in range(len(x_mm)):
            points_outer.add(adsk.core.Point3D.create(x_mm[i], r_mm[i] + small_offset, 0))
        
        topSpline_outer = sketch.sketchCurves.sketchFittedSplines.add(points_outer)
        
        inner_r_mm = [r - thickness for r in r_mm]
        points_inner = adsk.core.ObjectCollection.create()
        for i in range(len(x_mm)):
            points_inner.add(adsk.core.Point3D.create(x_mm[i], inner_r_mm[i] + small_offset, 0))
        
        topSpline_inner = sketch.sketchCurves.sketchFittedSplines.add(points_inner)
        
        leftEnd_outer = sketch.sketchCurves.sketchLines.addByTwoPoints(
            adsk.core.Point3D.create(x_mm[0], 0, 0),
            adsk.core.Point3D.create(x_mm[0], r_mm[0] + small_offset, 0)
        )
        rightEnd_outer = sketch.sketchCurves.sketchLines.addByTwoPoints(
            adsk.core.Point3D.create(x_mm[-1], 0, 0),
            adsk.core.Point3D.create(x_mm[-1], r_mm[-1] + small_offset, 0)
        )
        
        leftEnd_inner = sketch.sketchCurves.sketchLines.addByTwoPoints(
            adsk.core.Point3D.create(x_mm[0], inner_r_mm[0] + small_offset, 0),
            adsk.core.Point3D.create(x_mm[0], 0, 0)
        )
        rightEnd_inner = sketch.sketchCurves.sketchLines.addByTwoPoints(
            adsk.core.Point3D.create(x_mm[-1], inner_r_mm[-1] + small_offset, 0),
            adsk.core.Point3D.create(x_mm[-1], 0, 0)
        )
        
        baseLine = sketch.sketchCurves.sketchLines.addByTwoPoints(
            adsk.core.Point3D.create(x_mm[0], 0, 0),
            adsk.core.Point3D.create(x_mm[-1], 0, 0)
        )
        
        if sketch.profiles.count < 2:
            ui.messageBox('Failed: No valid profile found for revolution.')
            return
        
        revolves = rootComp.features.revolveFeatures
        outer_prof = sketch.profiles.item(0)
        inner_prof = sketch.profiles.item(1)
        xAxis = rootComp.xConstructionAxis
        
        revolveInput_outer = revolves.createInput(outer_prof, xAxis, adsk.fusion.FeatureOperations.NewBodyFeatureOperation)
        revolveInput_outer.setAngleExtent(False, adsk.core.ValueInput.createByReal(math.pi * 2))
        revolves.add(revolveInput_outer)
        
        revolveInput_inner = revolves.createInput(inner_prof, xAxis, adsk.fusion.FeatureOperations.CutFeatureOperation)
        revolveInput_inner.setAngleExtent(False, adsk.core.ValueInput.createByReal(math.pi * 2))
        revolves.add(revolveInput_inner)
        
        ui.messageBox(f'Nozzle created successfully!\nThroat Radius: {R_throat*1000*scale_factor:.2f} mm')
        
    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

def area_ratio(M, gamma=1.4):
    """Calculate isentropic area ratio A/A* for given Mach number M"""
    term = (2.0/(gamma+1)) * (1 + (gamma-1)/2 * M**2)
    exp = (gamma+1)/(2*(gamma-1))
    return (1.0/M) * (term ** exp)


def mach_from_area_ratio(AR_value, subsonic=False, gamma=1.4):
    """Calculate Mach number from area ratio using bisection method"""
    if abs(AR_value - 1.0) < 1e-6:
        return 0.99 if subsonic else 1.01
    
    M_low = 0.01 if subsonic else 1.01
    M_high = 0.99 if subsonic else 10.0
    
    for _ in range(50):
        M_mid = (M_low + M_high) / 2.0
        AR_mid = area_ratio(M_mid, gamma)
        
        if abs(AR_mid - AR_value) < 1e-6:
            return M_mid
        
        if subsonic:
            if AR_mid > AR_value:
                M_low = M_mid
            else:
                M_high = M_mid
        else:
            if AR_mid > AR_value:
                M_high = M_mid
            else:
                M_low = M_mid
    
    return (M_low + M_high) / 2.0


def create_nozzle_profile(Me, R_throat, gamma=1.4, L_con=None, L_div=None, R_chamber=None):
    """Generate nozzle profile coordinates (x, radius) without using NumPy"""
    # Default values based on throat radius
    L_con = L_con or 2.5 * R_throat
    L_div = L_div or 6.0 * R_throat
    R_chamber = R_chamber or 3.0 * R_throat
    
    # Calculate exit radius from area ratio
    AR_exit = area_ratio(Me, gamma)
    R_exit = R_throat * math.sqrt(AR_exit)
    
    # Convergent section using cubic polynomial
    a = 2.0 * (R_chamber - R_throat) / (L_con**3)
    b = 3.0 * (R_chamber - R_throat) / (L_con**2)
    
    # Generate points
    num_con_points = 100
    num_div_points = 200
    
    # Create lists to store points
    x_con = []
    r_con = []
    x_div = []
    r_div = []
    
    # Generate convergent section points
    for i in range(num_con_points):
        x = -L_con + i * (L_con / (num_con_points - 1))
        r = a * x**3 + b * x**2 + R_throat
        x_con.append(x)
        r_con.append(r)
    
    # Generate divergent section points
    for i in range(num_div_points):
        x = i * (L_div / (num_div_points - 1))
        r = R_throat + (R_exit - R_throat) * math.sin((math.pi/2) * (x / L_div))**2
        x_div.append(x)
        r_div.append(r)
    
    # Combine points (exclude duplicated throat point)
    x_full = x_con + x_div[1:]
    r_full = r_con + r_div[1:]
    
    return x_full, r_full


def calculate_mach_distribution(x, r, R_throat, gamma=1.4):
    """Calculate Mach number distribution along the nozzle"""
    M_vals = []
    for i in range(len(x)):
        AR_local = (r[i] / R_throat)**2
        if x[i] < 0:
            M_vals.append(mach_from_area_ratio(AR_local, subsonic=True, gamma=gamma))
        elif abs(x[i]) < 1e-6:
            M_vals.append(1.0)
        else:
            M_vals.append(mach_from_area_ratio(AR_local, subsonic=False, gamma=gamma))
    return M_vals