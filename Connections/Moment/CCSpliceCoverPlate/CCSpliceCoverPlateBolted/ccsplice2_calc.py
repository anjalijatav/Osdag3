'''
# Created on 30-Oct-2017
# Revised on 5-March-2018
# Revised on 13-April-2018
# Revised on 15-June-2018 (experts suggestions)
# Revised on 18-June-2018 (experts suggestions)
# Revised on 25-June-2018 (After launch)
@author: Kumari Anjali.
'''
'''
ASCII diagram- Column-Column Bolted Splice Connection with Cover Plates
'''

import sys

from Connections.Moment.CCSpliceCoverPlate.CCSpliceCoverPlateBolted.model import *
from utilities.is800_2007 import *
from utilities.other_standards import IS1363_part_1_2002, IS1363_part_3_2002, IS1367_Part3_2002
from Connections.connection_calculations import ConnectionCalculations
import math
import logging
import numpy as np
import sys

flag = 1
logger = None


def module_setup():
    global logger
    logger = logging.getLogger("osdag.ccsplice_calc")
    module_setup()
    #######################################################################
    #######################################################################
    ##check for section plastic or elastic###


def flange_force(column_d, column_f_t, column_b, column_area, factored_axial_force, moment_load):
    """
    Args:
       Column_d: Overall depth of the column section in mm (float)
       Column_b: width of the column section in mm (float)
       Column_f_t: Thickness of flange in mm (float)
       axial_force: Factored axial force in kN (float)
       moment_load: Factored bending moment in kN-m (float)
    Returns:
        Force in flange in kN (float)
    """

    area_f = column_b * column_f_t
    axial_force_f = (area_f * factored_axial_force / column_area)
    f_f = ((moment_load * 1000) / (column_d - column_f_t)) + (axial_force_f / 100)
    #print(f_f)
    return (f_f)
def web_force(column_d, column_f_t, column_t_w, factored_axial_force, column_area):
    """
    Args:
       column_d: Overall depth of the column section in mm (float)
       column_f_t: Thickness of flange in mm (float)
       column_t_w: Thickness of flange in mm (float)
       axial_force: Factored axial force in kN (float)

    Returns:
        Force in flange in kN (float)
    """
    axial_force_w = int((column_d - 2 *(column_f_t )) * column_t_w * factored_axial_force) / column_area
    return round(axial_force_w)
def thk_flange_plate(column_d, column_f_t, bolt_diameter, column_area, axial_force, moment_load, column_b,
                     column_fy, dia_hole):
    """
    Args:
        column_d: Overall depth of the column section in mm (float)
        column_f_t: Thickness of flange in mm (float)
        axial_force: Factored axial force in kN (float)
        moment_load: Factored bending moment in kN-m (float)
        column_b: Width of flange in mm (float)
        column_fy: Characteristic yield stress in N/mm^2 (float)
        pitch = gauge
    Returns:
    """

    edge_dist = 1.5 * dia_hole
    gauge = 2.5 * bolt_diameter
    gauge_req = (column_b - (2 * edge_dist)) / gauge  # number of gauge dist required along flange width
    n = gauge_req + 1  # number of bolts along flange width
    gamma_m0 = 1.10  # Partial safety factor against yield stress and buckling = 1.10 (float)
    ff = float(flange_force(column_d, column_f_t, column_b, column_area, axial_force, moment_load))
    flangeplatethickness = ff / ((column_b - n * dia_hole) * (column_fy / (gamma_m0 * 1000)))

    return round(flangeplatethickness, 2)  # mm

## Minimum thickness of flange splice plate
def flange_plate_min_t(column_t_w):
    """
    Args:
        column_t_w: thickness of web in mm(float)
    Returns: Maximum thickness of flange splice plate in mm (float)
    """
    flangeplatemint = min(column_t_w / 2, 10)
    return float(flangeplatemint)

# Thickness of inner flange splice plate [Reference: N. Subramanian (Page 428), M.L. Gambhir (Page 10.84)]
def func_flangeinnerplatet(column_d, column_f_t, column_area, bolt_diameter, axial_force, moment_load, column_b,
                           column_fy, dia_hole):
    """
    Args:
        column_d: Overall depth of the column section in mm (float)
        column_f_t: Thickness of flange in mm (float)
        axial_force: Factored axial force in kN (float)
        moment_load: Factored bending moment in kN-m (float)
        column_b: Width of flange in mm (float)
        column_fy: Characteristic yield stress in N/mm^2 (float)
        pitch = gauge
    Returns:
    """

    edge_dist = 1.5 * dia_hole
    gauge = 2.5 * bolt_diameter
    gauge_req = (column_b - (2 * edge_dist)) / gauge  # number of gauge dist required along flange width
    n = gauge_req + 1  # number of bolts along flange width
    gamma_m0 = 1.10  # Partial safety factor against yield stress and buckling = 1.10 (float)
    ff = float(flange_force(column_d, column_f_t, column_b, column_area, axial_force, moment_load))
    flangeinnerplatethickness = float(ff / ((column_b - n * dia_hole) * (column_fy / (gamma_m0 * 1000))))

    return round(flangeinnerplatethickness, 2)  # mm
## Minimum width of flange splice plate

def flange_plate_width(edge_dist, number_of_column_flange, flange_gauge):
    flangeplatewidthmin = float((edge_dist + (number_of_column_flange -1)*flange_gauge)*2)
    return round(flangeplatewidthmin, 2)

## Maximum width of flange splice plate
def flange_plate_w_max(column_b):
    flangeplatewidthmax = column_b
    return round(flangeplatewidthmax, 2)
## Maximum Height of flange splice plate
def flange_plate_l_max(end_dist_max_f,number_of_row_flange, pitch_dist_max_f):
    flangeplateheightmax = float((end_dist_max_f + (number_of_row_flange -1)* pitch_dist_max_f)*2)
    return round((flangeplateheightmax) ,2)
## Maximum Height of flange splice plate
def flange_plate_l_min(end_dist_f,number_of_row_flange, flange_pitch):
    flangeplateheightmax = float((end_dist_f + (number_of_row_flange -1)* flange_pitch)*2)
    return round((flangeplateheightmax) ,2)

#######################################################################
## Minimum thickness of web splice plate [Reference: N. subramanian, page 373]
def web_min_t( shear_load, column_d, gamma_m0, column_fy):
    """
    Args:
        column_t_w: thickness of web in mm(float)
    Returns: Maximum thickness of web splice plate in mm (float)
    """

    minwebt = float(math.sqrt(3) * (gamma_m0 * shear_load) / (0.6 * column_d * column_fy))
    return minwebt

## Maximum thickness of web splice plate [Reference: Handbook on structural steel detailing, INSDAG - Chapter 5, section 5.2.3 page 5.7]
def web_max_t(bolt_diameter):
    """
    Args:
        bolt_diameter: Nominal bolt diameter in mm (int)
    Returns: Maximum thickness of web splice plate in mm (float)
    """
    max_web_t = round((0.5 * bolt_diameter), 2)
    return max_web_t

# Maximum width of flange splice plate
def web_max_w(column_d, column_f_t, column_r1):
    """
    Args:
        column_d: Overall depth of supported beam (float) in mm
        column_f_t: Thickness of flange in mm (float)
        column_r1: Root radius of the beam section in mm (float)
    Returns: Maximum width of web splice plate in mm (float)
    """
    maxwebwidth = float(round((column_d - 2 * column_f_t - 2 * column_r1 - 2 * 5), 2))
    return maxwebwidth

    ## Width of web splice plate
# Minimum width of web splice plate [Reference: Steel Designer`s Manual - SCI - 6th edition, page 754]
def web_min_w(column_d):
    """
    Args:
    column_d: Overall depth of supported column (float) in mm
    Returns: Minimum width of web splice plate (float)
    """
    webminw = round((0.5 * column_d), 2)
    return webminw

## Height of Web splice plate
def web_plate_l_max(end_dist_max_w, number_of_row_web, pitch_dist_max_w):
    webplateheight = float((2*end_dist_max_w + (number_of_row_web - 1) * pitch_dist_max_w) * 2)
    return round((webplateheight), 2)

## Minimum thickness of web splice plate
def web_plate_l_min(end_dist_w, number_of_row_web, web_pitch):
    webplateheightmin = (float((2*end_dist_w + (number_of_row_web - 1) * web_pitch) * 2))
    return round(webplateheightmin,2)



    #print(hor_force, ver_force)

    # if shearresbolt > web_bolt_capacity:
    #     design_status = False
    #     logger.error(": Number of bolts is not sufficient")
    #     logger.warning(": shear_res_bolt should be less than web_bolt_capacity of web")
    #     logger.info(": Increase number of bolts and spacing between bolts")
    # else:
    #     pass
    #
  #
           #######################################################################
    # Start of Main Program


def coverplateboltedconnection(uiObj):
    global logger
    logger = logging.getLogger("osdag.ccsplice_calc")
    global design_status
    design_status = True

    connectivity = uiObj["Member"]["Connectivity"]
    column_section = uiObj["Member"]["ColumnSection"]
    column_fu = float(uiObj["Member"]["fu (MPa)"])
    column_fy = float(uiObj["Member"]["fy (MPa)"])

    axial_force = float(uiObj["Load"]["AxialForce"])
    moment_load = float(uiObj["Load"]["Moment (kNm)"])
    shear_load = float(uiObj["Load"]["ShearForce (kN)"])
    if shear_load == '':
        shear_load = 0
    else:
        shear_load = float(uiObj["Load"]["ShearForce (kN)"])

    bolt_diameter = int(uiObj["Bolt"]["Diameter (mm)"])
    bolt_grade = float(uiObj["Bolt"]["Grade"])
    bolt_type = (uiObj["Bolt"]["Type"])
    flange_plate_preference = uiObj['FlangePlate']['Preferences']
    gap = float(uiObj["detailing"]["gap(mm)"])  # gap between web  plate and column flange

    mu_f = float(uiObj["bolt"]["slip_factor"])
    dp_bolt_hole_type = str(uiObj["bolt"]["bolt_hole_type"])
    #dp_bolt_hole_type = ["standard","over_size", "short_slot","long_slot"]
    #dia_hole = IS800_2007.cl_10_2_1_bolt_hole_size(d =bolt_diameter, bolt_hole_type='standard')

    dia_hole = int(uiObj["bolt"]["bolt_hole_clrnce"]) + bolt_diameter
    bolt_fu = float(uiObj["bolt"]["bolt_fu"])
    # bolt_fy = float(uiObj["bolt"]["bolt_fy"])

    type_edge = str(uiObj["detailing"]["typeof_edge"])

    if uiObj["detailing"]["typeof_edge"] == "a - Sheared or hand flame cut":
        type_edge = 'hand_flame_cut'
    else:   # "b - Rolled, machine-flame cut, sawn and planed"
        type_edge = 'machine_flame_cut'
    flange_plate_t = float(uiObj["FlangePlate"]["Thickness (mm)"])
    flange_plate_w = float(uiObj["FlangePlate"]["Width (mm)"])
    flange_plate_l = str(uiObj["FlangePlate"]["Height (mm)"])
    flange_plate_l = float(flange_plate_l)
    # flange_plate_fu = float(uiObj["Member"]["fu (Mpa)"])
    # flange_plate_fy = float(uiObj["Member"]["fy (MPa)"])
    flange_plate_fu = float(column_fu)
    flange_plate_fy = float(column_fy)
    web_plate_t = float(uiObj['WebPlate']['Thickness (mm)'])
    web_plate_l = str(uiObj["WebPlate"]["Height (mm)"])
    web_plate_l = float(web_plate_l)
    web_plate_w = str(uiObj["WebPlate"]["Width (mm)"])
    web_plate_w = float(web_plate_w)
    web_plate_l = float(web_plate_l)
    web_plate_fu = float(column_fu)
    web_plate_fy = float(column_fy)

    old_column_section = get_oldcolumncombolist()

    if column_section in old_column_section:
        logger.warning(" : You are using a section (in red color) that is not available in latest version of IS 808")
    if column_fu < 410 or column_fy < 230 or flange_plate_fu or flange_plate_fy or web_plate_fu or web_plate_fy:
        logger.warning(" : You are using a section of grade that is not available in latest version of IS 2062")

    ########################################################################################################################
    # Read input values from Column  database

    dictcolumndata = get_columndata(column_section)

    column_t_w = float(dictcolumndata["tw"])
    column_f_t = float(dictcolumndata["T"])
    column_d = float(dictcolumndata["D"])

    column_r1 = float(dictcolumndata["R1"])
    column_b = float(dictcolumndata["B"])
    column_area = float(dictcolumndata["Area"])
    column_Zz = float(dictcolumndata["Zz"])
    column_Zy = float(dictcolumndata["Zy"])
    # Minimum Design Action (Cl. 10.7, IS 800:2007)
    # Axial Capacity #kN
    gamma_m0 = 1.1
    Axial_capacity = 0.3 * column_area *100 * column_fy / (gamma_m0*1000)
    factored_axial_force = max(Axial_capacity, axial_force)

    # Shear Capacity   # kN

    A_v = column_area  #todo
    design_shear_capacity = (A_v *100* column_fy) / (math.sqrt(3) * gamma_m0 * 1000)  # kN # A_v: Total cross sectional area in shear in mm^2 (float)
    if shear_load >= design_shear_capacity:
        design_status = False

    #######################################################################
    # Calculation of Spacing (Min values rounded to next multiple of 5)

    #  Minimum pitch & gauge distance  of flange and web plate(mm)
    pitch_dist_min_f = IS800_2007.cl_10_2_2_min_spacing(bolt_diameter)
    pitch_dist_min_w = IS800_2007.cl_10_2_2_min_spacing(bolt_diameter)

    platethickness_f = [column_f_t, flange_plate_t]
    platethickness_w = [column_t_w, web_plate_t]

    web_t_thinner = (min(column_t_w, web_plate_t))
    flange_t_thinner = (min(column_f_t, flange_plate_t))

    # Maximum pitch and gauge distance for Flange splice plate
    pitch_dist_max_f = IS800_2007.cl_10_2_3_1_max_spacing(plate_thicknesses =   platethickness_f)
    gauge_dist_max_f = pitch_dist_max_f

    # Maximim pitch and gauge distance for Web splice plate
    pitch_dist_max_w = IS800_2007.cl_10_2_3_1_max_spacing( plate_thicknesses =  platethickness_w )
    gauge_dist_max_w = pitch_dist_max_w
    flange_pitch = round(pitch_dist_min_f + 10)
    flange_gauge = flange_pitch
    web_pitch = round(pitch_dist_min_w + 10)
    web_gauge = round(pitch_dist_min_w + 10)
    # min_end_distance & max_end_distance = Minimum and Maximum end distance
    #       [Cl. 10.2.4.2 & Cl. 10.2.4.3, IS 800:2007]

    # min_end_distance and max end distance for flange plate
    corrosive_influences = False
    if uiObj['detailing']['is_env_corrosive'] == "Yes":
        corrosive_influences = True

    [bolt_shank_area, bolt_net_area] = IS1367_Part3_2002.bolt_area(bolt_diameter)

    end_dist_min_f = IS800_2007.cl_10_2_4_2_min_edge_end_dist(
        d=bolt_diameter, bolt_hole_type= dp_bolt_hole_type, edge_type=type_edge)
    end_dist_max_f = IS800_2007.cl_10_2_4_3_max_edge_dist(
        plate_thicknesses=platethickness_f, f_y=column_fy, corrosive_influences=corrosive_influences)

    # min_end_distance and max end distance for web plate
    end_dist_min_w = IS800_2007.cl_10_2_4_2_min_edge_end_dist(
        d=bolt_diameter, bolt_hole_type=dp_bolt_hole_type, edge_type=type_edge)
    end_dist_max_w = IS800_2007.cl_10_2_4_3_max_edge_dist(
        plate_thicknesses = platethickness_w, f_y=column_fy, corrosive_influences=corrosive_influences)

    edge_dist_min_f = end_dist_min_f
    edge_dist_max_f= end_dist_max_w

    end_dist_f = round(end_dist_min_f + 5)
    end_dist_w = round(end_dist_min_w + 5)
    edge_dist_f = end_dist_f
    edge_dist_w = end_dist_w
    flange_gauge2 = ((2 * edge_dist_f) + column_t_w + (2 * column_r1))
    #print(edge_dist_f)
    #print(flange_pitch)
    # Bolts arrrangement in flange plate #todo
    if (flange_plate_w - column_t_w - (2 * column_r1)) > ((2 * edge_dist_f) + (4 * flange_gauge) + flange_gauge2):
        number_of_column_flange = 6
    elif (flange_plate_w - column_t_w - (2 * column_r1)) > ((2 * edge_dist_f) + (2 * flange_gauge) + flange_gauge2):
        number_of_column_flange = 4
    else:
        number_of_column_flange = 2

#    print(number_of_row_flange)

        # Bolts arrrangement in web plate #todo
    if web_plate_w > ((2 * edge_dist_w) + (3 * web_gauge)):
        number_of_column_web = 4
    elif web_plate_w > ((2 * edge_dist_w) + (2 * web_gauge)):
        number_of_column_web = 3
    else:
        number_of_column_web = 2


    #######################################################################
    ###### Calculate bolt capacities ###
    # Calculation of kb for flange
    kbChk1 = end_dist_min_f / float(3 * dia_hole)
    kbChk2 = pitch_dist_min_f / float(3 * dia_hole) - 0.25
    kbChk3 = bolt_fu / float(column_fu)
    kbChk4 = 1
    kb = float(min(kbChk1, kbChk2, kbChk3, kbChk4))
    kb = round(kb, 3)

    # Bolt capacity calculation for flange splice
    if flange_plate_preference == "Outside":
        flange_t_thinner = min(column_f_t, flange_plate_t)
    else:
        flange_t_thinner = min(column_f_t, (2 * flange_plate_t))

    flange_bolt_planes = 1
    number_of_bolts = 1

    if bolt_type == "Bearing Bolt":
        flange_bolt_shear_capacity = round(ConnectionCalculations.bolt_shear( bolt_diameter = bolt_diameter,number_of_bolts = flange_bolt_planes,  bolt_fu = bolt_fu), 2)
        flange_bolt_bearing_capacity = round(ConnectionCalculations.bolt_bearing(bolt_diameter, number_of_bolts =1 ,thickness_plate= flange_t_thinner, \
                                                                                 k_b= kb, plate_fu =flange_plate_fu), 2)
        flange_bolt_capacity = min(flange_bolt_shear_capacity,flange_bolt_bearing_capacity)

    elif bolt_type == "Friction Grip Bolt":
        muf = mu_f
        bolt_hole_type = dp_bolt_hole_type  # 1 for standard, 0 for oversize hole
        n_e = 2  # number of effective surfaces offering frictional resistance
        flange_bolt_shear_capacity = round(ConnectionCalculations.bolt_shear_friction_grip_bolt(bolt_diameter = bolt_diameter, bolt_fu = bolt_fu, mu_f =mu_f, n_e=n_e, bolt_hole_type =dp_bolt_hole_type), 2)
        flange_bolt_bearing_capacity = 'N/A'
        flange_bolt_capacity = flange_bolt_shear_capacity

            #print(flange_bolt_bearing_capacity, flange_bolt_shear_capacity, flange_bolt_capacity)
    else:
        pass

    # Bolt capacity calculation for web splice
    web_bolt_planes = 2
    number_of_bolts = 1

    if bolt_type == "Bearing Bolt":
        web_bolt_shear_capacity = round(
            ConnectionCalculations.bolt_shear(bolt_diameter,  number_of_bolts ,bolt_fu = bolt_fu), 2)
        web_bolt_bearing_capacity = round(
            ConnectionCalculations.bolt_bearing(bolt_diameter, number_of_bolts, thickness_plate=web_t_thinner, \
                                                k_b = kb,plate_fu = web_plate_fu), 2)
        web_bolt_capacity = min(web_bolt_shear_capacity, web_bolt_bearing_capacity)

    elif bolt_type == "Friction Grip Bolt":
        #muf = mu_f
        bolt_hole_type = dp_bolt_hole_type  # 1 for standard, 0 for oversize hole
        n_e = 2  # number of effective surfaces offering frictional resistance
        web_bolt_shear_capacity = round(
            ConnectionCalculations.bolt_shear_friction_grip_bolt(bolt_diameter,bolt_fu = bolt_fu, mu_f = mu_f, n_e=2, bolt_hole_type=dp_bolt_hole_type), 2)
        web_bolt_bearing_capacity = 'N/A'
        web_bolt_capacity = flange_bolt_shear_capacity

        #print(flange_bolt_bearing_capacity, flange_bolt_shear_capacity, flange_bolt_capacity)
    else:
        pass

    ####check for long joints and large grip length in flange
    lj_flange = (number_of_column_flange - 1) * flange_gauge # lONG jOINT flange
    lj_web = (number_of_column_web  - 1) * web_gauge # long joint web
    if flange_plate_preference == "Outside + Inside": # large grip flange
        lg_flange = column_f_t + 2 * flange_plate_t
    else:
        lg_flange = column_f_t + flange_plate_t

    lg_web = column_t_w + (2* web_plate_t)
    beta_lj_flange = IS800_2007.cl_10_3_3_1_bolt_long_joint(d=bolt_diameter, l_j = lj_flange)
    beta_lg_flange =IS800_2007.cl_10_3_3_2_bolt_large_grip(d= bolt_diameter, l_g = lg_flange, l_j=0)
    beta_lj_web = IS800_2007.cl_10_3_3_1_bolt_long_joint(d=bolt_diameter, l_j = lj_web)
    beta_lg_web = IS800_2007.cl_10_3_3_2_bolt_large_grip(d=bolt_diameter, l_g=lg_web, l_j=0)

    flange_bolt_capacity_red = flange_bolt_capacity *beta_lj_flange *  beta_lg_flange
    web_bolt_capacity_red = web_bolt_capacity * beta_lj_web * beta_lg_web
    #print(flange_bolt_capacity_red, web_bolt_capacity_red)


    #####################################################################################################################
    # PACKING PLATES #TODO
    #####################################################################################################################


 # Calculation of number of bolts required for web splice plate
    if shear_load != 0:
        web_bolts_required = int(math.ceil(shear_load / web_bolt_capacity_red))
    else:
        web_bolts_required = 2

    total_web_plate_bolts = web_bolts_required * 2
    web_bolt_group_capacity = total_web_plate_bolts * web_bolt_capacity_red # Calculation of bolt group capacity for web splice plate

    # Bolts arrrangement in web plate #todo
    # if web_plate_w  > ((2 * edge_dist_w) + (3 * web_gauge)):
    #     number_of_column_web = 4
    # elif web_plate_w  > ((2 * edge_dist_w) + (2 * web_gauge)):
    #     number_of_column_web = 3
    # else:
    #     number_of_column_web = 2

    number_of_row_web = int(web_bolts_required / number_of_column_web)

    # Calculation of number of bolts required for flange splice plate

    #### 1. No. of bolts along column length; connecting each column = 1.05 * Force in flange [Reference: Annex F, Clause 10.6.1, F-2.2, Page 130]
    ff = flange_force(column_d, column_f_t, column_b, column_area, factored_axial_force, moment_load)
    if ff != 0:
        flange_bolts_required = int(math.ceil(1.05 * (ff/ ( flange_bolt_capacity_red) )))
    else:
        flange_bolts_required = 2


    # Number of bolts in even number (for design of flange splice plate)
    if flange_bolts_required % 2 == 0:
        flange_bolts_required = flange_bolts_required
    else:
        flange_bolts_required = flange_bolts_required + 1

    number_of_row_flange = int(flange_bolts_required / number_of_column_flange)
    # print(number_of_row_flange)
    total_flange_plate_bolts = int(flange_bolts_required * 4)
    flange_bolt_group_capacity = total_flange_plate_bolts * flange_bolt_capacity_red # Calculation of bolt group capacity for flange splice plate

    # # Bolts arrrangement in flange plate #todo
    # if (flange_plate_w - column_t_w - (2 * column_r1)) > ((2 * edge_dist_f) + (4 * flange_gauge) + flange_gauge2):
    #     number_of_column_flange = 6
    # elif (flange_plate_w - column_t_w - (2 * column_r1)) > ((2 * edge_dist_f) + (2 * flange_gauge) + flange_gauge2):
    #     number_of_column_flange = 4
    # else:
    #     number_of_column_flange = 2


    # ####check for long joints and large grip length in flange
    # lj_flange = (number_of_column_flange - 1) * flange_gauge # lONG jOINT flange
    # lj_web = (number_of_column_web  - 1) * web_gauge # long joint web
    # if flange_plate_preference == "Outside + Inside": # large grip flange
    #     lg_flange = column_f_t + 2 * flange_plate_t
    # else:
    #     lg_flange = column_f_t + flange_plate_t
    #
    # lg_web = column_t_w + (2* web_plate_t)
    # beta_lj_flange = IS800_2007.cl_10_3_3_1_bolt_long_joint(d=bolt_diameter, l_j = lj_flange)
    # beta_lg_flange =IS800_2007.cl_10_3_3_2_bolt_large_grip(d= bolt_diameter, l_g = lg_flange, l_j=0)
    # beta_lj_web = IS800_2007.cl_10_3_3_1_bolt_long_joint(d=bolt_diameter, l_j = lj_web)
    # beta_lg_web = IS800_2007.cl_10_3_3_2_bolt_large_grip(d=bolt_diameter, l_g=lg_web, l_j=0)
    #
    #
    #
    # flange_bolt_capacity_red = flange_bolt_capacity *beta_lj_flange *  beta_lg_flange
    # web_bolt_capacity_red = web_bolt_capacity * beta_lj_web * beta_lg_web
    #

    #print(flange_bolt_capacity_red, web_bolt_capacity_red)












    #number_of_row_flange = total_flange_plate_bolts / number_of_column_flange
    print(flange_bolts_required,web_bolts_required)
    print(total_flange_plate_bolts, total_web_plate_bolts)
    print(flange_bolt_group_capacity,web_bolt_group_capacity)

    #####################################################################################################################
    #####Calculation of resultant force on bolts in flange
    if ff > flange_bolt_group_capacity:
        design_status = False
        logger.error(": flange_bolt_group_capacity  is less than force in flange ff")
        logger.warning(": Minimum flange_bolt_capacity_red required is %2.2f kN" % (ff))
        logger.info(": Increase the size of column section")
    else:
        pass
    ###outside flange plate thickness
    min_thk_flange_plate = flange_plate_min_t(column_t_w)
    if flange_plate_t < min_thk_flange_plate:
        design_status = False
        logger.error(": flange_plate_t is less than min_thk_flange_plate:")
        logger.warning(": Minimum flange_plate_t required is %2.2f mm" % (min_thk_flange_plate))
        logger.info(": Increase the thickness of flange splice plate")
    else:
        pass
    opt_thk_flange_plate = thk_flange_plate (column_d, column_f_t, bolt_diameter, column_area, factored_axial_force, moment_load, column_b,
                     column_fy, dia_hole)
    if flange_plate_t < opt_thk_flange_plate:
        design_status = False
        logger.error(": flange_plate_t is less than opt_thk_flange_plate:")
        logger.warning(": Minimum flange_plate_t required is %2.2f mm" % (opt_thk_flange_plate))
        logger.info(": Increase the thickness of flange splice plate")
    else:
        flange_plate_t = float(uiObj["FlangePlate"]["Thickness (mm)"])


    ###outside flange plate width
    opt_flange_plate_width = flange_plate_width(edge_dist_f, number_of_column_flange, flange_gauge)
    if flange_plate_w < opt_flange_plate_width:
        design_status = False
        logger.error(": flange_plate_w is less than opt_flange_plate_width:")
        logger.warning(": Minimum flange_plate_w required is %2.2f mm" % (opt_flange_plate_width))
        logger.info(": Increase the width of flange splice plate")
    else:
        flange_plate_w = float(uiObj["FlangePlate"]["Width (mm)"])


    max_flange_width =flange_plate_w_max(column_b)
    if flange_plate_w > max_flange_width:
        design_status = False
        logger.error(": flange_plate_w is less than opt_flange_plate_width:")
        logger.warning(": Maximum flange_plate_w required is %2.2f mm" % (column_b))
        logger.info(": Decrease the width of flange splice plate")
    else:
        flange_plate_w = float(uiObj["FlangePlate"]["Width (mm)"])

    flangeheightmin = (2 * flange_plate_w)
    if flange_plate_l < flangeheightmin:
        design_status = False
        logger.error(": flange_plate_l is less than flangeheightmin:")
        logger.warning(": Maximum flange_plate_l required is %2.2f mm" % (flangeheightmin))
        logger.info(": Increase the height of flange splice plate")
    else:
        pass

    #number_of_row_flange = 1 #todo
    opt_flange_plate_l =  flange_plate_l_max(end_dist_f, number_of_row_flange, flange_pitch)
    # print(opt_flange_plate_l)
    if flange_plate_l < opt_flange_plate_l:
        design_status = False
        logger.error(": flange_plate_l is less than opt_flange_plate_l:")
        logger.warning(": Maximum flange_plate_l required is %2.2f mm" % (opt_flange_plate_l))
        logger.info(": Increase the height of flange splice plate")
    else:
        flange_plate_l = str(uiObj["FlangePlate"]["Height (mm)"])


    ###inner flange plate  thickness,w,l#todo
    opt_inner_thk_flange_plate = flange_plate_t

    #####################################################################################################################

    # Check for web plate thickness
    minwebt = web_min_t( shear_load, column_d, gamma_m0, column_fy)
    if web_plate_t < minwebt:
        design_status = False
        logger.error(": Chosen web splice plate thickness is not sufficient")
        logger.warning(": Minimum required thickness of web splice plate is %2.2f mm" % (minwebt))
        logger.info(": Increase the thickness of web splice plate")
    else:
        pass

    maxwebt = web_max_t(bolt_diameter)
    if web_plate_t > maxwebt:
        design_status = False
        logger.error(": Chosen web splice plate thickness is not sufficient")
        logger.warning(": Maximum required thickness of web splice plate is in between %2.2f mm" % (maxwebt))
        logger.info(": Increase the thickness of web splice plate")
    else:
        pass
    # Check for web plate width
    # Web splice plate width input and check for maximum and minimum values
    if flange_plate_preference == "Outside":
        webmaxw = web_max_w(column_d, column_f_t, column_r1)
    elif flange_plate_preference == "Outside + Inside":
        webmaxw = web_max_w(column_d, column_f_t, column_r1)- flange_plate_t



    if web_plate_w > webmaxw:
        design_status = False
        logger.error(": Width of web splice plate is greater than the clear depth of column")
        logger.warning(": Maximum required web splice plate width  is %2.2f mm" % (webmaxw))
        logger.info(": Reduce the width of web splice plate")
    else:
        pass

    webminw = web_min_w(column_d)
    if web_plate_w < webminw:
        design_status = False
        logger.error(": Width of web splice plate is less than minimum web width required")
        logger.warning(": Minimum required web splice plate width  is %2.2f mm" % (webminw))
        logger.info(": Increase the width of web splice plate")
    else:
        pass

    # Check for web plate height
    # Web splice plate height input and check for maximum and minimum values
    webplatelmin= web_plate_l_min(end_dist_w, number_of_row_web, web_pitch)
    if web_plate_l < webplatelmin:
        design_status = False
        logger.error(": web_plate_l is less than web_plate_l_min:")
        logger.warning(": Minimum web_plate_l required is %2.2f mm" % (webplatelmin))
        logger.info(": Increase the height of web splice plate")
    else:
        pass

    webplatelmax =  web_plate_l_max(end_dist_max_w, number_of_row_web, pitch_dist_max_w)
    if web_plate_l > webplatelmax:
        design_status = False
        logger.error(": web_plate_l is greater than web_plate_l_max:")
        logger.warning(": Maximun web_plate_l required is %2.2f mm" % (webplatelmax))
        logger.info(": Reduce the height of web splice plate")
    else:
        pass
    #####################################################################################################################

    ####Calculation of resultant force on bolts in web

    extreme_bolt_dist_z = number_of_column_web * web_gauge / 2 + web_gauge / 2
    extreme_bolt_dist_y = number_of_row_web * web_pitch / 2 + web_pitch / 2
    for i in np.arange(0, number_of_column_web / 2):
        sigma_z = (number_of_row_web - (2 * i - 1)) * web_gauge / 2
        sigma_y = (number_of_column_web - (2 * i - 1)) * web_pitch / 2
        ecc = ((number_of_row_web * web_pitch) + (end_dist_w)) / 2
        r_i = math.sqrt(sigma_z ** 2 + sigma_y ** 2)
        Z_w = float((column_t_w  * (column_d - 2 * (column_f_t))**2 )/4)
        mu_wz = moment_load  * Z_w / column_Zz
        mu_wy = moment_load  * Z_w / column_Zy

        hor_shear_force_bolts = (( mu_wz + shear_load * ecc) * extreme_bolt_dist_y) / r_i  # horizontal shear force acting on each bolt due to moment developed by eccentricity.
        # extreme bolt distance in X direction
        ver_shear_force_bolts = ((mu_wy + shear_load * ecc) * extreme_bolt_dist_z) / r_i  # vertical shear force acting on each bolt due to moment developed by eccentricity
        # extreme bolt distance in Y direction
        hor_force = shear_load / web_bolts_required
        area_web = (column_d - 2 * column_f_t) * column_t_w
        axial_force_w = int((column_d - 2 * (column_f_t)) * column_t_w * factored_axial_force) / column_area
        #axial_force_web = (factored_axial_force * area_web) / column_area  # horizontal force acting on each bolt (assuming uniform shear distribution)
        ver_force = axial_force_w / web_bolts_required  # vertical force acting on each bolt (assuming uniform axial distribution)
        shearresbolt = math.sqrt((hor_shear_force_bolts + hor_force) ** 2 + (ver_shear_force_bolts + ver_force) ** 2)

    if shearresbolt > web_bolt_capacity:
        design_status = False
        logger.error(": Number of bolts is not sufficient")
        logger.warning(": shear_res_bolt should be less than web_bolt_capacity of web")
        logger.info(": Increase number of bolts and spacing between bolts")
    else:
        pass
####################################################################################################################
    ################################ CAPACITY CHECK #####################################################################################

    ####Capacity of flange cover plate for bolted Outside #
    if flange_plate_preference == "Outside":
        net_eff_area_fp = (column_b - number_of_column_flange * dia_hole) * flange_plate_t
        gross_area_fp = flange_plate_t  * flange_plate_w
        Tdg_flange_plate =IS800_2007.tension_member_design_due_to_yielding_of_gross_section(A_g = gross_area_fp,
                                                                                            F_y = flange_plate_fy)
        Tdn_flange_plate = IS800_2007.tension_member_design_due_to_rupture_of_critical_section(A_n = net_eff_area_fp ,
                                                                                            F_u = flange_plate_fu)

        #  Block shear strength for flange
        Avg = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch) * column_f_t
        Avn = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch - (number_of_row_flange - 0.5) * dia_hole) * column_f_t
        Atg = (column_b - (number_of_column_flange - 1) * flange_gauge) * column_f_t
        Atn = (column_b - ((number_of_column_flange - 1) * flange_gauge) - (
                    number_of_column_flange - 1) * dia_hole) * column_f_t
        Tdb_flange_shear = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                    f_u=flange_plate_fu, f_y=flange_plate_fy)

        if Tdb_flange_shear < axial_force:
            design_status = False
        else:
            pass
        ##Block shear strength for flange cover plate

        Avg = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch) * flange_plate_t
        Avn = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch - (
                    number_of_row_flange - 0.5) * dia_hole) * flange_plate_t
        Atg = ((number_of_column_flange - 1) * flange_gauge) * flange_plate_t
        Atn = (((number_of_column_flange - 1) * flange_gauge) - (number_of_column_flange - 1) * dia_hole) * flange_plate_t

        Tdb_flange_plate = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                    f_u=flange_plate_fu, f_y=flange_plate_fy)
        if Tdb_flange_plate < factored_axial_force:
            design_status = False
        else:
            pass

        Td_flange_plate = min(Tdg_flange_plate ,Tdn_flange_plate ,Tdb_flange_plate, Tdb_flange_shear)

        ####Capacity of flange cover plate for bolted Outside + Inside #

    elif flange_plate_preference == "Outside + Inside" :
        net_eff_area_fp = (column_b - number_of_column_flange * dia_hole) * flange_plate_t
        gross_area_fp = flange_plate_t * flange_plate_w
        net_eff_area_ifp = (((column_b - number_of_column_flange * dia_hole) * flange_plate_t) - column_t_w - 2*column_r1) /2
        gross_area_ifp = flange_plate_t * (flange_plate_w - column_t_w- 2*column_r1)

        Tdg_flange_plate = IS800_2007.tension_member_design_due_to_yielding_of_gross_section(A_g=gross_area_fp,
                                                                                             F_y=flange_plate_fy)
        Tdn_flange_plate = IS800_2007.tension_member_design_due_to_rupture_of_critical_section(A_n=net_eff_area_fp,
                                                                                            F_u=flange_plate_fu)
        Tdg_flange_plate_i = IS800_2007.tension_member_design_due_to_yielding_of_gross_section(A_g=gross_area_ifp,
                                                                                             F_y=flange_plate_fy)
        Tdn_flange_plate_i = IS800_2007.tension_member_design_due_to_rupture_of_critical_section(A_n=net_eff_area_ifp,
                                                                                            F_u=flange_plate_fu)
        #  Block shear strength for flange under shear
        Avg = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch) * column_f_t
        Avn = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch - (number_of_row_flange - 0.5) * dia_hole) * column_f_t
        Atg = (column_b - (number_of_column_flange - 1) * flange_gauge) * column_f_t
        Atn = (column_b - ((number_of_column_flange - 1) * flange_gauge) - (
                number_of_column_flange - 1) * dia_hole) * column_f_t
        Tdb_flange_shear = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                    f_u=flange_plate_fu, f_y=flange_plate_fy)
        if Tdb_flange_shear < factored_axial_force:
            design_status = False
        else:
            pass
        #  Block shear strength for  outside flange plate
        Avg = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch) * flange_plate_t
        Avn = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch - (
                    number_of_row_flange - 0.5) * dia_hole) * flange_plate_t
        Atg = ((number_of_column_flange - 1) * flange_gauge) * flange_plate_t
        Atn = (((number_of_column_flange - 1) * flange_gauge) - (number_of_column_flange - 1) * dia_hole) * flange_plate_t
        Tdb_flange_plate_outside = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                    f_u=flange_plate_fu, f_y=flange_plate_fy)
        if Tdb_flange_plate_outside < factored_axial_force:
            design_status = False
        else:
            pass

        Avg = (end_dist_f + (number_of_row_flange - 1) * flange_pitch) * flange_plate_t
        Avn = (end_dist_f + (number_of_row_flange - 1) * flange_pitch - (number_of_row_flange - 0.5) * dia_hole) * flange_plate_t
        Atg = (edge_dist_f + (number_of_column_flange - 1) * flange_gauge) * flange_plate_t
        Atn = (edge_dist_f + (number_of_column_flange - 1) * flange_gauge - (
                    number_of_column_flange - 0.5) * dia_hole) * flange_plate_t
        Tdb_flange_plate_inside_shear = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                    f_u=flange_plate_fu, f_y=flange_plate_fy)
        if Tdb_flange_plate_inside_shear < factored_axial_force:
            design_status = False
        else:
            pass

        Avg = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch) * flange_plate_t
        Avn = 2 * (end_dist_f + (number_of_row_flange - 1) * flange_pitch- (
                    number_of_row_flange - 0.5) * dia_hole) * flange_plate_t
        Atg = ((number_of_column_flange - 1) * flange_gauge) * flange_plate_t
        Atn = (((number_of_column_flange - 1) * flange_gauge) - (number_of_column_flange - 1) * dia_hole) * flange_plate_t
        Tdb_flange_plate_inside_axial = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                                 f_u=flange_plate_fu,
                                                                                 f_y=flange_plate_fy)
        if Tdb_flange_plate_inside_axial < factored_axial_force:
            design_status = False
        else:
            pass

        Td_flange_plate = min(Tdg_flange_plate, Tdn_flange_plate, Tdg_flange_plate_i,Tdn_flange_plate_i ,
                Tdb_flange_shear, Tdb_flange_plate_outside,Tdb_flange_plate_inside_shear,Tdb_flange_plate_inside_axial)

    else:
        pass

    #print(Td_flange_plate)
    if ff > Td_flange_plate:
        design_status = False
        logger.error(": flange_force ff is greater than design strength of flange cover plate")
        logger.warning(": Maximum flange_force required is %2.2f kN" % (Tdb_flange_plate))
        logger.info(": decrease the size of column section")
    else:
        pass

    ####Capacity of web #
    ##Capacity of web cover plate for bolted
    if flange_plate_preference == "Outside":
        net_eff_area_wp = (column_d - column_f_t - number_of_column_web * dia_hole) * web_plate_t
        gross_area_wp = web_plate_t * web_plate_w
        Tdg_web_plate = IS800_2007.tension_member_design_due_to_yielding_of_gross_section(A_g=gross_area_wp,
                                                                                             F_y=web_plate_fy)
        Tdn_web_plate = IS800_2007.tension_member_design_due_to_rupture_of_critical_section(A_n=net_eff_area_wp,
                                                                                               F_u=web_plate_fu)
        Avg = 2 * (end_dist_w + (number_of_row_web - 1) * web_pitch) * column_t_w
        Avn = 2 * (end_dist_w + (number_of_row_web - 1) * web_pitch - (number_of_column_web - 0.5) * dia_hole) * column_t_w
        Atg = (column_d - column_f_t - ((number_of_column_web - 1) * web_gauge)) * column_t_w
        Atn = (column_d - column_f_t - ((number_of_column_web - 1) * web_gauge) - (
                    number_of_column_web - 0.5) * dia_hole) * column_t_w
        Tdb_web_axial = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                    f_u=web_plate_fu, f_y=web_plate_fy)

        Avg = (edge_dist_w + (number_of_column_web - 1) * web_gauge) * web_plate_t
        Avn = (edge_dist_w + (number_of_column_web - 1) * web_gauge - (
                    number_of_column_web - 0.5) * dia_hole) * web_plate_t
        Atg = ((number_of_row_web - 1) * web_pitch + end_dist_w) * web_plate_t
        Atn = ((number_of_row_web - 1) * web_pitch + (number_of_row_web - 1) * dia_hole + end_dist_w) * web_plate_t
        Tdb_web_shear = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                 f_u=web_plate_fu, f_y=web_plate_fy)
        Td_web_plate = min(Tdg_web_plate, Tdn_web_plate, Tdb_web_axial, Tdb_web_shear )

    elif flange_plate_preference == "Outside + Inside":
        net_eff_area_wp = (column_d - column_f_t - 2*column_r1 - number_of_column_web * dia_hole) * web_plate_t
        gross_area_wp = web_plate_t * web_plate_w
        net_eff_area_iwp = (column_d - column_f_t - flange_plate_t- number_of_column_web * dia_hole) * web_plate_t
        gross_area_iwp = web_plate_t * web_plate_w
        Tdg_web_plate = IS800_2007.tension_member_design_due_to_yielding_of_gross_section(A_g=gross_area_wp,
                                                                                          F_y=web_plate_fy)
        Tdn_web_plate = IS800_2007.tension_member_design_due_to_rupture_of_critical_section(A_n=net_eff_area_wp,
                                                                                            F_u=web_plate_fu)

        Tdg_web_plate_i = IS800_2007.tension_member_design_due_to_yielding_of_gross_section(A_g=gross_area_iwp,
                                                                                          F_y=web_plate_fy)
        Tdn_web_plate_i = IS800_2007.tension_member_design_due_to_rupture_of_critical_section(A_n=net_eff_area_iwp,
                                                                                            F_u=web_plate_fu)

        #web shear when flange plate is provided outside and inside
        Avg = 2 * (end_dist_w + (number_of_row_web - 1) * web_pitch) * column_t_w
        Avn = 2 * (end_dist_w + (number_of_row_web - 1) * web_pitch - (
                    number_of_column_web - 0.5) * dia_hole) * column_t_w
        Atg = (column_d - column_f_t - ((number_of_column_web - 1) * web_gauge)) * column_t_w
        Atn = (column_d - column_f_t - ((number_of_column_web - 1) * web_gauge) - (
                number_of_column_web - 0.5) * dia_hole) * column_t_w
        Tdb_web_axial = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                 f_u=web_plate_fu, f_y=web_plate_fy)

        Avg = (edge_dist_w + (number_of_column_web - 1) * web_gauge) * web_plate_t
        Avn = (edge_dist_w + (number_of_column_web - 1) * web_gauge - (
                number_of_column_web - 0.5) * dia_hole) * web_plate_t
        Atg = ((number_of_row_web - 1) * web_pitch + end_dist_f) * web_plate_t
        Atn = ((number_of_row_web - 1) * web_pitch + (number_of_row_web - 1) * dia_hole + end_dist_f) * web_plate_t
        Tdb_web_shear = IS800_2007.cl_6_4_1_block_shear_strength(A_vg=Avg, A_vn=Avn, A_tg=Atg, A_tn=Atn,
                                                                 f_u=web_plate_fu, f_y=web_plate_fy)

        Td_web_plate = min(Tdg_web_plate, Tdn_web_plate,Tdg_web_plate_i,Tdn_web_plate_i, Tdb_web_axial, Tdb_web_shear)
    else:
        pass

    wf = web_force(column_d, column_f_t, column_t_w, factored_axial_force, column_area)
    #print(wf)
    if wf > Td_web_plate:
        design_status = False
        logger.error(": web_force wf is greater than design strength of web cover plate")
        logger.warning(": Maximum web_force required is %2.2f kN" % (Td_web_plate))
        logger.info(": decrease the size of column section")
    else:
        pass
    #print(Td_web_plate)
    #####################################################################################################################
    # End of calculation for design of web splice plate
    # Output
    outputObj = {}
    outputObj["Bolt"] = {}
    outputObj["Bolt"]["status"] = design_status

    outputObj["WebBolt"] = {}
    outputObj["WebBolt"]["ShearCapacity"] = web_bolt_shear_capacity
    outputObj["WebBolt"]["BearingCapacity"] =  web_bolt_bearing_capacity
    outputObj["WebBolt"]["CapacityBolt"] = web_bolt_capacity_red
    outputObj["WebBolt"]["BoltsRequired"] = web_bolts_required
    outputObj["WebBolt"]["TotalBoltsRequired"] = total_web_plate_bolts

    outputObj["WebBolt"]["Pitch"] = web_pitch
    outputObj["WebBolt"]["End"] = end_dist_w
    outputObj["WebBolt"]["Edge"] = edge_dist_w
    outputObj["WebBolt"]["WebPlateHeight"] = web_plate_l
    outputObj["WebBolt"]["WebGauge"] = web_gauge
    outputObj["WebBolt"]["WebGaugeMax"] =  gauge_dist_max_w
    outputObj["WebBolt"]["webPlateDemand"] = web_bolt_group_capacity  ## capacity of web group bolt
    outputObj["WebBolt"]["WebPlateWidth"] = web_plate_w
    outputObj["WebBolt"]["WebPlateCapacity"] = Td_web_plate

    print(web_bolt_shear_capacity)
    print(web_bolt_bearing_capacity)
    print(web_bolt_capacity_red)
    print(web_bolts_required)
    print(total_web_plate_bolts)
    print(web_pitch)
    print(end_dist_w)
    print(edge_dist_w)
    print(web_plate_l)
    print(web_gauge)
    print(gauge_dist_max_w)
    print(web_bolt_group_capacity)
    print(web_plate_w)
    print(Td_web_plate)
    # End of calculation for design of web splice plate

    outputObj["FlangeBolt"] = {}
    outputObj["FlangeBolt"]["ShearCapacityF"] = flange_bolt_shear_capacity
    outputObj["FlangeBolt"]["BearingCapacityF"] =  flange_bolt_bearing_capacity
    outputObj["FlangeBolt"]["CapacityBoltF"] = flange_bolt_capacity_red
    outputObj["FlangeBolt"]["BoltsRequiredF"] = flange_bolts_required
    # Note: This outputs number of bolts required in one side of splice
    outputObj["FlangeBolt"]["TotalBoltsRequiredF"] = total_flange_plate_bolts
    outputObj["FlangeBolt"]["NumberBoltColFlange"] = number_of_row_flange
    outputObj["FlangeBolt"]["PitchF"] = flange_pitch
    outputObj["FlangeBolt"]["EndF"] =  end_dist_f
    outputObj["FlangeBolt"]["EdgeF"] =  edge_dist_f
    outputObj["FlangeBolt"]["FlangePlateHeight"] = flange_plate_l
    outputObj["FlangeBolt"]["FlangePlateWidth"] = flange_plate_w
    outputObj["FlangeBolt"]["ThicknessFlangePlate"] = flange_plate_t
    outputObj["FlangeBolt"]["FlangeGauge"] = flange_gauge
    outputObj["FlangeBolt"]["FlangePlateDemand"] = flange_bolt_group_capacity##### capacity of flange group bolt
    outputObj["FlangeBolt"]["FlangeCapacity"] = Td_flange_plate
    #outputObj["FlangeBolt"]["edge_dist_gauge"] = edge_dist  # For 3D model


    if design_status == True:

        logger.info(": Overall bolted cover plate splice connection design is safe \n")
        logger.debug(" :=========End Of design===========")
    else:
        logger.error(": Design is not safe \n ")
        logger.debug(" :=========End Of design===========")

    return  outputObj