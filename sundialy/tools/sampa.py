"""
SAMPA (Solar and Moon Position Algorithm) (https://midcdmz.nrel.gov/sampa)

This algorithm calculates the solar and lunar zenith and azimuth angles in the period from the year -2000 to 6000, with
uncertainties of +/- 0.0003 degrees for the Sun and +/- 0.003 degrees for the Moon, based on the date, time, and
location on Earth. The algorithm can be used for solar eclipse monitoring and estimating the reduction in solar
irradiance for many applications, such as smart grid, solar energy, etc.

Reda, I. (2010). Solar Eclipse Monitoring for Solar Energy Applications Using the Solar and Moon Position Algorithms.
35 pp.; NREL Report No. TP-3B0-47681.:
http://www.nrel.gov/docs/fy10osti/47681.pdf
"""

import math
from . import spa
from . import bird
from .constants import DT
from typing import Union, Optional, List, Tuple

NUMBER_TYPE = Union[int, float]
INT_NONE_TYPE = Optional[int]
NUMBER_NONE_TYPE = Union[int, float, None]
DATA_TYPE = List[List[NUMBER_TYPE]]
DATA_INT_TYPE = List[List[int]]
DATA_INT_MISSING_TYPE = List[List[INT_NONE_TYPE]]
DATA_MISSING_TYPE = List[List[NUMBER_NONE_TYPE]]
SAMPA_RETURN_TYPE = Tuple[Tuple[float, float, float, float, float, float, float, float],
                          Tuple[float, float, NUMBER_TYPE, float, Tuple[NUMBER_TYPE, NUMBER_TYPE, NUMBER_TYPE,
                                                                        NUMBER_TYPE, NUMBER_TYPE, NUMBER_TYPE], str]]


def _sun_distance(year: int, month: int, day: NUMBER_TYPE, hour: NUMBER_TYPE, minute: NUMBER_TYPE, second: NUMBER_TYPE,
                  microsecond: NUMBER_TYPE, dt: NUMBER_TYPE = DT) -> float:
    """
    Calculates the Sun's distance from the center of the Earth in astronomical units for use inside the sampa function.

    :param year: Year
    :param month: Month
    :param day: Day
    :param hour: Hour
    :param minute: Minute
    :param second: Second
    :param microsecond: Microsecond
    :param dt: The difference between the Earth rotation time and the Terrestrial Time (TT)
    :return: The Sun's distance from the center of the Earth in astronomical units
    """

    earth_periodic_terms_r0: DATA_TYPE = [[100013989, 0, 0], [1670700, 3.0984635, 6283.07585],
                                          [13956, 3.05525, 12566.1517], [3084, 5.1985, 77713.7715],
                                          [1628, 1.1739, 5753.3849], [1576, 2.8469, 7860.4194],
                                          [925, 5.453, 11506.77], [542, 4.564, 3930.21],
                                          [472, 3.661, 5884.927], [346, 0.964, 5507.553],
                                          [329, 5.9, 5223.694], [307, 0.299, 5573.143],
                                          [243, 4.273, 11790.629], [212, 5.847, 1577.344],
                                          [186, 5.022, 10977.079], [175, 3.012, 18849.228],
                                          [110, 5.055, 5486.778], [98, 0.89, 6069.78],
                                          [86, 5.69, 15720.84], [86, 1.27, 161000.69],
                                          [65, 0.27, 17260.15], [63, 0.92, 529.69],
                                          [57, 2.01, 83996.85], [56, 5.24, 71430.7],
                                          [49, 3.25, 2544.31], [47, 2.58, 775.52],
                                          [45, 5.54, 9437.76], [43, 6.01, 6275.96],
                                          [39, 5.36, 4694], [38, 2.39, 8827.39],
                                          [37, 0.83, 19651.05], [37, 4.9, 12139.55],
                                          [36, 1.67, 12036.46], [35, 1.84, 2942.46],
                                          [33, 0.24, 7084.9], [32, 0.18, 5088.63],
                                          [32, 1.78, 398.15], [28, 1.21, 6286.6],
                                          [28, 1.9, 6279.55], [26, 4.59, 10447.39]]

    earth_periodic_terms_r1: DATA_TYPE = [[103019, 1.10749, 6283.07585], [1721, 1.0644, 12566.1517], [702, 3.142, 0],
                                          [32, 1.02, 18849.23], [31, 2.84, 5507.55], [25, 1.32, 5223.69],
                                          [18, 1.42, 1577.34], [10, 5.91, 10977.08], [9, 1.42, 6275.96],
                                          [9, 0.27, 5486.78]]

    earth_periodic_terms_r2: DATA_TYPE = [[4359, 5.7846, 6283.0758], [124, 5.579, 12566.152], [12, 3.14, 0],
                                          [9, 3.63, 77713.77], [6, 1.87, 5573.14], [3, 5.47, 18849.23]]

    earth_periodic_terms_r3: DATA_TYPE = [[145, 4.273, 6283.076], [7, 3.92, 12566.15]]

    earth_periodic_terms_r4: DATA_TYPE = [[4, 2.56, 6283.08]]

    ept_r0 = earth_periodic_terms_r0
    ept_r1 = earth_periodic_terms_r1
    ept_r2 = earth_periodic_terms_r2
    ept_r3 = earth_periodic_terms_r3
    ept_r4 = earth_periodic_terms_r4

    # b_for_jd is 0 for julian calendar and (2 - a + int(a / 4)) for gregorian calendar
    if month <= 2:
        year -= 1
        month += 12
    second += microsecond / 1000000
    minute += second / 60
    hour += minute / 60
    day += hour / 24
    a = int(year / 100)
    b_for_jd = 2 - a + int(a / 4)
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + b_for_jd - 1524.5
    jde = jd + (dt / 86400)
    jce = (jde - 2451545) / 36525
    jme = jce / 10

    r0_i = []
    for i in range(len(ept_r0)):
        r0_i.append(ept_r0[i][0] * math.cos(ept_r0[i][1] + ept_r0[i][2] * jme))
    r0 = sum(r0_i)
    r1_i = []
    for i in range(len(ept_r1)):
        r1_i.append(ept_r1[i][0] * math.cos(ept_r1[i][1] + ept_r1[i][2] * jme))
    r1 = sum(r1_i)
    r2_i = []
    for i in range(len(ept_r2)):
        r2_i.append(ept_r2[i][0] * math.cos(ept_r2[i][1] + ept_r2[i][2] * jme))
    r2 = sum(r2_i)
    r3_i = []
    for i in range(len(ept_r3)):
        r3_i.append(ept_r3[i][0] * math.cos(ept_r3[i][1] + ept_r3[i][2] * jme))
    r3 = sum(r3_i)
    r4_i = []
    for i in range(len(ept_r4)):
        r4_i.append(ept_r4[i][0] * math.cos(ept_r4[i][1] + ept_r4[i][2] * jme))
    r4 = sum(r4_i)
    r = (r0 + r1 * jme + r2 * jme ** 2 + r3 * jme ** 3 + r4 * jme ** 4) / 100000000
    return r


def sampa(year: int, month: int, day: NUMBER_TYPE, hour: NUMBER_TYPE, minute: NUMBER_TYPE, second: NUMBER_TYPE,
          microsecond: NUMBER_TYPE, latitude: NUMBER_TYPE, longitude: NUMBER_TYPE, elevation: NUMBER_TYPE,
          pressure: NUMBER_TYPE, temperature: NUMBER_TYPE, ozone: NUMBER_TYPE = 0.3, water: NUMBER_TYPE = 1.5,
          aerosol: NUMBER_TYPE = 0.04, albedo: NUMBER_TYPE = 0.2, ba: NUMBER_TYPE = 0.85, k1: NUMBER_TYPE = 0.1,
          dt: NUMBER_TYPE = DT) -> SAMPA_RETURN_TYPE:
    """
    SAMPA (Solar and Moon Position Algorithm).

    This algorithm calculates the solar and lunar zenith and azimuth angles in the period from the year -2000 to 6000, with
    uncertainties of +/- 0.0003 degrees for the Sun and +/- 0.003 degrees for the Moon, based on the date, time, and
    location on Earth. The algorithm can be used for solar eclipse monitoring and estimating the reduction in solar
    irradiance for many applications, such as smart grid, solar energy, etc.

    :param year: Year
    :param month: Month
    :param day: Day
    :param hour: Hour
    :param minute: Minute
    :param second: Second
    :param microsecond: Microsecond
    :param latitude: Latitude
    :param longitude: Longitude
    :param elevation: Observer elevation (in meters)
    :param pressure: Annual average local pressure (in millibars)
    :param temperature: Annual average local temperature (in Celsius)
    :param ozone: Amount of ozone in a vertical column from surface (cm)
    :param water: Amount of precipitable water in a vertical column from surface (cm)
    :param aerosol: Broadband aerosol optical depth from surface in a vertical path (broadband turbidity)
    :param albedo: Ground albedo
    :param ba: Ratio of the forward-scattered irradiance to the total scattered irradiance due to aerosols
    :param k1: Constant used in Bird model associated with aerosol absorptance
    :param dt: The difference between the Earth rotation time and the Terrestrial Time (TT)
    :return: ((sun's topocentric zenith angle, sun's topocentric azimuth angle, moon's topocentric zenith angle,
        moon's topocentric azimuth angle, topocentric moon’s declination, topocentric moon's local hour angle,
        moon’s topocentric right ascension, true obliquity of the ecliptic), (radius of the sun’s disk,
        radius of the moon’s disk, area of SUL, percentage area of the SUL with respect to the area of the sun’s disk,
        (Estimated direct normal irradiance [W/m^2], Estimated direct normal irradiance from sun's unshaded lune [W/m^2],
        Estimated global horiz irradiance [W/m^2], Estimated global horiz irradiance from sun's unshaded lune [W/m^2],
        Estimated diffuse horiz irradiance [W/m^2], Estimated duffuse horiz irradiance from sun's unshaded lune [W/m^2]),
        solar eclipse notes))
    """
    # In the paper for irradiance it suggests ba=0.84 but in the c code it uses ba=0.85

    # Main Report
    # A lot of code is copied from spa.py. Some revisions were copied even if they were not in the SAMPA paper.
    original_time = (year, month, day, hour, minute, second, microsecond)

    moon_periodic_terms: DATA_INT_MISSING_TYPE = [[0, 0, 1, 0, 6288774, -20905355], [2, 0, -1, 0, 1274027, -3699111],
                                                  [2, 0, 0, 0, 658314, -2955968], [0, 0, 2, 0, 213618, -569925],
                                                  [0, 1, 0, 0, -185116, 48888], [0, 0, 0, 2, -114332, -3149],
                                                  [2, 0, -2, 0, 58793, 246158], [2, -1, -1, 0, 57066, -152138],
                                                  [2, 0, 1, 0, 53322, -170733], [2, -1, 0, 0, 45758, -204586],
                                                  [0, 1, -1, 0, -40923, -129620], [1, 0, 0, 0, -34720, 108743],
                                                  [0, 1, 1, 0, -30383, 104755], [2, 0, 0, -2, 15327, 10321],
                                                  [0, 0, 1, 2, -12528, None], [0, 0, 1, -2, 10980, 79661],
                                                  [4, 0, -1, 0, 10675, -34782], [0, 0, 3, 0, 10034, -23210],
                                                  [4, 0, -2, 0, 8548, -21636], [2, 1, -1, 0, -7888, 24208],
                                                  [2, 1, 0, 0, -6766, 30824], [1, 0, -1, 0, -5163, -8379],
                                                  [1, 1, 0, 0, 4987, -16675], [2, -1, 1, 0, 4036, -12831],
                                                  [2, 0, 2, 0, 3994, -10445], [4, 0, 0, 0, 3861, -11650],
                                                  [2, 0, -3, 0, 3665, 14403], [0, 1, -2, 0, -2689, -7003],
                                                  [2, 0, -1, 2, -2602, None], [2, -1, -2, 0, 2390, 10056],
                                                  [1, 0, 1, 0, -2348, 6322], [2, -2, 0, 0, 2236, -9884],
                                                  [0, 1, 2, 0, -2120, 5751], [0, 2, 0, 0, -2069, None],
                                                  [2, -2, -1, 0, 2048, -4950], [2, 0, 1, -2, -1773, 4130],
                                                  [2, 0, 0, 2, -1595, None], [4, -1, -1, 0, 1215, -3958],
                                                  [0, 0, 2, 2, -1110, None], [3, 0, -1, 0, -892, 3258],
                                                  [2, 1, 1, 0, -810, 2616], [4, -1, -2, 0, 759, -1897],
                                                  [0, 2, -1, 0, -713, -2117], [2, 2, -1, 0, -700, 2354],
                                                  [2, 1, -2, 0, 691, None], [2, -1, 0, -2, 596, None],
                                                  [4, 0, 1, 0, 549, -1423], [0, 0, 4, 0, 537, -1117],
                                                  [4, -1, 0, 0, 520, -1571], [1, 0, -2, 0, -487, -1739],
                                                  [2, 1, 0, -2, -399, None], [0, 0, 2, -2, -381, -4421],
                                                  [1, 1, 1, 0, 351, None], [3, 0, -2, 0, -340, None],
                                                  [4, 0, -3, 0, 330, None], [2, -1, 2, 0, 327, None],
                                                  [0, 2, 1, 0, -323, 1165], [1, 1, -1, 0, 299, None],
                                                  [2, 0, 3, 0, 294, None], [2, 0, -1, -2, None, 8752]]

    periodic_terms_for_moon_latitude: DATA_INT_TYPE = [[0, 0, 0, 1, 5128122], [0, 0, 1, 1, 280602],
                                                       [0, 0, 1, -1, 277693], [2, 0, 0, -1, 173237],
                                                       [2, 0, -1, 1, 55413], [2, 0, -1, -1, 46271],
                                                       [2, 0, 0, 1, 32573], [0, 0, 2, 1, 17198],
                                                       [2, 0, 1, -1, 9266], [0, 0, 2, -1, 8822],
                                                       [2, -1, 0, -1, 8216], [2, 0, -2, -1, 4324],
                                                       [2, 0, 1, 1, 4200], [2, 1, 0, -1, -3359],
                                                       [2, -1, -1, 1, 2463], [2, -1, 0, 1, 2211],
                                                       [2, -1, -1, -1, 2065], [0, 1, -1, -1, -1870],
                                                       [4, 0, -1, -1, 1828], [0, 1, 0, 1, -1794],
                                                       [0, 0, 0, 3, -1749], [0, 1, -1, 1, -1565],
                                                       [1, 0, 0, 1, -1491], [0, 1, 1, 1, -1475],
                                                       [0, 1, 1, -1, -1410], [0, 1, 0, -1, -1344],
                                                       [1, 0, 0, -1, -1335], [0, 0, 3, 1, 1107],
                                                       [4, 0, 0, -1, 1021], [4, 0, -1, 1, 833],
                                                       [0, 0, 1, -3, 777], [4, 0, -2, 1, 671],
                                                       [2, 0, 0, -3, 607], [2, 0, 2, -1, 596],
                                                       [2, -1, 1, -1, 491], [2, 0, -2, 1, -451],
                                                       [0, 0, 3, -1, 439], [2, 0, 2, 1, 422],
                                                       [2, 0, -3, -1, 421], [2, 1, -1, 1, -366],
                                                       [2, 1, 0, 1, -351], [4, 0, 0, 1, 331],
                                                       [2, -1, 1, 1, 315], [2, -2, 0, -1, 302],
                                                       [0, 0, 1, 3, -283], [2, 1, 1, -1, -229],
                                                       [1, 1, 0, -1, 223], [1, 1, 0, 1, 223],
                                                       [0, 1, -2, -1, -220], [2, 1, -1, -1, -220],
                                                       [1, 0, 1, 1, -185], [2, -1, -2, -1, 181],
                                                       [0, 1, 2, 1, -177], [4, 0, -2, -1, 176],
                                                       [4, -1, -1, -1, 166], [1, 0, 1, -1, -164],
                                                       [4, 0, 1, -1, 132], [1, 0, -1, -1, -119],
                                                       [4, -1, 0, -1, 115], [2, -2, 0, 1, 107]]

    mpt: DATA_INT_TYPE = []  # moon_periodic_terms
    ptml = periodic_terms_for_moon_latitude

    for item in moon_periodic_terms:
        mpt.append(list(map(lambda number: 0 if number is None else number, item)))

    coefficients_for_sin_terms: DATA_INT_TYPE = [[0, 0, 0, 0, 1], [-2, 0, 0, 2, 2], [0, 0, 0, 2, 2], [0, 0, 0, 0, 2],
                                                 [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [-2, 1, 0, 2, 2], [0, 0, 0, 2, 1],
                                                 [0, 0, 1, 2, 2], [-2, -1, 0, 2, 2], [-2, 0, 1, 0, 0], [-2, 0, 0, 2, 1],
                                                 [0, 0, -1, 2, 2], [2, 0, 0, 0, 0], [0, 0, 1, 0, 1], [2, 0, -1, 2, 2],
                                                 [0, 0, -1, 0, 1], [0, 0, 1, 2, 1], [-2, 0, 2, 0, 0], [0, 0, -2, 2, 1],
                                                 [2, 0, 0, 2, 2], [0, 0, 2, 2, 2], [0, 0, 2, 0, 0], [-2, 0, 1, 2, 2],
                                                 [0, 0, 0, 2, 0], [-2, 0, 0, 2, 0], [0, 0, -1, 2, 1], [0, 2, 0, 0, 0],
                                                 [2, 0, -1, 0, 1], [-2, 2, 0, 2, 2], [0, 1, 0, 0, 1], [-2, 0, 1, 0, 1],
                                                 [0, -1, 0, 0, 1], [0, 0, 2, -2, 0], [2, 0, -1, 2, 1], [2, 0, 1, 2, 2],
                                                 [0, 1, 0, 2, 2], [-2, 1, 1, 0, 0], [0, -1, 0, 2, 2], [2, 0, 0, 2, 1],
                                                 [2, 0, 1, 0, 0], [-2, 0, 2, 2, 2], [-2, 0, 1, 2, 1], [2, 0, -2, 0, 1],
                                                 [2, 0, 0, 0, 1], [0, -1, 1, 0, 0], [-2, -1, 0, 2, 1], [-2, 0, 0, 0, 1],
                                                 [0, 0, 2, 2, 1], [-2, 0, 2, 0, 1], [-2, 1, 0, 2, 1], [0, 0, 1, -2, 0],
                                                 [-1, 0, 1, 0, 0], [-2, 1, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 1, 2, 0],
                                                 [0, 0, -2, 2, 2], [-1, -1, 1, 0, 0], [0, 1, 1, 0, 0], [0, -1, 1, 2, 2],
                                                 [2, -1, -1, 2, 2], [0, 0, 3, 2, 2], [2, -1, 0, 2, 2]]

    coefficients_for_dy: DATA_MISSING_TYPE = [[-171996, -174.2], [-13187, -1.6], [-2274, -0.2], [2062, 0.2],
                                              [1426, -3.4], [712, 0.1], [-517, 1.2], [-386, -0.4],
                                              [-301, None], [217, -0.5], [-158, None], [129, 0.1],
                                              [123, None], [63, None], [63, 0.1], [-59, None],
                                              [-58, -0.1], [-51, None], [48, None], [46, None],
                                              [-38, None], [-31, None], [29, None], [29, None],
                                              [26, None], [-22, None], [21, None], [17, -0.1],
                                              [16, None], [-16, 0.1], [-15, None], [-13, None],
                                              [-12, None], [11, None], [-10, None], [-8, None],
                                              [7, None], [-7, None], [-7, None], [-7, None],
                                              [6, None], [6, None], [6, None], [-6, None],
                                              [-6, None], [5, None], [-5, None], [-5, None],
                                              [-5, None], [4, None], [4, None], [4, None],
                                              [-4, None], [-4, None], [-4, None], [3, None],
                                              [-3, None], [-3, None], [-3, None], [-3, None],
                                              [-3, None], [-3, None], [-3, None]]

    coefficients_for_de: DATA_MISSING_TYPE = [[92025, 8.9], [5736, -3.1], [977, -0.5], [-895, 0.5],
                                              [54, -0.1], [-7, None], [224, -0.6], [200, None],
                                              [129, -0.1], [-95, 0.3], [None, None], [-70, None],
                                              [-53, None], [None, None], [-33, None], [26, None],
                                              [32, None], [27, None], [None, None], [-24, None],
                                              [16, None], [13, None], [None, None], [-12, None],
                                              [None, None], [None, None], [-10, None], [None, None],
                                              [-8, None], [7, None], [9, None], [7, None],
                                              [6, None], [None, None], [5, None], [3, None],
                                              [-3, None], [None, None], [3, None], [3, None],
                                              [None, None], [-3, None], [-3, None], [3, None],
                                              [3, None], [None, None], [3, None], [3, None],
                                              [3, None], [None, None], [None, None], [None, None],
                                              [None, None], [None, None], [None, None], [None, None],
                                              [None, None], [None, None], [None, None], [None, None],
                                              [None, None], [None, None], [None, None]]

    cfst = coefficients_for_sin_terms
    cfdy: DATA_TYPE = []  # coefficients_for_dy
    cfde: DATA_TYPE = []  # coefficients_for_de
    # Replace None with 0
    for row in coefficients_for_dy:
        cfdy.append(list(map(lambda value: 0 if value is None else value, row)))
    for row in coefficients_for_de:
        cfde.append(list(map(lambda value: 0 if value is None else value, row)))

    # b_for_jd is 0 for julian calendar and (2 - a + int(a / 4)) for gregorian calendar
    if month <= 2:
        year -= 1
        month += 12
    second += microsecond / 1000000
    minute += second / 60
    hour += minute / 60
    day += hour / 24
    a = int(year / 100)
    b_for_jd = 2 - a + int(a / 4)
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + b_for_jd - 1524.5
    jde = jd + (dt / 86400)
    jc = (jd - 2451545) / 36525
    jce = (jde - 2451545) / 36525
    jme = jce / 10
    # T in the paper is JCE
    l_prime = 218.3164477 + 481267.88123421 * jce - 0.0015786 * jce ** 2 + jce ** 3 / 538841 - jce ** 4 / 65194000
    d = 297.8501921 + 445267.1114034 * jce - 0.0018819 * jce ** 2 + jce ** 3 / 545868 - jce ** 4 / 113065000
    m = 357.5291092 + 35999.0502909 * jce - 0.0001536 * jce ** 2 + jce ** 3 / 24490000
    m_prime = 134.9633964 + 477198.8675055 * jce + 0.0087414 * jce ** 2 + jce ** 3 / 69699 - jce ** 4 / 14712000
    f = 93.2720950 + 483202.0175233 * jce - 0.0036539 * jce ** 2 - jce ** 3 / 3526000 + jce ** 4 / 863310000

    e = 1 - 0.002516 * jce - 0.0000074 * jce ** 2
    d = d % 360
    m = m % 360
    m_prime = m_prime % 360
    f = f % 360
    l_prime = l_prime % 360

    l_terms = []
    for i in range(len(mpt)):
        multiply_l_i = 1.
        if abs(mpt[i][1]) == 1:
            multiply_l_i = e
        elif abs(mpt[i][1]) == 2:
            multiply_l_i = e ** 2
        l_terms.append(mpt[i][4] * multiply_l_i * math.sin(
            math.radians(mpt[i][0] * d + mpt[i][1] * m + mpt[i][2] * m_prime + mpt[i][3] * f)))
    l = sum(l_terms)

    r_terms = []
    for i in range(len(mpt)):
        multiply_l_i = 1.
        if abs(mpt[i][1]) == 1:
            multiply_l_i = e
        elif abs(mpt[i][1]) == 2:
            multiply_l_i = e ** 2
        r_terms.append(mpt[i][5] * multiply_l_i * math.cos(
            math.radians(mpt[i][0] * d + mpt[i][1] * m + mpt[i][2] * m_prime + mpt[i][3] * f)))
    r = sum(r_terms)

    b_terms = []
    for i in range(len(ptml)):
        multiply_l_i = 1.
        if abs(ptml[i][1]) == 1:
            multiply_l_i = e
        elif abs(ptml[i][1]) == 2:
            multiply_l_i = e ** 2
        b_terms.append(ptml[i][4] * multiply_l_i * math.sin(
            math.radians(ptml[i][0] * d + ptml[i][1] * m + ptml[i][2] * m_prime + ptml[i][3] * f)))
    b = sum(b_terms)

    a1 = 119.75 + 131.849 * jce
    a2 = 53.09 + 479264.29 * jce
    a3 = 313.45 + 481266.484 * jce

    dl = 3958 * math.sin(math.radians(a1)) + 1962 * math.sin(math.radians(l_prime - f)) + 318 * math.sin(
        math.radians(a2))
    db = -2235 * math.sin(math.radians(l_prime)) + 382 * math.sin(math.radians(a3)) + 175 * math.sin(
        math.radians(a1 - f)) + 175 * math.sin(math.radians(a1 + f)) + 127 * math.sin(
        math.radians(l_prime - m_prime)) - 115 * math.sin(math.radians(l_prime + m_prime))
    l_small_prime = l_prime + (l + dl) / 1000000
    beta = (b + db) / 1000000
    l_small_prime = l_small_prime % 360
    beta = beta % 360
    delta = 385000.56 + r / 1000
    p = math.degrees(math.asin(
        6378.14 / delta))  # In the paper they wrote asin while they usually wrote arcsin but the sampa.c code uses asin

    x = [297.85036 + 445267.111480 * jce - 0.0019142 * jce ** 2 + jce ** 3 / 189474.,
         357.52772 + 35999.050340 * jce - 0.0001603 * jce ** 2 - jce ** 3 / 300000.,
         134.96298 + 477198.867398 * jce + 0.0086972 * jce ** 2 + jce ** 3 / 56250.,
         93.27191 + 483202.017538 * jce - 0.0036825 * jce ** 2 + jce ** 3 / 327270.,
         125.04452 - 1934.136261 * jce + 0.0020708 * jce ** 2 + jce ** 3 / 450000.]

    x = list(map(math.radians, x))

    dy_i: List[float] = []
    for i in range(63):
        x_i = []
        for j in range(5):
            x_i.append(x[j] * cfst[i][j])
        dy_i.append((cfdy[i][0] + cfdy[i][1] * jce) * math.sin(sum(x_i)))
    de_i: List[float] = []
    for i in range(63):
        x_i = []
        for j in range(5):
            x_i.append(x[j] * cfst[i][j])
        de_i.append((cfde[i][0] + cfde[i][1] * jce) * math.cos(sum(x_i)))
    dy = sum(dy_i) / 36000000
    de = sum(de_i) / 36000000
    u = jme / 10
    e0 = 84381.448 - 4680.93 * u - 1.55 * u ** 2 + 1999.25 * u ** 3 - 51.38 * u ** 4 - 249.67 * u ** 5 \
        - 39.05 * u ** 6 + 7.12 * u ** 7 + 27.87 * u ** 8 + 5.79 * u ** 9 + 2.45 * u ** 10
    e = e0 / 3600 + de
    l_small = l_small_prime + dy
    n0 = 280.46061837 + 360.98564736629 * (jd - 2451545) + 0.000387933 * jc ** 2 - jc ** 3 / 38710000
    n = n0 + dy * math.cos(math.radians(e))
    n = n % 360
    a2 = math.degrees(math.atan2(
        math.sin(math.radians(l_small)) * math.cos(math.radians(e)) - math.tan(math.radians(beta)) * math.sin(
            math.radians(e)), math.cos(math.radians(l_small))))
    a2 = a2 % 360
    d = math.degrees(math.asin(
        math.sin(math.radians(beta)) * math.cos(math.radians(e)) + math.cos(math.radians(beta)) * math.sin(
            math.radians(e)) * math.sin(math.radians(l_small))))
    h = n + longitude - a2
    h = h % 360
    u2 = math.degrees(math.atan(0.99664719 * math.tan(math.radians(latitude))))
    x2 = math.cos(math.radians(u2)) + elevation / 6378140 * math.cos(math.radians(latitude))
    y = 0.99664719 * math.sin(math.radians(u2)) + elevation / 6378140 * math.sin(math.radians(latitude))
    da = math.degrees(math.atan2(-x2 * math.sin(math.radians(p)) * math.sin(math.radians(h)),
                                 math.cos(math.radians(d)) - x2 * math.sin(math.radians(p)) * math.cos(
                                     math.radians(h))))
    a_prime = a2 + da
    # On the SAMPA paper they haven't changed y to x2 as on the SPA paper where they revised it in 2008.
    # It is changed here.
    d_prime = math.degrees(
        math.atan2((math.sin(math.radians(d)) - y * math.sin(math.radians(p))) * math.cos(math.radians(da)),
                   math.cos(math.radians(d)) - x2 * math.sin(math.radians(p)) * math.cos(math.radians(h))))
    h_prime = h - da
    e0 = math.degrees(math.asin(math.sin(math.radians(latitude)) * math.sin(math.radians(d_prime)) + math.cos(
        math.radians(latitude)) * math.cos(math.radians(d_prime)) * math.cos(math.radians(h_prime))))
    # In the SPA paper de2 = 0 when the sun is below the horizon.
    # The SAMPA paper doesn't mention this but the code seems to be using it.
    if e0 >= -1 * (0.26667 + 0.5667):
        de2 = (pressure / 1010) * (283 / (273 + temperature)) * (
                1.02 / (60 * math.tan(math.radians(e0 + 10.3 / (e0 + 5.11)))))
    else:
        de2 = 0
    e = e0 + de2
    theta_m = 90 - e
    gamma = math.degrees(math.atan2(math.sin(math.radians(h_prime)),
                                    math.cos(math.radians(h_prime)) * math.sin(math.radians(latitude)) - math.tan(
                                        math.radians(d_prime)) * math.cos(math.radians(latitude))))
    gamma = gamma % 360
    f_m = gamma + 180
    f_m = f_m % 360

    # Section 5, 6 and Appendix A.2
    year, month, day, hour, minute, second, microsecond = original_time
    spa_results = spa.spa(year, month, day, hour, minute, second, microsecond, latitude, longitude, elevation, pressure,
                          temperature, 0, 0, dt=dt)
    theta_s = spa_results[0][1]
    f_s = spa_results[0][2]
    ems = math.degrees(math.acos(
        math.cos(math.radians(theta_s)) * math.cos(math.radians(theta_m)) + math.sin(math.radians(theta_s)) * math.sin(
            math.radians(theta_m)) * math.cos(math.radians(f_s - f_m))))
    r_capital_s = _sun_distance(year, month, day, hour, minute, second, microsecond, dt=dt)
    r_s = 959.63 / (3600 * r_capital_s)
    r_m = (358473400 * (1 + math.sin(math.radians(e)) * math.sin(math.radians(p)))) / (3600 * delta)
    if ems > (r_m + r_s):  # No eclipse
        solar_eclipse = 'No Eclipse'
    elif ems == (r_m + r_s):  # Start and End of Eclipse
        solar_eclipse = 'Start or End of Eclipse'
    else:  # elif ems < (r_m + r_s):  # Start and End of Eclipse
        solar_eclipse = 'Solar Eclipse'

    if ems < r_m + r_s:
        if ems <= abs(r_m - r_s):
            solar_eclipse = 'Total Solar Eclipse'
            a_i = math.pi * r_m ** 2
        else:
            solar_eclipse = 'Partial Solar Eclipse'
            s = (ems ** 2 + r_s ** 2 - r_m ** 2) / (2 * ems)
            m = (ems ** 2 - r_s ** 2 + r_m ** 2) / (2 * ems)
            h = math.sqrt(4 * ems ** 2 * r_s ** 2 - (ems ** 2 + r_s ** 2 - r_m ** 2) ** 2) / (2 * ems)
            t_s = h * s
            t_m = h * m
            a_s = r_s ** 2 * math.acos(s / r_s)
            a_capital_1 = a_s - t_s
            a_m = r_m ** 2 * math.acos(m / r_m)
            a_capital_2 = a_m - t_m
            a_i = a_capital_1 + a_capital_2
    else:
        a_i = 0
    a_sul = math.pi * r_s ** 2 - a_i
    if a_sul < 0:
        a_sul = 0
    a_sul_percent = (a_sul * 100) / (math.pi * r_s ** 2)
    irradiance = bird.bird(r_capital_s, theta_s, pressure, ozone, water, aerosol, albedo, a_sul_percent / 100, ba, k1)
    I_e = irradiance[1], irradiance[4], irradiance[2], irradiance[5], irradiance[3], irradiance[6]
    # This code is not needed because this is calculated in the bird model by using the dni_mod argument
    # for irradiance in irradiances:
    # for i in range(1, 7):
    #     I_e.append((irradiance[i] * a_sul_percent) / 100)
    return (theta_s, f_s, theta_m, f_m, d_prime, h_prime, a_prime, e), (r_s, r_m, a_sul, a_sul_percent, I_e, solar_eclipse)
