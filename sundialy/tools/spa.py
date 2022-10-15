"""
SPA (Solar Position Algorithm) (https://midcdmz.nrel.gov/spa)

This algorithm calculates the solar zenith and azimuth angles in the period from the year -2000 to 6000, with
uncertainties of +/- 0.0003 degrees based on the date, time, and location on Earth. (Reference: Reda, I.; Andreas, A.,
Solar Position Algorithm for Solar Radiation Applications, Solar Energy. Vol. 76(5), 2004; pp. 577-589).

Reda, I.; Andreas, A. (2003). Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No.
TP-560-34302, Revised January 2008.:
http://www.nrel.gov/docs/fy08osti/34302.pdf
"""

from .constants import DT
import math
from datetime import datetime, timedelta
from typing import Union, Optional, List, Tuple

NUMBER_TYPE = Union[int, float]
INT_NONE_TYPE = Optional[int]
NUMBER_NONE_TYPE = Union[int, float, None]
DATA_TYPE = List[List[NUMBER_TYPE]]
DATA_INT_TYPE = List[List[int]]
DATA_INT_MISSING_TYPE = List[List[INT_NONE_TYPE]]
DATA_MISSING_TYPE = List[List[NUMBER_NONE_TYPE]]
SPA_RETURN_TYPE = Tuple[Tuple[float, float, float, float, float, float, float], Tuple[float, float, float, float, str],
                        Tuple[int, int, float]]


def _sideral_time(year: int, month: int, day: NUMBER_TYPE, hour: NUMBER_TYPE, minute: NUMBER_TYPE, second: NUMBER_TYPE,
                  microsecond: NUMBER_TYPE, dt: NUMBER_TYPE = DT) -> float:
    """
    Calculates the apparent sidereal time at Greenwich in degrees for use inside the spa function.

    :param year: Year
    :param month: Month
    :param day: Day
    :param hour: Hour
    :param minute: Minute
    :param second: Second
    :param microsecond: Microsecond
    :param dt: The difference between the Earth rotation time and the Terrestrial Time (TT)
    :return: The apparent sidereal time at Greenwich in degrees
    """

    # Main Report
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

    x = [297.85036 + 445267.111480 * jce - 0.0019142 * jce ** 2 + jce ** 3 / 189474.,
         357.52772 + 35999.050340 * jce - 0.0001603 * jce ** 2 - jce ** 3 / 300000.,
         134.96298 + 477198.867398 * jce + 0.0086972 * jce ** 2 + jce ** 3 / 56250.,
         93.27191 + 483202.017538 * jce - 0.0036825 * jce ** 2 + jce ** 3 / 327270.,
         125.04452 - 1934.136261 * jce + 0.0020708 * jce ** 2 + jce ** 3 / 450000.]

    x = list(map(math.radians, x))

    dy_i = []
    for i in range(63):
        x_i = []
        for j in range(5):
            x_i.append(x[j] * cfst[i][j])
        dy_i.append((cfdy[i][0] + cfdy[i][1] * jce) * math.sin(sum(x_i)))
    de_i = []
    for i in range(63):
        x_i = []
        for j in range(5):
            x_i.append(x[j] * cfst[i][j])
        de_i.append((cfde[i][0] + cfde[i][1] * jce) * math.cos(sum(x_i)))
    dy = sum(dy_i) / 36000000
    de = sum(de_i) / 36000000
    u = jme / 10
    e0_2 = 84381.448 - 4680.93 * u - 1.55 * u ** 2 + 1999.25 * u ** 3 - 51.38 * u ** 4 - 249.67 * u ** 5 \
        - 39.05 * u ** 6 + 7.12 * u ** 7 + 27.87 * u ** 8 + 5.79 * u ** 9 + 2.45 * u ** 10
    e_2 = e0_2 / 3600 + de
    n0 = 280.46061837 + 360.98564736629 * (jd - 2451545) + 0.000387933 * jc ** 2 - jc ** 3 / 38710000
    n0 = n0 % 360
    n = n0 + dy * math.cos(math.radians(e_2))
    return n


def _geocentric_ra_and_d(year: int, month: int, day: NUMBER_TYPE, hour: NUMBER_TYPE, minute: NUMBER_TYPE,
                         second: NUMBER_TYPE, microsecond: NUMBER_TYPE, dt: NUMBER_TYPE = DT) -> Tuple[float, float]:
    """
    Calculates the geocentric sun right ascension and sun declination in degrees for use inside the spa function.

    :param year: Year
    :param month: Month
    :param day: Day
    :param hour: Hour
    :param minute: Minute
    :param second: Second
    :param microsecond: Microsecond
    :param dt: The difference between the Earth rotation time and the Terrestrial Time (TT)
    :return: The geocentric sun right ascension and sun declination in degrees
    """

    # Main Report
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

    earth_periodic_terms_l0: DATA_TYPE = [[175347046, 0, 0], [3341656, 4.6692568, 6283.07585],
                                          [34894, 4.6261, 12566.1517], [3497, 2.7441, 5753.3849],
                                          [3418, 2.8289, 3.5231], [3136, 3.6277, 77713.7715],
                                          [2676, 4.4181, 7860.4194], [2343, 6.1352, 3930.2097],
                                          [1324, 0.7425, 11506.7698], [1273, 2.0371, 529.691],
                                          [1199, 1.1096, 1577.3435], [990, 5.233, 5884.927],
                                          [902, 2.045, 26.298], [857, 3.508, 398.149],
                                          [780, 1.179, 5223.694], [753, 2.533, 5507.553],
                                          [505, 4.583, 18849.228], [492, 4.205, 775.523],
                                          [357, 2.92, 0.067], [317, 5.849, 11790.629],
                                          [284, 1.899, 796.298], [271, 0.315, 10977.079],
                                          [243, 0.345, 5486.778], [206, 4.806, 2544.314],
                                          [205, 1.869, 5573.143], [202, 2.458, 6069.777],
                                          [156, 0.833, 213.299], [132, 3.411, 2942.463],
                                          [126, 1.083, 20.775], [115, 0.645, 0.98],
                                          [103, 0.636, 4694.003], [102, 0.976, 15720.839],
                                          [102, 4.267, 7.114], [99, 6.21, 2146.17],
                                          [98, 0.68, 155.42], [86, 5.98, 161000.69],
                                          [85, 1.3, 6275.96], [85, 3.67, 71430.7],
                                          [80, 1.81, 17260.15], [79, 3.04, 12036.46],
                                          [75, 1.76, 5088.63], [74, 3.5, 3154.69],
                                          [74, 4.68, 801.82], [70, 0.83, 9437.76],
                                          [62, 3.98, 8827.39], [61, 1.82, 7084.9],
                                          [57, 2.78, 6286.6], [56, 4.39, 14143.5],
                                          [56, 3.47, 6279.55], [52, 0.19, 12139.55],
                                          [52, 1.33, 1748.02], [51, 0.28, 5856.48],
                                          [49, 0.49, 1194.45], [41, 5.37, 8429.24],
                                          [41, 2.4, 19651.05], [39, 6.17, 10447.39],
                                          [37, 6.04, 10213.29], [37, 2.57, 1059.38],
                                          [36, 1.71, 2352.87], [36, 1.78, 6812.77],
                                          [33, 0.59, 17789.85], [30, 0.44, 83996.85],
                                          [30, 2.74, 1349.87], [25, 3.16, 4690.48]]

    earth_periodic_terms_l1: DATA_TYPE = [[628331966747, 0, 0], [206059, 2.678235, 6283.07585],
                                          [4303, 2.6351, 12566.1517], [425, 1.59, 3.523],
                                          [119, 5.796, 26.298], [109, 2.966, 1577.344],
                                          [93, 2.59, 18849.23], [72, 1.14, 529.69],
                                          [68, 1.87, 398.15], [67, 4.41, 5507.55],
                                          [59, 2.89, 5223.69], [56, 2.17, 155.42],
                                          [45, 0.4, 796.3], [36, 0.47, 775.52],
                                          [29, 2.65, 7.11], [21, 5.34, 0.98],
                                          [19, 1.85, 5486.78], [19, 4.97, 213.3],
                                          [17, 2.99, 6275.96], [16, 0.03, 2544.31],
                                          [16, 1.43, 2146.17], [15, 1.21, 10977.08],
                                          [12, 2.83, 1748.02], [12, 3.26, 5088.63],
                                          [12, 5.27, 1194.45], [12, 2.08, 4694],
                                          [11, 0.77, 553.57], [10, 1.3, 6286.6],
                                          [10, 4.24, 1349.87], [9, 2.7, 242.73],
                                          [9, 5.64, 951.72], [8, 5.3, 2352.87],
                                          [6, 2.65, 9437.76], [6, 4.67, 4690.48]]

    earth_periodic_terms_l2: DATA_TYPE = [[52919, 0, 0], [8720, 1.0721, 6283.0758], [309, 0.867, 12566.152],
                                          [27, 0.05, 3.52], [16, 5.19, 26.3], [16, 3.68, 155.42],
                                          [10, 0.76, 18849.23], [9, 2.06, 77713.77], [7, 0.83, 775.52],
                                          [5, 4.66, 1577.34], [4, 1.03, 7.11], [4, 3.44, 5573.14],
                                          [3, 5.14, 796.3], [3, 6.05, 5507.55], [3, 1.19, 242.73],
                                          [3, 6.12, 529.69], [3, 0.31, 398.15], [3, 2.28, 553.57],
                                          [2, 4.38, 5223.69], [2, 3.75, 0.98]]

    earth_periodic_terms_l3: DATA_TYPE = [[289, 5.844, 6283.076], [35, 0, 0], [17, 5.49, 12566.15], [3, 5.2, 155.42],
                                          [1, 4.72, 3.52], [1, 5.3, 18849.23], [1, 5.97, 242.73]]

    earth_periodic_terms_l4: DATA_TYPE = [[114, 3.142, 0], [8, 4.13, 6283.08], [1, 3.84, 12566.15]]

    earth_periodic_terms_l5: DATA_TYPE = [[1, 3.14, 0]]

    earth_periodic_terms_b0: DATA_TYPE = [[280, 3.199, 84334.662], [102, 5.422, 5507.553], [80, 3.88, 5223.69],
                                          [44, 3.7, 2352.87], [32, 4, 1577.34]]

    earth_periodic_terms_b1: DATA_TYPE = [[9, 3.9, 5507.55], [6, 1.73, 5223.69]]

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

    ept_l0 = earth_periodic_terms_l0
    ept_l1 = earth_periodic_terms_l1
    ept_l2 = earth_periodic_terms_l2
    ept_l3 = earth_periodic_terms_l3
    ept_l4 = earth_periodic_terms_l4
    ept_l5 = earth_periodic_terms_l5
    ept_b0 = earth_periodic_terms_b0
    ept_b1 = earth_periodic_terms_b1
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

    x = [297.85036 + 445267.111480 * jce - 0.0019142 * jce ** 2 + jce ** 3 / 189474.,
         357.52772 + 35999.050340 * jce - 0.0001603 * jce ** 2 - jce ** 3 / 300000.,
         134.96298 + 477198.867398 * jce + 0.0086972 * jce ** 2 + jce ** 3 / 56250.,
         93.27191 + 483202.017538 * jce - 0.0036825 * jce ** 2 + jce ** 3 / 327270.,
         125.04452 - 1934.136261 * jce + 0.0020708 * jce ** 2 + jce ** 3 / 450000.]

    x = list(map(math.radians, x))

    l0_i = []
    for i in range(len(ept_l0)):
        l0_i.append(ept_l0[i][0] * math.cos(ept_l0[i][1] + ept_l0[i][2] * jme))
    l0 = sum(l0_i)
    l1_i = []
    for i in range(len(ept_l1)):
        l1_i.append(ept_l1[i][0] * math.cos(ept_l1[i][1] + ept_l1[i][2] * jme))
    l1 = sum(l1_i)
    l2_i = []
    for i in range(len(ept_l2)):
        l2_i.append(ept_l2[i][0] * math.cos(ept_l2[i][1] + ept_l2[i][2] * jme))
    l2 = sum(l2_i)
    l3_i = []
    for i in range(len(ept_l3)):
        l3_i.append(ept_l3[i][0] * math.cos(ept_l3[i][1] + ept_l3[i][2] * jme))
    l3 = sum(l3_i)
    l4_i = []
    for i in range(len(ept_l4)):
        l4_i.append(ept_l4[i][0] * math.cos(ept_l4[i][1] + ept_l4[i][2] * jme))
    l4 = sum(l4_i)
    l5_i = []
    for i in range(len(ept_l5)):
        l5_i.append(ept_l5[i][0] * math.cos(ept_l5[i][1] + ept_l5[i][2] * jme))
    l5 = sum(l5_i)
    b0_i = []
    for i in range(len(ept_b0)):
        b0_i.append(ept_b0[i][0] * math.cos(ept_b0[i][1] + ept_b0[i][2] * jme))
    b0 = sum(b0_i)
    b1_i = []
    for i in range(len(ept_b1)):
        b1_i.append(ept_b1[i][0] * math.cos(ept_b1[i][1] + ept_b1[i][2] * jme))
    b1 = sum(b1_i)
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
    l = (l0 + l1 * jme + l2 * jme ** 2 + l3 * jme ** 3 + l4 * jme ** 4 + l5 * jme ** 5) / 100000000
    B = (b0 + b1 * jme) / 100000000
    B = math.degrees(B)
    r = (r0 + r1 * jme + r2 * jme ** 2 + r3 * jme ** 3 + r4 * jme ** 4) / 100000000
    l = math.degrees(l)
    l = l % 360
    b = -B
    theta2 = l + 180
    theta2 = theta2 % 360
    dy_i = []
    for i in range(63):
        x_i = []
        for j in range(5):
            x_i.append(x[j] * cfst[i][j])
        dy_i.append((cfdy[i][0] + cfdy[i][1] * jce) * math.sin(sum(x_i)))
    de_i = []
    for i in range(63):
        x_i = []
        for j in range(5):
            x_i.append(x[j] * cfst[i][j])
        de_i.append((cfde[i][0] + cfde[i][1] * jce) * math.cos(sum(x_i)))
    dy = sum(dy_i) / 36000000
    de = sum(de_i) / 36000000
    u = jme / 10
    e0_2 = 84381.448 - 4680.93 * u - 1.55 * u ** 2 + 1999.25 * u ** 3 - 51.38 * u ** 4 - 249.67 * u ** 5 \
        - 39.05 * u ** 6 + 7.12 * u ** 7 + 27.87 * u ** 8 + 5.79 * u ** 9 + 2.45 * u ** 10
    e_2 = e0_2 / 3600 + de
    dt_2 = -(20.4898 / (3600 * r))
    l = theta2 + dy + dt_2
    a2 = math.degrees(math.atan2(math.sin(math.radians(l)) * math.cos(math.radians(e_2)) - math.tan(math.radians(
        b)) * math.sin(math.radians(e_2)), math.cos(math.radians(l))))
    a2 = a2 % 360
    d = math.degrees(math.asin(math.sin(math.radians(b)) * math.cos(math.radians(e_2)) + math.cos(math.radians(
        b)) * math.sin(math.radians(e_2)) * math.sin(math.radians(l))))
    return d, a2


def spa(year: int, month: int, day: NUMBER_TYPE, hour: NUMBER_TYPE, minute: NUMBER_TYPE, second: NUMBER_TYPE,
        microsecond: NUMBER_TYPE, latitude: NUMBER_TYPE, longitude: NUMBER_TYPE, elevation: NUMBER_TYPE,
        pressure: NUMBER_TYPE, temperature: NUMBER_TYPE,
        omega: NUMBER_TYPE, gamma: NUMBER_TYPE, dt: NUMBER_TYPE = DT) -> SPA_RETURN_TYPE:
    """
    SPA (Solar Position Algorithm).

    This algorithm calculates the solar zenith and azimuth angles in the period from the year -2000 to 6000, with
    uncertainties of +/- 0.0003 degrees based on the date, time, and location on Earth.

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
    :param omega: The slope of the surface measured from the horizontal plane
    :param gamma: The surface azimuth rotation angle
    :param dt: The difference between the Earth rotation time and the Terrestrial Time (TT)
    :return: ((incidence angle, topocentric zenith angle, topocentric azimuth angle, topocentric sun declination,
        topocentric local hour angle, topocentric sun right ascension, e topocentric elevation angle), (Equation of Time,
        sun transit, sunrise, sunset, note), (year, month, day))
    """

    # Main Report
    original_time = (year, month, day, hour, minute, second, microsecond)
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

    earth_periodic_terms_l0: DATA_TYPE = [[175347046, 0, 0], [3341656, 4.6692568, 6283.07585],
                                          [34894, 4.6261, 12566.1517], [3497, 2.7441, 5753.3849],
                                          [3418, 2.8289, 3.5231], [3136, 3.6277, 77713.7715],
                                          [2676, 4.4181, 7860.4194], [2343, 6.1352, 3930.2097],
                                          [1324, 0.7425, 11506.7698], [1273, 2.0371, 529.691],
                                          [1199, 1.1096, 1577.3435], [990, 5.233, 5884.927],
                                          [902, 2.045, 26.298], [857, 3.508, 398.149],
                                          [780, 1.179, 5223.694], [753, 2.533, 5507.553],
                                          [505, 4.583, 18849.228], [492, 4.205, 775.523],
                                          [357, 2.92, 0.067], [317, 5.849, 11790.629],
                                          [284, 1.899, 796.298], [271, 0.315, 10977.079],
                                          [243, 0.345, 5486.778], [206, 4.806, 2544.314],
                                          [205, 1.869, 5573.143], [202, 2.458, 6069.777],
                                          [156, 0.833, 213.299], [132, 3.411, 2942.463],
                                          [126, 1.083, 20.775], [115, 0.645, 0.98],
                                          [103, 0.636, 4694.003], [102, 0.976, 15720.839],
                                          [102, 4.267, 7.114], [99, 6.21, 2146.17],
                                          [98, 0.68, 155.42], [86, 5.98, 161000.69],
                                          [85, 1.3, 6275.96], [85, 3.67, 71430.7],
                                          [80, 1.81, 17260.15], [79, 3.04, 12036.46],
                                          [75, 1.76, 5088.63], [74, 3.5, 3154.69],
                                          [74, 4.68, 801.82], [70, 0.83, 9437.76],
                                          [62, 3.98, 8827.39], [61, 1.82, 7084.9],
                                          [57, 2.78, 6286.6], [56, 4.39, 14143.5],
                                          [56, 3.47, 6279.55], [52, 0.19, 12139.55],
                                          [52, 1.33, 1748.02], [51, 0.28, 5856.48],
                                          [49, 0.49, 1194.45], [41, 5.37, 8429.24],
                                          [41, 2.4, 19651.05], [39, 6.17, 10447.39],
                                          [37, 6.04, 10213.29], [37, 2.57, 1059.38],
                                          [36, 1.71, 2352.87], [36, 1.78, 6812.77],
                                          [33, 0.59, 17789.85], [30, 0.44, 83996.85],
                                          [30, 2.74, 1349.87], [25, 3.16, 4690.48]]

    earth_periodic_terms_l1: DATA_TYPE = [[628331966747, 0, 0], [206059, 2.678235, 6283.07585],
                                          [4303, 2.6351, 12566.1517], [425, 1.59, 3.523],
                                          [119, 5.796, 26.298], [109, 2.966, 1577.344],
                                          [93, 2.59, 18849.23], [72, 1.14, 529.69],
                                          [68, 1.87, 398.15], [67, 4.41, 5507.55],
                                          [59, 2.89, 5223.69], [56, 2.17, 155.42],
                                          [45, 0.4, 796.3], [36, 0.47, 775.52],
                                          [29, 2.65, 7.11], [21, 5.34, 0.98],
                                          [19, 1.85, 5486.78], [19, 4.97, 213.3],
                                          [17, 2.99, 6275.96], [16, 0.03, 2544.31],
                                          [16, 1.43, 2146.17], [15, 1.21, 10977.08],
                                          [12, 2.83, 1748.02], [12, 3.26, 5088.63],
                                          [12, 5.27, 1194.45], [12, 2.08, 4694],
                                          [11, 0.77, 553.57], [10, 1.3, 6286.6],
                                          [10, 4.24, 1349.87], [9, 2.7, 242.73],
                                          [9, 5.64, 951.72], [8, 5.3, 2352.87],
                                          [6, 2.65, 9437.76], [6, 4.67, 4690.48]]

    earth_periodic_terms_l2: DATA_TYPE = [[52919, 0, 0], [8720, 1.0721, 6283.0758], [309, 0.867, 12566.152],
                                          [27, 0.05, 3.52], [16, 5.19, 26.3], [16, 3.68, 155.42],
                                          [10, 0.76, 18849.23], [9, 2.06, 77713.77], [7, 0.83, 775.52],
                                          [5, 4.66, 1577.34], [4, 1.03, 7.11], [4, 3.44, 5573.14],
                                          [3, 5.14, 796.3], [3, 6.05, 5507.55], [3, 1.19, 242.73],
                                          [3, 6.12, 529.69], [3, 0.31, 398.15], [3, 2.28, 553.57],
                                          [2, 4.38, 5223.69], [2, 3.75, 0.98]]

    earth_periodic_terms_l3: DATA_TYPE = [[289, 5.844, 6283.076], [35, 0, 0], [17, 5.49, 12566.15], [3, 5.2, 155.42],
                                          [1, 4.72, 3.52], [1, 5.3, 18849.23], [1, 5.97, 242.73]]

    earth_periodic_terms_l4: DATA_TYPE = [[114, 3.142, 0], [8, 4.13, 6283.08], [1, 3.84, 12566.15]]

    earth_periodic_terms_l5: DATA_TYPE = [[1, 3.14, 0]]

    earth_periodic_terms_b0: DATA_TYPE = [[280, 3.199, 84334.662], [102, 5.422, 5507.553], [80, 3.88, 5223.69],
                                          [44, 3.7, 2352.87], [32, 4, 1577.34]]

    earth_periodic_terms_b1: DATA_TYPE = [[9, 3.9, 5507.55], [6, 1.73, 5223.69]]

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

    ept_l0 = earth_periodic_terms_l0
    ept_l1 = earth_periodic_terms_l1
    ept_l2 = earth_periodic_terms_l2
    ept_l3 = earth_periodic_terms_l3
    ept_l4 = earth_periodic_terms_l4
    ept_l5 = earth_periodic_terms_l5
    ept_b0 = earth_periodic_terms_b0
    ept_b1 = earth_periodic_terms_b1
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
    jc = (jd - 2451545) / 36525
    jce = (jde - 2451545) / 36525
    jme = jce / 10

    x = [297.85036 + 445267.111480 * jce - 0.0019142 * jce ** 2 + jce ** 3 / 189474.,
         357.52772 + 35999.050340 * jce - 0.0001603 * jce ** 2 - jce ** 3 / 300000.,
         134.96298 + 477198.867398 * jce + 0.0086972 * jce ** 2 + jce ** 3 / 56250.,
         93.27191 + 483202.017538 * jce - 0.0036825 * jce ** 2 + jce ** 3 / 327270.,
         125.04452 - 1934.136261 * jce + 0.0020708 * jce ** 2 + jce ** 3 / 450000.]

    x = list(map(math.radians, x))

    l0_i = []
    for row_index in range(len(ept_l0)):
        l0_i.append(ept_l0[row_index][0] * math.cos(ept_l0[row_index][1] + ept_l0[row_index][2] * jme))
    l0 = sum(l0_i)
    l1_i = []
    for row_index in range(len(ept_l1)):
        l1_i.append(ept_l1[row_index][0] * math.cos(ept_l1[row_index][1] + ept_l1[row_index][2] * jme))
    l1 = sum(l1_i)
    l2_i = []
    for row_index in range(len(ept_l2)):
        l2_i.append(ept_l2[row_index][0] * math.cos(ept_l2[row_index][1] + ept_l2[row_index][2] * jme))
    l2 = sum(l2_i)
    l3_i = []
    for row_index in range(len(ept_l3)):
        l3_i.append(ept_l3[row_index][0] * math.cos(ept_l3[row_index][1] + ept_l3[row_index][2] * jme))
    l3 = sum(l3_i)
    l4_i = []
    for row_index in range(len(ept_l4)):
        l4_i.append(ept_l4[row_index][0] * math.cos(ept_l4[row_index][1] + ept_l4[row_index][2] * jme))
    l4 = sum(l4_i)
    l5_i = []
    for row_index in range(len(ept_l5)):
        l5_i.append(ept_l5[row_index][0] * math.cos(ept_l5[row_index][1] + ept_l5[row_index][2] * jme))
    l5 = sum(l5_i)
    b0_i = []
    for row_index in range(len(ept_b0)):
        b0_i.append(ept_b0[row_index][0] * math.cos(ept_b0[row_index][1] + ept_b0[row_index][2] * jme))
    b0 = sum(b0_i)
    b1_i = []
    for row_index in range(len(ept_b1)):
        b1_i.append(ept_b1[row_index][0] * math.cos(ept_b1[row_index][1] + ept_b1[row_index][2] * jme))
    b1 = sum(b1_i)
    r0_i = []
    for row_index in range(len(ept_r0)):
        r0_i.append(ept_r0[row_index][0] * math.cos(ept_r0[row_index][1] + ept_r0[row_index][2] * jme))
    r0 = sum(r0_i)
    r1_i = []
    for row_index in range(len(ept_r1)):
        r1_i.append(ept_r1[row_index][0] * math.cos(ept_r1[row_index][1] + ept_r1[row_index][2] * jme))
    r1 = sum(r1_i)
    r2_i = []
    for row_index in range(len(ept_r2)):
        r2_i.append(ept_r2[row_index][0] * math.cos(ept_r2[row_index][1] + ept_r2[row_index][2] * jme))
    r2 = sum(r2_i)
    r3_i = []
    for row_index in range(len(ept_r3)):
        r3_i.append(ept_r3[row_index][0] * math.cos(ept_r3[row_index][1] + ept_r3[row_index][2] * jme))
    r3 = sum(r3_i)
    r4_i = []
    for row_index in range(len(ept_r4)):
        r4_i.append(ept_r4[row_index][0] * math.cos(ept_r4[row_index][1] + ept_r4[row_index][2] * jme))
    r4 = sum(r4_i)
    l = (l0 + l1 * jme + l2 * jme ** 2 + l3 * jme ** 3 + l4 * jme ** 4 + l5 * jme ** 5) / 100000000
    B = (b0 + b1 * jme) / 100000000
    B = math.degrees(B)
    r = (r0 + r1 * jme + r2 * jme ** 2 + r3 * jme ** 3 + r4 * jme ** 4) / 100000000
    l = math.degrees(l)
    l = l % 360
    b = -B
    theta2 = l + 180
    theta2 = theta2 % 360
    dy_i = []
    for row_index in range(63):
        x_i = []
        for j in range(5):
            x_i.append(x[j] * cfst[row_index][j])
        dy_i.append((cfdy[row_index][0] + cfdy[row_index][1] * jce) * math.sin(sum(x_i)))
    de_i = []
    for row_index in range(63):
        x_i = []
        for j in range(5):
            x_i.append(x[j] * cfst[row_index][j])
        de_i.append((cfde[row_index][0] + cfde[row_index][1] * jce) * math.cos(sum(x_i)))
    dy = sum(dy_i) / 36000000
    de = sum(de_i) / 36000000
    u = jme / 10
    e0_2 = 84381.448 - 4680.93 * u - 1.55 * u ** 2 + 1999.25 * u ** 3 - 51.38 * u ** 4 - 249.67 * u ** 5 \
        - 39.05 * u ** 6 + 7.12 * u ** 7 + 27.87 * u ** 8 + 5.79 * u ** 9 + 2.45 * u ** 10
    e_2 = e0_2 / 3600 + de
    dt_2 = -(20.4898 / (3600 * r))
    l = theta2 + dy + dt_2
    n0 = 280.46061837 + 360.98564736629 * (jd - 2451545) + 0.000387933 * jc ** 2 - jc ** 3 / 38710000
    n0 = n0 % 360
    n = n0 + dy * math.cos(math.radians(e_2))
    a2 = math.degrees(math.atan2(
        math.sin(math.radians(l)) * math.cos(math.radians(e_2)) - math.tan(math.radians(b)) * math.sin(
            math.radians(e_2)), math.cos(math.radians(l))))
    a2 = a2 % 360
    d = math.degrees(math.asin(
        math.sin(math.radians(b)) * math.cos(math.radians(e_2)) + math.cos(math.radians(b)) * math.sin(
            math.radians(e_2)) * math.sin(math.radians(l))))
    h = n + longitude - a2
    h = h % 360
    xi = 8.794 / (3600 * r)
    u2 = math.degrees(math.atan(0.99664719 * math.tan(math.radians(latitude))))
    x2 = math.cos(math.radians(u2)) + elevation / 6378140 * math.cos(math.radians(latitude))
    y = 0.99664719 * math.sin(math.radians(u2)) + elevation / 6378140 * math.sin(math.radians(latitude))
    da = math.degrees(math.atan2(-x2 * math.sin(math.radians(xi)) * math.sin(math.radians(h)),
                                 math.cos(math.radians(d)) - x2 * math.sin(math.radians(xi)) * math.cos(
                                     math.radians(h))))
    a_prime = a2 + da
    d_prime = math.degrees(
        math.atan2((math.sin(math.radians(d)) - y * math.sin(math.radians(xi))) * math.cos(math.radians(da)),
                   math.cos(math.radians(d)) - x2 * math.sin(math.radians(xi)) * math.cos(math.radians(h))))
    h_prime = h - da
    e0 = math.degrees(math.asin(math.sin(math.radians(latitude)) * math.sin(math.radians(d_prime)) + math.cos(
        math.radians(latitude)) * math.cos(math.radians(d_prime)) * math.cos(math.radians(h_prime))))
    if e0 >= -1 * (0.26667 + 0.5667):
        de2 = (pressure / 1010) * (283 / (273 + temperature)) * (1.02 / (60 * math.tan(math.radians(e0 + 10.3 /
                                                                                                    (e0 + 5.11)))))
    else:
        de2 = 0
    e = e0 + de2
    theta = 90 - e
    gamma2 = math.degrees(math.atan2(math.sin(math.radians(h_prime)),
                                     math.cos(math.radians(h_prime)) * math.sin(math.radians(latitude)) - math.tan(
                                         math.radians(d_prime)) * math.cos(math.radians(latitude))))
    gamma2 = gamma2 % 360
    f = gamma2 + 180
    f = f % 360
    i = math.degrees(math.acos(
        math.cos(math.radians(theta)) * math.cos(math.radians(omega)) + math.sin(math.radians(omega)) * math.sin(
            math.radians(theta)) * math.cos(math.radians(gamma2 - gamma))))

    # Appendix
    m = 280.4664567 + 360007.6982779 * jme + 0.03032028 * jme ** 2 + jme ** 3 / 49931 - jme ** 4 / 15300 \
        - jme ** 5 / 2000000
    m = m % 360
    eot = m - 0.0057183 - a2 + dy * math.cos(math.radians(e_2))
    eot_minutes = eot * 4
    if eot_minutes > 20:
        eot_minutes -= 1440
    elif eot_minutes < -20:
        eot_minutes += 1440
    year, month, day, hour, minute, second, microsecond = original_time
    h_prime_0 = -1 * (0.26667 + 0.5667)
    n = _sideral_time(year, month, int(day), 0, 0, 0, 0, dt=dt)
    day_m1 = datetime(year, month, int(day)) - timedelta(days=1) - timedelta(seconds=dt)
    day_0 = datetime(year, month, int(day)) - timedelta(seconds=dt)
    day_p1 = datetime(year, month, int(day)) + timedelta(days=1) - timedelta(seconds=dt)
    d_m1, a_m1 = _geocentric_ra_and_d(day_m1.year, day_m1.month, day_m1.day, day_m1.hour, day_m1.minute, day_m1.second,
                                      day_m1.microsecond, dt=dt)
    d_0, a_0 = _geocentric_ra_and_d(day_0.year, day_0.month, day_0.day, day_0.hour, day_0.minute, day_0.second,
                                    day_0.microsecond, dt=dt)
    d_p1, a_p1 = _geocentric_ra_and_d(day_p1.year, day_p1.month, day_p1.day, day_p1.hour, day_p1.minute, day_p1.second,
                                      day_p1.microsecond, dt=dt)
    m_0 = (a_0 - longitude - n) / 360
    acos_arguement = (math.sin(math.radians(h_prime_0)) - math.sin(math.radians(latitude)) * math.sin(
        math.radians(d_0))) / (math.cos(math.radians(latitude)) * math.cos(math.radians(d_0)))
    if abs(acos_arguement) <= 1:
        H_0 = math.degrees(math.acos(acos_arguement))
        H_0 = H_0 % 180
        sun_time_notes = ''
    else:
        H_0 = math.degrees(math.acos(acos_arguement))
        H_0 = H_0 % 180
        # In the c code if the abs of the arccos is not smaller or equal to 1 it has H_0 = -99999,
        # but the paper didn't mention anything, so this isn't done here.
        # H_0 = -99999
        sun_time_notes = 'Sun is always above or below the horizon for that day'
    m_1 = m_0 - (H_0 / 360)
    m_2 = m_0 + (H_0 / 360)
    m_0 = m_0 % 1
    m_1 = m_1 % 1
    m_2 = m_2 % 1
    n_0 = n + 360.985647 * m_0
    n_1 = n + 360.985647 * m_1
    n_2 = n + 360.985647 * m_2
    n_english_0 = m_0 + dt / 86400
    n_english_1 = m_1 + dt / 86400
    n_english_2 = m_2 + dt / 86400
    a_parameter = a_0 - a_m1
    a_prime_parameter = d_0 - d_m1
    b_parameter = a_p1 - a_0
    b_prime_parameter = d_p1 - d_0
    # In the paper it says if it is greater than 2 but in the c code it is if it greater or equal to 2.
    if abs(a_parameter) > 2:
        a_parameter = a_parameter % 1
    if abs(a_prime_parameter) > 2:
        a_prime_parameter = a_prime_parameter % 1
    if abs(b_parameter) > 2:
        b_parameter = b_parameter % 1
    if abs(b_prime_parameter) > 2:
        b_prime_parameter = b_prime_parameter % 1
    c_parameter = b_parameter - a_parameter
    c_prime_parameter = b_prime_parameter - a_prime_parameter
    a_prime_0 = a_0 + (n_english_0 * (a_parameter + b_parameter + c_parameter * n_english_0)) / 2
    a_prime_1 = a_0 + (n_english_1 * (a_parameter + b_parameter + c_parameter * n_english_1)) / 2
    a_prime_2 = a_0 + (n_english_2 * (a_parameter + b_parameter + c_parameter * n_english_2)) / 2
    d_prime_0 = d_0 + (n_english_0 * (a_prime_parameter + b_prime_parameter + c_prime_parameter * n_english_0)) / 2
    d_prime_1 = d_0 + (n_english_1 * (a_prime_parameter + b_prime_parameter + c_prime_parameter * n_english_1)) / 2
    d_prime_2 = d_0 + (n_english_2 * (a_prime_parameter + b_prime_parameter + c_prime_parameter * n_english_2)) / 2
    H_prime_0 = n_0 + longitude - a_prime_0
    H_prime_1 = n_1 + longitude - a_prime_1
    H_prime_2 = n_2 + longitude - a_prime_2
    H_prime_0 = H_prime_0 / 360
    H_prime_0 = 360 * (H_prime_0 - math.floor(H_prime_0))
    H_prime_1 = H_prime_1 / 360
    H_prime_1 = 360 * (H_prime_1 - math.floor(H_prime_1))
    H_prime_2 = H_prime_2 / 360
    H_prime_2 = 360 * (H_prime_2 - math.floor(H_prime_2))
    h_0 = math.degrees(math.asin(math.sin(math.radians(latitude)) * math.sin(math.radians(d_prime_0)) + math.cos(
        math.radians(latitude)) * math.cos(math.radians(d_prime_0)) * math.cos(math.radians(H_prime_0))))
    h_1 = math.degrees(math.asin(math.sin(math.radians(latitude)) * math.sin(math.radians(d_prime_1)) + math.cos(
        math.radians(latitude)) * math.cos(math.radians(d_prime_1)) * math.cos(math.radians(H_prime_1))))
    h_2 = math.degrees(math.asin(math.sin(math.radians(latitude)) * math.sin(math.radians(d_prime_2)) + math.cos(
        math.radians(latitude)) * math.cos(math.radians(d_prime_2)) * math.cos(math.radians(H_prime_2))))
    T = m_0 - (H_prime_0 / 360)
    R = m_1 + (h_1 - h_prime_0) / (360 * math.cos(math.radians(d_prime_1)) * math.cos(math.radians(
        latitude)) * math.sin(math.radians(H_prime_1)))
    S = m_2 + (h_2 - h_prime_0) / (360 * math.cos(math.radians(d_prime_2)) * math.cos(math.radians(
        latitude)) * math.sin(math.radians(H_prime_2)))
    T = T % 1
    R = R % 1
    S = S % 1
    T *= 24
    R *= 24
    S *= 24

    # Appendix section A.3
    jd_plus_5 = jd + 0.5
    z_capital = int(jd_plus_5)
    f_capital = jd_plus_5 - z_capital
    if z_capital < 2299161:
        a_capital = z_capital
    else:
        b_capital = int((z_capital - 1867216.25) / 36524.25)
        a_capital = z_capital + 1 + b_capital - int(b_capital / 4)
    c_capital = a_capital + 1524
    d_capital = int((c_capital - 122.1) / 365.25)
    g_capital = int(365.25 * d_capital)
    i_capital = int((c_capital - g_capital) / 30.6001)
    day_from_jd = c_capital - g_capital - int(30.6001 * i_capital) + f_capital
    if i_capital < 14:
        month_from_jd = i_capital - 1
    else:
        month_from_jd = i_capital - 13
    if month_from_jd > 2:
        year_from_jd = d_capital - 4716
    else:
        year_from_jd = d_capital - 4715
    return ((i, theta, f, d_prime, h_prime, a_prime, e), (eot_minutes, T, R, S, sun_time_notes),
            (year_from_jd, month_from_jd, day_from_jd))
