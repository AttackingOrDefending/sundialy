"""
Bird (Bird Clear Sky Model) (https://www.nrel.gov/grid/solar-resource/clear-sky.html)

The Bird Clear Sky Model, authored by Richard Bird, is a broadband algorithm that produces estimates of clear sky direct
beam, hemispherical diffuse, and total hemispherical solar radiation on a horizontal surface.

The model is based on comparisons with results from rigorous radiative transfer codes. It is composed of simple
algebraic expressions with 10 user-provided inputs. Model results should be expected to agree within ±10% with rigorous
radiative transfer codes. The model computes hourly average solar radiation for every hour of the year, based on the 10
user input parameters; however, variable atmospheric parameters such as aerosol optical depth, ozone, and water vapor
are fixed for the entire year.

Simplified Clear Sky Model for Direct and Diffuse Insolation on Horizontal Surfaces:
https://www.nrel.gov/grid/solar-resource/assets/data/tr-642-761.pdf
"""

import math
from typing import Union, Tuple

NUMBER_TYPE = Union[int, float]
BIRD_RETURN_TYPE = Tuple[NUMBER_TYPE, NUMBER_TYPE, NUMBER_TYPE, NUMBER_TYPE, NUMBER_TYPE, NUMBER_TYPE, NUMBER_TYPE]


def bird(r: NUMBER_TYPE, zenith: NUMBER_TYPE, pressure: NUMBER_TYPE, ozone: NUMBER_TYPE, water: NUMBER_TYPE,
         aerosol: NUMBER_TYPE, albedo: NUMBER_TYPE, dni_mod: NUMBER_TYPE, ba: NUMBER_TYPE = 0.85,
         k1: NUMBER_TYPE = 0.1) -> BIRD_RETURN_TYPE:
    """
    Bird (Bird Clear Sky Model).

    The Bird Clear Sky Model, authored by Richard Bird, is a broadband algorithm that produces estimates of clear sky direct
    beam, hemispherical diffuse, and total hemispherical solar radiation on a horizontal surface.

    :param r: The Sun’s distance from the center of the Earth in astronomical units
    :param zenith: The zenith
    :param pressure: The pressure
    :param ozone: Amount of ozone in a vertical column from surface (cm)
    :param water: Amount of precipitable water in a vertical column from surface (cm)
    :param aerosol: Broadband aerosol optical depth from surface in a vertical path (broadband turbidity)
    :param albedo: Ground albedo
    :param dni_mod: Area of the SUL With Respect to the Area of the Sun’s Disk
    :param ba: Ratio of the forward-scattered irradiance to the total scattered irradiance due to aerosols
    :param k1: Constant used in Bird model associated with aerosol absorptance
    :return: (Relative optical airmass (not pressure corrected), Estimated direct normal irradiance [W/m^2],
        Estimated global horiz irradiance [W/m^2], Estimated diffuse horiz irradiance [W/m^2],
        Estimated direct normal irradiance from sun's unshaded lune [W/m^2],
        Estimated global horiz irradiance from sun's unshaded lune [W/m^2],
        Estimated duffuse horiz irradiance from sun's unshaded lune [W/m^2])
    """
    # In the paper for it is ba=0.84 but in the c code it is ba=0.85
    if 0 <= zenith < 90 and r > 0:
        # In the paper it is 0.15 but in the c code it is 0.50572
        # In the paper it is 93.885 but in the c code it is 96.07995
        # In the paper it is -1.25 but in the c code it is -1.6364
        m = (math.cos(math.radians(zenith)) + 0.50572 * (96.07995 - zenith) ** -1.6364) ** -1
        m_prime = m * pressure / 1013
        xo = ozone * m
        xw = water * m
        tr = math.exp(-0.0903 * m_prime ** 0.84 * (1 + m_prime - m_prime ** 1.01))
        # In the paper it is -0.3035 but in the c code it is -0.3034
        to = 1 - 0.1611 * xo * (1 + 139.48 * xo) ** -0.3034 - 0.002715 * xo * (1 + 0.044 * xo + 0.0003 * xo ** 2) ** -1
        tum = math.exp(-0.0127 * m_prime ** 0.26)
        tw = 1 - 2.4959 * xw * ((1 + 79.034 * xw) ** 0.6828 + 6.385 * xw) ** -1
        ta = math.exp(-(aerosol ** 0.873) * (1 + aerosol - aerosol ** 0.7088) * m ** 0.9108)
        taa = 1 - k1 * (1 - m + m ** 1.06) * (1 - ta)
        tas = ta / taa
        rs = 0.0685 + (1 - ba) * (1 - tas)
        # In the paper it is 1353 but in the bird.c it is 1367. This number is most likely the solar constant.
        io = 1367 / r ** 2
        id2 = io * math.cos(math.radians(zenith)) * 0.9662 * tr * to * tum * tw * ta
        ias = io * math.cos(math.radians(zenith)) * 0.79 * to * tw * tum * taa * (0.5 * (1 - tr) + ba * (1 - tas)) / (
                1 - m + m ** 1.02)
        it = (id2 + ias) / (1 - albedo * rs)
        direct_normal = 0.9662 * io * tr * to * tum * tw * ta
        if direct_normal < 0:
            direct_normal = 0
        direct_horizontal = direct_normal * math.cos(math.radians(zenith))
        diffuse_horizontal = it - direct_horizontal
        direct_normal_mod, it_mod, diffuse_horizontal_mod = 0., 0, 0.
        if dni_mod >= 0:
            direct_normal_mod = direct_normal * dni_mod
            direct_horizontal_mod = direct_normal_mod * math.cos(math.radians(zenith))
            it_mod = (direct_horizontal_mod + ias) / (1 - albedo * rs)
            diffuse_horizontal_mod = it_mod - direct_horizontal_mod
        return m, direct_normal, it, diffuse_horizontal, direct_normal_mod, it_mod, diffuse_horizontal_mod
    else:
        return 0, 0, 0, 0, 0, 0, 0
