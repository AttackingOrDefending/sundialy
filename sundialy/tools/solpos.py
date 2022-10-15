"""
SOLPOS (Solar Position and Intensity) (https://www.nrel.gov/grid/solar-resource/solpos.html)

References
Astronomical Solar Position
Michalsky, J. 1988. The Astronomical Almanac's algorithm for approximate solar position (1950-2050). Solar Energy 40
(3), 227-235.

Michalsky, J. 1988. ERRATA: The astronomical almanac's algorithm for approximate solar position (1950-2050). Solar
Energy 41 (1), 113.

Distance from Sun to Earth
Spencer, J. W. 1971. Fourier series representation of the position of the sun. Search 2 (5), 172.
NOTE: This paper gives solar position algorithms as well, but the Michalsky/Almanac algorithm above is more accurate.

Atmospheric Refraction Correction
Zimmerman, John C. 1981. Sun-pointing programs and their accuracy. SAND81-0761, Experimental Systems Operation Division\
4721, Sandia National Laboratories, Albuquerque, NM.

Shadow Band Correction Factor
Drummond, A. J. 1956. A contribution to absolute pyrheliometry. Q. J. R. Meteorol.2 Soc. 82, 481-493.

Relative Optical Air Mass
Kasten, F. and Young, A. 1989. Revised optical air mass tables and approximation formula. Applied Optics 28 (22),
4735-4738.

Renormalization of KT (“PRIME”)
Perez, R., P. Ineichen, Seals, R., & Zelenka, A. 1990. Making full use of the clearness index for parameterizing hourly
insolation conditions. Solar Energy 45 (2), 111-114.

Solar Position Relative to Earth
Iqbal, M. 1983. An Introduction to Solar Radiation. Academic Press, NY.
"""

import math
from typing import Union, Tuple

NUMBER_TYPE = Union[int, float]
SOLPOS_RETURN_TYPE = Tuple[float, float, float, float, float, float, float, float, float, float, float, float, float,
                           float, float, float, float, float]


def solpos(year: int, month: int, day: int, hour: NUMBER_TYPE, minute: NUMBER_TYPE, second: NUMBER_TYPE,
           timezone: NUMBER_TYPE, latitude: NUMBER_TYPE, longitude: NUMBER_TYPE, pressure: NUMBER_TYPE,
           temperature: NUMBER_TYPE, aspect: NUMBER_TYPE = 180, tilt: NUMBER_TYPE = 0, sb_width: NUMBER_TYPE = 7.6,
           sb_radius: NUMBER_TYPE = 31.7, sb_sky: NUMBER_TYPE = 0.04, interval: NUMBER_TYPE = 0) -> SOLPOS_RETURN_TYPE:
    """
    SOLPOS (Solar Position and Intensity).

    :param year: Year
    :param month: Month
    :param day: Day
    :param hour: Hour
    :param minute: Minute
    :param second: Second
    :param timezone: Timezone
    :param latitude: Latitude
    :param longitude: Longitude
    :param pressure: Pressure (in millibars)
    :param temperature: Temperature (in Celsius)
    :param aspect: Where the surface faces. Default: 180 (south)
    :param tilt: The tilt of the surface. Default: 0 (horizontal)
    :param sb_width: Eppley shadow band width
    :param sb_radius: Eppley shadow band radius
    :param sb_sky: Drummond factor for partly cloudy skies
    :param interval: Instantaneous measurement interval
    :return: (airmass, pressure corrected airmass, zenith (refracted), azimuth angle, elevation (ETR),
        refracted solar elevation angle, ETR global, ETR direct, ETR tilt, sunrise time, sunset time,
        shadowband correction factor, prime, unprime, day angle, declination, equation of time, right ascension)
    """
    month_days = [[0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],
                  [0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]]
    daynum = day + month_days[0][month]
    if ((year % 4) == 0) and (((year % 100) != 0) or ((year % 400) == 0)) and (month > 2):
        daynum += 1
    if ((year % 4) == 0) and (((year % 100) != 0) or ((year % 400) == 0)):
        leap = 1
    else:
        leap = 0
    month = 12
    while daynum <= month_days[leap][month]:
        month -= 1
    # day = daynum - month_days[leap][month]
    day_angle = 360. * (daynum - 1) / 365.
    sd = math.sin(math.radians(day_angle))
    cd = math.cos(math.radians(day_angle))
    d2 = 2 * day_angle
    c2 = math.cos(math.radians(d2))
    s2 = math.sin(math.radians(d2))
    erv = 1.000110 + 0.034221 * cd + 0.001280 * sd
    erv += 0.000719 * c2 + 0.000077 * s2
    utime = hour * 3600 + minute * 60 + second - interval / 2
    utime = utime / 3600. - timezone
    delta = year - 1949
    leap = int(delta / 4)
    julday = 32916.5 + (delta * 365.) + leap + daynum + utime / 24.
    ectime = julday - 51545
    mnlong = 280.460 + 0.9856474 * ectime
    mnlong = mnlong % 360
    mnanom = 357.528 + 0.9856003 * ectime
    mnanom = mnanom % 360
    eclong = mnlong + 1.915 * math.sin(math.radians(mnanom)) + 0.020 * math.sin(math.radians(2.0 * mnanom))
    eclong = eclong % 360
    ecobli = 23.439 - 4.0e-07 * ectime
    declin = math.degrees(math.asin(math.sin(math.radians(ecobli)) * math.sin(math.radians(eclong))))
    top = math.cos(math.radians(ecobli)) * math.sin(math.radians(eclong))
    bottom = math.cos(math.radians(eclong))
    rascen = math.degrees(math.atan2(top, bottom))
    rascen = rascen % 360
    gmst = 6.697375 + 0.0657098242 * ectime + utime
    gmst = gmst % 24
    lmst = gmst * 15. + longitude
    lmst = lmst % 360
    hrang = lmst - rascen
    if hrang < -180:
        hrang += 360
    elif hrang > 180:
        hrang -= 360
    tdat_cd = math.cos(math.radians(declin))
    tdat_ch = math.cos(math.radians(hrang))
    tdat_cl = math.cos(math.radians(latitude))
    tdat_sd = math.sin(math.radians(declin))
    tdat_sl = math.sin(math.radians(latitude))
    cz = tdat_sd * tdat_sl + tdat_cd * tdat_cl * tdat_ch
    if abs(cz) > 1:
        if cz >= 0:
            cz = 1
        else:
            cz = -1
    zenetr = math.degrees(math.acos(cz))
    if zenetr > 99:
        zenetr = 99
    elevetr = 90. - zenetr
    cdcl = tdat_cd * tdat_cl
    if abs(cdcl) >= 0.001:
        cssha = -tdat_sl * tdat_sd / cdcl
        if cssha < -1:
            ssha = 180.
        elif cssha > 1:
            ssha = 0.
        else:
            ssha = math.degrees(math.acos(cssha))
    elif ((declin >= 0.) and (latitude > 0.)) or ((declin < 0.) and (latitude < 0.)):
        ssha = 180.
    else:
        ssha = 0.
    p = 0.6366198 * sb_width / sb_radius * pow(tdat_cd, 3)
    t1 = math.radians(tdat_sl * tdat_sd * ssha)
    t2 = tdat_cl * tdat_cd * math.sin(math.radians(ssha))
    sbcf = sb_sky + 1 / (1 - p * (t1 + t2))
    tst = (180 + hrang) * 4
    tstfix = tst - hour * 60 - minute - second / 60 + interval / 120
    while tstfix > 720:
        tstfix -= 1440
    while tstfix < -720:
        tstfix += 1440
    eqntim = tstfix + 60 * timezone - 4 * longitude
    if ssha <= 1.:
        sretr = 2999.
        ssetr = -2999.
    elif ssha >= 179:
        sretr = -2999.
        ssetr = 2999.
    else:
        sretr = 720. - 4 * ssha - tstfix
        ssetr = 720. + 4 * ssha - tstfix
    ce = math.cos(math.radians(elevetr))
    se = math.sin(math.radians(elevetr))
    azim = 180.
    cecl = ce * tdat_cl
    if abs(cecl) >= 0.001:
        ca = (se * tdat_sl - tdat_sd) / cecl
        if ca > 1:
            ca = 1.
        elif ca < -1:
            ca = -1.0
        azim = 180.0 - math.degrees(math.acos(ca))
        if hrang > 0:
            azim = 360 - azim
    if elevetr > 85:
        refcor = 0.
    else:
        tanelev = math.tan(math.radians(elevetr))
        if elevetr >= 5:
            refcor = 58.1 / tanelev - 0.07 / (pow(tanelev, 3)) + 0.000086 / (pow(tanelev, 5))
        elif elevetr >= -0.575:
            refcor = 1735 + elevetr * (-518.2 + elevetr * (103.4 + elevetr * (-12.79 + elevetr * 0.711)))
        else:
            refcor = -20.774 / tanelev
        prestemp = (pressure * 283) / (1013 * (273 + temperature))
        refcor *= prestemp / 3600.0
    elevref = elevetr + refcor
    if elevref < -9:
        elevref = -9
    zenref = 90 - elevref
    coszen = math.cos(math.radians(zenref))
    if zenref > 93:
        amass = -1
        ampress = -1.
    else:
        amass = 1 / (math.cos(math.radians(zenref)) + 0.50572 * pow((96.07995 - zenref), -1.6364))
        ampress = amass * pressure / 1013
    unprime = 1.031 * math.exp(-1.4 / (0.9 + 9.4 / amass)) + 0.1
    prime = 1 / unprime
    solcon = 1367
    if coszen > 0:
        etrn = solcon * erv
        etr = etrn * coszen
    else:
        etrn = 0
        etr = 0
    ca = math.cos(math.radians(azim))
    cp = math.cos(math.radians(aspect))
    ct = math.cos(math.radians(tilt))
    sa = math.sin(math.radians(azim))
    sp = math.sin(math.radians(aspect))
    st = math.sin(math.radians(tilt))
    sz = math.sin(math.radians(zenref))
    cosinc = coszen * ct + sz * st * (ca * cp + sa * sp)
    if cosinc > 0:
        etrtilt = etrn * cosinc
    else:
        etrtilt = 0
    return (amass, ampress, zenref, azim, elevetr, elevref, etr, etrn, etrtilt, sretr, ssetr, sbcf, prime, unprime,
            day_angle, declin, eqntim, rascen)
