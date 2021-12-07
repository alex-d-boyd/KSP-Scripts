import math
from math import pi
from typing import Iterator, Tuple, Union, Optional
from dataclasses import dataclass

from kerbal_utils import Body, bodies


KERBIN_RADIUS: float = 600000

@dataclass
class Scanner:
    min_alt: int
    best_alt: int
    max_alt: int
    fov: float

SCANNERS = {'R-3B Radar Altimeter': Scanner(5,70,250,1.5),
            'R-EO-1 Radar Antenna': Scanner(50,100,500,3.5),
            'SAR-X Antenna': Scanner(70,250,500,1.5),
            'SAR-C Antenna': Scanner(500,700,750,3.0),
            'SAR-L Antenna': Scanner(250,500,1000,4.0),
            'SCAN-R Resource Mapper': Scanner(20,70,250,1.0),
            'M4435 Narrow-Band Scanner': Scanner(10,150,500,2.0),
            'SCAN-R2 Advanced Resource Mapper': Scanner(70,250,500,2.5),
            'M700 Survey Scanner': Scanner(15,500,7500,3.0),
            'SCAN-RX Hyperspectral Resource Mapper': Scanner(100,500,750,3.0),
            'VS-1 High Resolution Imager': Scanner(20,70,250,1.5),
            'VS-11 Classified Reconnaissance Imager': Scanner(100,200,1000,4.0),
            'VS-3 Advanced High Resolution Imager': Scanner(70,250,500,2.5),
            'MS-1 Multispectral Scanner': Scanner(20,70,250,3.0),
            'MS-R Enhanced Multispectral Scanner': Scanner(70,300,400,1.5),
            'MS-2A Advanced Multispectral Scanner': Scanner(100,500,750,4.0),
            }


##@dataclass
##class Body:
##    mu: float
##    radius: float
##    sidereal_day: float
##    safe_altitude: float


def coprimes_of(n: int, limit: int = -1) -> Iterator[int]:
    """
    Create a generator that returns the coprimes of n upto a limit, or forever if none provided

    :param n: The number to find coprimes of
    :param limit: The cutoff point (inclusive). -1 if no cutoff
    :return: Coprimes of n
    """
    k: int = 1
    while k <= limit or limit == -1:
        if math.gcd(n, k) == 1:
            yield k
        k += 1


def cuberoot(n: float) -> float:
    if n < 0:
        return -cuberoot(-n)
    return n**(1/3)


def cubic_roots(a: float, b: float, c: float, d: float) -> Union[float, Tuple[float, float, float]]:
    p= (3*a*c - b**2)/(3 * a**2)
    q = (2*b**3 - 9*a*b*c + 27*d*a**2)/(27 * a**3)

    discriminant = (p**3)/27 + (q**2)/4

    xt = b/(3*a)
    if discriminant > 0:
        disc_root = math.sqrt(discriminant)

        t_plus = -(q/2) + disc_root
        t_minus = -(q/2) - disc_root
        t = cuberoot(t_plus) + cuberoot(t_minus)

        x = t - xt
        return x

    if discriminant == 0:
        if p == 0:
            return -xt, -xt, -xt
        t1 = 3*q/p
        t2 = -3*q/(2*p)

        x1 = t1 - xt
        x2 = t2 - xt

        return x1, x2, x2

    imag_part = math.sqrt(-discriminant)
    real_part = -q/2

    magnitude = math.sqrt(real_part**2 + imag_part**2)
    angle = math.atan2(imag_part, real_part)

    magnitude = cuberoot(magnitude)
    angle = angle/3

    t1 = 2 * magnitude*math.cos(angle)
    t2 = 2 * magnitude*math.cos(angle + 2*pi/3)
    t3 = 2 * magnitude*math.cos(angle - 2*pi/3)

    return t1-xt, t2-xt, t3-xt

def ap_peri(sma, ecc, body):
    linnear_ecc = ecc * sma
    apoapsis = linnear_ecc + sma
    periapsis = 2 * sma - apoapsis
    slr = sma * (1 - ecc**2)
    return (apoapsis-body.radius, periapsis-body.radius, slr-body.radius)

def find_eccentricity_bounds(body: Body,
                             semi_major_axis: float,
                             track_altitude: float,
                             min_altitude: float,
                             max_altitude: float) -> Optional[Tuple[float, float]]:

    ecc_polar = math.sqrt(max(1 - (body.radius + min_altitude)/semi_major_axis, 0))
    ecc_periapsis = 1 - (body.radius + body.safe_altitude)/semi_major_axis
    ecc_apoapsis = (body.radius + max_altitude)/semi_major_axis - 1

    ecc_max = min(ecc_polar, ecc_periapsis, ecc_apoapsis)

    ecc_min = (body.radius + track_altitude - semi_major_axis)/semi_major_axis

    if ecc_min >= ecc_max:
        return None

    ecc_fixed = math.sqrt((semi_major_axis - body.radius)/semi_major_axis)
    
    min_past_limit = ecc_min > track_altitude/(2*track_altitude + body.radius)
    max_past_limit = track_altitude > body.radius * ecc_fixed
    if min_past_limit or max_past_limit:
        a = 4 * track_altitude * semi_major_axis
        b = body.radius ** 2
        c = 2 * track_altitude * (body.radius - 2 * semi_major_axis)
        d = track_altitude ** 2
        roots = cubic_roots(a, b, c, d)
        if type(roots) is not tuple:
            return None

        roots = sorted(roots)
        if min_past_limit:
            ecc_min = roots[1]
        if max_past_limit:  # the condition can only look at the boundry of when ecc_fixed becomes incorrect, others may still be valid
            ecc_max = min(ecc_max, roots[2])

    return max(0, ecc_min), ecc_max

def find_eccentricity_bounds_reduced(body: Body,
                             semi_major_axis: float,
                             track_altitude: float,
                             min_altitude: float,
                             max_altitude: float) -> Optional[Tuple[float, float]]:

    ecc_polar = math.sqrt(max(1 - (body.radius + min_altitude)/semi_major_axis, 0))
    ecc_periapsis = 1 - (body.radius + body.safe_altitude)/semi_major_axis
    ecc_apoapsis = (body.radius + max_altitude)/semi_major_axis - 1

    ecc_max = min(ecc_polar, ecc_periapsis, ecc_apoapsis)

    ecc_min = (body.radius + track_altitude - semi_major_axis)/semi_major_axis

    if ecc_min >= ecc_max:
        return None

    return max(0, ecc_min), ecc_max

def find_orbits(body: Body, field_of_view: float, min_altitude: float, best_altitude: float, max_altitude: float,
                loose: bool, ecc_bounder=find_eccentricity_bounds):
    """
    Calculates the orbit with the smallest semi-major axis capable of scanning the body in the shortest time possible


    :param body: the body being orbited
    :param field_of_view: the base field of view for the satellite
    :param min_altitude: the lowest working altitude for the satellite
    :param best_altitude: the altitude at which the working fov is the base fov at kerbin
    :param max_altitude: the maximum working altitude of the satellite
    :return: The orbital properties for the found orbit as a dict {sma, min_ecc, max_ecc, tracks, skip}.<br>
        <bl>
            <li>sma is the semi-major axis</li>
            <li>min_ecc is the lower bound on eccentricity</li>
            <li>max_ecc us the upper bound on eccentricity</li>
            <li>tracks is the number of descending tracks made while scanning</li>
            <li>skip is the number of tracks between consecutively scanned tracks, also equates to the number of
            sidereal days (for the body) required to complete the scan</li>
        </bl>
    """
    # scale the field of view to the planet's radius as done in ScanSat
    scaled_fov: float = field_of_view
    if body.radius < KERBIN_RADIUS:
        scaled_fov *= math.sqrt(KERBIN_RADIUS / body.radius)

    # calculate the minimum number of ground tracks for 100% coverage at the equator
    min_ground_tracks: int = math.ceil(180 / min(scaled_fov, 20))

    # calculate orbital limits
    max_radius_apoapsis: float = body.radius + max_altitude  # cannot scan at apoapsis if above best altitude
    min_radius_periapsis: float = body.radius + body.safe_altitude  # if periapsis below safe altitude we will crash
    min_radius_polar: float = body.radius + min_altitude  # cannot scan at poles if radius below minimum altitude

    best_track_skip: int = -1
    best_orbit = None
    for n_tracks in range(min_ground_tracks, 2*min_ground_tracks):
        track_angle: float = 180 / n_tracks
        track_altitude = track_angle * best_altitude / scaled_fov

        for track_skip in coprimes_of(n_tracks, best_track_skip):

            orbit_period: float = track_skip * body.sidereal_day / n_tracks
            semi_major_axis: float = (body.mu * orbit_period**2 / (4 * pi**2)) ** (1/3)

            if semi_major_axis > max_radius_apoapsis:
                break  # orbit is too big. increase n_tracks and start again
            if semi_major_axis < min_radius_periapsis or semi_major_axis < min_radius_polar:
                continue  # axis is too small, advance to increase

            eccentricity = ecc_bounder(body, semi_major_axis, track_altitude, min_altitude, max_altitude)

            if eccentricity is None:
                continue  # no valid orbits, continue

            # no need to check if orbit is smaller as new period will always be lower as n_tracks increases
            # and once set best_track skip can only stay ths same or decrease
            best_orbit = {
                "sma": semi_major_axis,
                "min_ecc": eccentricity[0],
                "max_ecc": eccentricity[1],
                "tracks": n_tracks,
                "skip": track_skip
            }
            best_track_skip = track_skip - int(loose)
            break

    return best_orbit


def main():
##    bodies = {
##        "kerbin": Body(3.5316e12, KERBIN_RADIUS, 21549.425, 70000),
##        "mun": Body(6.5138398e10, 200000, 138984.38, 10000),
##        "minmus": Body(1.7658e9, 60000, 40400, 6000),
##        "gilly": Body(8289449.8, 13000, 28255, 8000),
##        "jool": Body(2.82528e14, 6000000, 36000, 200000)
##    }
    body = input("input target body:")
    while body not in bodies:
        print(f"body ${body} not found, expected one of", bodies.keys)
        body = input("input target body:")

    fov = float(input("fov: "))
    min_alt = float(input("min alt: "))
    best_alt = float(input("best alt: "))
    max_alt = float(input("max alt: "))
    loose = input("loose? (y/n)").lower() == "y"

    best_orbit = find_orbits(bodies[body], fov, min_alt, best_alt, max_alt, loose)
    print(best_orbit)

def test():
    for body_name, body in bodies.items():
        for scanner_name, scanner in SCANNERS.items():
            for loose_val in (True, False):
                orig = find_orbits(body,
                                   scanner.fov,
                                   scanner.min_alt*1_000,
                                   scanner.best_alt*1_000,
                                   scanner.max_alt*1_000,
                                   loose_val,
                                   )
                new =  find_orbits(body,
                                   scanner.fov,
                                   scanner.min_alt*1_000,
                                   scanner.best_alt*1_000,
                                   scanner.max_alt*1_000,
                                   loose_val,
                                   find_eccentricity_bounds_reduced,
                                   )
            if orig != new:
                print(body_name, scanner_name, loose_val)
                print(orig)
                print(new)
                print()

if __name__ == '__main__':
    test()
