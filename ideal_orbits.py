import math
from math import pi
from typing import Iterator, Tuple, Union, Optional


class Body:
    KERBIN_RADIUS: float = 600000

    def __init__(self, mu: float, radius: float, sidereal_day: float, safe_altitude: float):
        self._mu: float = mu
        self._radius: float = radius
        self._sidereal_day: float = sidereal_day
        self._safe_altitude: float = safe_altitude

    @property
    def mu(self) -> float:
        return self._mu

    @property
    def radius(self) -> float:
        return self._radius

    @property
    def sidereal_day(self) -> float:
        return self._sidereal_day

    @property
    def safe_altitude(self) -> float:
        return self._safe_altitude


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


def find_eccentricity_bounds(body: Body,
                             semi_major_axis: float,
                             track_altitude: float,
                             min_altitude: float,
                             max_altitude: float) -> Optional[Tuple[float, float]]:

    ecc_polar = math.sqrt(max(1 - (body.radius + min_altitude)/semi_major_axis, 0))
    ecc_periapsis = 1 - (body.radius + body.safe_altitude)/semi_major_axis
    ecc_apoapsis = (body.radius + max_altitude)/semi_major_axis - 1

    ecc_max = min(ecc_polar, ecc_periapsis, ecc_apoapsis)

    a = semi_major_axis
    b = -(track_altitude + body.radius)
    c = track_altitude + body.radius - semi_major_axis
    ecc_min = (- b - math.sqrt(b**2 - 4*a*c))/(2*a)

    if ecc_min >= ecc_max:
        return None

    ecc_warning = track_altitude/(2*track_altitude + body.radius)
    if ecc_min > ecc_warning:
        # pretty sure I proved this can't happen, but warning just in case
        print("\x1b[31mWarning: minimum eccentricity above divergent bound!\x1b[0m")

    return ecc_min, ecc_max


def find_orbits(body: Body, field_of_view: float, min_altitude: float, best_altitude: float, max_altitude: float, loose: bool):
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
    if body.radius < Body.KERBIN_RADIUS:
        scaled_fov *= math.sqrt(Body.KERBIN_RADIUS / body.radius)

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

            eccentricity = find_eccentricity_bounds(body, semi_major_axis, track_altitude, min_altitude, max_altitude)

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
    bodies = {
        "kerbin": Body(3.5316e12, Body.KERBIN_RADIUS, 21549.425, 70000),
        "mun": Body(6.5138398e10, 200000, 138984.38, 10000),
        "minmus": Body(1.7658e9, 60000, 40400, 6000),
        "gilly": Body(8289449.8, 13000, 28255, 8000),
        "jool": Body(2.82528e14, 6000000, 36000, 200000)
    }
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


if __name__ == '__main__':
    main()
