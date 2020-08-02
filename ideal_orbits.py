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
    k: int = 1
    while k <= limit or limit == -1:
        if math.gcd(n, k) == 1:
            yield k
        k += 1


def get_fov_at_radius(base_fov: float, best_altitude: float, radius: float, body: Body) -> float:
    scaled_fov = base_fov * (radius - body.radius)/best_altitude
    if body.radius < Body.KERBIN_RADIUS:
        scaled_fov *= math.sqrt(Body.KERBIN_RADIUS/body.radius)
    return scaled_fov


def calculate_eccentricity_lower_bound(semi_major_axis: float, body: Body, k: float) -> float:
    a = semi_major_axis
    b = -(k + body.radius)
    c = k + body.radius - semi_major_axis

    return (- b - math.sqrt(b**2 - 4*a*c))/(2*a)


def find_root_at_lower_bound(semi_major_axis: float, body: Body, bound: float, k: float) -> float:
    a = -bound * k
    b = k - bound * body.radius
    c = body.radius - semi_major_axis*(1 - bound**2)

    return (-b + math.sqrt(b**2 - 4*a*c))/(2*a)


def calculate_eccentricity_upper_bound(semi_major_axis: float, body: Body) -> float:
    return math.sqrt((semi_major_axis - body.radius) / semi_major_axis)


def find_root_at_upper_bound(semi_major_axis: float, body: Body, bound: float, k: float) -> float:
    a = -bound * k
    b = k - bound * body.radius
    c = body.radius - semi_major_axis*(1 - bound**2)

    return (-b - math.sqrt(b**2 - 4*a*c))/(2*a)


def rect_to_polar(re: float, im: float) -> Tuple[float, float]:
    r: float = math.sqrt(re**2 + im**2)
    a: float = math.atan2(im, re)
    return r, a


def polar_to_rect(r: float, a: float) -> Tuple[float, float]:
    re: float = r*math.cos(a)
    im: float = r*math.sin(a)
    return re, im


def complex_power(re: float, im: float, p: float) -> Tuple[float, float]:
    r, a = rect_to_polar(re, im)
    r **= p
    a *= p
    return polar_to_rect(r, a)


def rotate(re: float, im: float, angle: float) -> Tuple[float, float]:
    r, a = rect_to_polar(re, im)
    return polar_to_rect(r, a+angle)


def solve_cubic(a: float, b: float, c: float, d: float) -> Union[float, Tuple[float, float, float]]:
    """
    Solves a cubic of the form ax^3 + bx^2 + cx + d using Cardano's formula by converting to a depressed cubic.
    does not currently handle the two root case, as it is not required by the problem
    :return: the found roots to the cubic
    """
    # convert to depressed cubic px^3 + qx + r
    p = -b/(3*a)
    q = p**3 + (b*c - 3*a*d)/(6 * a**2)
    r = c/(3*a)

    discriminant = q**2 + (r - p**2)**3
    if discriminant >= 0:  # TODO, two roots when discriminant = 0, not needed for this problem
        root = math.sqrt(discriminant)
        left = q + root
        right = q - root
        return left**(1/3) + right**(1/3) + p

    re = q
    im = math.sqrt(-discriminant)

    # find all three cube roots by rotation 120°
    re, im = complex_power(re, im, 1/3)
    # both cube roots are complex conjugates, therefore we need only find one and use twice the real part
    root1 = 2*re + p
    re, im = rotate(re, im, 2*pi/3)
    root2 = 2*re + p
    re, im = rotate(re, im, 2*pi/3)
    root3 = 2*re + p

    # order the roots smallest to largest
    if root1 > root2:
        root1, root2 = root2, root1
    if root2 > root3:
        root2, root3 = root3, root2
    if root1 > root2:
        root1, root2 = root2, root1

    return root1, root2, root3


def calculate_eccentricity_bounds(body: Body, semi_major_axis: float, track_angle: float, fov_at_sma: float)\
        -> Optional[Tuple[float, float]]:
    k = track_angle * (semi_major_axis - body.radius) / fov_at_sma

    # calculates the eccentricity bounds where the fov from true anomaly graph touches in one place
    # lower bound touches at equator, upper bound touches at poles. done by solving the quadratic formula for
    # roots = 1, 0 respectively
    lower_bound = calculate_eccentricity_lower_bound(semi_major_axis, body, k)
    upper_bound = calculate_eccentricity_upper_bound(semi_major_axis, body)

    root_at_lower = find_root_at_lower_bound(semi_major_axis, body, lower_bound, k)
    root_at_upper = find_root_at_upper_bound(semi_major_axis, body, upper_bound, k)

    # epsilon used to compensate for floating point errors
    epsilon: float = 1e-5

    # checks additive root matches lower bound and subtractive matches upper bound. If not, the mismatch
    # needs to be replaced by the roots of the underlying cubic.
    if root_at_lower < (1-epsilon) or root_at_upper > epsilon:
        roots = solve_cubic(4*k, body.radius**2, 2*k*body.radius - 4*k*semi_major_axis, k**2)
        if type(roots) is not tuple:
            return None

        if root_at_lower < (1-epsilon):
            lower_bound = roots[1]
        if root_at_upper > epsilon:
            upper_bound = roots[2]
    return max(lower_bound, 0), min(upper_bound, 1)


def find_orbits(body: Body, field_of_view: float, min_altitude: float, best_altitude: float, max_altitude: float):

    # scale the field of view to the planet's radius as done in ScanSat
    scaled_fov: float = field_of_view
    if body.radius < Body.KERBIN_RADIUS:
        scaled_fov *= math.sqrt(Body.KERBIN_RADIUS / body.radius)
    scaled_fov = min(scaled_fov, 20)

    # calculate the minimum number of ground tracks for 100% coverage at the equator
    min_ground_tracks: int = math.ceil(180 / scaled_fov)

    # calculate orbital limits
    max_radius_apoapsis: float = body.radius + max_altitude  # cannot scan at apoapsis if above best altitude
    min_radius_periapsis: float = body.radius + body.safe_altitude  # if periapsis below safe altitude we will crash
    min_radius_polar: float = body.radius + min_altitude  # cannot scan at poles if radius below minimum altitude

    best_track_skip: int = -1
    best_orbit = None
    for n_tracks in range(min_ground_tracks, 2*min_ground_tracks):
        track_angle: float = 180 / n_tracks

        for track_skip in coprimes_of(n_tracks, best_track_skip):

            orbit_period: float = track_skip * body.sidereal_day / n_tracks
            semi_major_axis: float = (body.mu * orbit_period**2 / (4 * pi**2)) ** (1/3)
            if semi_major_axis < min_radius_periapsis or semi_major_axis < min_radius_polar:
                continue  # axis is too small, advance to increase
            if semi_major_axis > max_radius_apoapsis:
                break

            fov_at_radius: float = get_fov_at_radius(field_of_view, best_altitude, semi_major_axis, body)
            eccentricity = calculate_eccentricity_bounds(body, semi_major_axis, track_angle, fov_at_radius)

            if eccentricity is None:
                continue  # no valid orbits, continue

            low_eccentricity, high_eccentricity = eccentricity

            max_ecc_poles = math.sqrt(1 - min_radius_polar/semi_major_axis)
            max_ecc_periapsis = 1 - min_radius_periapsis/semi_major_axis
            max_ecc_apoapsis = max_radius_apoapsis/semi_major_axis - 1

            high_eccentricity = min(high_eccentricity, max_ecc_poles, max_ecc_periapsis, max_ecc_apoapsis)

            if high_eccentricity < low_eccentricity:
                continue

            best_orbit = {
                "sma": semi_major_axis,
                "min_ecc": low_eccentricity,
                "max_ecc": high_eccentricity,
                "tracks": n_tracks,
                "skip": track_skip
            }
            best_track_skip = track_skip
            break

    return best_orbit


def main():
    bodies = {
        "kerbin": Body(3.5316e12, Body.KERBIN_RADIUS, 21549.425, 87000),
        "mun": Body(6.5138398e10, 200000, 138984.38, 15000),
        "minmus": Body(1.7658e9, 60000, 40400, 8000)
    }
    body = input("input target body:")
    while body not in bodies:
        print(f"body ${body} not found, expected one of", bodies.keys)
        body = input("input target body:")

    fov = float(input("fov: "))
    minAlt = float(input("min alt: "))
    bestAlt = float(input("best alt: "))
    maxAlt = float(input("max alt: "))

    best_orbit = find_orbits(bodies[body], fov, minAlt, bestAlt, maxAlt)
    print(best_orbit)


if __name__ == '__main__':
    main()