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


########################################################################################################################
# Notes on the derivation of the following section:
#
# Starting with the equation for radius from the true anomaly we have:
#
#   r(a, e, v) = a * (1 - e^2)/(1 + e*cos(v))       (1)
#
# where
#   r is the distance from the focus point
#   a is the semi-major axis
#   e is the eccentricity of the orbit
#   v is the true anomaly
#
# We know that in ScanSat, the size of the body as
#
#   f(R) = f0 * Rk / R      (2)
#
# and with the size of the body as
#
#   f(A, R) = f(R) * A / A0     (3)
#
# where
#   f is the field of view --- the width in degrees scanned by the satellite at the equator
#   f0 is the base field of view for the scanner
#   A0 is the altitude at which the base field of view applies
#   R is the radius of the body being scanned
#   Rk is the radius of the planet Kerbin
#
# given that altitude can be rewritten as
#
#   A(a, e, v, R) = r(a, e, v) - R        (4)
#
# our field of view is
#
#   f(a, e, v, R) = f(R) * (r(a, e, v) - R) / A0       (5)
#
# We also know that the circumference, C, of each latitude, l, changes as
#
#   C(l) = C(0) * cos(l)        (6)
#
# when the apoapsis is above the equator, we are interested in 'v' in [90°, 180°]. In this case we also get
#
#   C(v) = C(0) * -cos(v)       (8)
#
# and can use the substitution `v = w + 90°` to make
#
#   C(w) = C(90°) * sin(w)      (9)
#
#   r(a, e, w) = a * (1 - e^2)/(1 - e*sin(w))       (10)
#
# If we have a tracks required to sweep an angle of k, we require that
#
#   k * C(w) <= f(a, e, w, R)       (11)
#
# for all values of w, given fixed a, e, R
#
# This can be rearranged into
#
#   -eS sin(w)^2 + (S - eR) sin(w) - a(1 - e^2) + R <= 0        (12)
#
# where
#   S = k * A0/f0
#
# This can be checked using the quadratic formula where
#   a = -eS
#   b = S - eR
#   c = R - a(1 - e^2)
# Remembering that as our equation is quadratic in sin, and our domain is [0, 90°] any roots outside [0,1] are not roots
# of our original problem
#
# To find bounds on the values of e, further steps are required.
# To create an initial lower bound, we need the solutions to the equation
#
#   (-b + sqrt(b^2 - 4ac))/(2a) = 1     (13)
#
# This becomes the following quadratic
#
#   a e^2 - (S+R) e + S + R - a = 0     (14)
#
# which is easily solved with the quadratic formula.
#
# similarly, we create the upper bound using the other form of the quadratic formula, looking for roots equal to zero.
#
#   (-b - sqrt(b^2 - 4ac))/(2a) = 0     (15)
#
# which comes to the incredibly nice result of
#
#   e = sqrt((a-R)/a)       (16)
#
# but is actually just the eccentricity ar which we would crash into the planet at v=90°.
#
# using these roots, we check how many solutions are inside our range of [0,1] at each. If either has more than 1 we
# must replace it with the corresponding root of the discriminant
#
#   b^2 - 4ac = 0       (17)
#
# this is the cubic equation
#
#   4Sa e^3 + R^2 e^2 + 2S(R-2a) e + S^2 = 0        (18)
#
# this can be solved using Cardano's formula.
# Once solved, if there are three real solutions, replace the previous bounds with it's respective root if  required for
# that bound (i.e. there were multiple roots there in the previous step).
#
# These are now your upper and lower bounds, but they don't include any information about minimum or maximum orbit
# restrictions. Reduce the upper bound until the orbit fits inside the maximum altitude, minimum safe altitude, and
# satisfies the minimum polar altitude
########################################################################################################################

def calculate_eccentricity_lower_bound(semi_major_axis: float, body: Body, k: float) -> float:
    """
    Finds a lower bound on the eccentricity of the orbit.
    :param semi_major_axis: the semi-major axis of the orbit
    :param body: the body being orbited
    :param k: the
    :return:
    """
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


def calculate_eccentricity_bounds(body: Body, semi_major_axis: float, track_altitude: float) \
        -> Optional[Tuple[float, float]]:

    # calculates the eccentricity bounds where the fov from true anomaly graph touches in one place
    # lower bound touches at equator, upper bound touches at poles. done by solving the quadratic formula for
    # roots = 1, 0 respectively
    lower_bound = calculate_eccentricity_lower_bound(semi_major_axis, body, track_altitude)
    upper_bound = calculate_eccentricity_upper_bound(semi_major_axis, body)

    root_at_lower = find_root_at_lower_bound(semi_major_axis, body, lower_bound, track_altitude)
    root_at_upper = find_root_at_upper_bound(semi_major_axis, body, upper_bound, track_altitude)

    # epsilon used to compensate for floating point errors
    epsilon: float = 1e-5

    # checks additive root matches lower bound and subtractive matches upper bound. If not, the mismatch
    # needs to be replaced by the roots of the underlying cubic.
    if root_at_lower < (1-epsilon) or root_at_upper > epsilon:
        roots = solve_cubic(4 * track_altitude,
                            body.radius ** 2,
                            2 * track_altitude * (body.radius - 2 * semi_major_axis),
                            track_altitude ** 2)
        if type(roots) is not tuple:
            return None

        if root_at_lower < (1-epsilon):
            lower_bound = roots[1]
        if root_at_upper > epsilon:
            upper_bound = roots[2]
    return max(lower_bound, 0), min(upper_bound, 1)


def find_orbits(body: Body, field_of_view: float, track_altitude: float, best_altitude: float, max_altitude: float):

    # scale the field of view to the planet's radius as done in ScanSat
    scaled_fov: float = field_of_view
    if body.radius < Body.KERBIN_RADIUS:
        scaled_fov *= math.sqrt(Body.KERBIN_RADIUS / body.radius)

    # calculate the minimum number of ground tracks for 100% coverage at the equator
    min_ground_tracks: int = math.ceil(180 / min(scaled_fov, 20))

    # calculate orbital limits
    max_radius_apoapsis: float = body.radius + max_altitude  # cannot scan at apoapsis if above best altitude
    min_radius_periapsis: float = body.radius + body.safe_altitude  # if periapsis below safe altitude we will crash
    min_radius_polar: float = body.radius + track_altitude  # cannot scan at poles if radius below minimum altitude

    best_track_skip: int = -1
    best_orbit = None
    for n_tracks in range(min_ground_tracks, 2*min_ground_tracks):
        track_angle: float = 180 / n_tracks

        track_altitude = track_angle * best_altitude / scaled_fov

        for track_skip in coprimes_of(n_tracks, best_track_skip):

            orbit_period: float = track_skip * body.sidereal_day / n_tracks
            semi_major_axis: float = (body.mu * orbit_period**2 / (4 * pi**2)) ** (1/3)
            if semi_major_axis > max_radius_apoapsis:
                break
            if semi_major_axis < min_radius_periapsis or semi_major_axis < min_radius_polar:
                continue  # axis is too small, advance to increase

            eccentricity = calculate_eccentricity_bounds(body, semi_major_axis, track_altitude)

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