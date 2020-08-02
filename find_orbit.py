from math import gcd, ceil, pi, sqrt

KERBIN_RADIUS = 600000

def coprimes_of(n):
    k = 1
    while True:
        if gcd(n, k) == 1:
            yield k
        k += 1

def get_scaled_fov(fov, body_radius):
    scaled_fov = fov * sqrt(max(KERBIN_RADIUS/body_radius, 1))
    return min(scaled_fov, 20)

def get_min_alt_from_sweep_fov(fov, best_alt, fov_base, body_radius):
    return (fov/fov_base)*sqrt(min(1, body_radius/KERBIN_RADIUS))*best_alt


def find_best_orbit(mu, day, body_radius, fov, min_alt, best_alt, max_alt):
    scaled_fov = get_scaled_fov(fov, body_radius)
    min_n_passes = ceil(180/scaled_fov)

    min_T = None
    sma = None
    np = None
    for n_passes in range(min_n_passes, min_n_passes*2):
        min_fov = 180/n_passes
        scaled_alt = get_min_alt_from_sweep_fov(min_fov, best_alt, fov, body_radius)
        scaled_alt = max(scaled_alt, min_alt)
        for co in coprimes_of(n_passes):
            if np is not None and co > np[0]:
                break

            t =  co * day / n_passes
            a = (mu * t**2 / (4 * pi**2)) ** (1/3)
            if a > max_alt + body_radius:
                break

            if a > scaled_alt + body_radius:
                T = co * day
                if min_T is None or T < min_T or (T == min_T and a < sma):
                    min_T = T
                    sma = a
                    np = (co, n_passes)
                break

    return min_T, sma, np

def main():
    bodies = {
        "kerbin": {
            "mu": 3.5316e12,
            "day": 21549.425,
            "radius": KERBIN_RADIUS
        },
        "mun": {
            "mu": 6.5138398e10,
            "day": 138984.38,
            "radius": 200000
        },
        "minmus": {
            "mu": 1.7658e9,
            "day": 40400,
            "radius": 60000
        }
    }

    body = input("input target body:")
    while body not in bodies:
        print(f"body ${body} not found, expected one of", bodies.keys)
        body = input("input target body:")

    fov = float(input("fov: "))
    minAlt = float(input("min alt: "))
    bestAlt = float(input("best alt: "))
    maxAlt = float(input("max alt: "))
    
    props = bodies[body]

    print(find_best_orbit(props["mu"], props["day"], props["radius"], fov, minAlt, bestAlt, maxAlt))


if __name__ == "__main__":
    main()