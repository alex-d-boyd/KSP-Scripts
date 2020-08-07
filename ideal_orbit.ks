@lazyGlobal off.

parameter __RUN__ is false.

function gcd {
    parameter a, b.
    if a = 0 {
        return b.
    }
    if b = 0 {
        return a.
    }
    return gcd(b, mod(a, b)).
}

function cuberoot {
    parameter a.
    if a < 0 {
        return -cuberoot(-a).
    }

    return a^(1/3).
}

function get_roots_cubic {
    parameter a, b, c, d.

    local p is (3*a*c - (b^2))/(3*(a^2)).
    local q is (2*(b^3) - 9*a*b*c + 27*d*(a^2))/(27*(a^3)).

    local discriminant is (p^3)/27 + (q^2)/4.
    local dxt is b/(3*a).

    if discriminant > 0 {
        local disc_root is sqrt(discriminant).
        local t is cuberoot(-(q/2) + disc_root) + cuberoot(-(q/2) - disc_root).

        return t - dxt.
    }
    if discriminant = 0 {
        if p = 0 {
            return list(-dxt, -dxt, -dxt).
        }

        local x1 is 3*q/p - dxt.
        local x2 is -3*q/(2*p) - dxt.
        if x2 < x1{
            return list(x2, x2, x1).
        }
        return list(x1, x2, x2).
    }

    // inside of cuberoot is the complex number zr +/- zi*j
    local zr is -q/2.
    local zi is sqrt(-discriminant).

    // convert to polar form
    local r is sqrt(zr^2 + zi^2).
    local theta is arcTan2(zi, zr).

    // cuberoot z = z^(1/3)
    // z^n = r^n * (cos(n*theta) + j*sin(n*theta))
    set r to cuberoot(r).
    set theta to theta/3.

    // z + z* = 2 r cos(theta)
    local x1 is 2*r*cos(theta) - dxt.
    // complex number has 3 cube roots, separated by 120° in the complex plane
    local x2 is 2*r*cos(theta + 120) - dxt.
    local x3 is 2*r*cos(theta - 120) - dxt.

    if x1 > x2 {
        local y is x1.
        set x1 to x2.
        set x2 to y.
    }
    if x2 > x3 {
        local y is x2.
        set x2 to x3.
        set x3 to y.
    }
    if x1 > x2 {
        local y is x1.
        set x1 to x2.
        set x2 to y.
    }
    
    return list(x1, x2, x3).
}

function eccentricity_bounds {
    parameter radius.
    parameter sma.
    parameter alt_track.
    parameter alt_min.
    parameter alt_max.
    parameter alt_safe.

    local ecc_polar is 1 - (radius + alt_min)/sma.
    if ecc_polar < 0 { // cannot scan the poles
        return 0.
    }
    set ecc_polar to sqrt(ecc_polar).
    
    local ecc_periapsis is 1 - (radius + alt_safe)/sma.
    local ecc_apoapsis is (radius + alt_max)/sma - 1.

    local max_eccentricity is min(ecc_periapsis, min(ecc_polar, ecc_apoapsis)).
    local min_eccentricity is (radius + alt_track - sma)/sma.

    if min_eccentricity >= max_eccentricity {
        return 0.
    }

    local divergent_lower is min_eccentricity > alt_track/(2*alt_track + radius).
    local divergent_upper is alt_track > radius * sqrt((sma - radius)/sma).

    if divergent_lower or divergent_upper {
        local a is 4*alt_track*sma.
        local b is radius^2.
        local c is 2*alt_track*(radius - 2*sma).
        local d is alt_track^2.

        local roots is get_roots_cubic(a, b, c, d).
        if not roots:istype("list") {
            return 0.
        }

        if divergent_lower {
            set min_eccentricity to roots[1].
        }
        if divergent_upper {
            set max_eccentricity to min(max_eccentricity, roots[2]).
        }
    }

    return list(max(0, min_eccentricity), max_eccentricity).
}

function find_orbit {
    parameter target_body.
    parameter fov.
    parameter alt_min.
    parameter alt_best.
    parameter alt_max.
    parameter alt_safe.
    parameter loose_orbit.

    local scaled_fov is fov.
    if target_body:radius < kerbin:radius {
        set scaled_fov to scaled_fov * sqrt(kerbin:radius/target_body:radius).
    }

    local min_n_tracks is ceiling(180 / min(scaled_fov, 20)).

    local max_sma is target_body:radius + alt_max.
    local min_sma is target_body:radius + max(alt_min, alt_safe).

    local best_skip is -1.
    local best_orbit is 0.

    from {local n_tracks is min_n_tracks.} until n_tracks > 2*min_n_tracks step {set n_tracks to n_tracks+1.} do {
        local alt_track is 180 * alt_best / (n_tracks * scaled_fov).

        from {local skip is 0.} until skip > best_skip and best_skip > 0 step {set skip to skip+1.} do {
            if gcd(n_tracks, skip) = 1 {
                local period is skip * target_body:rotationPeriod / n_tracks.
                local sma is cuberoot(target_body:mu * period^2 / (4*constant:pi^2)).

                if sma > max_sma {
                    break.
                } else if sma > min_sma {
                    local ecc is eccentricity_bounds(target_body:radius, sma, alt_track, alt_min, alt_max, alt_safe).
                    if ecc <> 0 {
                        set best_orbit to lex(
                            "sma", sma,
                            "period", period,
                            "ecc_min", ecc[0],
                            "ecc_max", ecc[1],
                            "n_tracks", n_tracks,
                            "skip", skip
                        ).
                        set best_skip to skip - (choose 1 if loose_orbit else 0).
                        break.
                    }
                }
            }
        }
    }

    return best_orbit.
}

function input {
    parameter col, row, decal is ">>".
    local ti is terminal:input.
    local value is "".
    until false {
        local dispStr is decal + value.
        print dispStr:padRight(terminal:width) at (col,row).
        local c is ti:getChar.
        if c = ti:backspace {
            set value to value:subString(0, max(0, value:length-1)).
        } else if c = ti:enter {
            return value.
        } else if unChar(c) > 32 and unChar(c) < 127 {
            set value to value + c.
        }
    }
}

function main {
    set terminal:width to 40.
    set terminal:height to 30.

    clearScreen.

    local row is 0.
    print "Input Target Body:" at (0, row).
    local target_body is input(0, row+1).
    until bodyExists(target_body) {
        print ("cannot find " + target_body + ". Check spelling!"):padRight(terminal:width) at (0, row+1).
        set target_body to input(0, row+1).
    }
    set target_body to body(target_body). 

    set row to row + 3.
    print "Input Minimum Safe Altitude: (meters)" at (0, row).
    local alt_safe is input(0, row+1).
    until alt_safe:matchesPattern("^\d+$") {
        print "Int expected.":padRight(terminal:width) at (0, row+1).
        set alt_safe to input(0, row+1).
    }
    set alt_safe to alt_safe:toNumber.

    set row to row + 3.
    print "Input FOV: (degrees)" at (0, row).
    local fov is input(0, row+1).
    until fov:matchesPattern("^\d+(\.\d+)?$") {
        print "Float expected.":padRight(terminal:width) at (0, row+1).
        set fov to input(0, row+1).
    }
    set fov to fov:toNumber.

    set row to row + 3.
    print "Input Minimum Scan Altitude: (meters)" at (0, row).
    local alt_min is input(0, row+1).
    until alt_min:matchesPattern("^\d+$") {
        print "Int expected.":padRight(terminal:width) at (0, row+1).
        set alt_min to input(0, row+1).
    }
    set alt_min to alt_min:toNumber.

    set row to row + 3.
    print "Input Best Altitude: (meters)" at (0, row).
    local alt_best is input(0, row+1).
    until alt_best:matchesPattern("^\d+$") {
        print "Int expected.":padRight(terminal:width) at (0, row+1).
        set alt_best to input(0, row+1).
    }
    set alt_best to alt_best:toNumber.

    set row to row + 3.
    print "Input Maximum Altitude: (meters)" at (0, row).
    local alt_max is input(0, row+1).
    until alt_max:matchesPattern("^\d+$") {
        print "Int expected.":padRight(terminal:width) at (0, row+1).
        set alt_max to input(0, row+1).
    }
    set alt_max to alt_max:toNumber.

    set row to row + 3.
    print "Keep orbit loose? (y/n)" at (0, row).
    local loose_orbit is input(0, row+1).
    until loose_orbit:matchesPattern("^([yn]|yes|no)$") {
        print "":padRight(terminal:width) at (0, row+1).
        set loose_orbit to input(0, row+1).
    }
    set loose_orbit to loose_orbit:startsWith("y").

    set row to row + 3.
    local best_orbit is find_orbit(target_body, fov, alt_min, alt_best, alt_max, alt_safe, loose_orbit).
    for key in best_orbit:keys {
        print key + " = " + best_orbit[key] at (0, row).
        set row to row+1.
    }

    return row.
}

if __RUN__ {
    local row is main().
    print "start again? (y/n)" at (0, row+1).
    until input(0, row+2):starsWith("n") {
        set row to main().
        print "start again? (y/n)" at (0, row+1).
    }
}