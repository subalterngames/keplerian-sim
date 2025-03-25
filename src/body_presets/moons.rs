//! This module contains presets for common moons, also known as
//! natural satellites.
//!
//! A natural satellite is, in the most common usage, an astronomical body
//! that orbits a planet, dwarf planet, or small Solar System body
//! (or sometimes another natural satellite).
//!
//! \- [Wikipedia](https://en.wikipedia.org/wiki/Natural_satellite)

use crate::{Body, Orbit, OrbitTrait};

/// Returns the Moon, the only natural satellite of Earth.
/// `include_orbit`: Whether to include the orbit of the Moon around the Earth.
pub fn luna(include_orbit: bool) -> Body {
    let orbit = if include_orbit {
        // Source: Wikipedia
        // ...except for the argument of periapsis and the longitude of ascending node,
        // where zero is used because Wikipedia doesn't provide the values,
        // as the Moon precesses too fast for them to be useful.
        Some(Orbit::with_apoapsis(
            405400.0,
            362600.0,
            5.145_f64.to_radians(),
            0.0,
            0.0,
            0.0,
            MASS_EARTH,
        ))
    } else {
        None
    };

    Body::new("Luna".to_string(), 7.342e22, 1.7371e6, orbit)
}

pub use luna as the_moon;

use super::{MASS_EARTH, MASS_ERIS, MASS_PLUTO, MASS_QUAOAR};

/// Returns (50000) Quaoar I, a.k.a. Weywot, the moon of the dwarf planet Quaoar.
/// `include_orbit`: Whether to include the orbit of Weywot around Quaoar.
pub fn weywot(include_orbit: bool) -> Body {
    let orbit = if include_orbit {
        // Source: Wikipedia
        Some(Orbit::new(
            0.056,
            1.3289e7,
            // I only found inclination to the ecliptic.
            // I couldn't find one in relation to Quaoar's equator.
            0.0,
            335_f64.to_radians(),
            1_f64.to_radians(),
            // I couldn't find the mean anomaly
            0.0,
            MASS_QUAOAR,
        ))
    } else {
        None
    };

    Body::new(
        "Weywot".to_string(),
        // Weywot's mass has not been measured.
        // I extrapolated from the mean density of Quaoar:
        // ~1.7 g/cm^3
        // Weywot's radius is about 1e5 meters
        // Therefore its volume is:
        // V = 4/3 * pi * (1e5)^3 in cubic meters
        // V = 4.18879e15 m^3
        // Its mass is then:
        // M = 1.7 g/cm^3 * 4.18879e15 m^3
        // M = 1.7e3 kg/m^3 * 4.18879e15 m^3
        // M = 7.12e18 kg
        // (approximately)
        7.12e18,
        1e5,
        orbit,
    )
}

/// Returns (134340) Pluto I, a.k.a., Charon, the largest moon orbiting Pluto.
/// `include_orbit`: Whether to include the orbit of Charon around Pluto.
pub fn charon(include_orbit: bool) -> Body {
    let orbit = if include_orbit {
        // Source: Wikipedia
        Some(Orbit::with_apoapsis(
            1.959892e7,
            1.959261e7,
            0.08_f64.to_radians(),
            // Could not find number for arg. of pe.
            0.0,
            223.046_f64.to_radians(),
            // Could not find number for mean anomaly
            0.0,
            MASS_PLUTO,
        ))
    } else {
        None
    };

    Body::new("Charon".to_string(), 1.5897e21, 6.06e5, orbit)
}

/// Returns (136199) Eris I Dysnomia, the moon of the dwarf planet Eris.
/// `include_orbit`: Whether to include the orbit of Dysnomia around Eris.
pub fn dysnomia(include_orbit: bool) -> Body {
    let orbit = if include_orbit {
        // Source: Wikipedia
        Some(Orbit::new(
            0.0062,
            3.7273e7,
            0.0,
            180.83_f64.to_radians(),
            126.17_f64.to_radians(),
            // Could not find mean anomaly number
            0.0,
            MASS_ERIS,
        ))
    } else {
        None
    };

    Body::new(
        "Dysnomia".to_string(),
        // Apparently not very precise;
        // it's plus or minus 5.7e19 kg
        8.2e19,
        6.15e5 / 2.0,
        orbit,
    )
}
