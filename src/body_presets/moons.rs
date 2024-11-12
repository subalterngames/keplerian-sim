use crate::{Body, Orbit};

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
            0.0
        ))
    } else { None };

    return Body::new(
        "Luna".to_string(),
        7.342e22,
        1.7371e6,
        orbit
    );
}

pub use luna as the_moon;