use crate::{Body, Orbit};

/// Returns Mercury, the closest planet to the Sun.
/// `include_orbit`: Whether to include the orbit of Mercury around the Sun.
pub fn mercury(include_orbit: bool) -> Body {
    let orbit = if include_orbit {
        // Source: Wikipedia
        Some(Orbit::new(
            6.982e7,
            4.6e7,
            7.005_f64.to_radians(),
            29.124_f64.to_radians(),
            48.331_f64.to_radians(),
            174.796_f64.to_radians()
        ))
    } else { None };

    return Body::new(
        "Mercury".to_string(),
        3.3011e23,
        2.4397e6,
        orbit
    );
}

/// Returns Venus, the second planet from the Sun.
/// `include_orbit`: Whether to include the orbit of Venus around the Sun.
pub fn venus(include_orbit: bool) -> Body {
    let orbit = if include_orbit {
        // Source: Wikipedia
        Some(Orbit::new(
            1.0894e8,
            1.0748e8,
            3.39458_f64.to_radians(),
            54.884_f64.to_radians(),
            76.680_f64.to_radians(),
            50.115_f64.to_radians()
        ))
    } else { None };

    return Body::new(
        "Venus".to_string(),
        4.8675e24,
        6.0518e6,
        orbit
    );
}

/// Returns Earth, the third planet from the Sun.
/// `include_orbit`: Whether to include the orbit of Earth around the Sun.
pub fn earth(include_orbit: bool) -> Body {
    let orbit = if include_orbit {
        // Source: Wikipedia
        Some(Orbit::new(
            1.52097597e11,
            1.47098450e11,
            0.00005_f64.to_radians(),
            114.20783_f64.to_radians(),
            -11.26064_f64.to_radians(),
            358.617_f64.to_radians()
        ))
    } else { None };

    return Body::new(
        "Earth".to_string(),
        5.972e24,
        6.371e6,
        orbit
    );
}