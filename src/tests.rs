#![cfg(test)]

use crate::{CompactOrbit, Orbit, OrbitTrait};
use std::f64::consts::{PI, TAU};

type Vec3 = (f64, f64, f64);
type Vec2 = (f64, f64);

const ALMOST_EQ_TOLERANCE: f64 = 1e-6;
const ORBIT_POLL_ANGLES: usize = 4096;

fn assert_almost_eq(a: f64, b: f64, what: &str) {
    let dist = (a - b).abs();
    let msg = format!(
        "Almost-eq assertion failed for '{what}'!\n\
        {a} and {b} has distance {dist}, which is more than max of {ALMOST_EQ_TOLERANCE}"
    );
    
    assert!(dist < ALMOST_EQ_TOLERANCE, "{msg}");
}

fn vec3_len(v: Vec3) -> f64 {
    return (
        v.0 * v.0 +
        v.1 * v.1 +
        v.2 * v.2
    ).sqrt();
}

fn assert_eq_vec3(a: Vec3, b: Vec3, what: &str) {
    assert_eq!(a.0.to_bits(), b.0.to_bits(), "X coord of {what}");
    assert_eq!(a.1.to_bits(), b.1.to_bits(), "Y coord of {what}");
    assert_eq!(a.2.to_bits(), b.2.to_bits(), "Z coord of {what}");
}

fn assert_eq_vec2(a: Vec2, b: Vec2, what: &str) {
    assert_eq!(a.0.to_bits(), b.0.to_bits(), "X coord of {what}");
    assert_eq!(a.1.to_bits(), b.1.to_bits(), "Y coord of {what}");
}

fn assert_almost_eq_vec3(a: Vec3, b: Vec3, what: &str) {
    assert_almost_eq(a.0, b.0, &("X coord of ".to_string() + what));
    assert_almost_eq(a.1, b.1, &("Y coord of ".to_string() + what));
    assert_almost_eq(a.2, b.2, &("Z coord of ".to_string() + what));
}

fn assert_almost_eq_vec2(a: Vec2, b: Vec2, what: &str) {
    assert_almost_eq(a.0, b.0, &("X coord of ".to_string() + what));
    assert_almost_eq(a.1, b.1, &("Y coord of ".to_string() + what));
}

fn assert_orbit_positions_3d(orbit: &impl OrbitTrait, tests: &[(&str, f64, Vec3)]) {
    for (what, angle, expected) in tests.iter() {
        let pos = orbit.get_position_at_angle(*angle);
        assert_almost_eq_vec3(pos, *expected, what);
    }
}

fn assert_orbit_positions_2d(orbit: &impl OrbitTrait, tests: &[(&str, f64, Vec2)]) {
    for (what, angle, expected) in tests.iter() {
        let pos = orbit.get_flat_position_at_angle(*angle);
        assert_almost_eq_vec2(pos, *expected, what);
    }
}

fn poll_orbit(orbit: &impl OrbitTrait) -> Vec<Vec3> {
    let mut vec: Vec<Vec3> = Vec::with_capacity(ORBIT_POLL_ANGLES);

    for i in 0..ORBIT_POLL_ANGLES {
        let angle = (i as f64) * 2.0 * PI / (ORBIT_POLL_ANGLES as f64);
        vec.push(orbit.get_position_at_angle(angle));
    }

    return vec;
}
fn poll_flat(orbit: &impl OrbitTrait) -> Vec<Vec2> {
    let mut vec: Vec<Vec2> = Vec::with_capacity(ORBIT_POLL_ANGLES);

    for i in 0..ORBIT_POLL_ANGLES {
        let angle = (i as f64) * 2.0 * PI / (ORBIT_POLL_ANGLES as f64);
        vec.push(orbit.get_flat_position_at_angle(angle));
    }

    return vec;
}
fn poll_transform(orbit: &impl OrbitTrait) -> Vec<Vec3> {
    let mut vec: Vec<Vec3> = Vec::with_capacity(ORBIT_POLL_ANGLES);

    for i in 0..ORBIT_POLL_ANGLES {
        let angle = (i as f64) * 2.0 * PI / (ORBIT_POLL_ANGLES as f64);
        vec.push(orbit.tilt_flat_position(1.0 * angle.cos(), 1.0 * angle.sin()));
    }

    return vec;
}
fn poll_eccentric_anomaly(orbit: &impl OrbitTrait) -> Vec<f64> {
    let mut vec: Vec<f64> = Vec::with_capacity(ORBIT_POLL_ANGLES);

    for i in 0..ORBIT_POLL_ANGLES {
        let angle = (i as f64) * 2.0 * PI / (ORBIT_POLL_ANGLES as f64);
        vec.push(orbit.get_eccentric_anomaly(angle));
    }

    return vec;
}

fn unit_orbit() -> Orbit {
    return Orbit::new(
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0
    );
}

#[test]
fn unit_orbit_angle_3d() {
    let orbit = unit_orbit();

    assert_orbit_positions_3d(&orbit, &[
        ("unit orbit 1", 0.0 * PI, (1.0, 0.0, 0.0)),
        ("unit orbit 2", 0.5 * PI, (0.0, 1.0, 0.0)),
        ("unit orbit 3", 1.0 * PI, (-1.0, 0.0, 0.0)),
        ("unit orbit 4", 1.5 * PI, (0.0, -1.0, 0.0)),
        ("unit orbit 5", 2.0 * PI, (1.0, 0.0, 0.0)),
    ]);
}

#[test]
fn unit_orbit_angle_2d() {
    let orbit = unit_orbit();

    assert_orbit_positions_2d(&orbit, &[
        ("unit orbit 1", 0.0 * PI, (1.0, 0.0)),
        ("unit orbit 2", 0.5 * PI, (0.0, 1.0)),
        ("unit orbit 3", 1.0 * PI, (-1.0, 0.0)),
        ("unit orbit 4", 1.5 * PI, (0.0, -1.0)),
        ("unit orbit 5", 2.0 * PI, (1.0, 0.0)),
    ]);
}

#[test]
fn unit_orbit_transformation() {
    // Test how the inclination and LAN tilts points in the orbit.
    // Since inclination is zero, it should not do anything.
    let orbit = unit_orbit();

    let tests = [
        (1.0, 1.0),
        (1.0, 0.0),
        (0.0, 1.0),
        (0.0, 0.0),
    ];

    for point in tests {
        let transformed = orbit.tilt_flat_position(point.0, point.1);

        assert_eq!(transformed.0, point.0);
        assert_eq!(transformed.1, point.1);
        assert_eq!(transformed.2, 0.0);
    }
}

#[test]
fn tilted_equidistant() {
    let orbit = Orbit::new(
        0.0,
        1.0,
        2.848915582093,
        1.9520945821,
        2.1834987325,
        0.69482153021
    );

    // Test for equidistance
    let points = poll_orbit(&orbit);

    for point in points {
        let distance = vec3_len(point);

        assert_almost_eq(distance, 1.0, "Distance");
    }
}

#[test]
fn tilted_90deg() {
    let orbit = Orbit::new(
        0.0,
        1.0,
        PI / 2.0,
        0.0,
        0.0,
        0.0
    );

    // Transform test
    let tests = [
        // Before and after transformation
        (("Vector 1"), (1.0, 0.0),  (1.0, 0.0, 0.0)),
        (("Vector 2"), (0.0, 1.0),  (0.0, 0.0, 1.0)),
        (("Vector 3"), (-1.0, 0.0), (-1.0, 0.0, 0.0)),
        (("Vector 4"), (0.0, -1.0), (0.0, 0.0, -1.0)),
    ];

    for (what, point, expected) in tests.iter() {
        let transformed = orbit.tilt_flat_position(point.0, point.1);

        assert_almost_eq_vec3(transformed, *expected, what);
    }
}

#[test]
fn apoapsis_of_two() {
    let orbit = Orbit::with_apoapsis(
        2.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0
    );
    
    let point_at_apoapsis = orbit.get_position_at_angle(PI);
    let point_at_periapsis = orbit.get_position_at_angle(0.0);

    assert_almost_eq_vec3(point_at_apoapsis, (-2.0, 0.0, 0.0), "Ap");
    assert_almost_eq_vec3(point_at_periapsis, (1.0, 0.0, 0.0), "Pe");
}

#[test]
fn huge_apoapsis() {
    let orbit = Orbit::with_apoapsis(
        10000.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0
    );

    let point_at_apoapsis = orbit.get_position_at_angle(PI);
    let point_at_periapsis = orbit.get_position_at_angle(0.0);

    assert_almost_eq_vec3(point_at_apoapsis, (-10000.0, 0.0, 0.0), "Ap");
    assert_almost_eq_vec3(point_at_periapsis, (1.0, 0.0, 0.0), "Pe");
}

const JUST_BELOW_ONE: f64 = 0.9999999999999999;

#[test]
fn almost_parabolic() {
    let orbit = Orbit::new(
        // The largest f64 below 1
        JUST_BELOW_ONE,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0
    );

    let eccentric_anomalies = poll_eccentric_anomaly(&orbit);

    for ecc in eccentric_anomalies {
        assert!(
            ecc.is_finite(),
            "Eccentric anomaly algorithm instability at near-parabolic edge case"
        );
    }

    let positions = poll_flat(&orbit);

    for pos in positions {
        assert!(
            pos.0.is_finite() && pos.1.is_finite(),
            "2D position algorithm instability at near-parabolic edge case"
        );
    }

    let positions = poll_orbit(&orbit);

    for pos in positions {
        assert!(
            pos.0.is_finite() && pos.1.is_finite() && pos.2.is_finite(),
            "3D position algorithm instability at near-parabolic edge case"
        );
    }

    let position_at_periapsis = orbit.get_position_at_angle(TAU);

    assert_almost_eq_vec3(
        position_at_periapsis,
        (1.0, 0.0, 0.0),
        "Periapsis"
    );

    let position_at_apoapsis = orbit.get_position_at_angle(PI);

    assert!(position_at_apoapsis.0.abs() > 1e12, "Apoapsis is not far enough");
}

#[test]
fn parabolic() {
    let orbit = Orbit::new(
        1.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0
    );

    let point_near_infinity = orbit.get_position_at_angle(PI - 1e-7);
    let point_at_periapsis = orbit.get_position_at_angle(0.0);

    assert!(vec3_len(point_near_infinity) > 1e9, "Point near infinity is not far enough");
    assert!(point_near_infinity.1.abs() > 0.0, "Y coord near infinity should move a little");
    assert_almost_eq(point_near_infinity.2, 0.0, "Point near infinity should be flat");
    assert_almost_eq_vec3(point_at_periapsis, (1.0, 0.0, 0.0), "Pe");

    let point_at_asymptote = orbit.get_position_at_angle(PI);

    assert!(point_at_asymptote.0.is_nan(), "X at asymptote should be undefined");
    assert!(point_at_asymptote.1.is_nan(), "Y at asymptote should be undefined");
    assert!(point_at_asymptote.2.is_nan(), "Z at asymptote should be undefined");
}

fn orbit_conversion_base_test(orbit: Orbit, what: &str) {
    let compact_orbit = CompactOrbit::from(orbit.clone());
    let reexpanded_orbit = Orbit::from(compact_orbit.clone());
    
    let compact_message = format!("Original / Compact ({what})");
    let reexpanded_message = format!("Compact /  Reexpanded ({what})");

    {
        let original_transforms = poll_transform(&orbit);
        let compact_transforms = poll_transform(&compact_orbit);
        let reexpanded_transforms = poll_transform(&reexpanded_orbit);

        let compact_message = format!("{compact_message} (transform)");
        let reexpanded_message = format!("{reexpanded_message} (transform)");

        for i in 0..original_transforms.len() {
            assert_eq_vec3(original_transforms[i], compact_transforms[i], &compact_message);
            assert_eq_vec3(original_transforms[i], reexpanded_transforms[i], &reexpanded_message);
        }
    }
    {
        let original_ecc = poll_eccentric_anomaly(&orbit);
        let compact_ecc = poll_eccentric_anomaly(&compact_orbit);
        let reexpanded_ecc = poll_eccentric_anomaly(&reexpanded_orbit);

        for i in 0..original_ecc.len() {
            assert_eq!(original_ecc[i], compact_ecc[i], "{compact_message} (eccentric anomaly)");
            assert_eq!(original_ecc[i], reexpanded_ecc[i], "{reexpanded_message} (eccentric anomaly)");
        }
    }
    {
        let original_true = {
            let mut vec: Vec<f64> = Vec::with_capacity(ORBIT_POLL_ANGLES);

            for i in 0..ORBIT_POLL_ANGLES {
                let angle = (i as f64) * 2.0 * PI / (ORBIT_POLL_ANGLES as f64);
                let true_anomaly = orbit.get_true_anomaly(angle);
                vec.push(true_anomaly);
            }

            vec
        };
        let compact_true = {
            let mut vec: Vec<f64> = Vec::with_capacity(ORBIT_POLL_ANGLES);

            for i in 0..ORBIT_POLL_ANGLES {
                let angle = (i as f64) * 2.0 * PI / (ORBIT_POLL_ANGLES as f64);
                let true_anomaly = compact_orbit.get_true_anomaly(angle);
                vec.push(true_anomaly);
            }

            vec
        };
        let reexpanded_true = {
            let mut vec: Vec<f64> = Vec::with_capacity(ORBIT_POLL_ANGLES);

            for i in 0..ORBIT_POLL_ANGLES {
                let angle = (i as f64) * 2.0 * PI / (ORBIT_POLL_ANGLES as f64);
                vec.push(reexpanded_orbit.get_true_anomaly(angle));
            }

            vec
        };

        for i in 0..original_true.len() {
            assert_eq!(original_true[i], compact_true[i], "{compact_message} (true anomaly) (i={i})");
            assert_eq!(original_true[i], reexpanded_true[i], "{reexpanded_message} (true anomaly) (i={i})");
        }
    }
    {
        let original_eccentricity = orbit.get_eccentricity();
        let compact_eccentricity = compact_orbit.eccentricity;
        let reexpanded_eccentricity = reexpanded_orbit.get_eccentricity();

        assert_eq!(original_eccentricity, compact_eccentricity, "{compact_message} (eccentricity)");
        assert_eq!(original_eccentricity, reexpanded_eccentricity, "{reexpanded_message} (eccentricity)");
    }
    {
        let original_flat = poll_flat(&orbit);
        let compact_flat = poll_flat(&compact_orbit);
        let reexpanded_flat = poll_flat(&reexpanded_orbit);

        let compact_message = format!("{compact_message} (flat)");
        let reexpanded_message = format!("{reexpanded_message} (flat)");

        for i in 0..original_flat.len() {
            assert_eq_vec2(original_flat[i], compact_flat[i], &compact_message);
            assert_eq_vec2(original_flat[i], reexpanded_flat[i], &reexpanded_message);
        }
    }
    {
        let original_positions = poll_orbit(&orbit);
        let compact_positions = poll_orbit(&compact_orbit);
        let reexpanded_positions = poll_orbit(&reexpanded_orbit);

        let compact_message = format!("{compact_message} (position)");
        let reexpanded_message = format!("{reexpanded_message} (position)");

        for i in 0..original_positions.len() {
            assert_eq_vec3(original_positions[i], compact_positions[i], &compact_message);
            assert_eq_vec3(original_positions[i], reexpanded_positions[i], &reexpanded_message);
        }
    }
}

#[test]
fn orbit_conversions() {
    let orbits = [
        (
            "Unit orbit",
            Orbit::new(
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0
            )
        ),
        (
            "Mildly eccentric orbit",
            Orbit::new(
                0.39,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0
            )
        ),
        (
            "Very eccentric orbit",
            Orbit::new(
                0.99,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0
            )
        ),
        (
            "Parabolic trajectory",
            Orbit::new(
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0
            )
        ),
        (
            "Barely hyperbolic trajectory",
            Orbit::new(
                1.01,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0
            )
        ),
        (
            "Very hyperbolic trajectory",
            Orbit::new(
                9.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0
            )
        ),
        (
            "Extremely hyperbolic trajectory",
            Orbit::new(
                100.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0
            )
        ),
        (
            "Tilted orbit",
            Orbit::new(
                0.0,
                1.0,
                2.848915582093,
                1.9520945821,
                2.1834987325,
                0.69482153021
            )
        ),
        (
            "Tilted eccentric",
            Orbit::new(
                0.39,
                1.0,
                2.848915582093,
                1.9520945821,
                2.1834987325,
                0.69482153021
            )
        ),
        (
            "Tilted parabolic",
            Orbit::new(
                1.0,
                1.0,
                2.848915582093,
                1.9520945821,
                2.1834987325,
                0.69482153021
            )
        ),
        (
            "Tilted hyperbolic",
            Orbit::new(
                1.8,
                1.0,
                2.848915582093,
                1.9520945821,
                2.1834987325,
                0.69482153021
            )
        ),
    ];

    for (what, orbit) in orbits.iter() {
        orbit_conversion_base_test(orbit.clone(), what);
    }
}