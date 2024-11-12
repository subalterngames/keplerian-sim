#![cfg(test)]

use crate::{CompactOrbit, Orbit, OrbitTrait};
use std::f64::consts::PI;
use std::f64::{
    INFINITY as INF,
    NEG_INFINITY as NINF
};

type Vec3 = (f64, f64, f64);
type Vec2 = (f64, f64);

const ALMOST_EQ_TOLERANCE: f64 = 1e-6;
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

fn assert_almost_eq_vec3(a: Vec3, b: Vec3) {
    assert_almost_eq(a.0, b.0, "X coord");
    assert_almost_eq(a.1, b.1, "Y coord");
    assert_almost_eq(a.2, b.2, "Z coord");
}

fn assert_almost_eq_vec2(a: Vec2, b: Vec2) {
    assert_almost_eq(a.0, b.0, "X coord");
    assert_almost_eq(a.1, b.1, "Y coord");
}

fn assert_orbit_positions_3d(orbit: &impl OrbitTrait, tests: &[(f64, Vec3)]) {
    for (angle, expected) in tests.iter() {
        let pos = orbit.get_position_at_angle(*angle);
        assert_almost_eq_vec3(pos, *expected);
    }
}

fn assert_orbit_positions_2d(orbit: &impl OrbitTrait, tests: &[(f64, Vec2)]) {
    for (angle, expected) in tests.iter() {
        let pos = orbit.get_flat_position_at_angle(*angle);
        assert_almost_eq_vec2(pos, *expected);
    }
}

const ORBIT_POLL_ANGLES: usize = 255;
fn poll_orbit(orbit: &impl OrbitTrait) -> Vec<Vec3> {
    let mut vec: Vec<Vec3> = Vec::with_capacity(ORBIT_POLL_ANGLES);

    for i in 0..ORBIT_POLL_ANGLES {
        let angle = (i as f64) * 2.0 * PI / (ORBIT_POLL_ANGLES as f64);
        vec.push(orbit.get_position_at_angle(angle));
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
        (0.0 * PI, (1.0, 0.0, 0.0)),
        (0.5 * PI, (0.0, 1.0, 0.0)),
        (1.0 * PI, (-1.0, 0.0, 0.0)),
        (1.5 * PI, (0.0, -1.0, 0.0)),
        (2.0 * PI, (1.0, 0.0, 0.0)),
    ]);
}

#[test]
fn unit_orbit_angle_2d() {
    let orbit = unit_orbit();

    assert_orbit_positions_2d(&orbit, &[
        (0.0 * PI, (1.0, 0.0)),
        (0.5 * PI, (0.0, 1.0)),
        (1.0 * PI, (-1.0, 0.0)),
        (1.5 * PI, (0.0, -1.0)),
        (2.0 * PI, (1.0, 0.0)),
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
        ((1.0, 0.0), (1.0, 0.0, 0.0)),
        ((0.0, 1.0), (0.0, 0.0, 1.0)),
        ((-1.0, 0.0), (-1.0, 0.0, 0.0)),
        ((0.0, -1.0), (0.0, 0.0, -1.0)),
    ];

    for (point, expected) in tests.iter() {
        let transformed = orbit.tilt_flat_position(point.0, point.1);

        assert_almost_eq_vec3(transformed, *expected);
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

    assert_almost_eq_vec3(point_at_apoapsis, (-2.0, 0.0, 0.0));
    assert_almost_eq_vec3(point_at_periapsis, (1.0, 0.0, 0.0));
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

    assert_almost_eq_vec3(point_at_apoapsis, (-10000.0, 0.0, 0.0));
    assert_almost_eq_vec3(point_at_periapsis, (1.0, 0.0, 0.0));
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
    assert_almost_eq(point_near_infinity.2, 0.0, "Point near infinity should not be tilted");
    assert_almost_eq_vec3(point_at_periapsis, (1.0, 0.0, 0.0));

    let point_at_asymptote = orbit.get_position_at_angle(PI);

    assert!(point_at_asymptote.0.is_nan(), "X at asymptote should be undefined");
    assert!(point_at_asymptote.1.is_nan(), "Y at asymptote should be undefined");
    assert!(point_at_asymptote.2.is_nan(), "Z at asymptote should be undefined");
}

#[test]
fn orbit_conversion() {
    let original_orbit = Orbit::new_default();
    let original_positions = poll_orbit(&original_orbit);

    let compact_orbit = CompactOrbit::from(original_orbit.clone());
    let compact_positions = poll_orbit(&compact_orbit);

    let reexpanded_orbit = Orbit::from(compact_orbit.clone());
    let reexpanded_positions = poll_orbit(&reexpanded_orbit);

    for i in 0..original_positions.len() {
        assert_almost_eq_vec3(original_positions[i], compact_positions[i]);
        assert_almost_eq_vec3(original_positions[i], reexpanded_positions[i]);
    }
}