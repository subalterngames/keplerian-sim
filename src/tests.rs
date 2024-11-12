#![cfg(test)]

use crate::{CompactOrbit, Orbit, OrbitTrait};
use std::f64::consts::PI;

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

const ORBIT_POLL_ANGLES: [f64; 5] = [
    0.0 * PI,
    0.5 * PI,
    1.0 * PI,
    1.5 * PI,
    2.0 * PI
];
fn poll_orbit(orbit: &impl OrbitTrait) -> Vec<Vec3> {
    let mut vec: Vec<Vec3> = Vec::with_capacity(ORBIT_POLL_ANGLES.len());

    for angle in ORBIT_POLL_ANGLES.iter() {
        vec.push(orbit.get_position_at_time(*angle));
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