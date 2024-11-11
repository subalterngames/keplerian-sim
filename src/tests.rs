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

#[test]
fn unit_orbit_3d() {
    let orbit = Orbit::new(
        1.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
    );

    assert_orbit_positions_3d(&orbit, &[
        (0.0 * PI, (1.0, 0.0, 0.0)),
        (0.5 * PI, (0.0, 1.0, 0.0)),
        (1.0 * PI, (-1.0, 0.0, 0.0)),
        (1.5 * PI, (0.0, -1.0, 0.0)),
        (2.0 * PI, (1.0, 0.0, 0.0)),
    ]);
}

#[test]
fn unit_orbit_2d() {
    let orbit = Orbit::new(
        1.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
    );

    assert_orbit_positions_2d(&orbit, &[
        (0.0 * PI, (1.0, 0.0)),
        (0.5 * PI, (0.0, 1.0)),
        (1.0 * PI, (-1.0, 0.0)),
        (1.5 * PI, (0.0, -1.0)),
        (2.0 * PI, (1.0, 0.0)),
    ]);
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