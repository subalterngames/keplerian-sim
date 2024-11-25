use keplerian_rust::{Orbit, OrbitTrait};
use std::hint::black_box;
use criterion::{criterion_group, criterion_main, Criterion};

const POLL_ANGLES: usize = 1024;

#[inline(always)]
fn normal(orbit: &Orbit) {
    let multiplier = std::f64::consts::TAU / POLL_ANGLES as f64;

    for i in 0..POLL_ANGLES {
        let angle = i as f64 * multiplier;
        black_box(orbit.get_eccentric_anomaly(black_box(angle)));
    }
}

#[inline(always)]
fn experimental(orbit: &Orbit) {
    let multiplier = std::f64::consts::TAU / POLL_ANGLES as f64;

    for i in 0..POLL_ANGLES {
        let angle = i as f64 * multiplier;
        black_box(orbit.get_eccentric_anomaly_hyperbolic_experimental(
            black_box(angle)));
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    let orbit = Orbit::new(
        2.86,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
    );

    c.bench_function("_EXPERIMENTAL hyp ecc normal", |b| b.iter(||
        normal(black_box(&orbit))
    ));

    c.bench_function("_EXPERIMENTAL hyp ecc exp", |b| b.iter(||
        experimental(black_box(&orbit))
    ));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);