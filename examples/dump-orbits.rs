use keplerian_sim::{OrbitTrait, Orbit};
use std::fs;

const SIMULATION_TICKS: usize = 10_000;
const CSV_PATH: &str = "out/output-orbit-dump.csv";

struct OrbitWrapper {
    name: &'static str,
    orbit: Orbit,
    time_mult: f64,
    time_offset: f64,
}

struct OrbitLog<'a> {
    name: &'a str,

    iter: usize,
    time: f64,
    mean_anom: f64,
    ecc_anom: f64,
    angle: f64,

    x: f64,
    y: f64,
    z: f64,
    
    flat_x: f64,
    flat_y: f64,

    altitude: f64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let orbits = [
        OrbitWrapper {
            name: "The Moon",
            orbit: Orbit::with_apoapsis(
                405400.0,
                362600.0,
                5.145_f64.to_radians(),
                0.0,
                0.0,
                0.0
            ),
            time_mult: 1.0,
            time_offset: 0.0,
        },
        OrbitWrapper {
            name: "Elliptic",
            orbit: Orbit::new(
                0.85,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 1.0,
            time_offset: 0.0,
        },
        OrbitWrapper {
            name: "Parabolic",
            orbit: Orbit::new(
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 10.0,
            time_offset: -5.0,
        },
        OrbitWrapper {
            name: "Barely Hyperbolic 1",
            orbit: Orbit::new(
                1.001,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 10.0,
            time_offset: -5.0,
        },
        OrbitWrapper {
            name: "Barely Hyperbolic 2",
            orbit: Orbit::new(
                1.002,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 10.0,
            time_offset: -5.0,
        },
        OrbitWrapper {
            name: "Barely Hyperbolic 3",
            orbit: Orbit::new(
                1.005,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 10.0,
            time_offset: -5.0,
        },
        OrbitWrapper {
            name: "Barely Hyperbolic 4",
            orbit: Orbit::new(
                1.01,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 10.0,
            time_offset: -5.0,
        },
        OrbitWrapper {
            name: "Barely Hyperbolic 5",
            orbit: Orbit::new(
                1.02,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 10.0,
            time_offset: -5.0,
        },
        OrbitWrapper {
            name: "Hyperbolic",
            orbit: Orbit::new(
                2.5,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 10.0,
            time_offset: -5.0,
        },
        OrbitWrapper {
            name: "Very Hyperbolic",
            orbit: Orbit::new(
                12.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 50.0,
            time_offset: -25.0,
        },
        OrbitWrapper {
            name: "Extreme Hyperbolic",
            orbit: Orbit::new(
                1000.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 1000.0,
            time_offset: -500.0,
        },
        OrbitWrapper {
            name: "Voyager",
            orbit: Orbit::new(
                2.0,
                1.496e11,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
            time_mult: 100.0,
            time_offset: -5.0,
        }
    ];

    let mut logs: Vec<OrbitLog> = Vec::with_capacity(orbits.len() * SIMULATION_TICKS);

    println!(">> Begin simulations.");

    for named_orbit in orbits {
        let (name, orbit) =
            (named_orbit.name, named_orbit.orbit);

        println!("Simulating orbit {name} for {SIMULATION_TICKS} ticks");

        for iter in 0..SIMULATION_TICKS {
            let time = iter as f64 / SIMULATION_TICKS as f64;
            let time = named_orbit.time_mult * time + named_orbit.time_offset;

            let mean_anom = orbit.get_mean_anomaly_at_time(time);
            let ecc_anom = orbit.get_eccentric_anomaly_at_time(time);
            let angle = orbit.get_true_anomaly_at_eccentric_anomaly(ecc_anom);

            let altitude = orbit.get_altitude_at_angle(angle);
            let (flat_x, flat_y) = orbit.get_flat_position_at_angle(angle);
            let (x, y, z) = orbit.get_position_at_angle(angle);

            logs.push(OrbitLog {
                name,
                iter,
                time,
                mean_anom,
                ecc_anom,
                angle,
                x,
                y,
                z,
                flat_x,
                flat_y,
                altitude,
            });
        }
    }

    println!(">> Creating CSV file.");

    let csv = create_csv(&logs);

    println!("CSV file generated, length: {} KiB", csv.len() / 1024);
    println!(">> Writing CSV file to '{CSV_PATH}'");

    let _ = fs::write(CSV_PATH, csv)
        .inspect_err(|err| eprintln!("Failed to write CSV file: {err}"));

    println!(">> All done!");

    return Ok(());
}

fn create_csv(logs: &Vec<OrbitLog>) -> String {
    let mut string = String::new();
    string += "orbit,iter,time,mean-anom,ecc-anom,angle,flat-x,flat-y,dfx,dfy,x,y,z,dx,dy,dz,altitude,d-altitude,speed,accel\n";

    let mut prev_flat_pos: Option<(f64, f64)> = None;
    let mut prev_pos: Option<(f64, f64, f64)> = None;
    let mut prev_altitude: Option<f64> = None;
    let mut prev_speed: Option<f64> = None;

    for log in logs {
        let (flat_x, flat_y) = (log.flat_x, log.flat_y);
        let (x, y, z) = (log.x, log.y, log.z);
        let altitude = log.altitude;

        let (prev_flat_x, prev_flat_y) = prev_flat_pos.unwrap_or((0.0, 0.0));
        let (prev_x, prev_y, prev_z) = prev_pos.unwrap_or((0.0, 0.0, 0.0));
        let prev_altitude_unwrapped = prev_altitude.unwrap_or(0.0);

        let dx = x - prev_x;
        let dy = y - prev_y;
        let dz = z - prev_z;

        let dfx = flat_x - prev_flat_x;
        let dfy = flat_y - prev_flat_y;

        let d_altitude = altitude - prev_altitude_unwrapped;

        let speed = (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();
        let accel = speed - prev_speed.unwrap_or(0.0);

        string += &format!(
            "{orbit},{iter},{time},{mean_anom},{ecc_anom},{angle},{flat_x},{flat_y},{dfx},{dfy},{x},{y},{z},{dx},{dy},{dz},{altitude},{d_altitude},{speed},{accel}\n",
            orbit = log.name,
            iter = log.iter,
            time = log.time,
            mean_anom = log.mean_anom,
            ecc_anom = log.ecc_anom,
            angle = log.angle,
        );

        prev_flat_pos = Some((flat_x, flat_y));
        prev_pos = Some((x, y, z));
        prev_altitude = Some(altitude);
        prev_speed = Some(speed);
    }

    string
}