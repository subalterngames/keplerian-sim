use keplerian_sim::{OrbitTrait, body_presets};
use std::{fs, path::PathBuf};

const SIMULATION_TICKS: u128 = 10_000;
const CSV_PATH: &str = "out/output-dump-earthmoon.csv";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let moon_orbit = body_presets::moons::the_moon(true)
        .orbit.unwrap();
    let mut moon_positions: Vec<(f64, f64, f64)> = Vec::with_capacity(SIMULATION_TICKS as usize);

    print!("Simulating Earth-Moon system for {SIMULATION_TICKS} ticks...");

    for i in 0..SIMULATION_TICKS {
        let time = i as f64 / SIMULATION_TICKS as f64;
        let pos = moon_orbit.get_position_at_time(time);
        moon_positions.push(pos);
    }

    println!("Done!");

    let path = PathBuf::from(CSV_PATH);
    print!("Saving CSV to {:#?}...", &path);

    fs::create_dir_all(
        &path.parent()
        .expect("Failed to get parent of CSV path")
    )?;

    let contents = create_csv(&moon_positions);
    fs::write(&path, contents)?;

    println!("Done!");

    return Ok(());
}

fn create_csv(positions: &Vec<(f64, f64, f64)>) -> String {
    let mut string = String::new();
    string += "time,x,y,z,dx,dy,dz,altitude,d-altitude,speed,d-speed\n";

    let mut prev_pos: Option<(f64, f64, f64)> = None;
    let mut prev_altitude: Option<f64> = None;
    let mut prev_speed: Option<f64> = None;

    for (time, &(x, y, z)) in positions.iter().enumerate() {
        let (dx, dy, dz) = match prev_pos {
            Some((px, py, pz)) => (x - px, y - py, z - pz),
            None => (0.0, 0.0, 0.0),
        };

        let altitude = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
        let speed = (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt();

        let d_altitude = match prev_altitude {
            Some(prev_altitude) => altitude - prev_altitude,
            None => 0.0,
        };

        let d_speed = match prev_speed {
            Some(prev_speed) => speed - prev_speed,
            None => 0.0,
        };

        string += &format!(
            "{time},{x},{y},{z},{dx},{dy},{dz},{altitude},{d_altitude},{speed},{d_speed}\n"
        );

        prev_pos = Some((x, y, z));
        prev_altitude = Some(altitude);
        prev_speed = Some(speed);
    }

    string
}