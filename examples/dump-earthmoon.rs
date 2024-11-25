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

    string += "time,x,y,z\n";

    let iterator = positions.iter().enumerate();
    
    for (time, (x, y, z)) in iterator {
        string += &format!("{time},{x},{y},{z}\n");
    }

    return string;
}