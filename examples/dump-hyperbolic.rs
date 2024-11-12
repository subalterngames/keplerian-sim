use keplerian_rust::{Orbit, OrbitTrait};
use std::{fs, path::PathBuf};

const SIMULATION_TICKS: u128 = 10_000;
const CSV_PATH: &str = "out/output-dump-hyperbolic.csv";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let orbit = Orbit::new(
        8.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
    );
    let mut positions: Vec<(f64, f64, f64)> = Vec::with_capacity(SIMULATION_TICKS as usize);

    print!("Simulating parabolic trajectory for {SIMULATION_TICKS} ticks...");

    for i in 0..SIMULATION_TICKS {
        let time = i as f64 / SIMULATION_TICKS as f64 - 0.5;
        let pos = orbit.get_position_at_time(time);
        positions.push(pos);
    }

    println!("Done!");

    let path = PathBuf::from(CSV_PATH);
    print!("Saving CSV to {:#?}...", &path);

    fs::create_dir_all(
        &path.parent()
        .expect("Failed to get parent of CSV path")
    )?;

    let contents = create_csv(&positions);
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