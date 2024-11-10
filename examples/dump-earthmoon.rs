use keplerian_rust::{
    Universe, body_presets
};

const SIMULATION_TICKS: u128 = 1_000_000;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut universe = create_earth_moon();
    let mut moon_positions: Vec<(f64, f64, f64)> = Vec::with_capacity(SIMULATION_TICKS as usize);

    print!("Simulating Earth-Moon system for {SIMULATION_TICKS} ticks...");

    for _ in 0..SIMULATION_TICKS {
        universe.tick();
        let moon = universe.get_body(1);
        moon_positions.push(moon.get_relative_position());
    }

    println!("Done!");
    println!("Saving CSV...");

    // TODO: save CSV

    return Ok(());
}

fn create_earth_moon() -> Universe {
    let mut universe = Universe::new_default();

    universe.time_step = 10.0;

    let earth = body_presets::planets::earth(false);
    let earth_idx = universe.add_body(earth, None);

    let moon = body_presets::moons::the_moon(true);
    universe.add_body(moon, Some(earth_idx));

    return universe;
}