use keplerian_rust::{
    Universe, body_presets
};

const SIMULATION_TICKS: u128 = 1_000_000;
fn main() {
    let mut universe = generate_solar_system();
    describe_universe(&universe);
    print!("Simulating {SIMULATION_TICKS} ticks...");
    universe.warp(SIMULATION_TICKS);
    println!(" done");
    print_all_body_positions(&universe);
}

fn generate_solar_system<'a>() -> Universe {
    let mut universe = Universe::new_default();
    
    let sun = body_presets::stars::the_sun(false);
    let sun_idx = universe.add_body(sun, None);

    let mercury = body_presets::planets::mercury(true);
    universe.add_body(mercury, Some(sun_idx));

    let venus = body_presets::planets::venus(true);
    universe.add_body(venus, Some(sun_idx));
    
    let earth = body_presets::planets::earth(true);
    let earth_idx = universe.add_body(earth, Some(sun_idx));

    let moon = body_presets::moons::the_moon(true);
    universe.add_body(moon, Some(earth_idx));

    let mars = body_presets::planets::mars(true);
    universe.add_body(mars, Some(sun_idx));

    return universe;
}

fn describe_universe(universe: &Universe) {
    println!("Simulation universe with {} bodies", universe.get_bodies().len());
    for (i, body) in universe.get_bodies().iter().enumerate() {
        println!("    {}: {:?}", i, body.name);
        println!("      Mass: {}", body.mass);
        println!("      Radius: {}", body.radius);
        if let Some(orbit) = &body.orbit {
            let location = universe.get_body_position(i);
            println!("      Orbit: {:?}", location);
            println!("        Semi-major axis: {}", orbit.get_semi_major_axis());
            println!("        Eccentricity: {}", orbit.get_eccentricity());
            println!("        Inclination: {}", orbit.get_inclination());
            println!("        Argument of periapsis: {}", orbit.get_arg_pe());
            println!("        Longitude of ascending node: {}", orbit.get_long_asc_node());
            println!("        Mean anomaly at epoch: {}", orbit.get_mean_anomaly());

        }
    }
}

fn print_all_body_positions(universe: &Universe) {
    for (i, body) in universe.get_bodies().iter().enumerate() {
        let location = universe.get_body_position(i);
        println!("{}: {:?}", body.name, location);
    }
}