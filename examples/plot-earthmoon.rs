use keplerian_rust::{
    Universe, body_presets
};
use plotters::prelude::*;
use plotters_backend::piston::PistonBackend;

const SIMULATION_TICKS: u128 = 1_000_000;
const OUTPUT_PATH: &str = "~earth_moon_orbits.png";

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
    println!("Showing plot...");

    // Create a window
    let mut window: PistonWindow = WindowSettings::new("3D Orbit", [800, 600])
    .exit_on_esc(true)
    .build()
    .unwrap();

    // Create a drawing area
    let root = PistonBackend::new(&mut window).into_drawing_area();
    root.fill(&WHITE).unwrap();

    // Create a 3D chart
    let mut chart = ChartBuilder::on(&root)
        .caption("Earth-Moon Orbit", ("sans-serif", 50))
        .build_cartesian_3d(-1.0..1.0, -1.0..1.0, -1.0..1.0)
        .unwrap();

    chart.configure_axes().draw().unwrap();

    // Plot the moon positions
    chart
        .draw_series(LineSeries::new(
            moon_positions.iter().map(|&(x, y, z)| (x, y, z)),
            &RED,
        ))
        .unwrap()
        .label("Moon Orbit")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart.configure_series_labels().draw().unwrap();

    // Render the window
    while let Some(event) = window.next() {
        window.draw_2d(&event, |c, g, _| {
            root.draw(&c, g).unwrap();
        });
    }

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