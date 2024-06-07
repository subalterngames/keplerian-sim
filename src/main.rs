mod celestial;
use celestial::{
    Body, Orbit, Universe, body_presets
};

fn main() -> Result<(), String> {
    let mut universe = generate_solar_system();
    dbg!(universe.clone());
    universe.warp(10_000_000)?;
    dbg!(universe.clone());
    return Ok(());
}

fn generate_solar_system<'a>() -> Universe {
    let mut universe = Universe::new_default();
    
    let sun = body_presets::stars::the_sun(false);
    let sun_idx = universe.add_body(sun, None);
    
    let earth = body_presets::planets::earth(true);
    universe.add_body(earth, Some(sun_idx));
    return universe;
}