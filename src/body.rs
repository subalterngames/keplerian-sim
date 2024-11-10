use crate::Orbit;
use std::f64::consts::TAU as TAU;

/// A struct representing a celestial body.
#[derive(Clone, Debug, PartialEq)]
pub struct Body {
    /// The name of the celestial body.
    pub name: String,

    /// The mass of the celestial body, in kilograms.
    pub mass: f64,

    /// The radius of the celestial body, in meters.
    pub radius: f64,

    /// The orbit of the celestial body, if it is orbiting one.
    pub orbit: Option<Orbit>,

    /// The orbit progress, between 0 and 1.
    pub progress: f64,
}

impl Body {
    pub fn new(
        name: String, mass: f64, radius: f64,
        orbit: Option<Orbit>
    ) -> Body {
        return Body {
            name, mass, radius, orbit,
            progress: 0.0
        };
    }
    pub fn new_default() -> Body {
        return Body {
            name: "Earth".to_string(),
            mass: 5.972e24,
            radius: 6.371e6,
            orbit: None,
            progress: 0.0,
        };
    }
    pub fn release_from_orbit(&mut self) {
        self.orbit = None;
        self.progress = 0.0;
    }
    pub fn get_orbital_period(&self, g: f64) -> Result<f64, String> {
        let orbit = self.orbit.as_ref().ok_or("Body is not in orbit")?;
        let mu = g * self.mass;

        let semi_major_axis = orbit.get_semi_major_axis();

        return Ok(TAU * (semi_major_axis / mu).sqrt());
    }
    pub fn progress_orbit(&mut self, dt: f64, g: f64) -> Result<(), String> {
        let period = self.get_orbital_period(g)?;
        let delta_progress = dt / period;
        self.progress += delta_progress;
        self.progress = self.progress.rem_euclid(1.0);

        return Ok(());
    }
    pub fn get_relative_position(&self) -> (f64, f64, f64) {
        let orbit = self.orbit.as_ref();

        if orbit.is_none() {
            return (0.0, 0.0, 0.0);
        }

        return orbit.unwrap().get_position_at_time(self.progress);
    }
}