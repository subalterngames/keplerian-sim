use crate::{body_presets::planets::earth, Orbit, OrbitTrait};
use core::f64::consts::TAU;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A struct representing a celestial body.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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
    /// Creates a new `Body` instance.
    ///
    /// # Arguments
    ///
    /// * `name` - The name of the celestial body.
    /// * `mass` - The mass of the celestial body, in kilograms.
    /// * `radius` - The radius of the celestial body, in meters.
    /// * `orbit` - An optional orbit for the celestial body.
    ///
    /// # Returns
    ///
    /// A new `Body` instance.
    pub fn new(name: String, mass: f64, radius: f64, orbit: Option<Orbit>) -> Self {
        Self {
            name,
            mass,
            radius,
            orbit,
            progress: 0.0,
        }
    }

    /// Releases the body from its orbit.
    pub fn release_from_orbit(&mut self) {
        self.orbit = None;
        self.progress = 0.0;
    }

    /// Get the amount of time it takes for the body to complete one orbit,
    /// given a gravitational constant.
    pub fn get_orbital_period(&self, g: f64) -> Option<f64> {
        let orbit = self.orbit.as_ref()?;
        let mu = orbit.get_mu(g);

        Some(if orbit.get_eccentricity() >= 1.0 {
            f64::INFINITY
        } else {
            TAU * (orbit.get_semi_major_axis() / mu).sqrt()
        })
    }
    /// Progresses this body's orbit, given a time step and the gravitational
    /// acceleration towards the parent body.
    pub fn progress_orbit(&mut self, dt: f64, g: f64) -> Result<(), String> {
        let orbit = self.orbit.as_ref().ok_or("Body is not in orbit")?;

        if orbit.get_eccentricity() >= 1.0 {
            self.progress += dt * g.sqrt();
        } else {
            let period = self.get_orbital_period(g).unwrap();
            self.progress += dt / period;
            self.progress = self.progress.rem_euclid(1.0);
        }

        Ok(())
    }
}

impl Default for Body {
    /// Creates a default `Body` instance.
    ///
    /// Currently, this function returns the Earth.  
    /// However, do not rely on this behavior, as it may change in the future.
    fn default() -> Self {
        earth(true)
    }
}
