//! # Keplerian Orbital Mechanics
//! This library crate contains logic for Keplerian orbits, similar to the ones
//! you'd find in a game like Kerbal Space Program.  
//! 
//! Keplerian orbits are special in that they are more stable and predictable than
//! Newtonian orbits. In fact, unlike Newtonian orbits, Keplerian orbits don't use
//! time steps to calculate the next position of an object. Keplerian orbits use
//! state vectors to determine the object's *full trajectory* at any given time.  
//! This way, you don't need to worry about lag destabilizing Keplerian orbits.  
//! 
//! However, Keplerian orbits are significantly more complex to calculate than
//! just using Newtonian physics. It's also a two-body simulation, meaning that
//! it doesn't account for external forces like gravity from other bodies or the
//! engines of a spacecraft.
//! 
//! The way Kerbal Space Program handles this is to have an "on-rails" physics
//! system utilizing Keplerian orbits, and an "active" physics system utilizing
//! Newtonian two-body physics.
//! 
//! ## Getting started
//! This crate provides four main structs:
//! - [`Orbit`]: A struct representing an orbit around a celestial body.
//!   Each instance of this struct has some cached data to speed up
//!   certain calculations, and has a larger memory footprint.
//! - [`CompactOrbit`]: A struct representing an orbit around a celestial body.
//!   This struct has a smaller memory footprint than the regular `Orbit` struct,
//!   but some calculations may take 2~10x slower because it doesn't save any
//!   cached calculations.
//! - [`Body`]: A struct representing a celestial body. This struct contains
//!   information about the body's mass, radius, and orbit.
//! - [`Universe`]: A struct representing the entire simulation. This struct
//!   contains a list of all the bodies in the simulation, and can calculate
//!   the absolute position of any body at any given time.
//!   To do this, it stores parent-child relationships between bodies.
//! 
//! We also provide a [`body_presets`] module, which contains some preset celestial
//! bodies to use in your simulation. It contains many celestial bodies, like
//! the Sun, the Moon, and all the planets in the Solar System.
//! 
//! ## Example
//! 
//! ```rust
//! use keplerian_rust::{Orbit, OrbitTrait};
//! 
//! # fn main() {
//! // Create a perfectly circular orbit with a radius of 1 meter
//! let orbit = Orbit::new_default();
//! assert_eq!(orbit.get_position_at_time(0.0), (1.0, 0.0, 0.0));
//! # }
//! #
//! ```

#![warn(missing_docs)]

mod cached_orbit;
mod compact_orbit;
mod body;
mod universe;
pub mod body_presets;

pub use cached_orbit::Orbit;
pub use compact_orbit::CompactOrbit;
pub use body::Body;
pub use universe::Universe;

/// A struct representing a 3x2 matrix.
/// 
/// This struct is used to store the transformation matrix
/// for transforming a 2D vector into a 3D vector.
/// 
/// Namely, it is used in the [`tilt_flat_position`][OrbitTrait::tilt_flat_position]
/// method to tilt a 2D position into 3D, using the orbital parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct Matrix3x2<T> {
    // Element XY
    e11: T, e12: T,
    e21: T, e22: T,
    e31: T, e32: T
}

impl Matrix3x2<f64> {
    /// Create a new Matrix3x2 instance where each
    /// element is initialized with the same value.
    fn filled_with<T: Copy>(element: T) -> Matrix3x2<T> {
        return Matrix3x2 {
            e11: element, e12: element,
            e21: element, e22: element,
            e31: element, e32: element,
        };
    }
}

type Vec3 = (f64, f64, f64);
type Vec2 = (f64, f64);

/// A trait that defines the methods that a Keplerian orbit must implement.
/// 
/// This trait is implemented by both [`Orbit`] and [`CompactOrbit`].
pub trait OrbitTrait {
    /// Gets the eccentric anomaly at a given mean anomaly in the orbit.
    /// 
    /// The method to get the eccentric anomaly often uses numerical
    /// methods like Newton's method, and so it is not very performant.  
    /// It is recommended to cache this value if you can.
    /// 
    /// When the orbit is open (has an eccentricity of at least 1),
    /// the [hyperbolic eccentric anomaly](https://space.stackexchange.com/questions/27602/what-is-hyperbolic-eccentric-anomaly-f)
    /// would be returned instead.
    /// 
    /// The eccentric anomaly is an angular parameter that defines the position
    /// of a body that is moving along an elliptic Kepler orbit.
    /// 
    /// \- [Wikipedia](https://en.wikipedia.org/wiki/Eccentric_anomaly)
    fn get_eccentric_anomaly(&self, mean_anomaly: f64) -> f64;

    /// Gets the true anomaly at a given mean anomaly in the orbit.
    /// 
    /// The true anomaly is derived from the eccentric anomaly, which
    /// uses numerical methods and so is not very performant.  
    /// It is recommended to cache this value if you can.
    /// 
    /// The true anomaly is the angle between the direction of periapsis
    /// and the current position of the body, as seen from the main focus
    /// of the ellipse.
    /// 
    /// \- [Wikipedia](https://en.wikipedia.org/wiki/True_anomaly)
    fn get_true_anomaly(&self, mean_anomaly: f64) -> f64;

    /// Gets the mean anomaly at a given time in the orbit.
    /// 
    /// The mean anomaly is the fraction of an elliptical orbit's period
    /// that has elapsed since the orbiting body passed periapsis,
    /// expressed as an angle which can be used in calculating the position
    /// of that body in the classical two-body problem.
    /// 
    /// \- [Wikipedia](https://en.wikipedia.org/wiki/Mean_anomaly)
    fn get_mean_anomaly_at_time(&self, t: f64) -> f64;

    /// Gets the eccentric anomaly at a given time in the orbit.
    /// 
    /// The method to get the eccentric anomaly often uses numerical
    /// methods like Newton's method, and so it is not very performant.  
    /// It is recommended to cache this value if you can.
    /// 
    /// When the orbit is open (has an eccentricity of at least 1),
    /// the [hyperbolic eccentric anomaly](https://space.stackexchange.com/questions/27602/what-is-hyperbolic-eccentric-anomaly-f)
    /// would be returned instead.
    /// 
    /// The eccentric anomaly is an angular parameter that defines the position
    /// of a body that is moving along an elliptic Kepler orbit.
    /// 
    /// \- [Wikipedia](https://en.wikipedia.org/wiki/Eccentric_anomaly)
    fn get_eccentric_anomaly_at_time(&self, t: f64) -> f64;

    /// Gets the true anomaly at a given time in the orbit.
    /// 
    /// The true anomaly is derived from the eccentric anomaly, which
    /// uses numerical methods and so is not very performant.  
    /// It is recommended to cache this value if you can.
    /// 
    /// The true anomaly is the angle between the direction of periapsis
    /// and the current position of the body, as seen from the main focus
    /// of the ellipse.
    /// 
    /// \- [Wikipedia](https://en.wikipedia.org/wiki/True_anomaly)
    fn get_true_anomaly_at_time(&self, t: f64) -> f64;

    /// Gets the 3D position at a given angle (true anomaly) in the orbit.
    /// 
    /// The angle is expressed in radians, and ranges from 0 to tau.  
    /// Anything out of range will get wrapped around.
    fn get_position_at_angle(&self, angle: f64) -> Vec3;

    /// Gets the 2D position at a given angle (true anomaly) in the orbit.
    /// 
    /// This ignores "orbital tilting" parameters, namely the inclination and
    /// the longitude of ascending node.
    /// 
    /// The angle is expressed in radians, and ranges from 0 to tau.  
    /// Anything out of range will get wrapped around.
    fn get_flat_position_at_angle(&self, angle: f64) -> Vec2;

    /// Gets the 3D position at a given time in the orbit.
    /// 
    /// This involves calculating the true anomaly at a given time,
    /// and so is not very performant.  
    /// It is recommended to cache this value when possible.
    /// 
    /// For closed orbits (with an eccentricity less than 1), the
    /// `t` (time) value ranges from 0 to 1.  
    /// Anything out of range will get wrapped around.
    /// 
    /// For open orbits (with an eccentricity of at least 1), the
    /// `t` (time) value is unbounded.  
    /// Note that due to floating-point imprecision, values of extreme
    /// magnitude may not be accurate.
    fn get_position_at_time(&self, t: f64) -> Vec3;

    /// Gets the 2D position at a given time in the orbit.
    /// 
    /// This involves calculating the true anomaly at a given time,
    /// and so is not very performant.
    /// It is recommended to cache this value when possible.
    /// 
    /// This ignores "orbital tilting" parameters, namely the inclination
    /// and longitude of ascending node.
    /// 
    /// For closed orbits (with an eccentricity less than 1), the
    /// `t` (time) value ranges from 0 to 1.  
    /// Anything out of range will get wrapped around.
    /// 
    /// For open orbits (with an eccentricity of at least 1), the
    /// `t` (time) value is unbounded.  
    /// Note that due to floating-point imprecision, values of extreme
    /// magnitude may not be accurate.
    fn get_flat_position_at_time(&self, t: f64) -> Vec2;

    /// Tilts a 2D position into 3D, using the orbital parameters.
    /// 
    /// This uses the "orbital tilting" parameters, namely the inclination
    /// and longitude of ascending node, to tilt that position into the same
    /// plane that the orbit resides in.
    /// 
    /// This function performs 10x faster in the cached version of the
    /// [`Orbit`] struct, as it doesn't need to recalculate the transformation
    /// matrix needed to transform 2D vector.
    fn tilt_flat_position(&self, x: f64, y: f64) -> Vec3;
}

/// An error to describe why setting the periapsis of an orbit failed.
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum ApoapsisSetterError {
    /// ### Attempt to set apoapsis to a value less than periapsis.
    /// By definition, an orbit's apoapsis is the highest point in the orbit, 
    /// and its periapsis is the lowest point in the orbit.  
    /// Therefore, it doesn't make sense for the apoapsis to be lower than the periapsis.
    ApoapsisLessThanPeriapsis,

    /// ### Attempt to set apoapsis to a negative value.
    /// By definition, the apoapsis is the highest point in the orbit.  
    /// You can't be a negative distance away from the center of mass of the parent body.  
    /// Therefore, it doesn't make sense for the apoapsis to be lower than zero.
    ApoapsisNegative
}

#[cfg(test)]
mod tests;


fn keplers_equation(mean_anomaly: f64, eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentric_anomaly - (eccentricity * eccentric_anomaly.sin()) - mean_anomaly;
}
fn keplers_equation_derivative(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return 1.0 - (eccentricity * eccentric_anomaly.cos());
}

fn keplers_equation_hyperbolic(mean_anomaly: f64, eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentricity * eccentric_anomaly.sinh() - eccentric_anomaly - mean_anomaly;
}

fn keplers_equation_hyperbolic_derivative(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentricity * eccentric_anomaly.cosh() - 1.0;
}

fn keplers_equation_hyperbolic_second_derivative(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentricity * eccentric_anomaly.sinh();
}