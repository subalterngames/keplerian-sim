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
/// 
/// Each element is named `eXY`, where `X` is the row and `Y` is the column.
/// 
/// # Example
/// ```
/// use keplerian_rust::Matrix3x2;
/// 
/// let matrix: Matrix3x2<f64> = Matrix3x2 {
///    e11: 1.0, e12: 0.0,
///    e21: 0.0, e22: 1.0,
///    e31: 0.0, e32: 0.0,
/// };
/// 
/// let vec = (1.0, 2.0);
/// 
/// let result = matrix.dot_vec(vec);
/// 
/// assert_eq!(result, (1.0, 2.0, 0.0));
/// ```
#[allow(missing_docs)]
#[derive(Clone, Debug, PartialEq)]
pub struct Matrix3x2<T> {
    // Element XY
    pub e11: T, pub e12: T,
    pub e21: T, pub e22: T,
    pub e31: T, pub e32: T
}

impl<T: Copy> Copy for Matrix3x2<T> {}
impl<T: Eq> Eq for Matrix3x2<T> {}

impl<T: Copy> Matrix3x2<T> {
    /// Create a new Matrix3x2 instance where each
    /// element is initialized with the same value.
    /// 
    /// # Example
    /// ```
    /// use keplerian_rust::Matrix3x2;
    /// 
    /// let matrix = Matrix3x2::filled_with(0.0);
    /// 
    /// assert_eq!(matrix, Matrix3x2 {
    ///    e11: 0.0, e12: 0.0,
    ///    e21: 0.0, e22: 0.0,
    ///    e31: 0.0, e32: 0.0,
    /// });
    /// ```
    pub fn filled_with(element: T) -> Matrix3x2<T> {
        return Matrix3x2 {
            e11: element, e12: element,
            e21: element, e22: element,
            e31: element, e32: element,
        };
    }
}

impl<T> Matrix3x2<T>
where
    T: Copy + std::ops::Mul<Output = T> + std::ops::Add<Output = T>
{
    /// Computes a dot product between this matrix and a 2D vector.
    /// 
    /// # Example
    /// ```
    /// use keplerian_rust::Matrix3x2;
    /// 
    /// let matrix: Matrix3x2<f64> = Matrix3x2 {
    ///     e11: 1.0, e12: 0.0,
    ///     e21: 0.0, e22: 1.0,
    ///     e31: 1.0, e32: 1.0,
    /// };
    /// 
    /// let vec = (1.0, 2.0);
    /// 
    /// let result = matrix.dot_vec(vec);
    /// 
    /// assert_eq!(result, (1.0, 2.0, 3.0));
    /// ```
    pub fn dot_vec(&self, vec: (T, T)) -> (T, T, T) {
        return (
            vec.0 * self.e11 + vec.1 * self.e12,
            vec.0 * self.e21 + vec.1 * self.e22,
            vec.0 * self.e31 + vec.1 * self.e32
        );
    }
}

type Vec3 = (f64, f64, f64);
type Vec2 = (f64, f64);

/// A trait that defines the methods that a Keplerian orbit must implement.
/// 
/// This trait is implemented by both [`Orbit`] and [`CompactOrbit`].
/// 
/// # Examples
/// ```
/// use keplerian_rust::{Orbit, OrbitTrait, CompactOrbit};
/// 
/// fn accepts_orbit(orbit: &impl OrbitTrait) {
///     println!("That's an orbit!");
/// }
/// 
/// fn main() {
///     let orbit = Orbit::new_default();
///     accepts_orbit(&orbit);
/// 
///     let compact = CompactOrbit::new_default();
///     accepts_orbit(&compact);
/// }
/// ```
/// 
/// ```compile_fail
/// # use keplerian_rust::{Orbit, OrbitTrait, CompactOrbit};
/// # 
/// # fn accepts_orbit(orbit: &impl OrbitTrait) {
/// #     println!("That's an orbit!");
/// # }
/// # 
/// # fn main() {
/// #     let orbit = Orbit::new_default();
/// #     accepts_orbit(&orbit);
/// #  
/// #     let compact = CompactOrbit::new_default();
/// #     accepts_orbit(&compact);
/// let not_orbit = (0.0, 1.0);
/// accepts_orbit(&not_orbit);
/// # }
/// ```
pub trait OrbitTrait {
    /// Gets the semi-major axis of the orbit.
    /// 
    /// In an elliptic orbit, the semi-major axis is the
    /// average of the apoapsis and periapsis.  
    /// This function uses a generalization which uses
    /// eccentricity instead.
    /// 
    /// Learn more: <https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes>
    /// 
    /// # Example
    /// ```
    /// use keplerian_rust::{Orbit, OrbitTrait};
    /// 
    /// let mut orbit = Orbit::new_default();
    /// orbit.set_periapsis(50.0);
    /// orbit.set_apoapsis_force(100.0);
    /// let sma = orbit.get_semi_major_axis();
    /// let expected = 75.0;
    /// assert!((sma - expected).abs() < 1e-6);
    /// ```
    fn get_semi_major_axis(&self) -> f64;

    /// Gets the semi-minor axis of the orbit.
    /// 
    /// In an elliptic orbit, the semi-minor axis is half of the "width" of the orbit.
    /// 
    /// Learn more: <https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes>
    fn get_semi_minor_axis(&self) -> f64;

    /// Gets the semi-latus rectum of the orbit.
    /// 
    /// Learn more: <https://en.wikipedia.org/wiki/Ellipse#Semi-latus_rectum>
    fn get_semi_latus_rectum(&self) -> f64;

    /// Gets the linear eccentricity of the orbit, in meters.
    /// 
    /// In an elliptic orbit, the linear eccentricity is the distance
    /// between its center and either of its two foci (focuses).
    /// 
    /// # Example
    /// ```
    /// use keplerian_rust::{Orbit, OrbitTrait};
    /// 
    /// let mut orbit = Orbit::new_default();
    /// orbit.set_periapsis(50.0);
    /// orbit.set_apoapsis_force(100.0);
    /// 
    /// // Let's say the periapsis is at x = -50.
    /// // The apoapsis would be at x = 100.
    /// // The midpoint would be at x = 25.
    /// // The parent body - one of its foci - is always at the origin (x = 0).
    /// // This means the linear eccentricity is 25.
    /// 
    /// let linear_eccentricity = orbit.get_linear_eccentricity();
    /// let expected = 25.0;
    /// 
    /// assert!((linear_eccentricity - expected).abs() < 1e-6);
    /// ```
    fn get_linear_eccentricity(&self) -> f64;

    /// Gets the apoapsis of the orbit.  
    /// Returns infinity for parabolic orbits.  
    /// Returns negative values for hyperbolic orbits.  
    /// 
    /// # Examples
    /// ```
    /// use keplerian_rust::{Orbit, OrbitTrait};
    /// 
    /// let mut orbit = Orbit::new_default();
    /// orbit.set_eccentricity(0.5); // Elliptic
    /// assert!(orbit.get_apoapsis() > 0.0);
    /// 
    /// orbit.set_eccentricity(1.0); // Parabolic
    /// assert!(orbit.get_apoapsis().is_infinite());
    /// 
    /// orbit.set_eccentricity(2.0); // Hyperbolic
    /// assert!(orbit.get_apoapsis() < 0.0);
    /// ```
    fn get_apoapsis(&self) -> f64;

    /// Sets the apoapsis of the orbit.  
    /// Errors when the apoapsis is less than the periapsis, or less than zero.  
    /// If you want a setter that does not error, use `set_apoapsis_force`, which will
    /// try its best to interpret what you might have meant, but may have
    /// undesirable behavior.
    /// 
    /// # Examples
    /// ```
    /// use keplerian_rust::{Orbit, OrbitTrait};
    /// 
    /// let mut orbit = Orbit::new_default();
    /// orbit.set_periapsis(50.0);
    /// 
    /// assert!(
    ///     orbit.set_apoapsis(100.0)
    ///         .is_ok()
    /// );
    /// 
    /// let result = orbit.set_apoapsis(25.0);
    /// assert!(result.is_err());
    /// assert!(
    ///     result.unwrap_err() ==
    ///     keplerian_rust::ApoapsisSetterError::ApoapsisLessThanPeriapsis
    /// );
    /// 
    /// let result = orbit.set_apoapsis(-25.0);
    /// assert!(result.is_err());
    /// assert!(
    ///     result.unwrap_err() ==
    ///     keplerian_rust::ApoapsisSetterError::ApoapsisNegative
    /// );
    /// ```
    fn set_apoapsis(&mut self, apoapsis: f64) -> Result<(), ApoapsisSetterError>;

    /// Sets the apoapsis of the orbit, with a best-effort attempt at interpreting
    /// possibly-invalid values.  
    /// This function will not error, but may have undesirable behavior:
    /// - If the given apoapsis is less than the periapsis but more than zero,
    ///   the orbit will be flipped and the periapsis will be set to the given apoapsis.
    /// - If the given apoapsis is less than zero, the orbit will be hyperbolic
    ///   instead.
    /// 
    /// If these behaviors are undesirable, consider creating a custom wrapper around
    /// `set_eccentricity` instead.
    fn set_apoapsis_force(&mut self, apoapsis: f64);

    /// Gets the transformation matrix needed to tilt a 2D vector into the
    /// tilted orbital plane.
    /// 
    /// # Example
    /// ```
    /// use keplerian_rust::{Orbit, OrbitTrait};
    /// 
    /// let orbit = Orbit::new_default();
    /// let matrix = orbit.get_transformation_matrix();
    /// 
    /// assert_eq!(matrix, keplerian_rust::Matrix3x2 {
    ///     e11: 1.0, e12: 0.0,
    ///     e21: 0.0, e22: 1.0,
    ///     e31: 0.0, e32: 0.0,
    /// });
    /// ```
    fn get_transformation_matrix(&self) -> Matrix3x2<f64>;

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

    /// Gets the true anomaly at a given eccentric anomaly in the orbit.
    /// 
    /// This function is faster than the function which takes mean anomaly as input,
    /// as the eccentric anomaly is hard to calculate.
    fn get_true_anomaly_at_eccentric_anomaly(&self, eccentric_anomaly: f64) -> f64;

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
    fn get_true_anomaly(&self, mean_anomaly: f64) -> f64 {
        self.get_true_anomaly_at_eccentric_anomaly(
            self.get_eccentric_anomaly(mean_anomaly)
        )
    }

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
    fn get_eccentric_anomaly_at_time(&self, t: f64) -> f64 {
        self.get_eccentric_anomaly(
            self.get_mean_anomaly_at_time(t)
        )
    }

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
    fn get_true_anomaly_at_time(&self, t: f64) -> f64 {
        self.get_true_anomaly(
            self.get_mean_anomaly_at_time(t)
        )
    }

    /// Gets the 3D position at a given angle (true anomaly) in the orbit.
    /// 
    /// The angle is expressed in radians, and ranges from 0 to tau.  
    /// Anything out of range will get wrapped around.
    /// 
    /// # Example
    /// ```
    /// use keplerian_rust::{Orbit, OrbitTrait};
    /// 
    /// let mut orbit = Orbit::new_default();
    /// orbit.set_periapsis(100.0);
    /// orbit.set_eccentricity(0.0);
    /// 
    /// let pos = orbit.get_position_at_angle(0.0);
    /// 
    /// assert_eq!(pos, (100.0, 0.0, 0.0));
    /// ```
    fn get_position_at_angle(&self, angle: f64) -> Vec3 {
        let (x, y) = self.get_flat_position_at_angle(angle);
        self.tilt_flat_position(x, y)
    }

    /// Gets the 2D position at a given angle (true anomaly) in the orbit.
    /// 
    /// This ignores "orbital tilting" parameters, namely the inclination and
    /// the longitude of ascending node.
    /// 
    /// The angle is expressed in radians, and ranges from 0 to tau.  
    /// Anything out of range will get wrapped around.
    /// 
    /// # Example
    /// ```
    /// use keplerian_rust::{Orbit, OrbitTrait};
    /// 
    /// let mut orbit = Orbit::new_default();
    /// orbit.set_periapsis(100.0);
    /// orbit.set_eccentricity(0.0);
    /// 
    /// let pos = orbit.get_flat_position_at_angle(0.0);
    /// 
    /// assert_eq!(pos, (100.0, 0.0));
    /// ```
    fn get_flat_position_at_angle(&self, angle: f64) -> Vec2 {
        let alt = self.get_altitude_at_angle(angle);
        return (
            alt * angle.cos(),
            alt * angle.sin()
        );
    }

    /// Gets the altitude of the body from its parent at a given angle (true anomaly) in the orbit.
    /// 
    /// # Example
    /// ```
    /// use keplerian_rust::{Orbit, OrbitTrait};
    /// 
    /// let mut orbit = Orbit::new_default();
    /// orbit.set_periapsis(100.0);
    /// orbit.set_eccentricity(0.0);
    /// 
    /// let altitude = orbit.get_altitude_at_angle(0.0);
    /// 
    /// assert_eq!(altitude, 100.0);
    /// ```
    fn get_altitude_at_angle(&self, angle: f64) -> f64;

    /// Gets the altitude of the body from its parent at a given time in the orbit.
    /// 
    /// This involves calculating the true anomaly at a given time, and so is not very performant.  
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
    fn get_altitude_at_time(&self, t: f64) -> f64 {
        self.get_altitude_at_angle(
            self.get_true_anomaly_at_time(t)
        )
    }

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
    fn get_position_at_time(&self, t: f64) -> Vec3 {
        self.get_position_at_angle(
            self.get_true_anomaly_at_time(t)
        )
    }

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
    fn get_flat_position_at_time(&self, t: f64) -> Vec2 {
        self.get_flat_position_at_angle(
            self.get_true_anomaly_at_time(t)
        )
    }

    /// Tilts a 2D position into 3D, using the orbital parameters.
    /// 
    /// This uses the "orbital tilting" parameters, namely the inclination
    /// and longitude of ascending node, to tilt that position into the same
    /// plane that the orbit resides in.
    /// 
    /// This function performs 10x faster in the cached version of the
    /// [`Orbit`] struct, as it doesn't need to recalculate the transformation
    /// matrix needed to transform 2D vector.
    fn tilt_flat_position(&self, x: f64, y: f64) -> Vec3 {
        self.get_transformation_matrix().dot_vec((x, y))
    }
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

#[inline]
fn keplers_equation(mean_anomaly: f64, eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentric_anomaly - (eccentricity * eccentric_anomaly.sin()) - mean_anomaly;
}
#[inline]
fn keplers_equation_derivative(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return 1.0 - (eccentricity * eccentric_anomaly.cos());
}
#[inline]
fn keplers_equation_second_derivative(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentricity * eccentric_anomaly.sin();
}

#[inline]
fn keplers_equation_hyperbolic(mean_anomaly: f64, eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentricity * eccentric_anomaly.sinh() - eccentric_anomaly - mean_anomaly;
}
#[inline]
fn keplers_equation_hyperbolic_derivative(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentricity * eccentric_anomaly.cosh() - 1.0;
}
#[inline]
fn keplers_equation_hyperbolic_second_derivative(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentricity * eccentric_anomaly.sinh();
}

/// Solve a cubic equation to get its real root.
/// 
/// The cubic equation is in the form of:
/// ax^3 + bx^2 + cx + d
/// 
/// The cubic equation is assumed to be monotone.  
/// If it isn't monotone (i.e., the discriminant
/// is negative), it may return an incorrect value
/// or NaN.
fn solve_monotone_cubic(a: f64, b: f64, c: f64, d: f64) -> f64 {
    // Normalize coefficients so that a = 1
    // ax^3 + bx^2 + cx + d
    // ...where b, c, d are the normalized coefficients,
    // and a = 1
    let b = b / a;
    let c = c / a;
    let d = d / a;

    // Depress the cubic equation
    // t^3 + pt + q = 0
    // ...where:
    // p = (3ac - b^2) / (3a^2)
    // q = (2b^3 - 9abc + 27da^2) / (27a^3)
    // ...since a = 1, we can simplify them to:
    // p = (3c - b^2) / 3
    // q = (2b^3 - 9bc + 27d) / 27
    let b_sq = b * b;

    let p = (3.0 * c - b_sq) / 3.0;
    let q = (2.0 * b_sq * b - 9.0 * b * c + 27.0 * d) / 27.0;

    let q_div_two = q / 2.0;
    let p_div_three = p / 3.0;
    let p_div_three_cubed = p_div_three * p_div_three * p_div_three;
    let discriminant =
        q_div_two * q_div_two +
        p_div_three_cubed;

    if discriminant < 0.0 {
        // Function is not monotone
        return f64::NAN;
    }
    
    let t = {
        let sqrt_discriminant = discriminant.sqrt();
        let neg_q_div_two = -q_div_two;
        let u = (neg_q_div_two + sqrt_discriminant).cbrt();
        let v = (neg_q_div_two - sqrt_discriminant).cbrt();
        u + v
    };

    // x_i = t_i - b / 3a
    // here, a = 1
    return t - b / 3.0;
}

mod generated_sinh_approximator;

#[cfg(test)]
mod sinh_approx_bench {
    use crate::generated_sinh_approximator::sinh_approx_lt5;
    use std::hint::black_box;
    use std::time::Instant;

    // Run this benchmark:
    // `RUST_TEST_NOCAPTURE=1 cargo test -r bench_sinh` (bash)
    #[test]
    fn bench_sinh() {
        let base_start_time = Instant::now();
        for i in 0..5000000 {
            let f = i as f64 * 1e-6;
            black_box(black_box(f).sinh());
        }
        let base_duration = base_start_time.elapsed();

        let approx_start_time = Instant::now();
        for i in 0..5000000 {
            let f = i as f64 * 1e-6;
            black_box(sinh_approx_lt5(black_box(f)));
        }
        let approx_duration = approx_start_time.elapsed();

        eprintln!("Real sinh call: {:?}", base_duration);
        eprintln!("Approximation: {:?}", approx_duration);
        let speed_factor = base_duration.as_nanos() as f64 / approx_duration.as_nanos() as f64;
        eprintln!(
            "\nThe approximation is {speed_factor}x the speed of the real sinh call.",
        );
        assert!(speed_factor > 1.0, "The approximation is slower than the real sinh call.");
    }
}