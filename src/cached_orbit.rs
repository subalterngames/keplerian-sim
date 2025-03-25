#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::{ApoapsisSetterError, CompactOrbit, Matrix3x2, OrbitTrait};

/// A struct representing a Keplerian orbit with some cached values.
///
/// This struct consumes significantly more memory because of the cache.  
/// However, this will speed up orbital calculations.  
/// If memory efficiency is your goal, you may consider using the `CompactOrbit` struct instead.  
///
/// # Example
/// ```
/// use keplerian_sim::{Orbit, OrbitTrait};
///
/// let orbit = Orbit::new(
///     // Initialize using eccentricity, periapsis, inclination,
///     // argument of periapsis, longitude of ascending node,
///     // and mean anomaly at epoch
///
///     // Eccentricity
///     0.0,
///
///     // Periapsis
///     1.0,
///
///     // Inclination
///     0.0,
///
///     // Argument of periapsis
///     0.0,
///
///     // Longitude of ascending node
///     0.0,
///
///     // Mean anomaly at epoch
///     0.0,
///
///     // Mass of the parent body
///     1.0
/// );
///
/// let orbit = Orbit::with_apoapsis(
///     // Initialize using apoapsis in place of eccentricity
///     
///     // Apoapsis
///     2.0,
///
///     // Periapsis
///     1.0,
///
///     // Inclination
///     0.0,
///
///     // Argument of periapsis
///     0.0,
///
///     // Longitude of ascending node
///     0.0,
///
///     // Mean anomaly at epoch
///     0.0,
///
///     // Mass of the parent body
///     1.0
/// );
/// ```
/// See [Orbit::new] and [Orbit::with_apoapsis] for more information.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Orbit {
    /// The eccentricity of the orbit.  
    /// e < 1: ellipse  
    /// e = 1: parabola  
    /// e > 1: hyperbola  
    eccentricity: f64,

    /// The periapsis of the orbit, in meters.
    periapsis: f64,

    /// The inclination of the orbit, in radians.
    inclination: f64,

    /// The argument of periapsis of the orbit, in radians.
    arg_pe: f64,

    /// The longitude of ascending node of the orbit, in radians.
    long_asc_node: f64,

    /// The mean anomaly at orbit epoch, in radians.
    mean_anomaly: f64,

    /// The mass of the parent body, in kilograms.
    parent_mass: f64,

    cache: OrbitCachedCalculations,
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
struct OrbitCachedCalculations {
    /// The semi-major axis of the orbit, in meters.
    semi_major_axis: f64,

    /// The semi-minor axis of the orbit, in meters.
    semi_minor_axis: f64,

    /// The linear eccentricity of the orbit, in meters.
    linear_eccentricity: f64,

    /// The transformation matrix to tilt the 2D planar orbit into 3D space.
    transformation_matrix: Matrix3x2,

    /// A value based on the orbit's eccentricity, used to calculate
    /// the true anomaly from the eccentric anomaly.  
    /// https://en.wikipedia.org/wiki/True_anomaly#From_the_eccentric_anomaly
    beta: f64,
}
// Initialization and cache management
impl Orbit {
    fn update_cache(&mut self) {
        self.cache = Self::get_cached_calculations(
            self.eccentricity,
            self.periapsis,
            self.inclination,
            self.arg_pe,
            self.long_asc_node,
        );
    }

    fn get_cached_calculations(
        eccentricity: f64,
        periapsis: f64,
        inclination: f64,
        arg_pe: f64,
        long_asc_node: f64,
    ) -> OrbitCachedCalculations {
        let semi_major_axis = periapsis / (1.0 - eccentricity);
        let semi_minor_axis = semi_major_axis * (1.0 - eccentricity * eccentricity).abs().sqrt();
        let linear_eccentricity = semi_major_axis * eccentricity;
        let transformation_matrix =
            Self::get_transformation_matrix(inclination, arg_pe, long_asc_node);
        let beta = eccentricity / (1.0 + (1.0 - eccentricity * eccentricity).sqrt());

        OrbitCachedCalculations {
            semi_major_axis,
            semi_minor_axis,
            linear_eccentricity,
            transformation_matrix,
            beta,
        }
    }

    fn get_transformation_matrix(inclination: f64, arg_pe: f64, long_asc_node: f64) -> Matrix3x2 {
        let mut matrix = Matrix3x2::default();

        let (sin_inc, cos_inc) = inclination.sin_cos();
        let (sin_arg_pe, cos_arg_pe) = arg_pe.sin_cos();
        let (sin_lan, cos_lan) = long_asc_node.sin_cos();

        // https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
        matrix.e11 = cos_arg_pe * cos_lan - sin_arg_pe * cos_inc * sin_lan;
        matrix.e12 = -(sin_arg_pe * cos_lan + cos_arg_pe * cos_inc * sin_lan);

        matrix.e21 = cos_arg_pe * sin_lan + sin_arg_pe * cos_inc * cos_lan;
        matrix.e22 = cos_arg_pe * cos_inc * cos_lan - sin_arg_pe * sin_lan;

        matrix.e31 = sin_arg_pe * sin_inc;
        matrix.e32 = cos_arg_pe * sin_inc;

        matrix
    }
}

impl OrbitTrait for Orbit {
    fn new(
        eccentricity: f64,
        periapsis: f64,
        inclination: f64,
        arg_pe: f64,
        long_asc_node: f64,
        mean_anomaly: f64,
        parent_mass: f64,
    ) -> Self {
        let cache = Self::get_cached_calculations(
            eccentricity,
            periapsis,
            inclination,
            arg_pe,
            long_asc_node,
        );
        Self {
            eccentricity,
            periapsis,
            inclination,
            arg_pe,
            long_asc_node,
            mean_anomaly,
            parent_mass,
            cache,
        }
    }

    fn with_apoapsis(
        apoapsis: f64,
        periapsis: f64,
        inclination: f64,
        arg_pe: f64,
        long_asc_node: f64,
        mean_anomaly: f64,
        parent_mass: f64,
    ) -> Self {
        let eccentricity = (apoapsis - periapsis) / (apoapsis + periapsis);
        Self::new(
            eccentricity,
            periapsis,
            inclination,
            arg_pe,
            long_asc_node,
            mean_anomaly,
            parent_mass,
        )
    }

    fn get_parent_mass(&self) -> f64 {
        self.parent_mass
    }

    fn get_semi_major_axis(&self) -> f64 {
        self.cache.semi_major_axis
    }

    fn get_semi_minor_axis(&self) -> f64 {
        self.cache.semi_minor_axis
    }

    fn get_linear_eccentricity(&self) -> f64 {
        self.cache.linear_eccentricity
    }

    fn set_apoapsis(&mut self, apoapsis: f64) -> Result<(), ApoapsisSetterError> {
        if apoapsis < 0.0 {
            Err(ApoapsisSetterError::ApoapsisNegative)
        } else if apoapsis < self.periapsis {
            Err(ApoapsisSetterError::ApoapsisLessThanPeriapsis)
        } else {
            self.eccentricity = (apoapsis - self.periapsis) / (apoapsis + self.periapsis);
            self.update_cache();

            Ok(())
        }
    }

    fn set_apoapsis_force(&mut self, apoapsis: f64) {
        let mut apoapsis = apoapsis;
        if apoapsis < self.periapsis && apoapsis > 0.0 {
            (apoapsis, self.periapsis) = (self.periapsis, apoapsis);
        }

        self.eccentricity = (apoapsis - self.periapsis) / (apoapsis + self.periapsis);
        self.update_cache();
    }

    fn get_transformation_matrix(&self) -> Matrix3x2 {
        self.cache.transformation_matrix
    }

    #[inline]
    fn get_eccentricity(&self) -> f64 {
        self.eccentricity
    }

    #[inline]
    fn get_periapsis(&self) -> f64 {
        self.periapsis
    }

    #[inline]
    fn get_inclination(&self) -> f64 {
        self.inclination
    }

    #[inline]
    fn get_arg_pe(&self) -> f64 {
        self.arg_pe
    }

    #[inline]
    fn get_long_asc_node(&self) -> f64 {
        self.long_asc_node
    }

    #[inline]
    fn get_mean_anomaly_at_epoch(&self) -> f64 {
        self.mean_anomaly
    }

    fn set_eccentricity(&mut self, value: f64) {
        self.eccentricity = value;
        self.update_cache();
    }
    fn set_periapsis(&mut self, value: f64) {
        self.periapsis = value;
        self.update_cache();
    }
    fn set_inclination(&mut self, value: f64) {
        self.inclination = value;
        self.update_cache();
    }
    fn set_arg_pe(&mut self, value: f64) {
        self.arg_pe = value;
        self.update_cache();
    }
    fn set_long_asc_node(&mut self, value: f64) {
        self.long_asc_node = value;
        self.update_cache();
    }
    fn set_mean_anomaly_at_epoch(&mut self, value: f64) {
        self.mean_anomaly = value;
        self.update_cache();
    }
}

impl From<CompactOrbit> for Orbit {
    fn from(compact: CompactOrbit) -> Self {
        Self::new(
            compact.eccentricity,
            compact.periapsis,
            compact.inclination,
            compact.arg_pe,
            compact.long_asc_node,
            compact.mean_anomaly,
            compact.parent_mass,
        )
    }
}

impl Default for Orbit {
    /// Creates a unit orbit.
    ///
    /// The unit orbit is a perfect circle of radius 1 and no "tilt".
    fn default() -> Orbit {
        Self::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    }
}
