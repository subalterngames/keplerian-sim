use crate::{
    ApoapsisSetterError,
    CompactOrbit,
    Matrix3x2,
    OrbitTrait,
    keplers_equation,
    keplers_equation_derivative,
    keplers_equation_hyperbolic,
    keplers_equation_hyperbolic_derivative,
    keplers_equation_hyperbolic_second_derivative
};
use std::f64::consts::{PI, TAU};

/// A struct representing a Keplerian orbit with some cached values.
/// 
/// This struct consumes significantly more memory because of the cache.  
/// However, this will speed up orbital calculations.  
/// If memory efficiency is your goal, you may consider using the `CompactOrbit` struct instead.  
#[derive(Clone, Debug, PartialEq)]
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
    cache: OrbitCachedCalculations,
}
#[derive(Clone, Debug, PartialEq)]
struct OrbitCachedCalculations {
    /// The semi-major axis of the orbit, in meters.
    semi_major_axis: f64,

    /// The semi-minor axis of the orbit, in meters.
    semi_minor_axis: f64,

    /// The linear eccentricity of the orbit, in meters.
    linear_eccentricity: f64,

    /// The transformation matrix to tilt the 2D planar orbit into 3D space.
    transformation_matrix: Matrix3x2<f64>,
}
// Initialization and cache management
impl Orbit {
    /// Creates a new orbit with the given parameters.
    /// 
    /// Note: This function uses eccentricity instead of apoapsis.  
    /// If you want to provide an apoapsis instead, consider using the
    /// [`Orbit::with_apoapsis`] function instead.
    /// 
    /// ### Parameters
    /// - `eccentricity`: The eccentricity of the orbit.
    /// - `periapsis`: The periapsis of the orbit, in meters.
    /// - `inclination`: The inclination of the orbit, in radians.
    /// - `arg_pe`: The argument of periapsis of the orbit, in radians.
    /// - `long_asc_node`: The longitude of ascending node of the orbit, in radians.
    /// - `mean_anomaly`: The mean anomaly of the orbit, in radians.
    pub fn new(
        eccentricity: f64, periapsis: f64,
        inclination: f64, arg_pe: f64, long_asc_node: f64,
        mean_anomaly: f64
    ) -> Orbit {
        let cache = Self::get_cached_calculations(
            eccentricity, periapsis,
            inclination, arg_pe, long_asc_node
        );
        return Orbit {
            eccentricity, periapsis,
            inclination, arg_pe, long_asc_node,
            mean_anomaly,
            cache
        };
    }

    /// Creates a new orbit with the given parameters.
    /// 
    /// Note: This function uses apoapsis instead of eccentricity.  
    /// Because of this, it's not recommended to initialize
    /// parabolic or hyperbolic 'orbits' with this function.  
    /// If you're looking to initialize a parabolic or hyperbolic
    /// trajectory, consider using the [`Orbit::new`] function instead.
    /// 
    /// ### Parameters
    /// - `apoapsis`: The apoapsis of the orbit, in meters.
    /// - `periapsis`: The periapsis of the orbit, in meters.
    /// - `inclination`: The inclination of the orbit, in radians.
    /// - `arg_pe`: The argument of periapsis of the orbit, in radians.
    /// - `long_asc_node`: The longitude of ascending node of the orbit, in radians.
    /// - `mean_anomaly`: The mean anomaly of the orbit, in radians.
    pub fn with_apoapsis(
        apoapsis: f64, periapsis: f64,
        inclination: f64, arg_pe: f64, long_asc_node: f64,
        mean_anomaly: f64
    ) -> Orbit {
        let eccentricity = (apoapsis - periapsis) / (apoapsis + periapsis);
        return Self::new(eccentricity, periapsis, inclination, arg_pe, long_asc_node, mean_anomaly);
    }

    /// Creates a unit orbit.
    /// 
    /// The unit orbit is a perfect circle of radius 1 and no "tilt".
    pub fn new_default() -> Orbit {
        return Self::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    }

    fn update_cache(&mut self) {
        self.cache = Self::get_cached_calculations(
            self.eccentricity,
            self.periapsis,
            self.inclination,
            self.arg_pe,
            self.long_asc_node
        );
    }

    fn get_cached_calculations(
        eccentricity: f64, periapsis: f64,
        inclination: f64, arg_pe: f64, long_asc_node: f64
    ) -> OrbitCachedCalculations {
        let semi_major_axis = periapsis / (1.0 - eccentricity);
        let semi_minor_axis =
            semi_major_axis * (1.0 - eccentricity * eccentricity).abs().sqrt();
        let linear_eccentricity = semi_major_axis * eccentricity;
        let transformation_matrix = Self::get_transformation_matrix(inclination, arg_pe, long_asc_node);        

        return OrbitCachedCalculations {
            semi_major_axis,
            semi_minor_axis,
            linear_eccentricity,
            transformation_matrix
        };
    }

    fn get_transformation_matrix(inclination: f64, arg_pe: f64, long_asc_node: f64) -> Matrix3x2<f64> {
        let mut matrix = Matrix3x2::<f64>::filled_with(0.0);
        {
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
        }
        return matrix;
    }
}

/// A constant used to get the initial seed for the eccentric anomaly.
/// 
/// It's very arbitrary, but according to some testing, a value just
/// below 1 works better than exactly 1.
/// 
/// Source:
/// "Two fast and accurate routines for solving the elliptic Kepler
/// equation for all values of the eccentricity and mean anomaly"
/// by Daniele Tommasini and David N. Olivieri,
/// section 2.1.2, 'The "rational seed"'
/// 
/// https://doi.org/10.1051/0004-6361/202141423
const B: f64 = 0.999999;

/// The maximum number of iterations for the Newton-Raphson method.
/// 
/// This is used to prevent infinite loops in case the method fails to converge.
const NEWTON_MAX_ITERS: u32 = 1000;

const PI_SQUARED: f64 = PI * PI;

/// The target accuracy for Newton's method.
/// 
/// Except not. It's more complicated than that.
/// 
/// "Two fast and accurate routines for solving the elliptic Kepler
/// equation for all values of the eccentricity and mean anomaly"
/// by Daniele Tommasini and David N. Olivieri,
/// section 2.1.1. 'The iteration stopping condition' says:  
/// "As we shall demonstrate in Sect. 4.2, Eq. (9) holds
/// whenever the accuracy is set to a level BIG_EPSILON ≲ 10−4 rad."
/// 
/// The paper represents this value as a fancy E.  
/// It looks like a big epsilon.
/// 
/// Because of this I set it to 1e-11.
/// 
/// https://doi.org/10.1051/0004-6361/202141423
const TARGET_ACCURACY: f64 = 1e-11;

/// The machine epsilon for f64.
/// 
/// Source:
/// "Two fast and accurate routines for solving the elliptic Kepler
/// equation for all values of the eccentricity and mean anomaly"
/// by Daniele Tommasini and David N. Olivieri,
/// section 2.1.1. 'The iteration stopping condition' says:  
/// "the machine epsilon ϵ has been introduced"
/// 
/// The paper represents this value as a lowercase epsilon.
const MACHINE_EPSILON: f64 = f64::EPSILON;

// Kepler solvers
impl Orbit {
    // "Two fast and accurate routines for solving
    // the elliptic Kepler equation for all values
    // of the eccentricity and mean anomaly" by
    // Daniele Tommasini and David N. Olivieri
    // 
    // https://doi.org/10.1051/0004-6361/202141423
    fn get_eccentric_anomaly_elliptic(&self, mean_anomaly: f64) -> f64 {
        // Use the symmetry and periodicity of the eccentric anomaly
        // Equation 2 of the aforementioned paper
        if mean_anomaly > PI {
            return self.get_eccentric_anomaly_elliptic(mean_anomaly - TAU);
        }
        if mean_anomaly < 0.0 {
            return -self.get_eccentric_anomaly_elliptic(-mean_anomaly);
        }
        
        // Starting guess
        // Section 2.1.2, 'The "rational seed"',
        // Equation 19, of the aforementioned paper
        //
        // E_0 = M + (4beM(pi - M)) / (8eM + 4e(e-pi) + pi^2)
        // where:
        // e = eccentricity
        // M = mean anomaly
        // pi = the constant PI
        // b = the constant B
        let mut eccentric_anomaly =
            mean_anomaly +
            (4.0 * self.eccentricity * B * mean_anomaly * (PI - mean_anomaly)) /
            (
                8.0 * self.eccentricity * mean_anomaly +
                4.0 * self.eccentricity * (self.eccentricity - PI) +
                PI_SQUARED
            );
        
        for _ in 0..NEWTON_MAX_ITERS {
            // NEWTON'S METHOD
            // x_n+1 = x_n - f(x_n)/f'(x_n)

            let next_value = 
                eccentric_anomaly - 
                (
                    keplers_equation(mean_anomaly, eccentric_anomaly, self.eccentricity) /
                    keplers_equation_derivative(eccentric_anomaly, self.eccentricity)
                );

            let diff = (eccentric_anomaly - next_value).abs();
            eccentric_anomaly = next_value;

            // Section 2.1.1, 'The iteration stopping condition',
            // Equation 9, of the aforementioned paper, says:
            // 
            // delta_n^2 < (2(1 - e cos E_n) * fancy_e) / (e + machine_epsilon)
            //
            // we can rearrange it to remove the slow division into a multiplication:
            //
            // delta_n^2 * (e + machine_epsilon) < 2(1 - e cos E_n) * fancy_e

            if
                diff * diff * (self.eccentricity + MACHINE_EPSILON) <
                2.0 * (1.0 - self.eccentricity * eccentric_anomaly.cos()) * TARGET_ACCURACY
            {
                break;
            }
        }

        return eccentric_anomaly;
    }

    // "Two fast and accurate routines for solving
    // the elliptic Kepler equation for all values
    // of the eccentricity and mean anomaly" by
    // Daniele Tommasini and David N. Olivieri
    // 
    // https://doi.org/10.1051/0004-6361/202141423
    #[doc(hidden)]
    #[cfg(debug_assertions)]
    pub fn get_eccentric_anomaly_elliptic_debug(&self, mean_anomaly: f64) -> (f64, u32) {
        // Use the symmetry and periodicity of the eccentric anomaly
        // Equation 2 of the aforementioned paper
        if mean_anomaly > PI {
            return self.get_eccentric_anomaly_elliptic_debug(mean_anomaly - TAU);
        }
        if mean_anomaly < 0.0 {
            // return -self.get_eccentric_anomaly_elliptic_debug(-mean_anomaly);
            let (res, iters) = self.get_eccentric_anomaly_elliptic_debug(-mean_anomaly);
            return (-res, iters);
        }
        
        let mut iterations: u32 = 0;

        // Starting guess
        // Section 2.1.2, 'The "rational seed"',
        // Equation 19, of the aforementioned paper
        //
        // E_0 = M + (4beM(pi - M)) / (8eM + 4e(e-pi) + pi^2)
        // where:
        // e = eccentricity
        // M = mean anomaly
        // pi = the constant PI
        // b = the constant B
        let mut eccentric_anomaly =
            mean_anomaly +
            (4.0 * self.eccentricity * B * mean_anomaly * (PI - mean_anomaly)) /
            (
                8.0 * self.eccentricity * mean_anomaly +
                4.0 * self.eccentricity * (self.eccentricity - PI) +
                PI_SQUARED
            );
        
        for _ in 0..NEWTON_MAX_ITERS {
            // NEWTON'S METHOD
            // x_n+1 = x_n - f(x_n)/f'(x_n)

            let next_value = 
                eccentric_anomaly - 
                (
                    keplers_equation(mean_anomaly, eccentric_anomaly, self.eccentricity) /
                    keplers_equation_derivative(eccentric_anomaly, self.eccentricity)
                );

            let diff = (eccentric_anomaly - next_value).abs();
            eccentric_anomaly = next_value;

            iterations += 1;

            // Section 2.1.1, 'The iteration stopping condition',
            // Equation 9, of the aforementioned paper, says:
            // 
            // delta_n^2 < (2(1 - e cos E_n) * fancy_e) / (e + machine_epsilon)
            //
            // we can rearrange it to remove the slow division into a multiplication:
            //
            // delta_n^2 * (e + machine_epsilon) < 2(1 - e cos E_n) * fancy_e

            if
                diff * diff * (self.eccentricity + MACHINE_EPSILON) <
                2.0 * (1.0 - self.eccentricity * eccentric_anomaly.cos()) * TARGET_ACCURACY
            {
                break;
            }
        }

        assert_eq!(
            eccentric_anomaly.to_bits(),
            self.get_eccentric_anomaly_elliptic(mean_anomaly).to_bits(),
            "Desync between debug and regular versions of get_eccentric_anomaly_elliptic!"
        );

        if iterations == NEWTON_MAX_ITERS {
            eprintln!(
                "Warning: get_eccentric_anomaly_elliptic_debug failed to converge after {iterations} iterations\n\
                ..With params:\n\
                ....mean_anomaly: {mean_anomaly}\n\
                ....eccentricity: {}", self.eccentricity,
            );
        }

        return (eccentric_anomaly, iterations);
    }

    fn get_eccentric_anomaly_hyperbolic(&self, mean_anomaly: f64) -> f64 {
        let target_accuracy = 1e-9;
        let max_iterations = 1000;
        let mut eccentric_anomaly =
            (2.0 * mean_anomaly.abs() / self.eccentricity).ln().max(0.01) *
            mean_anomaly.signum();

        for _ in 0..max_iterations {
            // HALLEY'S METHOD
            // https://en.wikipedia.org/wiki/Halley%27s_method
            // x_n+1 = x_n - f(x_n) / (f'(x_n) - f(x_n) * f''(x_n) / (2 * f'(x_n)))

            let f = keplers_equation_hyperbolic(mean_anomaly, eccentric_anomaly, self.eccentricity);
            let fp = keplers_equation_hyperbolic_derivative(eccentric_anomaly, self.eccentricity);
            let fpp = keplers_equation_hyperbolic_second_derivative(eccentric_anomaly, self.eccentricity);

            let denominator = fp - f * fpp / (2.0 * fp);

            if denominator.abs() < 1e-30 || !denominator.is_finite() {
                // dangerously close to div-by-zero, break out
                break;
            }

            let next_value = eccentric_anomaly - f / denominator;

            let diff = (eccentric_anomaly - next_value).abs();
            eccentric_anomaly = next_value;

            if diff < target_accuracy {
                break;
            }
        }

        return eccentric_anomaly;
    }
}

// The actual orbit position calculations
impl OrbitTrait for Orbit {
    fn get_semi_major_axis(&self) -> f64 {
        return self.cache.semi_major_axis;
    }

    fn get_semi_minor_axis(&self) -> f64 {
        return self.cache.semi_minor_axis;
    }

    fn get_linear_eccentricity(&self) -> f64 {
        return self.cache.linear_eccentricity;
    }

    fn get_apoapsis(&self) -> f64 {
        if self.eccentricity == 1.0 {
            return f64::INFINITY;
        } else {
            return self.cache.semi_major_axis * (1.0 + self.eccentricity);
        }
    }

    fn set_apoapsis(&mut self, apoapsis: f64) -> Result<(), ApoapsisSetterError> {
        if apoapsis < 0.0 {
            return Err(ApoapsisSetterError::ApoapsisNegative);
        } else if apoapsis < self.periapsis {
            return Err(ApoapsisSetterError::ApoapsisLessThanPeriapsis);
        }

        self.eccentricity = (apoapsis - self.periapsis) / (apoapsis + self.periapsis);
        self.update_cache();

        return Ok(());
    }

    fn set_apoapsis_force(&mut self, apoapsis: f64) {
        let mut apoapsis = apoapsis;
        if apoapsis < self.periapsis && apoapsis > 0.0 {
            (apoapsis, self.periapsis) = (self.periapsis, apoapsis);
        }

        self.eccentricity = (apoapsis - self.periapsis) / (apoapsis + self.periapsis);
        self.update_cache();
    }

    fn get_transformation_matrix(&self) -> Matrix3x2<f64> {
        return self.cache.transformation_matrix;
    }

    fn get_eccentric_anomaly(&self, mean_anomaly: f64) -> f64 {
        if self.eccentricity < 1.0 {
            self.get_eccentric_anomaly_elliptic(mean_anomaly)
        } else {
            self.get_eccentric_anomaly_hyperbolic(mean_anomaly)
        }
    }

    fn get_true_anomaly(&self, mean_anomaly: f64) -> f64 {
        let eccentric_anomaly = self.get_eccentric_anomaly(mean_anomaly);

        if self.eccentricity < 1.0 {
            // https://en.wikipedia.org/wiki/True_anomaly#From_the_eccentric_anomaly
            let eccentricity = self.eccentricity;
            let beta = eccentricity / (1.0 + (1.0 - eccentricity * eccentricity).sqrt());
    
            return eccentric_anomaly + 2.0 * 
                (beta * eccentric_anomaly.sin() / (1.0 - beta * eccentric_anomaly.cos()))
                .atan();
        } else {
            // idk copilot got this
            return 2.0 *
                ((self.eccentricity + 1.0) / (self.eccentricity - 1.0)).sqrt().atan() *
                eccentric_anomaly.tanh();
        }
    }

    fn get_flat_position_at_angle(&self, true_anomaly: f64) -> (f64, f64) {
        // Polar equation for conic section:
        // r = p / (1 + e*cos(theta))
        // ...where:
        // r = radius from focus
        // p = semi-latus rectum (periapsis * (1 + eccentricity))
        // e = eccentricity
        // theta = true anomaly

        let semi_latus_rectum = self.periapsis * (1.0 + self.eccentricity);
        let radius = semi_latus_rectum / (1.0 + self.eccentricity * true_anomaly.cos());

        // We then convert it into Cartesion (X, Y) coordinates:
        return (
            radius * true_anomaly.cos(),
            radius * true_anomaly.sin()
        );
    }

    fn get_mean_anomaly_at_time(&self, t: f64) -> f64 {
        return t * TAU + self.mean_anomaly;
    }
}

// Getters and setters (boring)
impl Orbit {
    // These functions are extremely intuitive and do not require documentation
    #![allow(missing_docs)]
    pub fn get_periapsis          (&self) -> f64 { self.periapsis }
    pub fn get_inclination        (&self) -> f64 { self.inclination }
    pub fn get_arg_pe             (&self) -> f64 { self.arg_pe }
    pub fn get_long_asc_node      (&self) -> f64 { self.long_asc_node }
    pub fn get_mean_anomaly       (&self) -> f64 { self.mean_anomaly }
    pub fn get_semi_major_axis    (&self) -> f64 { self.cache.semi_major_axis }
    pub fn get_semi_minor_axis    (&self) -> f64 { self.cache.semi_minor_axis }
    pub fn get_linear_eccentricity(&self) -> f64 { self.cache.linear_eccentricity }
    pub fn get_eccentricity       (&self) -> f64 { self.eccentricity }

    pub fn set_eccentricity (&mut self, value: f64) { self.eccentricity  = value; self.update_cache(); }
    pub fn set_periapsis    (&mut self, value: f64) { self.periapsis     = value; self.update_cache(); }
    pub fn set_inclination  (&mut self, value: f64) { self.inclination   = value; self.update_cache(); }
    pub fn set_arg_pe       (&mut self, value: f64) { self.arg_pe        = value; self.update_cache(); }
    pub fn set_long_asc_node(&mut self, value: f64) { self.long_asc_node = value; self.update_cache(); }
    pub fn set_mean_anomaly (&mut self, value: f64) { self.mean_anomaly  = value; self.update_cache(); }
}

impl From<CompactOrbit> for Orbit {
    fn from(compact: CompactOrbit) -> Self {
        return Self::new(
            compact.eccentricity,
            compact.periapsis,
            compact.inclination,
            compact.arg_pe,
            compact.long_asc_node,
            compact.mean_anomaly
        );
    }
}

impl CompactOrbit {
    /// Expand the compact orbit into a cached orbit to increase calculation speed
    /// while sacrificing memory efficiency.
    pub fn expand(self) -> Orbit {
        Orbit::from(self)
    }
}