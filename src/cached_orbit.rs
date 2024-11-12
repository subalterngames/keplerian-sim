// Shoutouts to Scott Anderson for these! I barely know anything about orbital mechanics LOL
// his repo on his implementation: https://github.com/ScottyRAnderson/Keplerian-Orbits
// his yt vid about this: https://www.youtube.com/watch?v=t89De819YMA (highly underrated btw, check him out)

// However his code is kinda incomplete and doesn't account for longitude of ascending node.
// I found an algorithm to account for it: https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

use crate::{
    ApoapsisSetterError,
    CompactOrbit,
    Matrix3x2,
    OrbitTrait,
    OrbitType,
    keplers_equation,
    keplers_equation_derivative,
    keplers_equation_hyperbolic,
    keplers_equation_hyperbolic_derivative
};
type Vec3 = (f64, f64, f64);
type Vec2 = (f64, f64);

/// A struct representing a Keplerian orbit with some cached values.  
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
    
    orbit_type: OrbitType,
}
// Initialization and cache management
impl Orbit {
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

    pub fn with_apoapsis(
        apoapsis: f64, periapsis: f64,
        inclination: f64, arg_pe: f64, long_asc_node: f64,
        mean_anomaly: f64
    ) -> Orbit {
        let eccentricity = (apoapsis - periapsis) / (apoapsis + periapsis);
        return Self::new(eccentricity, periapsis, inclination, arg_pe, long_asc_node, mean_anomaly);
    }

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
        let orbit_type = if eccentricity < 1.0 {
            OrbitType::Elliptic
        } else {
            OrbitType::Hyperbolic
        };

        let semi_major_axis = periapsis / (1.0 - eccentricity);
        let semi_minor_axis = match orbit_type {
            OrbitType::Elliptic => (semi_major_axis * (1.0 - eccentricity * eccentricity)).sqrt(),
            OrbitType::Hyperbolic => (semi_major_axis * (eccentricity * eccentricity - 1.0)).sqrt(),
        };
        let linear_eccentricity = semi_major_axis * eccentricity;
        let transformation_matrix = Self::get_transformation_matrix(inclination, arg_pe, long_asc_node);        

        return OrbitCachedCalculations {
            semi_major_axis,
            semi_minor_axis,
            linear_eccentricity,
            transformation_matrix,
            orbit_type
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

    fn get_eccentric_anomaly_elliptic(&self, mean_anomaly: f64) -> f64 {
        let target_accuracy = 1e-9;
        let max_iterations = 1000;

        // Starting guess
        let mut eccentric_anomaly =
            if self.eccentricity > 0.8 { 3.14 }
            else { self.eccentricity };
        
        for _ in 0..max_iterations {
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

            if diff < target_accuracy {
                break;
            }
        }

        return eccentric_anomaly;
    }

    fn get_eccentric_anomaly_hyperbolic(&self, mean_anomaly: f64) -> f64 {
        let target_accuracy = 1e-9;
        let max_iterations = 1000;
        let mut eccentric_anomaly = mean_anomaly;
        
        for _ in 0..max_iterations {
            // NEWTON'S METHOD
            // x_n+1 = x_n - f(x_n)/f'(x_n)

            let next_value =
                eccentric_anomaly - (
                    keplers_equation_hyperbolic(mean_anomaly, eccentric_anomaly, self.eccentricity) /
                    keplers_equation_hyperbolic_derivative(eccentric_anomaly, self.eccentricity)
                );

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
    /// Numerically approaches the eccentric anomaly using Newton's method. Not performant; cache the result if you can!
    fn get_eccentric_anomaly(&self, mean_anomaly: f64) -> f64 {
        match self.cache.orbit_type {
            OrbitType::Elliptic => self.get_eccentric_anomaly_elliptic(mean_anomaly),
            OrbitType::Hyperbolic => self.get_eccentric_anomaly_hyperbolic(mean_anomaly)
        }
    }

    /// Gets the true anomaly from the mean anomaly. Not performant; cache the result if you can!
    fn get_true_anomaly(&self, mean_anomaly: f64) -> f64 {
        let eccentric_anomaly = self.get_eccentric_anomaly(mean_anomaly);

        // https://en.wikipedia.org/wiki/True_anomaly#From_the_eccentric_anomaly
        let eccentricity = self.eccentricity;
        let beta = eccentricity / (1.0 + (1.0 - eccentricity * eccentricity).sqrt());

        return eccentric_anomaly + 2.0 * (beta * eccentric_anomaly.sin() / (1.0 - beta * eccentric_anomaly.cos())).atan();
    }

    /// Multiplies the input 2D vector with the 2x3 transformation matrix used to tilt the flat orbit into 3D space.
    fn tilt_flat_position(&self, x: f64, y: f64) -> Vec3 {
        let matrix = &self.cache.transformation_matrix;
        return (
            x * matrix.e11 + y * matrix.e12,
            x * matrix.e21 + y * matrix.e22,
            x * matrix.e31 + y * matrix.e32
        );
    }

    /// Gets the 2D position at a certain angle. True anomaly ranges from 0 to tau; anything out of range will wrap around.
    fn get_flat_position_at_angle(&self, true_anomaly: f64) -> Vec2 {
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

    /// Gets the 3D position at a certain angle. True anomaly ranges from 0 to tau; anything out of range will wrap around.
    fn get_position_at_angle(&self, true_anomaly: f64) -> Vec3 {
        let (x, y) = self.get_flat_position_at_angle(true_anomaly);
        return self.tilt_flat_position(x, y);
    }

    /// Gets the mean anomaly at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    fn get_mean_anomaly_at_time(&self, t: f64) -> f64 {
        return t * std::f64::consts::TAU + self.mean_anomaly;
    }

    /// Gets the eccentric anomaly at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    fn get_eccentric_anomaly_at_time(&self, t: f64) -> f64 {
        return self.get_eccentric_anomaly(self.get_mean_anomaly_at_time(t));
    }

    /// Gets the true anomaly at a certain time. t ranges form 0 to 1; anything out of range will wrap around.
    fn get_true_anomaly_at_time(&self, t: f64) -> f64 {
        return self.get_true_anomaly(self.get_mean_anomaly_at_time(t));
    }

    /// Gets the 2D position at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    fn get_flat_position_at_time(&self, t: f64) -> Vec2 {
        let true_anomaly = self.get_true_anomaly_at_time(t);
        return self.get_flat_position_at_angle(true_anomaly);
    }

    /// Gets the 3D position at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    fn get_position_at_time(&self, t: f64) -> Vec3 {
        let true_anomaly = self.get_true_anomaly_at_time(t);
        return self.get_position_at_angle(true_anomaly);
    }
}

// Getters and setters (boring)
impl Orbit {
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

// Apoapsis handling
impl Orbit {
    /// Gets the apoapsis of the orbit.  
    /// Returns infinity for parabolic orbits.  
    /// Returns negative values for hyperbolic orbits.  
    pub fn get_apoapsis(&self) -> f64 {
        if self.eccentricity == 1.0 {
            return f64::INFINITY;
        } else {
            return self.cache.semi_major_axis * (1.0 + self.eccentricity);
        }
    }

    /// Sets the apoapsis of the orbit.  
    /// Errors when the apoapsis is less than the periapsis, or less than zero.  
    /// If you want a setter that does not error, use `set_apoapsis_force`, which will
    /// try its best to interpret what you might have meant, but may have
    /// undesirable behavior.
    pub fn set_apoapsis(&mut self, apoapsis: f64) -> Result<(), ApoapsisSetterError> {
        if apoapsis < 0.0 {
            return Err(ApoapsisSetterError::ApoapsisNegative);
        } else if apoapsis < self.periapsis {
            return Err(ApoapsisSetterError::ApoapsisLessThanPeriapsis);
        }

        self.eccentricity = (apoapsis - self.periapsis) / (apoapsis + self.periapsis);
        self.update_cache();

        return Ok(());
    }

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
    pub fn set_apoapsis_force(&mut self, apoapsis: f64) {
        let mut apoapsis = apoapsis;
        if apoapsis < self.periapsis && apoapsis > 0.0 {
            (apoapsis, self.periapsis) = (self.periapsis, apoapsis);
        }

        self.eccentricity = (apoapsis - self.periapsis) / (apoapsis + self.periapsis);
        self.update_cache();
    }
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