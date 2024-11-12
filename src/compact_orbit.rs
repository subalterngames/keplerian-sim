// Shoutouts to Scott Anderson for these! I barely know anything about orbital mechanics LOL
// his repo on his implementation: https://github.com/ScottyRAnderson/Keplerian-Orbits
// his yt vid about this: https://www.youtube.com/watch?v=t89De819YMA (highly underrated btw, check him out)

// However his code is kinda incomplete and doesn't account for longitude of ascending node.
// I found an algorithm to account for it: https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

use crate::{
    keplers_equation, keplers_equation_derivative, keplers_equation_hyperbolic, keplers_equation_hyperbolic_derivative, ApoapsisSetterError, Matrix3x2, Orbit, OrbitTrait
};

/// A struct representing a Keplerian orbit.  
/// This struct minimizes memory footprint by not caching variables.  
/// Because of this, calculations can be slower than caching those variables.  
/// For this reason, you might consider using the `Orbit` struct instead.
#[derive(Clone, Debug)]
pub struct CompactOrbit {
    /// The eccentricity of the orbit.  
    /// e < 1: ellipse  
    /// e = 1: parabola  
    /// e > 1: hyperbola  
    pub eccentricity: f64,

    /// The periapsis of the orbit, in meters.
    pub periapsis: f64,

    /// The inclination of the orbit, in radians.
    pub inclination: f64,

    /// The argument of periapsis of the orbit, in radians.
    pub arg_pe: f64,

    /// The longitude of ascending node of the orbit, in radians.
    pub long_asc_node: f64,

    /// The mean anomaly at orbit epoch, in radians.
    pub mean_anomaly: f64,
}

// Initialization and cache management
impl CompactOrbit {
    pub fn new(
        eccentricity: f64, periapsis: f64,
        inclination: f64, arg_pe: f64, long_asc_node: f64,
        mean_anomaly: f64
    ) -> CompactOrbit {
        return CompactOrbit {
            eccentricity, periapsis,
            inclination, arg_pe, long_asc_node,
            mean_anomaly,
        };
    }

    pub fn with_apoapsis(
        apoapsis: f64, periapsis: f64,
        inclination: f64, arg_pe: f64, long_asc_node: f64,
        mean_anomaly: f64
    ) -> CompactOrbit {
        let eccentricity = (apoapsis - periapsis ) / (apoapsis + periapsis);
        return CompactOrbit::new(eccentricity, periapsis, inclination, arg_pe, long_asc_node, mean_anomaly);
    }

    pub fn new_default() -> CompactOrbit {
        return Self::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    }
}

// The actual orbit position calculations
impl CompactOrbit {
    pub fn get_semi_major_axis(&self) -> f64 {
        return self.periapsis / (1.0 - self.eccentricity);
    }
    
    pub fn get_semi_minor_axis(&self) -> f64 {
        let semi_major_axis = self.get_semi_major_axis();
        let eccentricity_squared = self.eccentricity * self.eccentricity;
        if self.eccentricity < 1.0 {
            return semi_major_axis * (1.0 - eccentricity_squared).sqrt();
        } else {
            return semi_major_axis * (eccentricity_squared - 1.0).sqrt();
        }
    }
    
    pub fn get_linear_eccentricity(&self) -> f64 {
        return self.get_semi_major_axis() - self.periapsis;
    }
    
    /// Gets the apoapsis of the orbit.  
    /// Returns infinity for parabolic orbits.  
    /// Returns negative values for hyperbolic orbits.  
    pub fn get_apoapsis(&self) -> f64 {
        if self.eccentricity == 1.0 {
            return f64::INFINITY;
        } else {
            return self.get_semi_major_axis() * (1.0 + self.eccentricity);
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
    }
}

impl CompactOrbit {
    fn get_transformation_matrix(&self) -> Matrix3x2<f64> {
        let mut matrix = Matrix3x2::<f64>::filled_with(0.0);
        {
            let (sin_inc, cos_inc) = self.inclination.sin_cos();
            let (sin_arg_pe, cos_arg_pe) = self.arg_pe.sin_cos();
            let (sin_lan, cos_lan) = self.long_asc_node.sin_cos();
    
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
        const MAX_ITERATIONS: u16 = 1000;
    
        // Starting guess
        let eccentricity = self.eccentricity;
        let mut eccentric_anomaly = {
            if eccentricity > 0.8 { std::f64::consts::PI }
            else { eccentricity }
        };
        
        for _ in 0..MAX_ITERATIONS {
            // NEWTON'S METHOD
            // x_n+1 = x_n - f(x_n)/f'(x_n)
    
            let next_value = 
                eccentric_anomaly - 
                (
                    keplers_equation(mean_anomaly, eccentric_anomaly, eccentricity) /
                    keplers_equation_derivative(eccentric_anomaly, eccentricity)
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

impl OrbitTrait for CompactOrbit {
    /// Numerically approaches the eccentric anomaly using Newton's method. Not performant; cache the result if you can!
    fn get_eccentric_anomaly(&self, mean_anomaly: f64) -> f64 {
        if self.eccentricity < 1.0 {
            return self.get_eccentric_anomaly_elliptic(mean_anomaly);
        } else {
            return self.get_eccentric_anomaly_hyperbolic(mean_anomaly);
        }
    }

    /// Gets the true anomaly from the mean anomaly. Not performant; cache the result if you can!
    fn get_true_anomaly(&self, mean_anomaly: f64) -> f64 {
        let eccentric_anomaly = self.get_eccentric_anomaly_elliptic(mean_anomaly);

        // https://en.wikipedia.org/wiki/True_anomaly#From_the_eccentric_anomaly
        let eccentricity = self.eccentricity;
        let beta = eccentricity / (1.0 + (1.0 - eccentricity * eccentricity).sqrt());

        return eccentric_anomaly + 2.0 * (beta * eccentric_anomaly.sin() / (1.0 - beta * eccentric_anomaly.cos())).atan();
    }

    /// Multiplies the input 2D vector with the 2x3 transformation matrix used to tilt the flat orbit into 3D space.
    fn tilt_flat_position(&self, x: f64, y: f64) -> (f64, f64, f64) {
        let matrix = &self.get_transformation_matrix();
        return (
            x * matrix.e11 + y * matrix.e12,
            x * matrix.e21 + y * matrix.e22,
            x * matrix.e31 + y * matrix.e32
        );
    }

    /// Gets the 2D position at a certain angle. True anomaly ranges from 0 to tau; anything out of range will wrap around.
    fn get_flat_position_at_angle(&self, true_anomaly: f64) -> (f64, f64) {
        return (
            self.get_semi_major_axis() * (true_anomaly.cos() - self.eccentricity),
            self.get_semi_minor_axis() * true_anomaly.sin()
        );
    }

    /// Gets the 3D position at a certain angle. True anomaly ranges from 0 to tau; anything out of range will wrap around.
    fn get_position_at_angle(&self, true_anomaly: f64) -> (f64, f64, f64) {
        let (x, y) = self.get_flat_position_at_angle(true_anomaly);
        return self.tilt_flat_position(x, y);
    }

    /// Gets the mean anomaly at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    fn get_mean_anomaly_at_time(&self, t: f64) -> f64 {
        let time = t.rem_euclid(1.0);
        return time * std::f64::consts::TAU + self.mean_anomaly;
    }

    /// Gets the eccentric anomaly at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    fn get_eccentric_anomaly_at_time(&self, t: f64) -> f64 {
        return self.get_eccentric_anomaly_elliptic(self.get_mean_anomaly_at_time(t));
    }

    /// Gets the true anomaly at a certain time. t ranges form 0 to 1; anything out of range will wrap around.
    fn get_true_anomaly_at_time(&self, t: f64) -> f64 {
        return self.get_true_anomaly(self.get_mean_anomaly_at_time(t));
    }

    /// Gets the 2D position at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    fn get_flat_position_at_time(&self, t: f64) -> (f64, f64) {
        let true_anomaly = self.get_true_anomaly_at_time(t);
        return self.get_flat_position_at_angle(true_anomaly);
    }

    /// Gets the 3D position at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    fn get_position_at_time(&self, t: f64) -> (f64, f64, f64) {
        let true_anomaly = self.get_true_anomaly_at_time(t);
        return self.get_position_at_angle(true_anomaly);
    }
}

impl From<Orbit> for CompactOrbit {
    fn from(cached: Orbit) -> Self {
        return Self {
            eccentricity: cached.get_eccentricity(),
            periapsis: cached.get_periapsis(),
            inclination: cached.get_inclination(),
            arg_pe: cached.get_arg_pe(),
            long_asc_node: cached.get_long_asc_node(),
            mean_anomaly: cached.get_mean_anomaly()
        };
    }
}

impl Orbit {
    /// Compactify the cached orbit into a compact orbit to increase memory efficiency
    /// while sacrificing calculation speed.
    pub fn compactify(self) -> CompactOrbit {
        CompactOrbit::from(self)
    }
}