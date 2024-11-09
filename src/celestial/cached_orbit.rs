// Shoutouts to Scott Anderson for these! I barely know anything about orbital mechanics LOL
// his repo on his implementation: https://github.com/ScottyRAnderson/Keplerian-Orbits
// his yt vid about this: https://www.youtube.com/watch?v=t89De819YMA (highly underrated btw, check him out)

// However his code is kinda incomplete and doesn't account for longitude of ascending node.
// I found an algorithm to account for it: https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

/// A struct representing a Keplerian orbit with some cached values.
#[derive(Clone, Debug)]
pub struct Orbit {
    /// The apoapsis of the orbit, in meters.
    apoapsis: f64,

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
#[derive(Clone, Debug)]
struct OrbitCachedCalculations {
    /// The semi-major axis of the orbit, in meters.
    semi_major_axis: f64,

    /// The semi-minor axis of the orbit, in meters.
    semi_minor_axis: f64,

    /// The linear eccentricity of the orbit, in meters.
    linear_eccentricity: f64,

    /// The eccentricity of the orbit, unitless.
    eccentricity: f64,

    /// The transformation matrix to tilt the 2D planar orbit into 3D space.
    transformation_matrix: Matrix3x2<f64>
}
#[derive(Clone, Debug)]
struct Matrix3x2<T> {
    // Element XY
    e11: T, e12: T,
    e21: T, e22: T,
    e31: T, e32: T
}

// Initialization and cache management
impl Orbit {
    pub fn new(
        apoapsis: f64, periapsis: f64,
        inclination: f64, arg_pe: f64, long_asc_node: f64,
        mean_anomaly: f64
    ) -> Orbit {
        let cache = Self::get_cached_calculations(apoapsis, periapsis, inclination, arg_pe, long_asc_node);
        return Orbit {
            apoapsis, periapsis,
            inclination, arg_pe, long_asc_node,
            mean_anomaly,
            cache
        };
    }

    pub fn new_default() -> Orbit {
        return Self::new(1.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    }

    fn update_cache(&mut self) {
        self.cache = Self::get_cached_calculations(
            self.apoapsis,
            self.periapsis,
            self.inclination,
            self.arg_pe,
            self.long_asc_node
        );
    }

    fn get_cached_calculations(
        apoapsis: f64, periapsis: f64,
        inclination: f64, arg_pe: f64, long_asc_node: f64
    ) -> OrbitCachedCalculations {
        let semi_major_axis = (apoapsis + periapsis) / 2.0;
        let semi_minor_axis = (apoapsis * periapsis).sqrt();
        let linear_eccentricity = semi_major_axis - periapsis;
        let eccentricity = linear_eccentricity / semi_major_axis;
        let transformation_matrix = Self::get_transformation_matrix(inclination, arg_pe, long_asc_node);        

        return OrbitCachedCalculations {
            semi_major_axis,
            semi_minor_axis,
            linear_eccentricity,
            eccentricity,
            transformation_matrix
        }
    }

    fn get_transformation_matrix(inclination: f64, arg_pe: f64, long_asc_node: f64) -> Matrix3x2<f64> {
        let mut matrix = Matrix3x2::<f64>::filled_with(0.0);
        {
            let (sin_inc, cos_inc) = inclination.sin_cos();
            let (sin_arg_pe, cos_arg_pe) = arg_pe.sin_cos();
            let (sin_lan, cos_lan) = long_asc_node.sin_cos();

            // https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
            matrix.e11 = sin_arg_pe * cos_lan - sin_arg_pe * cos_inc * sin_lan;
            matrix.e12 = -(sin_arg_pe * cos_lan + cos_arg_pe * cos_inc * sin_lan);
            
            matrix.e21 = cos_arg_pe * sin_lan + sin_arg_pe * cos_inc * cos_lan;
            matrix.e22 = cos_arg_pe * cos_inc * cos_lan - sin_arg_pe * sin_lan;

            matrix.e31 = sin_arg_pe * sin_inc;
            matrix.e32 = cos_arg_pe * sin_inc;
        }
        return matrix;
    }
}

// The actual orbit position calculations
impl Orbit {
    /// Numerically approaches the eccentric anomaly using Newton's method. Not performant; cache the result if you can!
    pub fn get_eccentric_anomaly(&self, mean_anomaly: f64) -> f64 {
        let target_accuracy = 1e-9;
        let max_iterations = 1000;

        // Starting guess
        let mut eccentric_anomaly =
            if self.cache.eccentricity > 0.8 { 3.14 }
            else { self.cache.eccentricity };
        
        for _ in 0..max_iterations {
            // NEWTON'S METHOD
            // x_n+1 = x_n - f(x_n)/f'(x_n)

            let next_value = 
                eccentric_anomaly - 
                (
                    keplers_equation(mean_anomaly, eccentric_anomaly, self.cache.eccentricity) /
                    keplers_equation_derivative(eccentric_anomaly, self.cache.eccentricity)
                );

            let diff = (eccentric_anomaly - next_value).abs();
            eccentric_anomaly = next_value;

            if diff < target_accuracy {
                break;
            }
        }

        return eccentric_anomaly;
    }

    /// Gets the true anomaly from the mean anomaly. Not performant; cache the result if you can!
    pub fn get_true_anomaly(&self, mean_anomaly: f64) -> f64 {
        let eccentric_anomaly = self.get_eccentric_anomaly(mean_anomaly);

        // https://en.wikipedia.org/wiki/True_anomaly#From_the_eccentric_anomaly
        let eccentricity = self.cache.eccentricity;
        let beta = eccentricity / (1.0 + (1.0 - eccentricity * eccentricity).sqrt());

        return eccentric_anomaly + 2.0 * (beta * eccentric_anomaly.sin() / (1.0 - beta * eccentric_anomaly.cos())).atan();
    }

    /// Multiplies the input 2D vector with the 2x3 transformation matrix used to tilt the flat orbit into 3D space.
    pub fn tilt_flat_position(&self, x: f64, y: f64) -> (f64, f64, f64) {
        let matrix = &self.cache.transformation_matrix;
        return (
            x * matrix.e11 + y * matrix.e12,
            x * matrix.e21 + y * matrix.e22,
            x * matrix.e31 + y * matrix.e32
        );
    }

    /// Gets the 2D position at a certain angle. True anomaly ranges from 0 to tau; anything out of range will wrap around.
    pub fn get_flat_position_at_angle(&self, true_anomaly: f64) -> (f64, f64) {
        return (
            self.cache.semi_major_axis * (true_anomaly.cos() - self.cache.eccentricity),
            self.cache.semi_minor_axis * true_anomaly.sin()
        );
    }

    /// Gets the 3D position at a certain angle. True anomaly ranges from 0 to tau; anything out of range will wrap around.
    pub fn get_position_at_angle(&self, true_anomaly: f64) -> (f64, f64, f64) {
        let (x, y) = self.get_flat_position_at_angle(true_anomaly);
        return self.tilt_flat_position(x, y);
    }

    /// Gets the mean anomaly at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    pub fn get_mean_anomaly_at_time(&self, t: f64) -> f64 {
        let time = t.rem_euclid(1.0);
        return time * std::f64::consts::TAU + self.mean_anomaly;
    }

    /// Gets the eccentric anomaly at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    pub fn get_eccentric_anomaly_at_time(&self, t: f64) -> f64 {
        return self.get_eccentric_anomaly(self.get_mean_anomaly_at_time(t));
    }

    /// Gets the true anomaly at a certain time. t ranges form 0 to 1; anything out of range will wrap around.
    pub fn get_true_anomaly_at_time(&self, t: f64) -> f64 {
        return self.get_true_anomaly(self.get_mean_anomaly_at_time(t));
    }

    /// Gets the 2D position at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    pub fn get_flat_position_at_time(&self, t: f64) -> (f64, f64) {
        let true_anomaly = self.get_true_anomaly_at_time(t);
        return self.get_flat_position_at_angle(true_anomaly);
    }

    /// Gets the 3D position at a certain time. t ranges from 0 to 1; anything out of range will wrap around.
    pub fn get_position_at_time(&self, t: f64) -> (f64, f64, f64) {
        let true_anomaly = self.get_true_anomaly_at_time(t);
        return self.get_position_at_angle(true_anomaly);
    }
}

// Getters and setters (boring)
impl Orbit {
    pub fn get_apoapsis           (&self) -> f64 { self.apoapsis }
    pub fn get_periapsis          (&self) -> f64 { self.periapsis }
    pub fn get_inclination        (&self) -> f64 { self.inclination }
    pub fn get_arg_pe             (&self) -> f64 { self.arg_pe }
    pub fn get_long_asc_node      (&self) -> f64 { self.long_asc_node }
    pub fn get_mean_anomaly       (&self) -> f64 { self.mean_anomaly }
    pub fn get_semi_major_axis    (&self) -> f64 { self.cache.semi_major_axis }
    pub fn get_semi_minor_axis    (&self) -> f64 { self.cache.semi_minor_axis }
    pub fn get_linear_eccentricity(&self) -> f64 { self.cache.linear_eccentricity }
    pub fn get_eccentricity       (&self) -> f64 { self.cache.eccentricity }

    pub fn set_apoapsis     (&mut self, value: f64) { self.apoapsis      = value; self.update_cache(); }
    pub fn set_periapsis    (&mut self, value: f64) { self.periapsis     = value; self.update_cache(); }
    pub fn set_inclination  (&mut self, value: f64) { self.inclination   = value; self.update_cache(); }
    pub fn set_arg_pe       (&mut self, value: f64) { self.arg_pe        = value; self.update_cache(); }
    pub fn set_long_asc_node(&mut self, value: f64) { self.long_asc_node = value; self.update_cache(); }
    pub fn set_mean_anomaly (&mut self, value: f64) { self.mean_anomaly  = value; self.update_cache(); }
}

impl Matrix3x2<f64> {
    fn filled_with<T: Copy>(element: T) -> Matrix3x2<T> {
        return Matrix3x2 {
            e11: element, e12: element,
            e21: element, e22: element,
            e31: element, e32: element,
        };
    }
}

fn keplers_equation(mean_anomaly: f64, eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentric_anomaly - (eccentricity * eccentric_anomaly.sin()) - mean_anomaly;
}
fn keplers_equation_derivative(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return 1.0 - (eccentricity * eccentric_anomaly.cos());
}