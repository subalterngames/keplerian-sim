// Shoutouts to Scott Anderson for these! I barely know anything about orbital mechanics LOL
// his repo on his implementation: https://github.com/ScottyRAnderson/Keplerian-Orbits
// his yt vid about this: https://www.youtube.com/watch?v=t89De819YMA (highly underrated btw, check him out)

// However his code is kinda incomplete and doesn't account for longitude of ascending node.
// I found an algorithm to account for it: https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf

use crate::Matrix3x2;

/// A struct representing a Keplerian orbit with some cached values.
#[derive(Clone, Debug)]
pub struct CompactOrbit {
    /// The apoapsis of the orbit, in meters.
    pub apoapsis: f64,

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
        apoapsis: f64, periapsis: f64,
        inclination: f64, arg_pe: f64, long_asc_node: f64,
        mean_anomaly: f64
    ) -> CompactOrbit {
        // let cache = Self::get_cached_calculations(apoapsis, periapsis, inclination, arg_pe, long_asc_node);
        return CompactOrbit {
            apoapsis, periapsis,
            inclination, arg_pe, long_asc_node,
            mean_anomaly,
        };
    }

    pub fn new_default() -> CompactOrbit {
        return Self::new(1.0, 1.0, 0.0, 0.0, 0.0, 0.0);
    }
}

// The actual orbit position calculations
impl CompactOrbit {
    pub fn get_semi_major_axis(&self) -> f64 {
        return (self.apoapsis + self.periapsis) / 2.0;
    }

    pub fn get_semi_minor_axis(&self) -> f64 {
        return (self.apoapsis * self.periapsis).sqrt();
    }

    pub fn get_linear_eccentricity(&self) -> f64 {
        return self.get_semi_major_axis() - self.periapsis;
    }

    pub fn get_eccentricity(&self) -> f64 {
        return self.get_linear_eccentricity() / self.get_semi_major_axis();
    }

    pub fn get_transformation_matrix(&self) -> Matrix3x2<f64> {
        let mut matrix = Matrix3x2::<f64>::filled_with(0.0);
        {
            let (sin_inc, cos_inc) = self.inclination.sin_cos();
            let (sin_arg_pe, cos_arg_pe) = self.arg_pe.sin_cos();
            let (sin_lan, cos_lan) = self.long_asc_node.sin_cos();

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

    /// Numerically approaches the eccentric anomaly using Newton's method. Not performant; cache the result if you can!
    pub fn get_eccentric_anomaly(&self, mean_anomaly: f64) -> f64 {
        let target_accuracy = 1e-9;
        const MAX_ITERATIONS: u16 = 1000;

        // Starting guess
        let eccentricity = self.get_eccentricity();
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

    /// Gets the true anomaly from the mean anomaly. Not performant; cache the result if you can!
    pub fn get_true_anomaly(&self, mean_anomaly: f64) -> f64 {
        let eccentric_anomaly = self.get_eccentric_anomaly(mean_anomaly);

        // https://en.wikipedia.org/wiki/True_anomaly#From_the_eccentric_anomaly
        let eccentricity = self.get_eccentricity();
        let beta = eccentricity / (1.0 + (1.0 - eccentricity * eccentricity).sqrt());

        return eccentric_anomaly + 2.0 * (beta * eccentric_anomaly.sin() / (1.0 - beta * eccentric_anomaly.cos())).atan();
    }

    /// Multiplies the input 2D vector with the 2x3 transformation matrix used to tilt the flat orbit into 3D space.
    pub fn tilt_flat_position(&self, x: f64, y: f64) -> (f64, f64, f64) {
        let matrix = &self.get_transformation_matrix();
        return (
            x * matrix.e11 + y * matrix.e12,
            x * matrix.e21 + y * matrix.e22,
            x * matrix.e31 + y * matrix.e32
        );
    }

    /// Gets the 2D position at a certain angle. True anomaly ranges from 0 to tau; anything out of range will wrap around.
    pub fn get_flat_position_at_angle(&self, true_anomaly: f64) -> (f64, f64) {
        return (
            self.get_semi_major_axis() * (true_anomaly.cos() - self.get_eccentricity()),
            self.get_semi_minor_axis() * true_anomaly.sin()
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


fn keplers_equation(mean_anomaly: f64, eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return eccentric_anomaly - (eccentricity * eccentric_anomaly.sin()) - mean_anomaly;
}
fn keplers_equation_derivative(eccentric_anomaly: f64, eccentricity: f64) -> f64 {
    return 1.0 - (eccentricity * eccentric_anomaly.cos());
}