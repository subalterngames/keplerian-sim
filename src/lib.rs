mod cached_orbit;
mod compact_orbit;
mod body;
mod universe;
pub mod body_presets;

pub use cached_orbit::Orbit;
pub use compact_orbit::CompactOrbit;
pub use body::Body;
pub use universe::Universe;

#[derive(Clone, Debug, PartialEq)]
pub struct Matrix3x2<T> {
    // Element XY
    e11: T, e12: T,
    e21: T, e22: T,
    e31: T, e32: T
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

type Vec3 = (f64, f64, f64);
type Vec2 = (f64, f64);

pub trait OrbitTrait {
    fn get_eccentric_anomaly(&self, mean_anomaly: f64) -> f64;
    fn get_true_anomaly(&self, mean_anomaly: f64) -> f64;
    fn get_mean_anomaly_at_time(&self, t: f64) -> f64;
    fn get_eccentric_anomaly_at_time(&self, t: f64) -> f64;
    fn get_true_anomaly_at_time(&self, t: f64) -> f64;
    fn get_position_at_angle(&self, angle: f64) -> Vec3;
    fn get_flat_position_at_angle(&self, angle: f64) -> Vec2;
    fn get_position_at_time(&self, t: f64) -> Vec3;
    fn get_flat_position_at_time(&self, t: f64) -> Vec2;
    fn tilt_flat_position(&self, x: f64, y: f64) -> Vec3;
}

#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum OrbitType {
    /// Eccentricity is less than 1.  
    /// Orbit is shaped mostly like an ellipse or circle.
    Elliptic,

    /// Eccentricity is more than or equal to 1.  
    /// Orbit is shaped mostly like a parabola or hyperbola.  
    /// The apoapsis does not exist.
    Hyperbolic
}

#[cfg(test)]
mod tests;