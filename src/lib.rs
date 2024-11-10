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