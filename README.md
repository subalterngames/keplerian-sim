# Keplerian-rust
A Rust program I ported from C# to simulate Keplerian orbits.

## Dependencies
You'll need Cargo and Rust to run this program. You can get it [at rustup.rs](https://rustup.rs/).

## Running
To run this project, simply clone the repository and run `cargo run` in the project directory.  
If you want to run it with more compiler optimizations, run `cargo run -r` instead to run it in release mode.

## Resources
I did not come up with the algorithms myself. Here are useful resources for converting Keplerian orbital parameters into Cartesian coordinates:  
- [Wikipedia: Keplerian elements](https://en.wikipedia.org/wiki/Orbital_elements)
- [Scott Anderson's Keplerian-Orbits repo](https://github.com/ScottyRAnderson/Keplerian-Orbits)
- [Scott Anderson's Keplerian Orbits video](https://www.youtube.com/watch?v=t89De819YMA) (highly underrated)
- [Converting Keplerian Orbital Elements to Cartesian State Vectors](https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf)