# Keplerian-rust
A Rust program I ported from C# to simulate Keplerian orbits.

## Dependencies
You'll need Cargo and Rust to run the example binary. You can get it [at rustup.rs](https://rustup.rs/).

## Running
This project is a library crate, so you can't really 'run' that, but you can run the example binary.  
Simply clone the repo using `git clone` and run `cargo run --example simulate` in the project's root directory.

## Resources
I did not come up with the algorithms myself. Here are useful resources for converting Keplerian orbital parameters into Cartesian coordinates:  
- [Wikipedia: Keplerian elements](https://en.wikipedia.org/wiki/Orbital_elements)
- [Scott Anderson's Keplerian-Orbits repo](https://github.com/ScottyRAnderson/Keplerian-Orbits)
- [Scott Anderson's Keplerian Orbits video](https://www.youtube.com/watch?v=t89De819YMA) (highly underrated)
- [Converting Keplerian Orbital Elements to Cartesian State Vectors](https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf)