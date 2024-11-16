use keplerian_rust::Orbit;
use std::{fs, path::PathBuf};

const CSV_PATH: &str = "out/output-debug-eccentric-anomaly.csv";
const MEAN_ANOMALY_TEST_COUNT: usize = 50;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut orbit = Orbit::new_default();

    let eccentricities = get_eccentricities();
    let mean_anomalies = get_mean_anomalies();

    let data_points = eccentricities.len() * mean_anomalies.len();
    let mut data: Vec<(f64, f64, u32, f64)> = Vec::with_capacity(
        data_points
    );

    println!("Running simulation to get {data_points} data points...");
    
    for e in eccentricities {
        orbit.set_eccentricity(e);

        for M in mean_anomalies.iter() {
            let (E, iters) = 
                orbit.get_eccentric_anomaly_elliptic_debug(*M);
            data.push((*M, e, iters, E));
        }
    }

    println!("Simulation complete\n");

    // Dump summary to stdout

    println!("Summary:");

    {
        let mut min_iters = (0.0, 0.0, u32::MAX, 0.0);
        let mut max_iters = (0.0, 0.0, 0, 0.0);
        
        let mut total_iters = 0;
        let mut diverge_count = 0;

        for (M, e, iters, E) in data.iter() {
            if iters < &min_iters.2 {
                min_iters = (*M, *e, *iters, *E);
            }

            if iters > &max_iters.2 {
                max_iters = (*M, *e, *iters, *E);
            }

            if iters == &1000 {
                diverge_count += 1;
            }

            total_iters += *iters;
        }

        let mean_iters = total_iters as f64 / data_points as f64;

        println!(
            "Minimum # of iterations: {}\n\
             ...with M = {}, e = {}\n\
             ...resulting in E = {}\n",
            min_iters.2, min_iters.0, min_iters.1, min_iters.3
        );
            
        println!(
            "Maximum # of iterations: {}\n\
             ...with M = {}, e = {}\n\
             ...resulting in E = {}\n",
            max_iters.2, max_iters.0, max_iters.1, max_iters.3
        );

        println!(
            "Out of {} data points, {} ({:.2?}%) did not converge\n",
            data_points, diverge_count, diverge_count as f64 / data_points as f64 * 100.0
        );

        println!(
            "Mean # of iterations: {:.3}\n",
            mean_iters
        );
    }


    let csv_string = create_csv(&data);

    let path = PathBuf::from(CSV_PATH);

    println!("Writing data to CSV file...");

    fs::create_dir_all(
        path.parent().expect("Failed to get parent of CSV file")
    )?;

    fs::write(&path, csv_string)?;

    println!("Data written to {path:?}");

    return Ok(());
}

fn create_csv(positions: &Vec<(f64, f64, u32, f64)>) -> String {
    let mut string = String::new();

    string += "M,e,iters,result\n";

    for (M, e, iters, result) in positions {
        string += &format!("{M},{e},{iters},{result}\n");
    }

    return string;
}

fn get_mean_anomalies() -> Vec<f64> {
    let mut vec: Vec<f64> = Vec::with_capacity(MEAN_ANOMALY_TEST_COUNT + 3);
    
    let multiplier = std::f64::consts::TAU / MEAN_ANOMALY_TEST_COUNT as f64;
    for i in 0..MEAN_ANOMALY_TEST_COUNT {
        vec.push(
            i as f64 * multiplier
        );
    }

    // Outstanding values
    vec.push(0.0001);
    vec.push(std::f64::consts::TAU);
    vec.push(999.0);

    return vec;
}

fn get_eccentricities() -> Vec<f64> {
    let mut vec: Vec<f64> = Vec::with_capacity(120);

    for i in 0..100 {
        vec.push(i as f64 / 100.0);
    }

    // Near-parabolic
    for i in 0..20 {
        vec.push(0.99 + i as f64 * 0.0005);
    }

    return vec;
}
