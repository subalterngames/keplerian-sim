#!/bin/sh

for example_file in examples/*.rs; do
    example_name=`basename $example_file .rs`
    echo "Running example $example_name"
    cargo run --example $example_name
done