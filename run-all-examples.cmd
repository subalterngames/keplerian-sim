@echo off
setlocal enabledelayedexpansion

for %%f in (examples\*.rs) do (
    set "filename=%%~nf"
    echo Running example: !filename!
    cargo run --example !filename!
)