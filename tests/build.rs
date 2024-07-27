use std::process::Command;
use std::env;
use std::path::{Path, PathBuf};

fn workspace_dir() -> PathBuf {
    let output = std::process::Command::new(env!("CARGO"))
        .arg("locate-project")
        .arg("--workspace")
        .arg("--message-format=plain")
        .output()
        .unwrap()
        .stdout;
    let cargo_path = Path::new(std::str::from_utf8(&output).unwrap().trim());
    cargo_path.parent().unwrap().to_path_buf()
}

fn main() {
    let workspace = workspace_dir();
    let str = workspace.to_str().unwrap();

    let status = Command::new("gcc").args(&["-o", "../target/release/ffi_test.o", "src/test.c",  &format!("-L{str}/target/release/"), "-lgbam_tools_cffi", "-lhts"])
                       .status().unwrap();

    if !status.success() {
        panic!("C code compilation failed");
    }
}