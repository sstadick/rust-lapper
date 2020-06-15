#[macro_use]
extern crate clap;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = clap_app!(bed_cov =>
        (version: crate_version!())
        (author: crate_authors!())
        (@arg file_a: -a +takes_value +required "Path to input file a")
        (@arg file_b: -b +takes_value +required "Path to input file b")
    )
    .get_matches();
    let file_a = matches.value_of("file_a").unwrap();
    let file_b = matches.value_of("file_b").unwrap();

    // Read in all of file 1 into hash / lapper structure
    let mut bed = HashMap::new();
    let file = File::open(file_a)?;
    let mut reader = BufReader::new(file);
    let mut buffer = String::new();
    while reader.read_line(&mut buffer).unwrap() > 0 {
        let mut iter = buffer[..buffer.len() - 1].split('\t');
        let chr = iter.next().unwrap();
        let start = iter.next().unwrap().parse::<u32>().unwrap();
        let stop = iter.next().unwrap().parse::<u32>().unwrap();
        if !bed.contains_key(chr) {
            bed.insert(chr.to_string(), vec![]);
        }
        bed.get_mut(chr).unwrap().push(Interval {
            start,
            stop,
            val: true,
        });
        buffer.clear();
    }

    // Convert to hash of lappers
    let mut lappers = HashMap::new();
    for (key, value) in bed.into_iter() {
        lappers.insert(key, Lapper::new(value));
    }

    // Iter over B and get the values as we go
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    let file = File::open(file_b)?;
    let mut reader = BufReader::new(file);
    buffer.clear();
    while reader.read_line(&mut buffer).unwrap() > 0 {
        let mut iter = buffer[..buffer.len() - 1].split('\t');
        let chr = iter.next().unwrap();
        if let Some(lapper) = lappers.get(chr) {
            let st0 = iter.next().unwrap().parse::<u32>().unwrap();
            let en0 = iter.next().unwrap().parse::<u32>().unwrap();
            let mut cov_st = 0;
            let mut cov_en = 0;
            let mut cov = 0;
            let mut n = 0;
            for iv in lapper.find(st0, en0) {
                n += 1;
                let st1 = if iv.start > st0 { iv.start } else { st0 };
                let en1 = if iv.stop < en0 { iv.stop } else { en0 };
                if st1 > cov_en {
                    cov += cov_en - cov_st;
                    cov_st = st1;
                    cov_en = en1;
                } else {
                    cov_en = if cov_en < en1 { en1 } else { cov_en };
                }
            }
            cov += cov_en - cov_st;
            writeln!(handle, "{}\t{}\t{}\t{}\t{}", chr, st0, en0, n, cov)?;
        } else {
            let start = iter.next().unwrap();
            let stop = iter.next().unwrap();
            // print the default stuff
            writeln!(handle, "{}\t{}\t{}\t0\t0", chr, start, stop)?;
        }
        buffer.clear();
    }
    Ok(())
}
