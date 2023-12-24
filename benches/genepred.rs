use std::{fs::read_to_string};

use atglib::models::TranscriptRead;
use criterion::{black_box, criterion_group, criterion_main, Criterion};


use atglib::genepredext;


fn read_from_bytes(data: &[u8]) {
    let transcripts = genepredext::Reader::new(data)
        .transcripts()
        .unwrap();
    assert_eq!(transcripts.len(), 1);    
}

fn read_genepred_bench(c: &mut Criterion) {

    c.bench_function("NM_001365057.2", |b| {
        let data = read_to_string("tests/data/NM_001365057.2.genepredext").unwrap();
        b.iter(|| read_from_bytes(black_box(data.as_bytes())))
    });

    c.bench_function("NM_001365408.1", |b| {
        let data = read_to_string("tests/data/NM_001365408.1.genepredext").unwrap();
        b.iter(|| read_from_bytes(black_box(data.as_bytes())))
    });

    c.bench_function("NM_001371720.1", |b| {
        let data = read_to_string("tests/data/NM_001371720.1.genepredext").unwrap();
        b.iter(|| read_from_bytes(black_box(data.as_bytes())))
    });

    c.bench_function("NM_201550.4", |b| {
        let data = read_to_string("tests/data/NM_201550.4.genepredext").unwrap();
        b.iter(|| read_from_bytes(black_box(data.as_bytes())))
    });
}

criterion_group!(
    genepred,
    read_genepred_bench,
);
criterion_main!(genepred);
