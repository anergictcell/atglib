use std::fs::read_to_string;

use atglib::models::TranscriptRead;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use atglib::gtf;

fn read_from_bytes(data: &[u8]) {
    let transcripts = gtf::Reader::new(data).transcripts().unwrap();
    assert!(!transcripts.is_empty());
}

fn read_gtf_bench(c: &mut Criterion) {
    c.bench_function("parse example GTF", |b| {
        let data = read_to_string("tests/data/example.gtf").unwrap();
        b.iter(|| read_from_bytes(black_box(data.as_bytes())))
    });
}

criterion_group!(gtf, read_gtf_bench,);
criterion_main!(gtf);
