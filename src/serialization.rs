use std::collections::HashMap;
use needletail::bitkmer::{bytes_to_bitmer, bitmer_to_bytes};

use filtering::{FilterParams, filter_sketch};
use minhashes::{KmerCount, hash_f};


#[allow(non_snake_case)]
#[derive(Serialize, Deserialize)]
pub struct SketchDistance {
    pub containment: f64,
    pub jaccard: f64,
    pub mashDistance: f64,
    pub commonHashes: u64,
    pub totalHashes: u64,
    pub query: String,
    pub reference: String,
}

#[allow(non_snake_case)]
#[derive(Serialize, Deserialize)]
pub struct JSONMultiSketch {
    pub kmer: u8,
    pub alphabet: String,
    pub preserveCase: bool,
    pub canonical: bool,
    pub sketchSize: u32,
    pub hashType: String,
    pub hashBits: u16,
    pub hashSeed: u64,
    pub sketches: Vec<JSONSketch>,
}

#[allow(non_snake_case)]
#[derive(Deserialize, Eq, PartialEq, Clone, Serialize)]
pub struct JSONSketch {
    pub name: String,
    pub seqLength: Option<u64>,
    pub comment: Option<String>,
    pub filters: Option<HashMap<String, String>>,
    hashes: Vec<String>,
    kmers: Option<Vec<String>>,
    counts: Option<Vec<u16>>,
}

impl JSONSketch {
    pub fn new(name: &str, length: u64, kmercounts: Vec<KmerCount>, filters: &HashMap<String, String>) -> Self {
        let mut hash_list = Vec::with_capacity(kmercounts.len());
        let mut kmer_list = Vec::with_capacity(kmercounts.len());
        let mut count_list = Vec::with_capacity(kmercounts.len());
        for hash in &kmercounts {
            hash_list.push(hash.hash.to_string());
            kmer_list.push(String::from_utf8(bitmer_to_bytes(hash.kmer)).unwrap());
            count_list.push(hash.count);
        }
        JSONSketch {
            name: String::from(name),
            seqLength: Some(length),
            comment: Some(String::from("")),
            filters: Some(filters.clone()),
            hashes: hash_list,
            kmers: Some(kmer_list),
            counts: Some(count_list),
        }
    }

    pub fn len(&self) -> usize {
        self.hashes.len()
    }

    pub fn get_kmers(&self) -> Option<Vec<KmerCount>> {
        let mut kmercount_list = Vec::with_capacity(self.hashes.len());
        for i in 0..self.hashes.len() {
            let hash;
            match self.hashes[i].parse::<usize>() {
                Ok(t) => hash = t,
                Err(_) => return None,
            }
            let kmer;
            match self.kmers {
                Some(ref v) => kmer = bytes_to_bitmer(v[i].as_bytes()),
                None => kmer = (0, 0),
            }
            let count;
            match self.counts {
                Some(ref v) => count = v[i],
                None => count = 1,
            }
            kmercount_list.push(KmerCount {
                hash: hash,
                kmer: kmer,
                count: count,
                extra_count: count / 2,
            });
        }
        Some(kmercount_list)
    }

    pub fn apply_filtering(&mut self, filters: &FilterParams) -> bool {
        let hashes: Vec<KmerCount>;
        match self.get_kmers() {
            Some(h) => hashes = h,
            None => return false,
        };
        let (filtered_hashes, filter_stats) = filter_sketch(&hashes, &filters);
        let mut hash_list = Vec::with_capacity(filtered_hashes.len());
        let mut kmer_list = Vec::with_capacity(filtered_hashes.len());
        let mut count_list = Vec::with_capacity(filtered_hashes.len());
        for hash in &filtered_hashes {
            hash_list.push(hash.hash.to_string());
            kmer_list.push(String::from_utf8(bitmer_to_bytes(hash.kmer)).unwrap());
            count_list.push(hash.count);
        }
        self.hashes = hash_list;
        self.kmers = Some(kmer_list);
        self.counts = Some(count_list);
        self.filters = Some(filter_stats);
        true
    }
}


#[repr(C, packed)]
pub struct BinarySketch {
    len: u32,
    kmer_size: u8,
    kmers: Box<[u64]>,
    counts: Box<[u16]>,
}

impl BinarySketch {
    pub fn new(kmercounts: Vec<KmerCount>) -> Self {
        let mut kmer_list = Vec::with_capacity(kmercounts.len());
        let mut count_list = Vec::with_capacity(kmercounts.len());
        for hash in &kmercounts {
            kmer_list.push(hash.kmer.0);
            count_list.push(hash.count);
        }
        BinarySketch {
            len: kmercounts.len() as u32,
            kmer_size: kmercounts[0].kmer.1,
            kmers: kmer_list.into_boxed_slice(), 
            counts: count_list.into_boxed_slice(),
        }
    }

    pub fn get_kmers(&self) -> Option<Vec<KmerCount>> {
        let mut kmercounts = Vec::with_capacity(self.len as usize);
        for i in 0..self.len {
            let bitmer = ((*self.kmers)[i as usize], self.kmer_size);
            let kmer = bitmer_to_bytes(bitmer);
            kmercounts.push(KmerCount {
                // there's an assumption here that the seed is 0
                hash: hash_f(&kmer, 0),
                kmer: bitmer,
                count: (*self.counts)[i as usize],
                extra_count: 0,
            });
        }
        Some(kmercounts)
    }
}
