use std::collections::HashMap;
use needletail::bitkmer::{str_to_bitmer, bitmer_to_str};

use minhashes::{KmerCount, hash_f};


#[derive(Deserialize, Eq, PartialEq, Serialize)]
pub struct JSONSketch {
    pub name: String,
    length: Option<u64>,
    comment: Option<String>,
    metadata: Option<HashMap<String, String>>,
    hashes: Vec<String>,
    kmers: Option<Vec<String>>,
    counts: Option<Vec<u16>>,
}

impl JSONSketch {
    pub fn new(name: &str, length: u64, kmercounts: Vec<KmerCount>) -> Self {
        let mut hash_list = Vec::with_capacity(kmercounts.len());
        let mut kmer_list = Vec::with_capacity(kmercounts.len());
        let mut count_list = Vec::with_capacity(kmercounts.len());
        for hash in &kmercounts {
            hash_list.push(hash.hash.to_string());
            kmer_list.push(bitmer_to_str(hash.kmer));
            count_list.push(hash.count);
        }
        JSONSketch {
            name: String::from(name),
            length: Some(length),
            comment: Some(String::from("")),
            metadata: Some(HashMap::new()),
            hashes: hash_list,
            kmers: Some(kmer_list),
            counts: Some(count_list),
        }
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
                Some(ref v) => kmer = str_to_bitmer(v[i].as_bytes()),
                None => return None,
            }
            let count;
            match self.counts {
                Some(ref v) => count = v[i],
                None => return None,
            }
            kmercount_list.push(KmerCount {
                hash: hash,
                kmer: kmer,
                count: count,
            });
        }
        Some(kmercount_list)
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
            kmer_list.push(hash.kmer.0 as u64);
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
            let bitmer = (*self.kmers)[i as usize];
            let kmer = &bitmer_to_str((bitmer, self.kmer_size));
            kmercounts.push(KmerCount {
                hash: hash_f(&kmer.as_bytes()),
                kmer: (bitmer, self.kmer_size),
                count: (*self.counts)[i as usize],
            });
        }
        Some(kmercounts)
    }
}