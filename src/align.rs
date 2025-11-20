// IMPORTANT: DO NOT ADD ANY COMMENTS

pub fn prefix_match(seq: &[u8], bc: &[u8], mismatches: usize) -> bool {
    if bc.len() > seq.len() { return false; }
    let mut mm = 0usize;
    for i in 0..bc.len() {
        let b = bc[i];
        if b == b'N' { continue; }
        if b != seq[i] { mm += 1; if mm > mismatches { return false; } }
    }
    true
}

pub fn suffix_match(seq: &[u8], bc: &[u8], mismatches: usize) -> bool {
    if bc.len() > seq.len() { return false; }
    let offset = seq.len() - bc.len();
    let mut mm = 0usize;
    for i in 0..bc.len() {
        let b = bc[i];
        if b == b'N' { continue; }
        if b != seq[offset + i] { mm += 1; if mm > mismatches { return false; } }
    }
    true
}

pub fn extract_umi_from_suffix(seq: &[u8], bc: &[u8]) -> Option<Vec<u8>> {
    if bc.len() > seq.len() { return None; }
    let offset = seq.len() - bc.len();
    let mut umi_positions = Vec::new();
    for i in 0..bc.len() { if bc[i] == b'N' { umi_positions.push(offset + i); } }
    if umi_positions.is_empty() { return Some(Vec::new()); }
    if *umi_positions.iter().min().unwrap() < offset { return None; }
    let mut umi = Vec::with_capacity(umi_positions.len());
    for &p in umi_positions.iter() { umi.push(seq[p]); }
    Some(umi)
}