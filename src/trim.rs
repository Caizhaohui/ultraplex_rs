pub fn quality_trim_index(qualities: &str, cutoff_front: i32, cutoff_back: i32, base: i32) -> (usize, usize) {
    let mut s = 0i32;
    let mut max_qual = 0i32;
    let mut start: usize = 0;
    let len = qualities.len();
    for (i, &b) in qualities.as_bytes().iter().enumerate() {
        s += cutoff_front - ((b as i32) - base);
        if s < 0 { break; }
        if s > max_qual { max_qual = s; start = i + 1; }
    }

    max_qual = 0;
    s = 0;
    let mut stop: usize = len;
    for i in (0..len).rev() {
        let b = qualities.as_bytes()[i] as i32;
        s += cutoff_back - (b - base);
        if s < 0 { break; }
        if s > max_qual { max_qual = s; stop = i; }
    }
    if start >= stop { return (0, 0); }
    (start, stop)
}

pub fn nextseq_trim_index(sequence: &str, qualities: &str, cutoff: i32, base: i32) -> usize {
    let mut s = 0i32;
    let mut max_qual = 0i32;
    let mut max_i: usize = qualities.len();
    let bases = sequence.as_bytes();
    let quals = qualities.as_bytes();
    for i in (0..quals.len()).rev() {
        let mut q = (quals[i] as i32) - base;
        if bases[i] == b'G' { q = cutoff - 1; }
        s += cutoff - q;
        if s < 0 { break; }
        if s > max_qual { max_qual = s; max_i = i; }
    }
    max_i
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_trim_no_trim() {
        let quals = "IIIIIIII"; // Phred 40 (73 ascii) if base=33
        let (start, stop) = quality_trim_index(quals, 30, 30, 33);
        assert!(start < stop);
    }

    #[test]
    fn test_quality_trim_all_trim() {
        let quals = "!!!!!!!!"; // Phred 0 (33 ascii)
        let (start, stop) = quality_trim_index(quals, 1, 1, 33);
        assert_eq!((start, stop), (0, 0));
    }

    #[test]
    fn test_nextseq_trim_basic() {
        let seq = "ACGTGGGG";
        let quals = "IIIIIIII";
        let idx = nextseq_trim_index(seq, quals, 20, 33);
        assert!(idx <= quals.len());
    }
}