//! A library for reversing and complementing array-like SAM tags for negative
//! facing alignments.
#![warn(missing_docs)]

use anyhow::Result;
use bio::alphabets::dna;
use log::*;
use proglog::ProgLogBuilder;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{Header, Read as BamRead, Reader, Record, Writer};
use std::error;
use std::path::PathBuf;

const CARGO_PKG_NAME: &str = env!("CARGO_PKG_NAME");
const CARGO_PKG_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Mutates a record by reversing and/or reverse complementing specified tags.
///
/// This function modifies the record in-place by:
/// - Reversing the order of array-like values in specified `rev` tags
/// - Reverse complementing array-like string values in specified `revcomp` tags
///
/// # Arguments
///
/// * `record` - The BAM record to mutate
/// * `rev` - SAM tags to reverse (e.g., base qualities) as 2-byte arrays
/// * `revcomp` - SAM tags to reverse complement (e.g., sequences) as 2-byte arrays
///
/// # Returns
///
/// Returns Ok(()) on success, or an error if tag manipulation fails.
///
fn reverse_tags_for(
    record: &mut Record,
    rev: &[[u8; 2]],
    revcomp: &[[u8; 2]],
) -> Result<(), Box<dyn error::Error>> {
    macro_rules! try_reverse_array {
        ($tag:expr, $variant:ident, $ty:ty) => {
            if let Ok(rust_htslib::bam::record::Aux::$variant(arr)) = record.aux($tag) {
                let mut values: Vec<$ty> = arr.iter().collect();
                values.reverse();
                record.remove_aux($tag)?;
                record.push_aux(
                    $tag,
                    rust_htslib::bam::record::Aux::$variant((&values[..]).into()),
                )?;
                continue;
            }
        };
    }

    for tag in rev {
        try_reverse_array!(tag, ArrayU8, u8);
        try_reverse_array!(tag, ArrayU16, u16);
        try_reverse_array!(tag, ArrayU32, u32);
        try_reverse_array!(tag, ArrayI8, i8);
        try_reverse_array!(tag, ArrayI16, i16);
        try_reverse_array!(tag, ArrayI32, i32);
        try_reverse_array!(tag, ArrayFloat, f32);

        if let Ok(rust_htslib::bam::record::Aux::String(s)) = record.aux(tag) {
            let reversed: String = s.chars().rev().collect();
            record.remove_aux(tag)?;
            record.push_aux(tag, rust_htslib::bam::record::Aux::String(&reversed))?;
        }
    }

    for tag in revcomp {
        if let Ok(rust_htslib::bam::record::Aux::String(s)) = record.aux(tag) {
            let revcomp_seq = dna::revcomp(s.as_bytes());
            let revcomp_str = String::from_utf8_lossy(&revcomp_seq).to_string();
            record.remove_aux(tag)?;
            record.push_aux(tag, rust_htslib::bam::record::Aux::String(&revcomp_str))?;
        } else if let Ok(rust_htslib::bam::record::Aux::ArrayU8(arr)) = record.aux(tag) {
            let values: Vec<u8> = arr.iter().collect();
            let revcomp_seq = dna::revcomp(&values);
            record.remove_aux(tag)?;
            record.push_aux(
                tag,
                rust_htslib::bam::record::Aux::ArrayU8((&revcomp_seq[..]).into()),
            )?;
        }
    }

    Ok(())
}

/// Validates and converts tag names to byte arrays.
///
/// # Arguments
///
/// * `tags` - Vector of tag name strings to validate and convert
///
/// # Returns
///
/// Returns a vector of 2-byte arrays on success, or an error if any tag name is not exactly 2 characters.
///
fn validate_tags(tags: &[String]) -> Result<Vec<[u8; 2]>, Box<dyn error::Error>> {
    let mut result = Vec::with_capacity(tags.len());
    for tag_name in tags {
        if tag_name.len() != 2 {
            return Err(format!("Tag name must be exactly 2 characters: {tag_name}").into());
        }
        let bytes = tag_name.as_bytes();
        result.push([bytes[0], bytes[1]]);
    }
    Ok(result)
}

/// Runs the tool `revtag` on an input SAM/BAM/CRAM file and writes the records to an output file.
///
/// For reverse strand alignments (flag 0x10 set), this function will:
/// - Reverse the order of array-like values in specified tags (--rev)
/// - Reverse complement array-like string values in specified tags (--revcomp)
///
/// # Arguments
///
/// * `input` - The input SAM/BAM/CRAM file path, or None/Some("-") for stdin
/// * `output` - The output SAM/BAM/CRAM file path, or None/Some("-") for stdout
/// * `rev` - SAM tags to reverse (e.g., base qualities)
/// * `revcomp` - SAM tags to reverse complement (e.g., sequences)
/// * `threads` - Extra threads for BAM/CRAM compression/decompression
///
/// # Returns
///
/// Returns the result of the execution with an integer exit code for success (0).
///
pub fn revtag(
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    rev: Vec<String>,
    revcomp: Vec<String>,
    threads: usize,
) -> Result<i32, Box<dyn error::Error>> {
    let rev_tags = validate_tags(&rev)?;
    let revcomp_tags = validate_tags(&revcomp)?;

    let mut reader = match &input {
        None => {
            info!("Input: stdin");
            Reader::from_stdin()?
        }
        Some(path) => {
            info!("Input: {path:?}");
            Reader::from_path(path)?
        }
    };

    if threads > 1 {
        reader.set_threads(threads - 1)?;
    }

    let mut header = Header::from_template(reader.header());

    header.push_record(
        HeaderRecord::new(b"PG")
            .push_tag(b"ID", CARGO_PKG_NAME)
            .push_tag(b"PN", CARGO_PKG_NAME)
            .push_tag(b"VN", CARGO_PKG_VERSION)
            .push_tag(b"CL", std::env::args().collect::<Vec<_>>().join(" ")),
    );

    let mut writer = match &output {
        None => {
            info!("Output: stdout");
            Writer::from_stdout(&header, rust_htslib::bam::Format::Sam)?
        }
        Some(path) => {
            info!("Output: {path:?}");
            let format = if path.to_str().map(|s| s.ends_with(".bam")).unwrap_or(false) {
                rust_htslib::bam::Format::Bam
            } else if path.to_str().map(|s| s.ends_with(".cram")).unwrap_or(false) {
                rust_htslib::bam::Format::Cram
            } else {
                rust_htslib::bam::Format::Sam
            };
            Writer::from_path(path, &header, format)?
        }
    };

    if threads > 1 {
        writer.set_threads(threads - 1)?;
    }

    let progress = ProgLogBuilder::new()
        .name("main")
        .verb("Processed")
        .noun("alignment records")
        .unit(100_000)
        .build();

    let mut record = Record::new();

    loop {
        match reader.read(&mut record) {
            Some(Ok(())) => {}
            None => break,
            Some(Err(e)) => return Err(Box::new(e)),
        }

        if record.is_reverse() {
            reverse_tags_for(&mut record, &rev_tags, &revcomp_tags)?;
        }

        writer.write(&record)?;
        progress.record();
    }

    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Aux;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Helper to create a minimal BAM record for testing
    fn create_test_record() -> Record {
        let mut record = Record::new();
        record.set_qname(b"test_read");
        record.set_pos(100);
        record.set_mapq(60);
        record.set_tid(0);
        record
    }

    /// Helper to convert tag strings to byte arrays for testing
    fn tags_to_bytes(tags: &[&str]) -> Vec<[u8; 2]> {
        tags.iter()
            .filter_map(|s| {
                if s.len() >= 2 {
                    let b = s.as_bytes();
                    Some([b[0], b[1]])
                } else {
                    None
                }
            })
            .collect()
    }

    #[test]
    fn test_reverse_u8_array() {
        let mut record = create_test_record();
        let values = vec![10u8, 20, 30, 40, 50];
        record
            .push_aux(b"QT", Aux::ArrayU8((&values[..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["QT"]), &[]).unwrap();

        if let Ok(Aux::ArrayU8(arr)) = record.aux(b"QT") {
            let result: Vec<u8> = arr.iter().collect();
            assert_eq!(result, vec![50, 40, 30, 20, 10]);
        } else {
            panic!("Expected ArrayU8");
        }
    }

    #[test]
    fn test_reverse_u16_array() {
        let mut record = create_test_record();
        let values = vec![100u16, 200, 300];
        record
            .push_aux(b"AB", Aux::ArrayU16((&values[..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["AB"]), &[]).unwrap();

        if let Ok(Aux::ArrayU16(arr)) = record.aux(b"AB") {
            let result: Vec<u16> = arr.iter().collect();
            assert_eq!(result, vec![300, 200, 100]);
        } else {
            panic!("Expected ArrayU16");
        }
    }

    #[test]
    fn test_reverse_u32_array() {
        let mut record = create_test_record();
        let values = vec![1000u32, 2000, 3000];
        record
            .push_aux(b"CD", Aux::ArrayU32((&values[..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["CD"]), &[]).unwrap();

        if let Ok(Aux::ArrayU32(arr)) = record.aux(b"CD") {
            let result: Vec<u32> = arr.iter().collect();
            assert_eq!(result, vec![3000, 2000, 1000]);
        } else {
            panic!("Expected ArrayU32");
        }
    }

    #[test]
    fn test_reverse_i8_array() {
        let mut record = create_test_record();
        let values = vec![-10i8, -5, 0, 5, 10];
        record
            .push_aux(b"EF", Aux::ArrayI8((&values[..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["EF"]), &[]).unwrap();

        if let Ok(Aux::ArrayI8(arr)) = record.aux(b"EF") {
            let result: Vec<i8> = arr.iter().collect();
            assert_eq!(result, vec![10, 5, 0, -5, -10]);
        } else {
            panic!("Expected ArrayI8");
        }
    }

    #[test]
    fn test_reverse_i16_array() {
        let mut record = create_test_record();
        let values = vec![-100i16, 0, 100];
        record
            .push_aux(b"GH", Aux::ArrayI16((&values[..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["GH"]), &[]).unwrap();

        if let Ok(Aux::ArrayI16(arr)) = record.aux(b"GH") {
            let result: Vec<i16> = arr.iter().collect();
            assert_eq!(result, vec![100, 0, -100]);
        } else {
            panic!("Expected ArrayI16");
        }
    }

    #[test]
    fn test_reverse_i32_array() {
        let mut record = create_test_record();
        let values = vec![-1000i32, 0, 1000];
        record
            .push_aux(b"IJ", Aux::ArrayI32((&values[..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["IJ"]), &[]).unwrap();

        if let Ok(Aux::ArrayI32(arr)) = record.aux(b"IJ") {
            let result: Vec<i32> = arr.iter().collect();
            assert_eq!(result, vec![1000, 0, -1000]);
        } else {
            panic!("Expected ArrayI32");
        }
    }

    #[test]
    fn test_reverse_float_array() {
        let mut record = create_test_record();
        let values = vec![1.5f32, 2.5, 3.5];
        record
            .push_aux(b"KL", Aux::ArrayFloat((&values[..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["KL"]), &[]).unwrap();

        if let Ok(Aux::ArrayFloat(arr)) = record.aux(b"KL") {
            let result: Vec<f32> = arr.iter().collect();
            assert_eq!(result, vec![3.5, 2.5, 1.5]);
        } else {
            panic!("Expected ArrayFloat");
        }
    }

    #[test]
    fn test_reverse_string() {
        let mut record = create_test_record();
        record.push_aux(b"MN", Aux::String("HELLO")).unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["MN"]), &[]).unwrap();

        if let Ok(Aux::String(s)) = record.aux(b"MN") {
            assert_eq!(s, "OLLEH");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_string_uppercase() {
        let mut record = create_test_record();
        record.push_aux(b"BC", Aux::String("ATCG")).unwrap();

        reverse_tags_for(&mut record, &[], &tags_to_bytes(&["BC"])).unwrap();

        if let Ok(Aux::String(s)) = record.aux(b"BC") {
            assert_eq!(s, "CGAT");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_string_lowercase() {
        let mut record = create_test_record();
        record.push_aux(b"BC", Aux::String("atcg")).unwrap();

        reverse_tags_for(&mut record, &[], &tags_to_bytes(&["BC"])).unwrap();

        if let Ok(Aux::String(s)) = record.aux(b"BC") {
            assert_eq!(s, "cgat");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_string_mixed_case() {
        let mut record = create_test_record();
        record.push_aux(b"BC", Aux::String("AtCg")).unwrap();

        reverse_tags_for(&mut record, &[], &tags_to_bytes(&["BC"])).unwrap();

        if let Ok(Aux::String(s)) = record.aux(b"BC") {
            assert_eq!(s, "cGaT");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_string_palindrome() {
        let mut record = create_test_record();
        record.push_aux(b"BC", Aux::String("TCGA")).unwrap();

        reverse_tags_for(&mut record, &[], &tags_to_bytes(&["BC"])).unwrap();

        if let Ok(Aux::String(s)) = record.aux(b"BC") {
            assert_eq!(s, "TCGA");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_revcomp_array_u8() {
        let mut record = create_test_record();
        let seq = b"ATCG";
        let values: Vec<u8> = seq.to_vec();
        record
            .push_aux(b"BC", Aux::ArrayU8((&values[..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &[], &tags_to_bytes(&["BC"])).unwrap();

        if let Ok(Aux::ArrayU8(arr)) = record.aux(b"BC") {
            let result: Vec<u8> = arr.iter().collect();
            assert_eq!(result, b"CGAT".to_vec());
        } else {
            panic!("Expected ArrayU8");
        }
    }

    #[test]
    fn test_revcomp_longer_sequence() {
        let mut record = create_test_record();
        record.push_aux(b"BC", Aux::String("ATCGATCGATCG")).unwrap();

        reverse_tags_for(&mut record, &[], &tags_to_bytes(&["BC"])).unwrap();

        if let Ok(Aux::String(s)) = record.aux(b"BC") {
            assert_eq!(s, "CGATCGATCGAT");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_multiple_rev_tags() {
        let mut record = create_test_record();
        record
            .push_aux(b"QT", Aux::ArrayU8((&[10u8, 20, 30][..]).into()))
            .unwrap();
        record
            .push_aux(b"AB", Aux::ArrayU8((&[1u8, 2, 3][..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["QT", "AB"]), &[]).unwrap();

        if let Ok(Aux::ArrayU8(arr)) = record.aux(b"QT") {
            let result: Vec<u8> = arr.iter().collect();
            assert_eq!(result, vec![30, 20, 10]);
        } else {
            panic!("Expected ArrayU8 for QT");
        }

        if let Ok(Aux::ArrayU8(arr)) = record.aux(b"AB") {
            let result: Vec<u8> = arr.iter().collect();
            assert_eq!(result, vec![3, 2, 1]);
        } else {
            panic!("Expected ArrayU8 for AB");
        }
    }

    #[test]
    fn test_multiple_revcomp_tags() {
        let mut record = create_test_record();
        record.push_aux(b"BC", Aux::String("ATCG")).unwrap();
        record.push_aux(b"XY", Aux::String("GGCC")).unwrap();

        reverse_tags_for(&mut record, &[], &tags_to_bytes(&["BC", "XY"])).unwrap();

        if let Ok(Aux::String(s)) = record.aux(b"BC") {
            assert_eq!(s, "CGAT");
        } else {
            panic!("Expected String for BC");
        }

        if let Ok(Aux::String(s)) = record.aux(b"XY") {
            assert_eq!(s, "GGCC");
        } else {
            panic!("Expected String for XY");
        }
    }

    #[test]
    fn test_both_rev_and_revcomp() {
        let mut record = create_test_record();
        record
            .push_aux(b"QT", Aux::ArrayU8((&[10u8, 20, 30][..]).into()))
            .unwrap();
        record.push_aux(b"BC", Aux::String("ATCG")).unwrap();

        reverse_tags_for(
            &mut record,
            &tags_to_bytes(&["QT"]),
            &tags_to_bytes(&["BC"]),
        )
        .unwrap();

        if let Ok(Aux::ArrayU8(arr)) = record.aux(b"QT") {
            let result: Vec<u8> = arr.iter().collect();
            assert_eq!(result, vec![30, 20, 10]);
        } else {
            panic!("Expected ArrayU8 for QT");
        }

        if let Ok(Aux::String(s)) = record.aux(b"BC") {
            assert_eq!(s, "CGAT");
        } else {
            panic!("Expected String for BC");
        }
    }

    #[test]
    fn test_nonexistent_tag() {
        let mut record = create_test_record();
        let result = reverse_tags_for(&mut record, &tags_to_bytes(&["ZZ"]), &[]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_empty_arrays() {
        let mut record = create_test_record();
        record
            .push_aux(b"QT", Aux::ArrayU8((&[][..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["QT"]), &[]).unwrap();

        if let Ok(Aux::ArrayU8(arr)) = record.aux(b"QT") {
            let result: Vec<u8> = arr.iter().collect();
            assert_eq!(result, Vec::<u8>::new());
        } else {
            panic!("Expected ArrayU8");
        }
    }

    #[test]
    fn test_single_element_array() {
        let mut record = create_test_record();
        record
            .push_aux(b"QT", Aux::ArrayU8((&[42u8][..]).into()))
            .unwrap();

        reverse_tags_for(&mut record, &tags_to_bytes(&["QT"]), &[]).unwrap();

        if let Ok(Aux::ArrayU8(arr)) = record.aux(b"QT") {
            let result: Vec<u8> = arr.iter().collect();
            assert_eq!(result, vec![42]);
        } else {
            panic!("Expected ArrayU8");
        }
    }

    #[test]
    fn test_empty_string() {
        let mut record = create_test_record();
        record.push_aux(b"BC", Aux::String("")).unwrap();

        reverse_tags_for(&mut record, &[], &tags_to_bytes(&["BC"])).unwrap();

        if let Ok(Aux::String(s)) = record.aux(b"BC") {
            assert_eq!(s, "");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_single_nucleotide() {
        let mut record = create_test_record();
        record.push_aux(b"BC", Aux::String("A")).unwrap();

        reverse_tags_for(&mut record, &[], &tags_to_bytes(&["BC"])).unwrap();

        if let Ok(Aux::String(s)) = record.aux(b"BC") {
            assert_eq!(s, "T");
        } else {
            panic!("Expected String");
        }
    }

    #[test]
    fn test_validate_tags_happy_path() {
        let tags = vec!["QT".to_string(), "BC".to_string()];
        let validated = validate_tags(&tags).expect("expected Ok for valid tags");
        assert_eq!(validated, vec!([b'Q', b'T'], [b'B', b'C']));
    }

    #[test]
    fn test_validate_tags_unhappy_path_invalid_length() {
        let tags = vec!["Q".to_string(), "ABC".to_string(), "BC".to_string()];
        let err = validate_tags(&tags).expect_err("expected Err for invalid tag lengths");
        let msg = err.to_string();
        assert!(
            msg.contains("exactly 2 characters"),
            "unexpected error: {}",
            msg
        );
    }

    // Helper to build a minimal SAM header
    fn sam_header() -> &'static str {
        "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:chr1\tLN:1000\n"
    }

    // Helper to create two SAM alignment lines: one forward and one reverse with tags
    fn sam_body_with_tags() -> String {
        // Fields: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL TAGS...
        // Reverse read flag: 16 (0x10), Forward read flag: 0
        // rev tags: QT (ArrayU8), MN (String)
        // revcomp tags: BC (String), SQ (ArrayU8 simulated sequence)
        let forward = "fwd\t0\tchr1\t1\t60\t10M\t*\t0\t0\tATCGATCGAA\tFFFFFFFFFF\tQT:B:C,10,20,30\tMN:Z:HELLO\tBC:Z:ATCG\tSQ:B:C,65,84,67,71"; // A T C G
        let reverse = "rev\t16\tchr1\t2\t60\t10M\t*\t0\t0\tGGCATCGTTA\tFFFFFFFFFF\tQT:B:C,1,2,3\tMN:Z:WORLD\tBC:Z:GATT\tSQ:B:C,71,65,84,84"; // G A T T
        format!("{forward}\n{reverse}\n")
    }

    // Parse SAM output produced by revtag, returning vector of (qname, tags_map)
    fn parse_sam_tags(sam: &str) -> Vec<(String, std::collections::HashMap<String, String>)> {
        sam.lines()
            .filter(|l| !l.starts_with('@') && !l.is_empty())
            .map(|line| {
                let parts: Vec<&str> = line.split('\t').collect();
                let qname = parts[0].to_string();
                let mut tags = std::collections::HashMap::new();
                for field in parts.iter().skip(11) {
                    if field.len() >= 5
                        && field.chars().nth(2) == Some(':')
                        && field.chars().nth(4) == Some(':')
                    {
                        // TAG:TYPE:VALUE simplistic parse
                        let tag = &field[0..2];
                        let value = &field[5..];
                        tags.insert(tag.to_string(), value.to_string());
                    } else if field.len() >= 5 && field.chars().nth(2) == Some(':') {
                        // TAG:B:TYPE,LIST  (we'll just store raw)
                        let tag = &field[0..2];
                        let value = &field[5..];
                        tags.insert(tag.to_string(), value.to_string());
                    }
                }
                (qname, tags)
            })
            .collect()
    }

    #[test]
    fn test_revtag_happy_path_forward_vs_reverse() {
        let mut infile = NamedTempFile::new().expect("temp sam input");
        write!(infile, "{}{}", sam_header(), sam_body_with_tags()).unwrap();

        let outfile = NamedTempFile::new().expect("temp sam output");

        let exit = revtag(
            Some(infile.path().to_path_buf()),
            Some(outfile.path().to_path_buf()),
            vec!["QT".into(), "MN".into()], // rev
            vec!["BC".into(), "SQ".into()], // revcomp
            1,
        )
        .expect("revtag should succeed");
        assert_eq!(exit, 0);

        let output_contents = std::fs::read_to_string(outfile.path()).unwrap();
        let parsed = parse_sam_tags(&output_contents);
        assert_eq!(parsed.len(), 2);

        // Find forward and reverse records
        let mut fwd = None;
        let mut rev_rec = None;
        for (q, t) in parsed {
            if q == "fwd" {
                fwd = Some(t);
            } else if q == "rev" {
                rev_rec = Some(t);
            }
        }
        let fwd = fwd.expect("forward record present");
        let rev_rec = rev_rec.expect("reverse record present");

        assert_eq!(fwd.get("MN").unwrap(), "HELLO");
        assert_eq!(fwd.get("BC").unwrap(), "ATCG");
        // MN reversed WORLD -> DLROW
        assert_eq!(rev_rec.get("MN").unwrap(), "DLROW");
        // BC revcomp of GATT -> AATC
        assert_eq!(rev_rec.get("BC").unwrap(), "AATC");
    }

    #[test]
    fn test_revtag_threads_two() {
        let mut infile = NamedTempFile::new().expect("temp sam input");
        write!(infile, "{}{}", sam_header(), sam_body_with_tags()).unwrap();
        let outfile = NamedTempFile::new().expect("temp sam output");

        let exit = revtag(
            Some(infile.path().to_path_buf()),
            Some(outfile.path().to_path_buf()),
            vec!["MN".into()],
            vec!["BC".into()],
            2,
        )
        .expect("revtag should succeed with threads=2");
        assert_eq!(exit, 0);
        let output_contents = std::fs::read_to_string(outfile.path()).unwrap();
        assert!(output_contents.contains("rev"));
        assert!(output_contents.contains("MN:Z:DLROW"));
    }

    #[test]
    fn test_revtag_empty_input() {
        let mut infile = NamedTempFile::new().expect("empty sam input");
        write!(infile, "{}", sam_header()).unwrap(); // header only
        let outfile = NamedTempFile::new().expect("empty sam output");

        let exit = revtag(
            Some(infile.path().to_path_buf()),
            Some(outfile.path().to_path_buf()),
            vec!["QT".into()],
            vec!["BC".into()],
            1,
        )
        .expect("revtag should succeed on empty input");
        assert_eq!(exit, 0);
        let output_contents = std::fs::read_to_string(outfile.path()).unwrap();
        assert!(output_contents.lines().all(|l| l.starts_with('@')));
    }

    #[test]
    fn test_revtag_output_bam_format() {
        let mut infile = NamedTempFile::new().expect("temp sam input");
        write!(infile, "{}{}", sam_header(), sam_body_with_tags()).unwrap();

        let tmpdir = tempfile::tempdir().unwrap();
        let bam_out = tmpdir.path().join("out.bam");

        let exit = revtag(
            Some(infile.path().to_path_buf()),
            Some(bam_out.clone()),
            vec!["QT".into(), "MN".into()],
            vec!["BC".into()],
            1,
        )
        .expect("revtag should succeed producing BAM");
        assert_eq!(exit, 0);

        let mut reader = Reader::from_path(&bam_out).expect("read BAM output");
        let mut rec = Record::new();
        let mut saw_fwd = false;
        let mut saw_rev = false;
        while let Some(Ok(())) = reader.read(&mut rec) {
            let qname = String::from_utf8_lossy(rec.qname()).to_string();
            if qname == "fwd" {
                saw_fwd = true;
                if let Ok(Aux::String(s)) = rec.aux(b"MN") {
                    assert_eq!(s, "HELLO");
                }
            } else if qname == "rev" {
                saw_rev = true;
                if let Ok(Aux::String(s)) = rec.aux(b"MN") {
                    assert_eq!(s, "DLROW");
                }
            }
        }
        assert!(saw_fwd && saw_rev);
    }

    #[test]
    fn test_revtag_output_cram_empty_input() {
        let mut infile = NamedTempFile::new().expect("empty sam input");
        write!(infile, "{}", sam_header()).unwrap();

        let tmpdir = tempfile::tempdir().unwrap();
        let cram_out = tmpdir.path().join("out.cram");

        let exit = revtag(
            Some(infile.path().to_path_buf()),
            Some(cram_out.clone()),
            vec![],
            vec![],
            1,
        )
        .expect("revtag should succeed for empty input CRAM");
        assert_eq!(exit, 0);
        let meta = std::fs::metadata(&cram_out).expect("cram file exists");
        assert!(meta.len() > 0);
    }
}
