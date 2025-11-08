#[cfg(test)]
mod tests {
    use assert_cmd::cmd::Command;
    use std::fs;
    use std::path::Path;
    use tempfile::NamedTempFile;

    fn count_sam_records(path: &Path) -> usize {
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('@'))
            .count()
    }

    fn get_tag_value(line: &str, tag: &str) -> Option<String> {
        for field in line.split('\t') {
            if field.starts_with(&format!("{}:", tag)) {
                return Some(field.to_string());
            }
        }
        None
    }

    #[test]
    fn test_reverse_qt_tag() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output_path = output.path();

        Command::cargo_bin(env!("CARGO_PKG_NAME"))?
            .arg("--input")
            .arg("tests/input.sam")
            .arg("--output")
            .arg(&output_path)
            .arg("--rev")
            .arg("QT")
            .assert()
            .success();

        assert!(output_path.exists());
        assert_eq!(count_sam_records(output_path), 4);

        let content = fs::read_to_string(output_path)?;
        let lines: Vec<&str> = content.lines().collect();

        let read2_line = lines.iter().find(|l| l.starts_with("read2\t")).unwrap();
        let qt_tag = get_tag_value(read2_line, "QT").unwrap();
        assert_eq!(qt_tag, "QT:Z:EFGH");

        let read4_line = lines.iter().find(|l| l.starts_with("read4\t")).unwrap();
        let qt_tag = get_tag_value(read4_line, "QT").unwrap();
        assert_eq!(qt_tag, "QT:Z:IHGF");

        Ok(())
    }

    #[test]
    fn test_reverse_bc_tag() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output_path = output.path();

        Command::cargo_bin(env!("CARGO_PKG_NAME"))?
            .arg("--input")
            .arg("tests/input.sam")
            .arg("--output")
            .arg(&output_path)
            .arg("--rev")
            .arg("BC")
            .assert()
            .success();

        let content = fs::read_to_string(output_path)?;
        let lines: Vec<&str> = content.lines().collect();

        // read2 should have BC reversed from GGAT to TAGG
        let read2_line = lines.iter().find(|l| l.starts_with("read2\t")).unwrap();
        let bc_tag = get_tag_value(read2_line, "BC").unwrap();
        assert_eq!(bc_tag, "BC:Z:TAGG");

        // read4 should have BC reversed from TCGA to AGCT
        let read4_line = lines.iter().find(|l| l.starts_with("read4\t")).unwrap();
        let bc_tag = get_tag_value(read4_line, "BC").unwrap();
        assert_eq!(bc_tag, "BC:Z:AGCT");

        Ok(())
    }

    #[test]
    fn test_revcomp_bc_tag() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output_path = output.path();

        Command::cargo_bin(env!("CARGO_PKG_NAME"))?
            .arg("--input")
            .arg("tests/input.sam")
            .arg("--output")
            .arg(&output_path)
            .arg("--revcomp")
            .arg("BC")
            .assert()
            .success();

        let content = fs::read_to_string(output_path)?;
        let lines: Vec<&str> = content.lines().collect();

        // read2 should have BC reverse complemented from GGAT to ATCC
        let read2_line = lines.iter().find(|l| l.starts_with("read2\t")).unwrap();
        let bc_tag = get_tag_value(read2_line, "BC").unwrap();
        assert_eq!(bc_tag, "BC:Z:ATCC");

        // read4 should have BC reverse complemented from TCGA to TCGA (palindrome!)
        let read4_line = lines.iter().find(|l| l.starts_with("read4\t")).unwrap();
        let bc_tag = get_tag_value(read4_line, "BC").unwrap();
        assert_eq!(bc_tag, "BC:Z:TCGA");

        Ok(())
    }

    #[test]
    fn test_multiple_tags() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output_path = output.path();

        Command::cargo_bin(env!("CARGO_PKG_NAME"))?
            .arg("--input")
            .arg("tests/input.sam")
            .arg("--output")
            .arg(&output_path)
            .arg("--rev")
            .arg("QT")
            .arg("--rev")
            .arg("BC")
            .assert()
            .success();

        let content = fs::read_to_string(output_path)?;
        let lines: Vec<&str> = content.lines().collect();
        let read2_line = lines.iter().find(|l| l.starts_with("read2\t")).unwrap();
        let qt_tag = get_tag_value(read2_line, "QT").unwrap();
        let bc_tag = get_tag_value(read2_line, "BC").unwrap();
        assert_eq!(qt_tag, "QT:Z:EFGH");
        assert_eq!(bc_tag, "BC:Z:TAGG");

        Ok(())
    }

    #[test]
    fn test_stdin_stdout() -> Result<(), Box<dyn std::error::Error>> {
        Command::cargo_bin(env!("CARGO_PKG_NAME"))?
            .arg("--rev")
            .arg("QT")
            .pipe_stdin("tests/input.sam")?
            .assert()
            .success();

        Ok(())
    }

    #[test]
    fn test_stdin_stdout_with_dash() -> Result<(), Box<dyn std::error::Error>> {
        Command::cargo_bin(env!("CARGO_PKG_NAME"))?
            .arg("--input")
            .arg("-")
            .arg("--output")
            .arg("-")
            .arg("--rev")
            .arg("QT")
            .pipe_stdin("tests/input.sam")?
            .assert()
            .success();

        Ok(())
    }

    #[test]
    fn test_no_tags_specified() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output_path = output.path();

        Command::cargo_bin(env!("CARGO_PKG_NAME"))?
            .arg("--input")
            .arg("tests/input.sam")
            .arg("--output")
            .arg(&output_path)
            .assert()
            .success();

        assert!(output_path.exists());
        assert_eq!(count_sam_records(output_path), 4);

        Ok(())
    }

    #[test]
    fn test_forward_reads_unchanged() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output_path = output.path();

        Command::cargo_bin(env!("CARGO_PKG_NAME"))?
            .arg("--input")
            .arg("tests/input.sam")
            .arg("--output")
            .arg(&output_path)
            .arg("--rev")
            .arg("QT")
            .assert()
            .success();

        let content = fs::read_to_string(output_path)?;
        let lines: Vec<&str> = content.lines().collect();

        // read1 and read3 are forward strand, should be unchanged
        let read1_line = lines.iter().find(|l| l.starts_with("read1\t")).unwrap();
        let qt_tag = get_tag_value(read1_line, "QT").unwrap();
        assert_eq!(qt_tag, "QT:Z:IJKL"); // Unchanged

        Ok(())
    }
}
