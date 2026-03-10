//! Section detection — language-aware structural boundary finder
//!
//! Detects logical sections within source files (functions, classes, modules,
//! headings) for section-level co-occurrence analysis.

/// A detected section within a file
#[derive(Debug, Clone)]
pub struct Section {
    /// Human-readable section name (e.g., "impl Foo", "## Chapter 3")
    pub name: String,
    /// Byte offset of section start within corpus
    pub start_offset: usize,
    /// Byte offset of section end within corpus
    pub end_offset: usize,
    /// Extracted identifier tokens in this section
    pub tokens: Vec<String>,
}

/// Detect sections within file content based on language heuristics
pub fn detect_sections(content: &str, path: &str) -> Vec<Section> {
    let detector = pick_detector(path);
    let raw_sections = detector(content);

    // If no sections detected, use sliding window fallback
    if raw_sections.is_empty() {
        return sliding_window_sections(content, 50, 25);
    }

    raw_sections
}

/// Pick the right section detector based on file extension
fn pick_detector(path: &str) -> fn(&str) -> Vec<Section> {
    let ext = path.rsplit('.').next().unwrap_or("");
    match ext {
        "rs" => detect_rust_sections,
        "py" => detect_python_sections,
        "erl" | "hrl" => detect_erlang_sections,
        "md" | "markdown" => detect_markdown_sections,
        "html" | "htm" => detect_html_sections,
        _ => detect_plaintext_sections,
    }
}

/// Extract identifier tokens from a text block
pub fn extract_tokens(text: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut chars = text.chars().peekable();

    while let Some(&c) = chars.peek() {
        if c.is_alphabetic() || c == '_' {
            let mut ident = String::new();
            while let Some(&ch) = chars.peek() {
                if ch.is_alphanumeric() || ch == '_' {
                    ident.push(ch);
                    chars.next();
                } else {
                    break;
                }
            }
            if ident.len() >= 2 {
                // Split camelCase
                let subs = split_camel_case(&ident);
                for s in &subs {
                    let lower = s.to_lowercase();
                    if lower.len() >= 2 {
                        tokens.push(lower);
                    }
                }
                // Also keep the original identifier
                let lower = ident.to_lowercase();
                if !tokens.contains(&lower) {
                    tokens.push(lower);
                }
            }
        } else {
            chars.next();
        }
    }

    tokens.sort();
    tokens.dedup();
    tokens
}

/// Split camelCase and PascalCase into subwords
fn split_camel_case(s: &str) -> Vec<String> {
    let mut parts = Vec::new();
    let mut current = String::new();

    for c in s.chars() {
        if c == '_' {
            if !current.is_empty() {
                parts.push(std::mem::take(&mut current));
            }
        } else if c.is_uppercase() && !current.is_empty() {
            parts.push(std::mem::take(&mut current));
            current.push(c);
        } else {
            current.push(c);
        }
    }
    if !current.is_empty() {
        parts.push(current);
    }
    parts
}

// ---- Language-specific detectors ----

/// Rust: mod, impl, fn at indent 0, // --- dividers
fn detect_rust_sections(content: &str) -> Vec<Section> {
    detect_by_line_patterns(content, &|line: &str| {
        let trimmed = line.trim();
        if trimmed.starts_with("mod ") || trimmed.starts_with("pub mod ") {
            return Some(trimmed.to_string());
        }
        if trimmed.starts_with("impl ") || trimmed.starts_with("impl<") {
            return Some(trimmed.split('{').next().unwrap_or(trimmed).trim().to_string());
        }
        if (trimmed.starts_with("fn ") || trimmed.starts_with("pub fn ")
            || trimmed.starts_with("pub(crate) fn "))
            && !line.starts_with(|c: char| c == ' ' || c == '\t')
        {
            return Some(trimmed.split('{').next().unwrap_or(trimmed).trim().to_string());
        }
        if trimmed.starts_with("// ---") || trimmed.starts_with("// ===") {
            return Some(trimmed.to_string());
        }
        None
    })
}

/// Python: class, def at indent 0, # --- dividers
fn detect_python_sections(content: &str) -> Vec<Section> {
    detect_by_line_patterns(content, &|line: &str| {
        let trimmed = line.trim();
        if (trimmed.starts_with("class ") || trimmed.starts_with("def ")
            || trimmed.starts_with("async def "))
            && !line.starts_with(|c: char| c == ' ' || c == '\t')
        {
            return Some(trimmed.split(':').next().unwrap_or(trimmed).trim().to_string());
        }
        if trimmed.starts_with("# ---") || trimmed.starts_with("# ===") {
            return Some(trimmed.to_string());
        }
        None
    })
}

/// Erlang: -module, -export, function definitions, %%--- dividers
fn detect_erlang_sections(content: &str) -> Vec<Section> {
    detect_by_line_patterns(content, &|line: &str| {
        let trimmed = line.trim();
        if trimmed.starts_with("-module(") {
            return Some(trimmed.to_string());
        }
        if trimmed.starts_with("%%---") || trimmed.starts_with("%% ---") {
            return Some(trimmed.to_string());
        }
        // Top-level function: starts with lowercase letter, not indented
        if !line.starts_with(|c: char| c == ' ' || c == '\t' || c == '%')
            && trimmed.len() > 1
            && trimmed.starts_with(|c: char| c.is_lowercase())
            && trimmed.contains('(')
        {
            let name = trimmed.split('(').next().unwrap_or(trimmed).trim();
            if !name.is_empty() && name.chars().all(|c| c.is_alphanumeric() || c == '_') {
                return Some(format!("{}()", name));
            }
        }
        None
    })
}

/// Markdown: # headings
fn detect_markdown_sections(content: &str) -> Vec<Section> {
    detect_by_line_patterns(content, &|line: &str| {
        let trimmed = line.trim();
        if trimmed.starts_with('#') {
            return Some(trimmed.to_string());
        }
        None
    })
}

/// HTML: h1-h3 tags
fn detect_html_sections(content: &str) -> Vec<Section> {
    detect_by_line_patterns(content, &|line: &str| {
        let trimmed = line.trim().to_lowercase();
        for tag in &["<h1", "<h2", "<h3"] {
            if trimmed.starts_with(tag) {
                // Extract text content
                let text = line
                    .trim()
                    .split('>')
                    .nth(1)
                    .unwrap_or("")
                    .split('<')
                    .next()
                    .unwrap_or("")
                    .trim();
                if !text.is_empty() {
                    return Some(text.to_string());
                }
                return Some(line.trim().to_string());
            }
        }
        None
    })
}

/// Plain text: blank-line separated paragraphs, ALL-CAPS headings
fn detect_plaintext_sections(content: &str) -> Vec<Section> {
    detect_by_line_patterns(content, &|line: &str| {
        let trimmed = line.trim();
        // ALL-CAPS heading (at least 3 chars, all uppercase letters/spaces)
        if trimmed.len() >= 3
            && trimmed.chars().all(|c| c.is_uppercase() || c == ' ' || c == '-')
            && trimmed.chars().any(|c| c.is_uppercase())
        {
            return Some(trimmed.to_string());
        }
        None
    })
}

/// Generic line-pattern based section detector
fn detect_by_line_patterns(
    content: &str,
    is_boundary: &dyn Fn(&str) -> Option<String>,
) -> Vec<Section> {
    let mut sections = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_start = 0;
    let mut offset = 0;

    for line in content.split('\n') {
        if let Some(name) = is_boundary(line) {
            // Close previous section
            if let Some(prev_name) = current_name.take() {
                let section_text = &content[current_start..offset];
                sections.push(Section {
                    name: prev_name,
                    start_offset: current_start,
                    end_offset: offset,
                    tokens: extract_tokens(section_text),
                });
            }
            current_name = Some(name);
            current_start = offset;
        }
        offset += line.len() + 1; // +1 for the newline
    }

    // Close final section
    if let Some(name) = current_name {
        let section_text = &content[current_start..content.len()];
        sections.push(Section {
            name: name,
            start_offset: current_start,
            end_offset: content.len(),
            tokens: extract_tokens(section_text),
        });
    }

    sections
}

/// Sliding window fallback when no structural markers found
fn sliding_window_sections(content: &str, window_lines: usize, overlap_lines: usize) -> Vec<Section> {
    let lines: Vec<&str> = content.split('\n').collect();
    let step = window_lines.saturating_sub(overlap_lines).max(1);
    let mut sections = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let end = (i + window_lines).min(lines.len());
        let window_text: String = lines[i..end].join("\n");

        // Calculate byte offsets
        let start_offset: usize = lines[..i].iter().map(|l| l.len() + 1).sum();
        let end_offset: usize = lines[..end].iter().map(|l| l.len() + 1).sum();

        sections.push(Section {
            name: format!("lines {}-{}", i + 1, end),
            start_offset,
            end_offset: end_offset.min(content.len()),
            tokens: extract_tokens(&window_text),
        });

        i += step;
    }

    sections
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_camel_case() {
        assert_eq!(split_camel_case("camelCase"), vec!["camel", "Case"]);
        assert_eq!(split_camel_case("PascalCase"), vec!["Pascal", "Case"]);
        assert_eq!(split_camel_case("snake_case"), vec!["snake", "case"]);
        assert_eq!(split_camel_case("simple"), vec!["simple"]);
        assert_eq!(
            split_camel_case("XMLParser"),
            vec!["X", "M", "L", "Parser"]
        );
    }

    #[test]
    fn test_extract_tokens() {
        let tokens = extract_tokens("fn handle_message(msg: &Message) -> Result<()>");
        assert!(tokens.contains(&"handle".to_string()));
        assert!(tokens.contains(&"message".to_string()));
        assert!(tokens.contains(&"handle_message".to_string()));
        assert!(tokens.contains(&"msg".to_string()));
        assert!(tokens.contains(&"result".to_string()));
    }

    #[test]
    fn test_detect_rust_sections() {
        let code = r#"mod networking {
    fn connect() {}
}

impl Server {
    pub fn start(&self) {}
}

fn main() {
    println!("hello");
}
"#;
        let sections = detect_rust_sections(code);
        assert!(sections.len() >= 3);
        assert!(sections[0].name.contains("mod networking"));
    }

    #[test]
    fn test_detect_erlang_sections() {
        let code = r#"-module(rabbit_channel).
-export([start/1]).

%%--- Public API ---

start(Config) ->
    ok.

handle_call(Req, State) ->
    {reply, ok, State}.
"#;
        let sections = detect_erlang_sections(code);
        assert!(sections.len() >= 2);
    }

    #[test]
    fn test_detect_markdown_sections() {
        let md = "# Introduction\nSome text\n## Methods\nMore text\n### Details\nFine text\n";
        let sections = detect_markdown_sections(md);
        assert_eq!(sections.len(), 3);
        assert!(sections[0].name.contains("Introduction"));
    }

    #[test]
    fn test_sliding_window_fallback() {
        let text = (0..100).map(|i| format!("line {}", i)).collect::<Vec<_>>().join("\n");
        let sections = sliding_window_sections(&text, 50, 25);
        assert!(sections.len() >= 3);
        // Each section should have tokens
        assert!(!sections[0].tokens.is_empty());
    }

    #[test]
    fn test_detect_sections_fallback() {
        // Unrecognized extension should still produce sections via sliding window
        let text = "some random content\nwith multiple lines\nand stuff\n".repeat(20);
        let sections = detect_sections(&text, "data.xyz");
        assert!(!sections.is_empty());
    }
}
