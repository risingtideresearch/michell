//! Minimal strict JSON parser (no dependencies). Supports objects, arrays,
//! strings (with escapes incl. \uXXXX), numbers, booleans, and null.

#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(kv) => kv.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Num(n) => Some(*n),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_arr(&self) -> Option<&[Json]> {
        match self {
            Json::Arr(a) => Some(a),
            _ => None,
        }
    }
}

pub fn parse(text: &str) -> Result<Json, String> {
    let bytes = text.as_bytes();
    let mut at = 0usize;
    let v = value(bytes, &mut at)?;
    skip_ws(bytes, &mut at);
    if at != bytes.len() {
        return Err(err_at(bytes, at, "trailing content after the JSON value"));
    }
    Ok(v)
}

fn err_at(bytes: &[u8], at: usize, msg: &str) -> String {
    let line = bytes[..at.min(bytes.len())]
        .iter()
        .filter(|&&b| b == b'\n')
        .count()
        + 1;
    format!("JSON line {line}: {msg}")
}

fn skip_ws(bytes: &[u8], at: &mut usize) {
    while *at < bytes.len() && matches!(bytes[*at], b' ' | b'\t' | b'\n' | b'\r') {
        *at += 1;
    }
}

fn value(bytes: &[u8], at: &mut usize) -> Result<Json, String> {
    skip_ws(bytes, at);
    match bytes.get(*at) {
        None => Err(err_at(bytes, *at, "unexpected end of input")),
        Some(b'{') => {
            *at += 1;
            let mut out = Vec::new();
            skip_ws(bytes, at);
            if bytes.get(*at) == Some(&b'}') {
                *at += 1;
                return Ok(Json::Obj(out));
            }
            loop {
                skip_ws(bytes, at);
                let Json::Str(key) = value(bytes, at)? else {
                    return Err(err_at(bytes, *at, "object keys must be strings"));
                };
                skip_ws(bytes, at);
                if bytes.get(*at) != Some(&b':') {
                    return Err(err_at(bytes, *at, "expected ':' after object key"));
                }
                *at += 1;
                let v = value(bytes, at)?;
                out.push((key, v));
                skip_ws(bytes, at);
                match bytes.get(*at) {
                    Some(b',') => *at += 1,
                    Some(b'}') => {
                        *at += 1;
                        return Ok(Json::Obj(out));
                    }
                    _ => return Err(err_at(bytes, *at, "expected ',' or '}' in object")),
                }
            }
        }
        Some(b'[') => {
            *at += 1;
            let mut out = Vec::new();
            skip_ws(bytes, at);
            if bytes.get(*at) == Some(&b']') {
                *at += 1;
                return Ok(Json::Arr(out));
            }
            loop {
                out.push(value(bytes, at)?);
                skip_ws(bytes, at);
                match bytes.get(*at) {
                    Some(b',') => *at += 1,
                    Some(b']') => {
                        *at += 1;
                        return Ok(Json::Arr(out));
                    }
                    _ => return Err(err_at(bytes, *at, "expected ',' or ']' in array")),
                }
            }
        }
        Some(b'"') => string(bytes, at).map(Json::Str),
        Some(b't') => literal(bytes, at, "true", Json::Bool(true)),
        Some(b'f') => literal(bytes, at, "false", Json::Bool(false)),
        Some(b'n') => literal(bytes, at, "null", Json::Null),
        Some(_) => number(bytes, at),
    }
}

fn literal(bytes: &[u8], at: &mut usize, word: &str, v: Json) -> Result<Json, String> {
    if bytes[*at..].starts_with(word.as_bytes()) {
        *at += word.len();
        Ok(v)
    } else {
        Err(err_at(bytes, *at, "invalid literal"))
    }
}

fn number(bytes: &[u8], at: &mut usize) -> Result<Json, String> {
    let start = *at;
    while *at < bytes.len()
        && matches!(bytes[*at], b'0'..=b'9' | b'-' | b'+' | b'.' | b'e' | b'E')
    {
        *at += 1;
    }
    std::str::from_utf8(&bytes[start..*at])
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .map(Json::Num)
        .ok_or_else(|| err_at(bytes, start, "invalid number"))
}

fn string(bytes: &[u8], at: &mut usize) -> Result<String, String> {
    *at += 1; // opening quote
    let mut out = String::new();
    loop {
        match bytes.get(*at) {
            None => return Err(err_at(bytes, *at, "unterminated string")),
            Some(b'"') => {
                *at += 1;
                return Ok(out);
            }
            Some(b'\\') => {
                *at += 1;
                match bytes.get(*at) {
                    Some(b'"') => out.push('"'),
                    Some(b'\\') => out.push('\\'),
                    Some(b'/') => out.push('/'),
                    Some(b'b') => out.push('\u{8}'),
                    Some(b'f') => out.push('\u{c}'),
                    Some(b'n') => out.push('\n'),
                    Some(b'r') => out.push('\r'),
                    Some(b't') => out.push('\t'),
                    Some(b'u') => {
                        let hex = bytes
                            .get(*at + 1..*at + 5)
                            .and_then(|h| std::str::from_utf8(h).ok())
                            .and_then(|h| u32::from_str_radix(h, 16).ok())
                            .ok_or_else(|| err_at(bytes, *at, "bad \\u escape"))?;
                        out.push(char::from_u32(hex).unwrap_or('\u{fffd}'));
                        *at += 4;
                    }
                    _ => return Err(err_at(bytes, *at, "bad escape")),
                }
                *at += 1;
            }
            Some(&b) if b < 0x80 => {
                out.push(b as char);
                *at += 1;
            }
            Some(_) => {
                // Multi-byte UTF-8: copy the full character.
                let s = std::str::from_utf8(&bytes[*at..])
                    .map_err(|_| err_at(bytes, *at, "invalid UTF-8"))?;
                let c = s.chars().next().unwrap();
                out.push(c);
                *at += c.len_utf8();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_typical_manifest_shapes() {
        let j = parse(
            r#"{ "name": "x", "n": -1.5e2, "ok": true, "none": null,
                 "arr": [1, 2.5, "s"], "obj": {"a": {"b": []}} }"#,
        )
        .unwrap();
        assert_eq!(j.get("name").unwrap().as_str().unwrap(), "x");
        assert_eq!(j.get("n").unwrap().as_f64().unwrap(), -150.0);
        assert_eq!(j.get("ok").unwrap(), &Json::Bool(true));
        assert_eq!(j.get("arr").unwrap().as_arr().unwrap().len(), 3);
        assert!(j.get("obj").unwrap().get("a").unwrap().get("b").is_some());
    }

    #[test]
    fn handles_escapes() {
        let j = parse(r#""a\"b\\c\ndA""#).unwrap();
        assert_eq!(j.as_str().unwrap(), "a\"b\\c\ndA");
    }

    #[test]
    fn rejects_garbage() {
        assert!(parse("{").is_err());
        assert!(parse("[1,]").is_err());
        assert!(parse("{\"a\" 1}").is_err());
        assert!(parse("12 34").is_err());
        assert!(parse("{'a': 1}").is_err());
    }
}
