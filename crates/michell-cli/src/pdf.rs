//! Minimal zero-dependency PDF writer: pages on a bottom-left, y-up point
//! grid; Helvetica / Helvetica-Bold text; vector strokes and fills;
//! `DeviceRGB` image XObjects (compressed with the same stored-deflate zlib
//! stream as the PNG writer); and internal *go-to* link annotations, so a
//! table of contents can hyperlink to per-row pages.
//!
//! Only the slice of PDF needed for the sweep report is implemented, but it
//! produces a fully conformant single-file document (header, body, xref,
//! trailer) that any reader opens.

use crate::png::zlib_stored;

/// An RGB8 raster embedded once and drawn by reference.
struct Image {
    w: usize,
    h: usize,
    rgb: Vec<u8>,
}

/// A clickable rectangle that jumps to another page (fit-to-page).
struct Link {
    rect: [f64; 4],
    target: usize,
}

/// One page: a content byte stream in PDF operators, plus the images it uses
/// and its link annotations.
pub struct Page {
    pub w: f64,
    pub h: f64,
    content: Vec<u8>,
    used_images: Vec<usize>,
    links: Vec<Link>,
}

impl Page {
    fn op(&mut self, s: &str) {
        self.content.extend_from_slice(s.as_bytes());
        self.content.push(b'\n');
    }

    /// Stroke a line in points, `rgb` in [0, 1].
    pub fn line(&mut self, x0: f64, y0: f64, x1: f64, y1: f64, width: f64, rgb: [f64; 3]) {
        self.op(&format!(
            "{} {} {} RG {width:.3} w {x0:.3} {y0:.3} m {x1:.3} {y1:.3} l S",
            rgb[0], rgb[1], rgb[2]
        ));
    }

    /// Stroke (and optionally fill) an axis-aligned rectangle.
    #[allow(clippy::too_many_arguments)]
    pub fn rect(
        &mut self,
        x: f64,
        y: f64,
        w: f64,
        h: f64,
        width: f64,
        stroke: Option<[f64; 3]>,
        fill: Option<[f64; 3]>,
    ) {
        if let Some(f) = fill {
            self.op(&format!("{} {} {} rg", f[0], f[1], f[2]));
        }
        if let Some(s) = stroke {
            self.op(&format!("{} {} {} RG {width:.3} w", s[0], s[1], s[2]));
        }
        self.op(&format!("{x:.3} {y:.3} {w:.3} {h:.3} re"));
        self.op(match (stroke.is_some(), fill.is_some()) {
            (true, true) => "B",
            (true, false) => "S",
            (false, true) => "f",
            (false, false) => "n",
        });
    }

    /// Stroke a polyline through `pts`; closes the path if `close`.
    pub fn polyline(&mut self, pts: &[(f64, f64)], width: f64, rgb: [f64; 3], close: bool) {
        if pts.len() < 2 {
            return;
        }
        self.op(&format!(
            "{} {} {} RG {width:.3} w",
            rgb[0], rgb[1], rgb[2]
        ));
        self.op(&format!("{:.3} {:.3} m", pts[0].0, pts[0].1));
        for p in &pts[1..] {
            self.op(&format!("{:.3} {:.3} l", p.0, p.1));
        }
        self.op(if close { "s" } else { "S" });
    }

    /// Fill a closed polygon.
    pub fn fill_poly(&mut self, pts: &[(f64, f64)], rgb: [f64; 3]) {
        if pts.len() < 3 {
            return;
        }
        self.op(&format!("{} {} {} rg", rgb[0], rgb[1], rgb[2]));
        self.op(&format!("{:.3} {:.3} m", pts[0].0, pts[0].1));
        for p in &pts[1..] {
            self.op(&format!("{:.3} {:.3} l", p.0, p.1));
        }
        self.op("f");
    }

    /// Left-aligned text with its baseline at `(x, y)`.
    pub fn text(&mut self, x: f64, y: f64, size: f64, rgb: [f64; 3], bold: bool, s: &str) {
        let font = if bold { "/F2" } else { "/F1" };
        self.op(&format!(
            "BT {font} {size:.2} Tf {} {} {} rg {x:.3} {y:.3} Td ({}) Tj ET",
            rgb[0],
            rgb[1],
            rgb[2],
            escape(s)
        ));
    }

    /// Text ending at `x` (right-aligned).
    pub fn text_right(&mut self, x: f64, y: f64, size: f64, rgb: [f64; 3], bold: bool, s: &str) {
        let w = text_width(s, size, bold);
        self.text(x - w, y, size, rgb, bold, s);
    }

    /// Text centered on `x`.
    pub fn text_center(&mut self, x: f64, y: f64, size: f64, rgb: [f64; 3], bold: bool, s: &str) {
        let w = text_width(s, size, bold);
        self.text(x - 0.5 * w, y, size, rgb, bold, s);
    }

    /// Draw image `id` into the rectangle `(x, y, w, h)` (points).
    pub fn image(&mut self, id: usize, x: f64, y: f64, w: f64, h: f64) {
        if !self.used_images.contains(&id) {
            self.used_images.push(id);
        }
        self.op(&format!(
            "q {w:.3} 0 0 {h:.3} {x:.3} {y:.3} cm /Im{id} Do Q"
        ));
    }

    /// Register a clickable rectangle jumping to `target` page index.
    pub fn link(&mut self, rect: [f64; 4], target: usize) {
        self.links.push(Link { rect, target });
    }
}

/// A PDF document under construction.
pub struct Pdf {
    pages: Vec<Page>,
    images: Vec<Image>,
}

impl Default for Pdf {
    fn default() -> Self {
        Pdf::new()
    }
}

impl Pdf {
    pub fn new() -> Pdf {
        Pdf {
            pages: Vec::new(),
            images: Vec::new(),
        }
    }

    /// Add a page of the given size (points) and return its index.
    pub fn add_page(&mut self, w: f64, h: f64) -> usize {
        self.pages.push(Page {
            w,
            h,
            content: Vec::new(),
            used_images: Vec::new(),
            links: Vec::new(),
        });
        self.pages.len() - 1
    }

    pub fn page(&mut self, idx: usize) -> &mut Page {
        &mut self.pages[idx]
    }

    /// Embed an RGB8 image (row-major, top row first) and return its id.
    pub fn add_image(&mut self, w: usize, h: usize, rgb: Vec<u8>) -> usize {
        assert_eq!(rgb.len(), 3 * w * h, "image buffer size");
        self.images.push(Image { w, h, rgb });
        self.images.len() - 1
    }

    /// Serialize the whole document.
    pub fn finish(&self) -> Vec<u8> {
        // Object numbering. 1 catalog, 2 pages tree, 3/4 fonts, then one per
        // image, then per page a page object, a content stream, and one
        // annotation object per link.
        let ni = self.images.len();
        let img_no = |k: usize| 5 + k;
        let mut next = 5 + ni;
        struct PagePlan {
            page_no: usize,
            content_no: usize,
            annot_nos: Vec<usize>,
        }
        let mut plans: Vec<PagePlan> = Vec::with_capacity(self.pages.len());
        for pg in &self.pages {
            let page_no = next;
            next += 1;
            let content_no = next;
            next += 1;
            let mut annot_nos = Vec::with_capacity(pg.links.len());
            for _ in &pg.links {
                annot_nos.push(next);
                next += 1;
            }
            plans.push(PagePlan {
                page_no,
                content_no,
                annot_nos,
            });
        }
        let total = next - 1;
        let mut objs: Vec<Vec<u8>> = vec![Vec::new(); total];
        let mut put = |no: usize, body: Vec<u8>| objs[no - 1] = body;

        // Catalog + page tree.
        put(1, b"<< /Type /Catalog /Pages 2 0 R >>".to_vec());
        let kids: String = plans
            .iter()
            .map(|p| format!("{} 0 R", p.page_no))
            .collect::<Vec<_>>()
            .join(" ");
        put(
            2,
            format!("<< /Type /Pages /Kids [{kids}] /Count {} >>", plans.len()).into_bytes(),
        );
        put(
            3,
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>"
                .to_vec(),
        );
        put(
            4,
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>"
                .to_vec(),
        );

        // Images.
        for (k, im) in self.images.iter().enumerate() {
            let z = zlib_stored(&im.rgb);
            let mut o = format!(
                "<< /Type /XObject /Subtype /Image /Width {} /Height {} /ColorSpace /DeviceRGB \
                 /BitsPerComponent 8 /Filter /FlateDecode /Length {} >>\nstream\n",
                im.w,
                im.h,
                z.len()
            )
            .into_bytes();
            o.extend_from_slice(&z);
            o.extend_from_slice(b"\nendstream");
            put(img_no(k), o);
        }

        // Pages, content streams, annotations.
        for (pg, plan) in self.pages.iter().zip(&plans) {
            let xobjects: String = pg
                .used_images
                .iter()
                .map(|&k| format!("/Im{k} {} 0 R", img_no(k)))
                .collect::<Vec<_>>()
                .join(" ");
            let annots: String = if plan.annot_nos.is_empty() {
                String::new()
            } else {
                let refs: Vec<String> =
                    plan.annot_nos.iter().map(|n| format!("{n} 0 R")).collect();
                format!(" /Annots [{}]", refs.join(" "))
            };
            put(
                plan.page_no,
                format!(
                    "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {:.2} {:.2}] \
                     /Resources << /Font << /F1 3 0 R /F2 4 0 R >> /XObject << {xobjects} >> >> \
                     /Contents {} 0 R{annots} >>",
                    pg.w, pg.h, plan.content_no
                )
                .into_bytes(),
            );
            let z = zlib_stored(&pg.content);
            let mut c = format!("<< /Length {} /Filter /FlateDecode >>\nstream\n", z.len())
                .into_bytes();
            c.extend_from_slice(&z);
            c.extend_from_slice(b"\nendstream");
            put(plan.content_no, c);

            for (link, &annot_no) in pg.links.iter().zip(&plan.annot_nos) {
                let [x0, y0, x1, y1] = link.rect;
                let target_no = plans[link.target].page_no;
                put(
                    annot_no,
                    format!(
                        "<< /Type /Annot /Subtype /Link /Rect [{x0:.2} {y0:.2} {x1:.2} {y1:.2}] \
                         /Border [0 0 0] /Dest [{target_no} 0 R /Fit] >>"
                    )
                    .into_bytes(),
                );
            }
        }

        // Assemble the file, recording byte offsets for the xref table.
        let mut out: Vec<u8> = Vec::new();
        out.extend_from_slice(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n");
        let mut offsets = vec![0usize; total + 1];
        for (i, body) in objs.iter().enumerate() {
            let no = i + 1;
            offsets[no] = out.len();
            out.extend_from_slice(format!("{no} 0 obj\n").as_bytes());
            out.extend_from_slice(body);
            out.extend_from_slice(b"\nendobj\n");
        }
        let xref_off = out.len();
        out.extend_from_slice(format!("xref\n0 {}\n", total + 1).as_bytes());
        out.extend_from_slice(b"0000000000 65535 f \n");
        for &off in &offsets[1..=total] {
            out.extend_from_slice(format!("{off:010} 00000 n \n").as_bytes());
        }
        out.extend_from_slice(
            format!(
                "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{xref_off}\n%%EOF\n",
                total + 1
            )
            .as_bytes(),
        );
        out
    }
}

/// Escape a string for a PDF literal `( ... )`, dropping non-WinAnsi bytes to
/// a question mark so the text tools never emit an invalid stream.
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 4);
    for ch in s.chars() {
        match ch {
            '(' => out.push_str("\\("),
            ')' => out.push_str("\\)"),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) >= 32 && (c as u32) < 127 => out.push(c),
            // Common typographic glyphs the report might carry.
            '\u{2013}' | '\u{2014}' => out.push('-'),
            '\u{00b0}' => out.push_str("\\260"), // degree sign (WinAnsi 0xB0)
            '\u{00b1}' => out.push_str("\\261"), // plus-minus (WinAnsi 0xB1)
            '\u{00b2}' => out.push_str("\\262"), // superscript two
            '\u{00b3}' => out.push_str("\\263"), // superscript three
            '\u{2207}' => out.push_str("del"),
            _ => out.push('?'),
        }
    }
    out
}

/// Width of `s` in points at the given font size, from the Adobe Core-14
/// Helvetica metrics (units of 1/1000 em). Non-printable/degree glyphs use
/// the space width as a safe estimate.
pub fn text_width(s: &str, size: f64, bold: bool) -> f64 {
    let table = if bold { &HELV_BOLD_W } else { &HELV_W };
    let mut units = 0u32;
    for ch in s.chars() {
        let c = ch as u32;
        units += if (32..127).contains(&c) {
            table[(c - 32) as usize] as u32
        } else {
            table[0] as u32
        };
    }
    units as f64 / 1000.0 * size
}

/// Helvetica glyph advance widths for ASCII 32..=126.
static HELV_W: [u16; 95] = [
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556,
    556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556, 1015, 667, 667, 722, 722, 667,
    611, 778, 722, 278, 500, 667, 556, 833, 722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667,
    667, 611, 278, 278, 278, 469, 556, 333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500,
    222, 833, 556, 556, 556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584,
];

/// Helvetica-Bold glyph advance widths for ASCII 32..=126.
static HELV_BOLD_W: [u16; 95] = [
    278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556,
    556, 556, 556, 556, 556, 556, 556, 333, 333, 584, 584, 584, 611, 975, 722, 722, 722, 722, 667,
    611, 778, 722, 278, 556, 722, 611, 833, 722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667,
    667, 611, 333, 278, 333, 584, 556, 333, 556, 611, 556, 611, 556, 333, 611, 611, 278, 278, 556,
    278, 889, 611, 611, 611, 611, 389, 556, 333, 611, 556, 778, 556, 556, 500, 389, 280, 389, 584,
];
