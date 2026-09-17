//! Minimal zero-dependency PDF 1.4 writer: multi-page documents with vector
//! line art, Helvetica text, embedded RGB raster images, and internal
//! (`GoTo`) links — just enough to build the sweep report, in the same
//! spirit as [`crate::png`]'s hand-rolled PNG encoder. No external crate: the
//! image streams reuse that encoder's zlib-with-stored-blocks trick (valid
//! `/FlateDecode` framing, no compression algorithm required), and text uses
//! the built-in Helvetica base-14 font, so nothing needs embedding.
//!
//! ## Coordinate convention
//!
//! PDF's default user space has **y up**, origin at the page's bottom-left,
//! in points (1/72 inch). Every drawing helper here takes and returns that
//! same convention; callers doing science plots with y-down data (image
//! rows, ship z-down) flip the sign themselves at the point of drawing, so
//! this module stays a dumb, general PDF layer with no plotting opinions.

use crate::png;

/// One page under construction: its content-stream operators, the
/// resources (images, link targets) it references, and its size.
pub struct Page {
    pub width: f64,
    pub height: f64,
    ops: String,
    images: Vec<(String, u32)>, // (resource name, object id)
    links: Vec<Link>,
}

struct Link {
    /// Rectangle in page points: [x0 y0 x1 y1].
    rect: [f64; 4],
    /// Target page index (0-based, into the document's page list).
    target: usize,
}

impl Page {
    pub fn new(width: f64, height: f64) -> Self {
        Page {
            width,
            height,
            ops: String::new(),
            images: Vec::new(),
            links: Vec::new(),
        }
    }

    /// A stroked polyline through `pts` (page points), `width` in points.
    pub fn polyline(&mut self, pts: &[[f64; 2]], rgb: [f64; 3], width: f64) {
        if pts.len() < 2 {
            return;
        }
        self.color_stroke(rgb);
        self.ops
            .push_str(&format!("{width:.2} w {} J {} j\n", 1, 1));
        self.ops.push_str(&format!("{:.2} {:.2} m\n", pts[0][0], pts[0][1]));
        for p in &pts[1..] {
            self.ops.push_str(&format!("{:.2} {:.2} l\n", p[0], p[1]));
        }
        self.ops.push_str("S\n");
    }

    /// A closed, filled polygon (e.g. a hull section).
    pub fn filled_polygon(&mut self, pts: &[[f64; 2]], fill: [f64; 3]) {
        if pts.len() < 3 {
            return;
        }
        self.color_fill(fill);
        self.ops.push_str(&format!("{:.2} {:.2} m\n", pts[0][0], pts[0][1]));
        for p in &pts[1..] {
            self.ops.push_str(&format!("{:.2} {:.2} l\n", p[0], p[1]));
        }
        self.ops.push_str("h f\n");
    }

    /// An unfilled rectangle outline, `[x, y, w, h]` (y up, bottom-left origin).
    pub fn rect_stroke(&mut self, rect: [f64; 4], rgb: [f64; 3], width: f64) {
        self.color_stroke(rgb);
        self.ops.push_str(&format!("{width:.2} w\n"));
        self.ops.push_str(&format!(
            "{:.2} {:.2} {:.2} {:.2} re S\n",
            rect[0], rect[1], rect[2], rect[3]
        ));
    }

    /// A filled rectangle, `[x, y, w, h]`.
    pub fn rect_fill(&mut self, rect: [f64; 4], rgb: [f64; 3]) {
        self.color_fill(rgb);
        self.ops.push_str(&format!(
            "{:.2} {:.2} {:.2} {:.2} re f\n",
            rect[0], rect[1], rect[2], rect[3]
        ));
    }

    /// Left-aligned Helvetica text with its baseline at `(x, y)`.
    pub fn text(&mut self, x: f64, y: f64, size: f64, rgb: [f64; 3], s: &str) {
        self.color_fill(rgb);
        self.ops.push_str("BT\n");
        self.ops.push_str(&format!("/F1 {size:.2} Tf\n"));
        self.ops.push_str(&format!("{x:.2} {y:.2} Td\n"));
        self.ops.push_str(&format!("({}) Tj\n", escape_pdf_string(s)));
        self.ops.push_str("ET\n");
    }

    /// Right-aligned Helvetica text ending at `x`.
    pub fn text_right(&mut self, x: f64, y: f64, size: f64, rgb: [f64; 3], s: &str) {
        let w = text_width(s, size);
        self.text(x - w, y, size, rgb, s);
    }

    /// Width in points of `s` set at `size` — for callers doing their own
    /// layout (column alignment, centring) before calling [`Page::text`].
    pub fn text_width(s: &str, size: f64) -> f64 {
        text_width(s, size)
    }

    /// Place an RGB image (registered via [`Document::image`]) filling
    /// `[x, y, w, h]` in page points.
    pub fn image(&mut self, id: u32, rect: [f64; 4]) {
        let name = format!("Im{id}");
        if !self.images.iter().any(|(n, _)| n == &name) {
            self.images.push((name.clone(), id));
        }
        self.ops.push_str("q\n");
        self.ops.push_str(&format!(
            "{:.2} 0 0 {:.2} {:.2} {:.2} cm\n",
            rect[2], rect[3], rect[0], rect[1]
        ));
        self.ops.push_str(&format!("/{name} Do\n"));
        self.ops.push_str("Q\n");
    }

    /// A clickable rectangle jumping to `target` (a page index into the
    /// document's page list, 0-based) when the report is opened in a
    /// PDF viewer.
    pub fn link_to_page(&mut self, rect: [f64; 4], target: usize) {
        self.links.push(Link { rect, target });
    }

    fn color_stroke(&mut self, rgb: [f64; 3]) {
        self.ops
            .push_str(&format!("{:.3} {:.3} {:.3} RG\n", rgb[0], rgb[1], rgb[2]));
    }

    fn color_fill(&mut self, rgb: [f64; 3]) {
        self.ops
            .push_str(&format!("{:.3} {:.3} {:.3} rg\n", rgb[0], rgb[1], rgb[2]));
    }
}

/// A document under construction: a list of pages, built up with
/// [`Document::page`], plus registered images. [`Document::finish`] renders
/// the whole object graph and returns the PDF bytes.
#[derive(Default)]
pub struct Document {
    pages: Vec<Page>,
    /// Registered images: `(object id, width, height, rgb bytes)`. IDs are
    /// pre-allocated by [`Document::image`] so a page can reference an image
    /// before the document is finished.
    images: Vec<(u32, usize, usize, Vec<u8>)>,
    next_id: u32,
}

impl Document {
    pub fn new() -> Self {
        Document {
            next_id: 1,
            ..Default::default()
        }
    }

    /// Register an RGB8 image (row-major, row 0 at the top — the usual
    /// raster convention; PDF's own image space also starts at the top-left
    /// regardless of the page's y-up user space, so no flip is needed here)
    /// and return its id for [`Page::image`].
    pub fn image(&mut self, width: usize, height: usize, rgb: Vec<u8>) -> u32 {
        assert_eq!(rgb.len(), 3 * width * height, "pixel buffer size");
        let id = self.alloc_id();
        self.images.push((id, width, height, rgb));
        id
    }

    fn alloc_id(&mut self) -> u32 {
        let id = self.next_id;
        self.next_id += 1;
        id
    }

    /// Append a finished page. Returns its index for [`Page::link_to_page`]
    /// targets computed ahead of time (e.g. an index page linking to the
    /// detail pages that follow it).
    pub fn push_page(&mut self, page: Page) -> usize {
        self.pages.push(page);
        self.pages.len() - 1
    }

    pub fn page_count(&self) -> usize {
        self.pages.len()
    }

    /// Render the whole document to PDF bytes.
    pub fn finish(mut self) -> Vec<u8> {
        // Object id plan: catalog, pages-root, font, one page-tree leaf +
        // one content stream per page (ids allocated here), and the image
        // XObjects (ids already allocated by `image()`).
        let catalog_id = self.alloc_id();
        let pages_root_id = self.alloc_id();
        let font_id = self.alloc_id();
        let page_ids: Vec<u32> = self.pages.iter().map(|_| 0).collect();
        let mut page_ids = page_ids;
        let mut content_ids = vec![0u32; self.pages.len()];
        for i in 0..self.pages.len() {
            page_ids[i] = self.alloc_id();
            content_ids[i] = self.alloc_id();
        }

        let mut objs: Vec<(u32, Vec<u8>)> = Vec::new();

        objs.push((
            catalog_id,
            format!("<< /Type /Catalog /Pages {pages_root_id} 0 R >>").into_bytes(),
        ));

        let kids: String = page_ids
            .iter()
            .map(|id| format!("{id} 0 R"))
            .collect::<Vec<_>>()
            .join(" ");
        objs.push((
            pages_root_id,
            format!(
                "<< /Type /Pages /Kids [{kids}] /Count {} >>",
                self.pages.len()
            )
            .into_bytes(),
        ));

        objs.push((
            font_id,
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>"
                .to_vec(),
        ));

        // Image XObjects: RGB8 raw, wrapped as a valid (uncompressed)
        // zlib/FlateDecode stream via png's stored-block helper.
        for (id, w, h, rgb) in &self.images {
            let z = png::zlib_store(rgb);
            let mut body = format!(
                "<< /Type /XObject /Subtype /Image /Width {w} /Height {h} \
                 /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode \
                 /Length {} >>\nstream\n",
                z.len()
            )
            .into_bytes();
            body.extend_from_slice(&z);
            body.extend_from_slice(b"\nendstream");
            objs.push((*id, body));
        }

        for (i, page) in self.pages.iter().enumerate() {
            let mut res = format!(
                "<< /Font << /F1 {font_id} 0 R >>"
            );
            if !page.images.is_empty() {
                res.push_str(" /XObject <<");
                for (name, id) in &page.images {
                    res.push_str(&format!(" /{name} {id} 0 R"));
                }
                res.push_str(" >>");
            }
            res.push_str(" >>");

            let annots = if page.links.is_empty() {
                String::new()
            } else {
                let mut s = String::from(" /Annots [");
                for link in &page.links {
                    let target_page = page_ids[link.target];
                    s.push_str(&format!(
                        "<< /Type /Annot /Subtype /Link /Rect [{:.2} {:.2} {:.2} {:.2}] \
                         /Border [0 0 0] /Dest [{target_page} 0 R /Fit] >> ",
                        link.rect[0], link.rect[1], link.rect[2], link.rect[3]
                    ));
                }
                s.push(']');
                s
            };

            objs.push((
                page_ids[i],
                format!(
                    "<< /Type /Page /Parent {pages_root_id} 0 R \
                     /MediaBox [0 0 {:.2} {:.2}] /Resources {res} \
                     /Contents {} 0 R{annots} >>",
                    page.width, page.height, content_ids[i]
                )
                .into_bytes(),
            ));

            let stream = page.ops.as_bytes();
            let mut body = format!("<< /Length {} >>\nstream\n", stream.len()).into_bytes();
            body.extend_from_slice(stream);
            body.extend_from_slice(b"\nendstream");
            objs.push((content_ids[i], body));
        }

        objs.sort_by_key(|(id, _)| *id);
        write_pdf(&objs, catalog_id)
    }
}

/// Serialize the finished object list (already `(id, body)` pairs, body
/// being the dictionary/stream content between `N 0 obj` and `endobj`) as a
/// complete PDF file: header, objects, xref table, trailer.
fn write_pdf(objs: &[(u32, Vec<u8>)], root_id: u32) -> Vec<u8> {
    let mut out = Vec::new();
    out.extend_from_slice(b"%PDF-1.4\n%\xE2\xE3\xCF\xD3\n");
    let max_id = objs.iter().map(|(id, _)| *id).max().unwrap_or(0);
    let mut offsets = vec![0u64; max_id as usize + 1];
    for (id, body) in objs {
        offsets[*id as usize] = out.len() as u64;
        out.extend_from_slice(format!("{id} 0 obj\n").as_bytes());
        out.extend_from_slice(body);
        out.extend_from_slice(b"\nendobj\n");
    }
    let xref_start = out.len();
    out.extend_from_slice(format!("xref\n0 {}\n", max_id + 1).as_bytes());
    out.extend_from_slice(b"0000000000 65535 f \n");
    for id in 1..=max_id {
        if offsets[id as usize] == 0 && !objs.iter().any(|(oid, _)| *oid == id) {
            out.extend_from_slice(b"0000000000 00000 f \n");
        } else {
            out.extend_from_slice(format!("{:010} 00000 n \n", offsets[id as usize]).as_bytes());
        }
    }
    out.extend_from_slice(
        format!(
            "trailer\n<< /Size {} /Root {root_id} 0 R >>\nstartxref\n{xref_start}\n%%EOF",
            max_id + 1
        )
        .as_bytes(),
    );
    out
}

fn escape_pdf_string(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '(' => out.push_str("\\("),
            ')' => out.push_str("\\)"),
            '\\' => out.push_str("\\\\"),
            // WinAnsiEncoding matches Latin-1 for the printable ASCII range
            // this tool ever emits (labels, numbers); anything outside it
            // is replaced rather than mis-rendered.
            c if (c as u32) < 128 => out.push(c),
            _ => out.push('?'),
        }
    }
    out
}

/// Approximate Helvetica advance widths (1/1000 em, the standard PDF Base-14
/// AFM metrics) for the printable ASCII range, indexed from `' '` (32).
/// Covers everything this tool ever sets: labels, numbers, punctuation.
#[rustfmt::skip]
const HELVETICA_WIDTHS: [u16; 95] = [
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278, // ' ' .. '/'
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556, // '0' .. '?'
    1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778, // '@' .. 'O'
    667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556, // 'P' .. '_'
    333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556, // '`' .. 'o'
    556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334, 260, 334, 584,      // 'p' .. '~'
];

fn text_width(s: &str, size: f64) -> f64 {
    let units: u32 = s
        .chars()
        .map(|c| {
            let i = c as u32;
            if (32..127).contains(&i) {
                HELVETICA_WIDTHS[(i - 32) as usize] as u32
            } else {
                556 // a plausible fallback advance for anything unmapped
            }
        })
        .sum();
    units as f64 * size / 1000.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_page_document_has_valid_structure() {
        let mut doc = Document::new();
        let mut page = Page::new(200.0, 100.0);
        page.text(10.0, 50.0, 12.0, [0.0, 0.0, 0.0], "hello");
        page.polyline(&[[0.0, 0.0], [100.0, 50.0]], [0.2, 0.4, 0.8], 1.0);
        doc.push_page(page);
        let bytes = doc.finish();
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.starts_with("%PDF-1.4"));
        assert!(s.contains("/Type /Catalog"));
        assert!(s.contains("/Type /Pages"));
        assert!(s.contains("/Type /Page"));
        assert!(s.contains("(hello) Tj"));
        assert!(s.contains("xref"));
        assert!(s.contains("trailer"));
        assert!(s.trim_end().ends_with("%%EOF"));
    }

    #[test]
    fn image_round_trips_through_zlib_store() {
        let mut doc = Document::new();
        let rgb = vec![255u8, 0, 0, 0, 255, 0, 0, 0, 255, 255, 255, 0]; // 2x2
        let id = doc.image(2, 2, rgb);
        let mut page = Page::new(100.0, 100.0);
        page.image(id, [0.0, 0.0, 100.0, 100.0]);
        doc.push_page(page);
        let bytes = doc.finish();
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.contains("/Subtype /Image"));
        assert!(s.contains("/FlateDecode"));
        assert!(s.contains(&format!("/Im{id} Do")));
    }

    #[test]
    fn index_page_link_targets_a_later_page() {
        let mut doc = Document::new();
        doc.push_page(Page::new(100.0, 100.0));
        let mut index = Page::new(100.0, 100.0);
        index.link_to_page([0.0, 0.0, 10.0, 10.0], 0);
        doc.push_page(index);
        let bytes = doc.finish();
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.contains("/Subtype /Link"));
        assert!(s.contains("/Dest ["));
    }

    #[test]
    fn text_width_is_positive_and_scales_with_size() {
        let w10 = text_width("Hello, World!", 10.0);
        let w20 = text_width("Hello, World!", 20.0);
        assert!(w10 > 0.0);
        assert!((w20 - 2.0 * w10).abs() < 1e-9);
    }
}
