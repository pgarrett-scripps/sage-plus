// =============================================================================
// FILL THIS IN. Everything project-specific about the manuscript lives here.
//
// This is the single source of truth for the manuscript's identity. paper.typ
// imports it for the PDF and the Word front matter, wordcount.typ counts the
// abstract out of it, and audio/config.py reads the title straight out of this
// file so the narration can never announce a title the paper no longer has.
//
// Nothing below this block should need editing to start a new paper.
// =============================================================================

#let paper-title = "Sage Plus Development and Benchmark Report"

// Short form used on the audiobook cover art. Keep it to a couple of words.
#let paper-wordmark = "sage plus"

// Shown under the wordmark on the cover. \n breaks the line.
#let paper-cover-subtitle = "Repository Technical Report"

#let paper-authors = ()

#let paper-keywords = (
  "proteomics",
  "database search",
  "memory efficiency",
  "entrapment",
  "benchmarking",
  "label-free quantification",
  "phosphorylation",
)

#let paper-date = "September 2026"

// Shown on the audiobook cover under the author line.
#let paper-institution = "Sage Plus"

// The bibliography style. Typst ships CSL styles by name, e.g.
// "american-chemical-society", "ieee", "nature", "apa".
#let paper-bib-style = "nature"

// The abstract. Kept as its own binding, rather than inline in the template
// call, because three separate consumers slice it out of this file by name: the
// Word export path, the word counter (journals cap the abstract separately), and
// the audiobook narrator.
#import "stats.typ": lit, s

#let paper-abstract = [
  Sage was developed by Michael R. Lazear and the upstream Sage contributors.
  Sage Plus is an independently versioned downstream distribution of that
  engine. Its extensions address input support, typed analytical output, memory
  use, execution controls, and modification-aware analysis. This report
  documents those changes and evaluates their consequences by comparing Sage
  #lit("v0.15.0-beta.2") with Sage Plus #lit("v0.1.0-beta.3") using matched
  public spectra, references, and search settings. The evaluation covered
  computational resources, worker scaling, identification agreement, peptide
  entrapment, synthetic phosphorylation, and label-free quantification. Repeated
  searches of selected human and mixed-species files used #s(
    "pilot.PXD001468.rss_reduction",
  ) and #s(
    "pilot.PXD028735.rss_reduction",
  ) percent less peak resident memory with Sage Plus. Runtime favored Sage Plus
  in the repeated public searches but favored Sage with entrapment-expanded
  references. Accepted peptide-spectrum match Jaccard indices ranged from #s(
    "pilot.overlap.min",
  ) to #s("pilot.overlap.max"). Disagreement included both changed assignments
  and identical assignments crossing different confidence thresholds. Mean
  peptide entrapment estimates were similar, with conditional difference
  intervals spanning zero. Matched quantification measured species-ratio
  accuracy, preparation variability, and coverage, while a human-only control
  revealed accepted foreign-species signal in both engines. Neither release
  produced a jointly accepted peptide set in the restricted synthetic
  phosphorylation challenge. The development achieved lower memory use while
  preserving substantial identification agreement on the tested workloads.
  Validation also exposed limits in peptide acceptance and recipient-file
  quantitative confidence, defining priorities for further development.
]

// -----------------------------------------------------------------------------
// Derived values. Nothing to edit below here.
// -----------------------------------------------------------------------------

// Unique affiliations in first-appearance order, so the Word front matter can
// number them the way the PDF template does. Deriving this rather than typing a
// second author line by hand removes the drift that a duplicated list invites.
#let paper-affiliations = {
  let seen = ()
  for a in paper-authors {
    if a.affiliation not in seen { seen.push(a.affiliation) }
  }
  seen
}

#let affiliation-number(affil) = (
  paper-affiliations.position(x => x == affil) + 1
)

// "Ada Lovelace^1, Grace Hopper^2" for the Word front matter, which has no
// template to build an author line for it. Derived rather than retyped, so the
// superscript markers cannot drift out of step with the PDF.
// Wrapped in a code block because a method chain broken across lines after
// `#let x =` would otherwise end at the first newline.
#let paper-author-line = {
  paper-authors
    .map(a => a.name + super(str(affiliation-number(a.affiliation))))
    .join(", ")
}

// Generational and post-nominal suffixes, so a surname lookup does not return
// "III" for "John R. Yates III". Compared case- and period-insensitively.
#let name-suffixes = (
  "jr",
  "sr",
  "ii",
  "iii",
  "iv",
  "v",
  "phd",
  "md",
  "dphil",
  "dsc",
  "esq",
)

// The family name: the last token that is not a suffix. "John R. Yates III"
// gives "Yates". Used for the audiobook artist tag and the cover art.
#let surname-of(full) = {
  let parts = full.split(" ").filter(p => p.trim() != "")
  let i = parts.len() - 1
  while i > 0 and lower(parts.at(i).replace(".", "")) in name-suffixes {
    i -= 1
  }
  parts.at(i)
}

// "Lovelace, Hopper" -- used as the audiobook artist tag and on the cover.
#let paper-surnames = paper-authors.map(a => surname-of(a.name))

// Restated on the Supporting Information title page.
#let si-authors = paper-authors.map(a => a.name).join(", ")
