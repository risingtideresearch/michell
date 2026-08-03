# Venue decision memo — Paper A

**Decision requested from Rob Story.** Policies checked 2 August 2026 against
the journals' current official author materials. The manuscript should not yet
be reframed or anonymized.

## A. Journal of Ship Research

**Fit.** This is the most direct subject-matter home. JSR explicitly seeks
highly technical applied research in hydrodynamics that advances ship and ocean
engineering while reaching applied mathematics and numerical analysis. The
paper can lead with the error-estimate honesty result: a reproducible gap
between reported and actual marcher error, the structural correction, and the
new endpoint-pair interference diagnostics. Its intellectual lineage is native
to JSR: Keller--Ahluwalia (1976), de Sendagorta--Grases (1988), and
Doctors--Beck (1987) all appeared there.

**Reframing effort and review culture.** Low. Preserve the hydrodynamic
motivation, elevate the published Wigley anchor, and make the operational value
of trustworthy error estimates explicit. Based on JSR's stated scope and
specialist editorial board, reviewers are likely to emphasize thin-ship
physics, physical interpretation, validation against ship-wave precedents, and
whether the engineering claim is fairly bounded. This is an inference, not a
published review rubric. The current SNAME template asks for a comprehensive
critical review, uncertainty in computations, a title of at most ten words,
an abstract of at most 250 words, and typically 3,000--8,000 words with no more
than ten illustrations.

**Blinding and AI.** SNAME's public JSR materials do not explicitly identify a
single- or double-blind model. The public template places author names and
affiliations in the manuscript, which suggests a named submission, but that is
not an explicit blinding rule. Neither the current journal template nor the
2024 SNAME publication-ethics statement states a generative-AI disclosure
policy. Both points should be confirmed with SNAME publications staff before
submission; substantive AI drafting should be disclosed regardless.

**Top risk.** A domain referee may regard the endpoint/Bickley composition as
too narrow unless the unobtained de Sendagorta--Grases full text is reviewed
and the distinction from earlier endpoint expansions is made airtight.

## B. Applied Numerical Mathematics

**Fit.** APNUM explicitly publishes high-quality computational mathematics and
applications in fluid dynamics and engineering. Here the paper should lead as
an oscillatory-quadrature method: exact B-spline endpoint reduction, analytic
continuation to Bickley kernels, and frequency-scaled Gaussian-contour
quadrature with cost independent of Froude frequency. Motygin (2017) supplies
an APNUM ship-wave precedent, and Huybrechs--Vandewalle supplies the numerical
steepest-descent lineage.

**Reframing effort and review culture.** Moderate to high. The introduction
would need to address a broad computational-mathematics audience; the paper
would benefit from separating the reusable quadrature theorem from Michell
specifics and comparing more directly with Filon, Levin, and numerical
steepest-descent alternatives. APNUM says both rigorous and heuristic full
papers are acceptable, but expects a complete, self-contained original
contribution. Its process uses an editor suitability screen followed by at
least two independent experts, so reviewers are likely to press hardest on
generality, error certification, conditioning on the imaginary axis, and
comparative numerical evidence.

**Blinding and AI.** APNUM explicitly uses **single-blind** review; author
anonymization is therefore not required. Elsevier requires a separate
generative-AI declaration immediately before the references, naming the tool,
purpose, human review, and author responsibility. AI-assisted code used in the
research must also be described in the methods. The policy was updated June
2026.

**Top risk.** The present evidence is one hull family within one physical
integral; numerical-analysis referees may require broader test classes and
stronger end-to-end error certification before accepting the method as a
general contribution.

## Recommendation

Choose **Journal of Ship Research**. It offers the sharper audience fit, the
smallest honest reframing, and the strongest venue-specific lineage for the
paper's validation and estimator story. Choose APNUM instead only if we are
prepared to add a broader quadrature test matrix and recast the paper around a
general computational-mathematics contribution.

Official sources: [SNAME journal scope](https://sname.org/journals-transactions),
[SNAME author resources and policies](https://sname.org/author-opportunities),
[SNAME journal template](https://www.sname.org/sites/default/files/2021-01/SNAME%20Journals%20Paper%20Template_0.pdf),
[SNAME publication-ethics statement](https://sname.org/sites/default/files/2024-04/SNAME%20Journals%20Publication%20Ethics%20and%20Publication%20Malpractice%20Statement_2024final.pdf),
[APNUM guide for authors](https://www.sciencedirect.com/journal/applied-numerical-mathematics/publish/guide-for-authors),
and [Elsevier generative-AI policy](https://www.elsevier.com/about/policies-and-standards/generative-ai-policies-for-journals).
