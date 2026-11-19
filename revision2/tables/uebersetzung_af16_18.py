# -*- coding: utf-8 -*-
"""Single source of the proposed changes to the results workbooks of the ChEMBL benchmark
(Additional files 16, 17, 18). Proposed 06.11.2026, used by translate_workbooks.py after approval."""

# 1. whole-cell replacements (sheet, column index from 0, old, new), applied to every cell equal to old
CELL = [
    ('Candidates', 11, 'dazwischen', 'intermediate'),
    ('Candidates', 12, 'Diversitaet weder kongener noch divers', 'diversity neither congeneric nor diverse'),
    ('Candidates', 12, 'Groesse ausserhalb 30-60', 'size outside 30-60'),
    ('Candidates', 12, 'Groesse >= 100 (nur Elternpool Groessenvariation)', 'size >= 100 (parent pool for the size subsets only)'),
    ('Selection', 3, 'Elternpool-Kandidat abgelehnt', 'rejected as parent pool candidate'),
]

# 2. replacements inside the decision text of Selection_log (column 2), applied in this order
LOG = [
    ('MUSS GELESEN WERDEN fuer ', 'MUST BE READ for '),
    ('SELECTED im 2. Durchgang (Zielwiederholung) ', 'SELECTED in pass 2 (repeated target allowed) '),
    ('SELECTED im 2. Durchgang ', 'SELECTED in pass 2 '),
    ('Publikation ausgeschlossen: review', 'publication excluded: review'),
    ('Publikation ausgeschlossen: keine_serie', 'publication excluded: values compiled from other publications'),
    ('Publikation ausgeschlossen: vorselektiert', 'publication excluded: compounds preselected by activity (screening hits)'),
    ('Publikation ausgeschlossen: methoden', 'publication excluded: method or assay paper without own series'),
    ("'Enzym'", "'Enzyme IC50'"), ("'Rezeptor'", "'Receptor Ki'"), ("'Zytotox'", "'Cytotoxicity'"),
    ("'kongener'", "'congeneric'"), ("'divers'", "'diverse'"),
    ('SELECTED Enzym / ', 'SELECTED Enzyme IC50 / '), ('SELECTED Rezeptor / ', 'SELECTED Receptor Ki / '),
    ('SELECTED Zytotox / ', 'SELECTED Cytotoxicity / '),
    (' / kongener / ', ' / congeneric / '), (' / divers / ', ' / diverse / '),
]

# 3. replacements inside single cells (sheet, old part, new part, expected number of cells per workbook)
PART = [
    ('Selection', '(Autoren: "Active compounds', '(authors: "Active compounds'),
    ('DS21', '(Autoren: "Active compounds', '(authors: "Active compounds'),
    ('DS22', '(Autoren: "Active compounds', '(authors: "Active compounds'),
    ('README', 'Honest performance estimation for small-sample QSAR (Journal of Cheminformatics, revision 2).',
               'Honest performance estimation for small-sample QSAR (Journal of Cheminformatics).'),
]

# 4. single sentences
EXCLUDED_NOTE = ('Each decision was taken from the article text and is documented with the source passage in the table of control and reference compounds per data set (Kontrollsubstanzen_je_Datensatz, revision archive).',
                 'Each decision was taken from the article text and is documented with the source passage in Additional file 14 (Table B1).')

# 5. sheet Deviations: the row 'open' is replaced
DEV_OPEN = ['open', 'Sensitivity analysis with the median rule not yet computed.', 'Requires a separate compute run with permutations, decision EW pending.',
            'Reported descriptively only, does not enter H1 to H4.']
DEV_MAIN = [
    ['25.10.2026', 'Addendum for the median-rule sensitivity analysis: treatment of values equal to the median and of censored values, seeds analysis code + 300. Frozen before any model with these classes was fitted.',
     'The pre-specification left these points open.', 'Descriptive only, not part of H1 to H4 (Additional file 14, section B8, and Additional file 18).'],
    ['30.10.2026', 'The seed of the bootstrap over data sets in the median-rule analysis was not specified. Used default_rng([93, 0, 308]).',
     'Gap in the addendum, set analogously to the main analysis (208).', 'Affects only the interval of the H4 effect estimate.'],
]
DEV_ORIGINAL = [
    ['note', 'The sensitivity analysis with the median rule was computed for the 16 series of the main analysis only (Additional file 18).',
     'As fixed in the addendum of 25.10.2026.', '-'],
]
